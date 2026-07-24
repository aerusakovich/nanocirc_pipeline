#!/usr/bin/env python3
"""
quant_chunk_remap.py

Tier1: chunked remap and classify for one sample against every reference
chunk built by quant_chunk_prep.py, with global cross-chunk exclusive
resolution. Runs each chunk in one sequential loop per sample, in a
single process.

Why minimap2 runs per chunk, not once against the whole catalog: one
minimap2 pass against a very large circle reference catalog undercounts
badly. Many similar circle references compete for the aligner's
attention, so a real locus's count can collapse (one case dropped from
5,064 to 481 reads). Running in a small chunk avoids this.

Cross-chunk exclusive resolution (bug fix): just summing each chunk's own
qualified read-locus assignments, with no cross-chunk read-ID dedup,
overcounts genomically nested or overlapping loci that a chunk split
happens to separate. The same physical read can qualify in each chunk
that holds one of the related loci, since a chunk cannot see any other
chunk. Confirmed up to 170x overcounting from this on real data. Fixing
only within one chunk does nothing, the inflation is entirely cross-chunk.
Both parts of the fix are needed:
  1. Track a score per qualifying (read, locus) pair, instead of a plain
     qualify-or-not boolean per chunk (quant_common.py's
     score_chunk_bam/score_chunk_psl, same thresholds as before, now
     keeping the total matched span).
  2. After all chunks are processed, resolve globally: each read keeps
     only its single highest-scoring locus across the whole catalog, not
     just within one chunk.
Since chunks are already processed serially in one loop, this last step
is one in-memory groupby at the end, no extra alignment passes needed.
quant_overlap_rescue.py (tier1.5) and quant_tier3_rescue.py (gene family)
add further, independent checks on top of this fix.

Usage:
    quant_chunk_remap.py \\
        --chunks_dir  chunks/ \\
        --genome_bam  sample1.genome_align.sorted.bam \\
        --sample      sample1 \\
        --minimap2 minimap2 --samtools samtools --pblat pblat \\
        --threads 4
"""
import argparse
from pathlib import Path
from collections import defaultdict

import pysam
import pandas as pd

from quant_common import (
    FLANK, score_chunk_bam, score_chunk_psl,
    run_minimap2_bam, run_pblat, parse_bam_hits, parse_psl_hits,
)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--chunks_dir', required=True, help='Directory of chunk_*.fa/_lengths.tsv (quant_chunk_prep.py)')
    p.add_argument('--genome_bam', required=True, help='This sample\'s genome-wide alignment BAM')
    p.add_argument('--sample', required=True)
    p.add_argument('--minimap2', default='minimap2')
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--pblat', default='pblat')
    p.add_argument('--threads', type=int, default=4)
    return p.parse_args()


def fetch_chunk_candidate_reads(bam_path, chunk_lengths_full: pd.DataFrame, flank=FLANK) -> dict:
    """Primary-alignment-only genomic-window overlap, returning full read
    sequences (needed to build this chunk's own small candidate FASTA)."""
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    bam_chroms = set(bam.references)
    reads = {}
    for row in chunk_lengths_full.itertuples(index=False):
        if row.chrom not in bam_chroms:
            continue
        span = row.end - row.start
        if span <= 2 * flank:
            windows = [(max(0, row.start), row.end)]
        else:
            windows = [(max(0, row.start - flank), row.start + flank),
                       (max(0, row.end - flank), row.end + flank)]
        for w_start, w_end in windows:
            for r in bam.fetch(row.chrom, max(0, w_start), w_end):
                if r.is_unmapped or r.is_secondary or r.is_supplementary:
                    continue
                if not r.query_sequence:
                    continue
                if r.query_name not in reads:
                    reads[r.query_name] = r.query_sequence
    return reads


def run_chunk(chunk_idx: int, chunks_dir: Path, genome_bam: Path, workdir: Path, out_dir: Path,
              minimap2_bin, samtools_bin, pblat_bin, threads=4) -> dict:
    """Returns {(read_id, bsj_id): score} for every (read, locus) pair that
    qualifies in this chunk. A locus belongs to exactly one chunk (chunks
    partition the catalog), so a given (read_id, bsj_id) pair can only ever
    come from one chunk. The caller accumulates these across all chunks
    and resolves per-read exclusivity globally at the end."""
    tag = f"chunk_{chunk_idx:04d}"
    chunk_fa = chunks_dir / f"{tag}.fa"
    chunk_lengths = pd.read_csv(chunks_dir / f"{tag}_lengths.tsv", sep="\t")

    reads = fetch_chunk_candidate_reads(genome_bam, chunk_lengths)

    join_pos_by_ref = dict(zip(chunk_lengths["safe_id"], chunk_lengths["join_pos"]))
    safe_to_bsj = dict(zip(chunk_lengths["safe_id"], chunk_lengths["bsj_id"]))

    stats_out = out_dir / f"{tag}_stats.tsv"
    if not reads:
        pd.DataFrame(columns=["bsj_id", "mm2_count", "total_count"]).to_csv(stats_out, sep="\t", index=False)
        return {}

    reads_fa = workdir / f"{tag}_reads.fa"
    with open(reads_fa, "w") as f:
        for qname, seq in reads.items():
            f.write(f">{qname}\n{seq}\n")

    bam_out = run_minimap2_bam(minimap2_bin, samtools_bin, chunk_fa, reads_fa,
                                workdir / f"{tag}.minimap2.bam", threads)
    mm2_scores = score_chunk_bam(parse_bam_hits(bam_out), join_pos_by_ref)

    psl_path = run_pblat(pblat_bin, chunk_fa, reads_fa, workdir / f"{tag}.psl", threads,
                          min_score=50, min_identity=90)
    psl_scores = score_chunk_psl(parse_psl_hits(psl_path), join_pos_by_ref)

    # A pair qualifying under both hit types keeps its best (larger) score.
    combined = dict(mm2_scores)
    for key, score in psl_scores.items():
        if key not in combined or score > combined[key]:
            combined[key] = score

    chunk_scores = {(qname, safe_to_bsj[rname]): score for (qname, rname), score in combined.items()}

    mm2_per_ref = defaultdict(set)
    total_per_ref = defaultdict(set)
    for qname, rname in mm2_scores:
        mm2_per_ref[rname].add(qname)
    for qname, rname in combined:
        total_per_ref[rname].add(qname)
    stats_rows = [{"bsj_id": safe_to_bsj[ref], "mm2_count": len(mm2_per_ref.get(ref, set())),
                   "total_count": len(total_per_ref.get(ref, set()))} for ref in join_pos_by_ref]
    pd.DataFrame(stats_rows).to_csv(stats_out, sep="\t", index=False)

    return chunk_scores


def resolve_global(all_scores: dict) -> pd.DataFrame:
    """all_scores: {(read_id, bsj_id): score} accumulated across every
    chunk. Keeps only each read's single highest-scoring locus (global
    cross-chunk exclusive resolution) before counting. See the module
    docstring for why this must be global, not per chunk."""
    best_per_read = {}
    for (read_id, bsj_id), score in all_scores.items():
        cur = best_per_read.get(read_id)
        if cur is None or score > cur[1]:
            best_per_read[read_id] = (bsj_id, score)

    counts = defaultdict(int)
    for read_id, (bsj_id, score) in best_per_read.items():
        counts[bsj_id] += 1

    return pd.DataFrame({"bsj_id": list(counts.keys()), "remap_count": list(counts.values())})


def main():
    args = parse_args()
    chunks_dir = Path(args.chunks_dir)
    n_chunks = len(list(chunks_dir.glob("chunk_*_lengths.tsv")))

    workdir = Path(f'{args.sample}_chunk_work')
    workdir.mkdir(parents=True, exist_ok=True)
    out_dir = Path(f'{args.sample}_chunk_results')
    out_dir.mkdir(parents=True, exist_ok=True)

    all_scores = {}
    for chunk_idx in range(n_chunks):
        chunk_scores = run_chunk(chunk_idx, chunks_dir, Path(args.genome_bam), workdir, out_dir,
                                  args.minimap2, args.samtools, args.pblat, args.threads)
        all_scores.update(chunk_scores)

    out_df = resolve_global(all_scores)
    out_path = Path(f'{args.sample}_tier1_counts.tsv')
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"[{args.sample}] tier1 (chunked remap, global cross-chunk exclusive): "
          f"{len(out_df)} loci with >=1 supporting read -> {out_path}")


if __name__ == '__main__':
    main()

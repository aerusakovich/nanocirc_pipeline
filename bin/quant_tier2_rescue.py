#!/usr/bin/env python3
"""
quant_tier2_rescue.py

Tier 2 (targeted rescue): runs BLAT and minimap2 on the full raw fastq
directly against a tiny per-sample reference built from
quant_cross_locus_similarity.py's tier2_candidate list. This skips the
genome-alignment candidate-pool prefilter entirely, since that prefilter
is exactly what starves these hard-to-place loci (they are flagged in
the first place because tier1 undercounted them).

Usage:
    quant_tier2_rescue.py \\
        --flagged_similarity sample1_flagged_loci_similarity.tsv \\
        --deduped_metadata   deduped_metadata.tsv \\
        --reads              sample1.fastq.gz \\
        --sample             sample1 \\
        --minimap2 minimap2 --samtools samtools --pblat pblat --threads 16
"""
import argparse
import gzip
from pathlib import Path

import pandas as pd

from quant_common import MIN_SEGMENT_MATCH, run_minimap2_bam, run_pblat, parse_bam_hits, parse_psl_hits, classify_exclusive
from quant_build_refs import build_circle_references


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--flagged_similarity', required=True)
    p.add_argument('--deduped_metadata', required=True,
                    help='needs bsj_id/chrom/start/strand/blockSizes/blockStarts')
    p.add_argument('--genome_fasta', required=True)
    p.add_argument('--reads', required=True, help='Full raw sample FASTQ')
    p.add_argument('--sample', required=True)
    p.add_argument('--minimap2', default='minimap2')
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--pblat', default='pblat')
    p.add_argument('--threads', type=int, default=16)
    return p.parse_args()


def _fastq_to_fasta(reads_fq, out_fa):
    opener = gzip.open if str(reads_fq).endswith(".gz") else open
    with opener(reads_fq, "rt") as f, open(out_fa, "w") as out:
        while True:
            name = f.readline()
            if not name:
                break
            seq = f.readline()
            f.readline()
            f.readline()
            rid = name[1:].strip().split()[0]
            out.write(f">{rid}\n{seq}")


def tier2_rescue(flagged_similarity: pd.DataFrame, deduped_metadata: pd.DataFrame, genome_fasta,
                  reads_fq: Path, sample: str, minimap2_bin, samtools_bin, pblat_bin, threads=16) -> pd.DataFrame:
    candidate_ids = flagged_similarity.loc[flagged_similarity["tier2_candidate"], "bsj_id"].tolist()
    workdir = Path(f'{sample}_tier2_work')
    workdir.mkdir(parents=True, exist_ok=True)

    if not candidate_ids:
        out = pd.DataFrame(columns=["bsj_id", "tier2_count"])
        out.to_csv(f'{sample}_tier2_counts.tsv', sep="\t", index=False)
        print(f"[{sample}] tier2: no candidates")
        return out

    candidates = deduped_metadata[deduped_metadata["bsj_id"].isin(candidate_ids)]
    ref_lengths = build_circle_references(genome_fasta, candidates, workdir)
    ref_fa = workdir / "circle_refs.fa"

    bam_path = run_minimap2_bam(minimap2_bin, samtools_bin, ref_fa, reads_fq,
                                 workdir / f"{sample}.tier2.minimap2.bam", threads)

    reads_fa = workdir / f"{sample}.reads.fa"
    _fastq_to_fasta(reads_fq, reads_fa)
    psl_path = run_pblat(pblat_bin, ref_fa, reads_fa, workdir / f"{sample}.tier2.psl", threads,
                          min_score=MIN_SEGMENT_MATCH, min_identity=90)

    join_pos_by_ref = dict(zip(ref_lengths["safe_id"], ref_lengths["join_pos"]))
    safe_to_bsj = dict(zip(ref_lengths["safe_id"], ref_lengths["bsj_id"]))
    support = classify_exclusive(parse_bam_hits(bam_path), parse_psl_hits(psl_path), join_pos_by_ref)

    rows = [{"bsj_id": safe_to_bsj[ref], "tier2_count": len(support.get(ref, set()))}
            for ref in join_pos_by_ref]
    out = pd.DataFrame(rows).sort_values("tier2_count", ascending=False)
    out_path = f'{sample}_tier2_counts.tsv'
    out.to_csv(out_path, sep="\t", index=False)
    print(f"[{sample}] tier2: {len(out)} candidates rescued -> {out_path}")
    return out


def main():
    args = parse_args()
    flagged_similarity = pd.read_csv(args.flagged_similarity, sep='\t')
    deduped_metadata = pd.read_csv(args.deduped_metadata, sep='\t')
    tier2_rescue(flagged_similarity, deduped_metadata, args.genome_fasta, Path(args.reads), args.sample,
                 args.minimap2, args.samtools, args.pblat, args.threads)


if __name__ == '__main__':
    main()

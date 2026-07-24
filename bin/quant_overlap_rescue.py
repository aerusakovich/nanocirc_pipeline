#!/usr/bin/env python3
"""
quant_overlap_rescue.py

Tier "1.5": re-arbitrates every genomically-overlapping locus cluster
(a looser criterion than quant_locus_dedup.py's near-duplicate boundary
tolerance: any span overlap, not just near-identical boundaries) that
tier1 (quant_chunk_remap.py) may have split across chunks and credited
independently. Reuses classify_exclusive (best-hit-wins) from
quant_common.py unchanged. Validated offline: +0.0388 mean Pearson r
(mouse), +0.0391 (human) vs tier1-only baseline, on real project
benchmark reruns.

For every locus cluster with more than 1 member: fetches candidate reads
from the sample's genome-wide BAM restricted to the cluster's combined
span (+/- FLANK bp), aligns with both minimap2 (--eqx, for CIGAR-based
edit counting) and pblat, and assigns each read to its single
best-scoring cluster member via classify_exclusive. Cluster members
always get the rescue count, even 0, overriding tier1. Non-cluster loci
keep tier1's count unchanged.

Usage:
    quant_overlap_rescue.py \\
        --catalog        deduped_metadata.tsv \\
        --tier1_counts   sample1_tier1_counts.tsv \\
        --genome_bam     sample1.genome_align.sorted.bam \\
        --refs_fasta     circle_refs.fa \\
        --refs_lengths   circle_refs_lengths.tsv \\
        --sample         sample1
"""
import argparse
from pathlib import Path
from collections import defaultdict

import pysam
import pandas as pd

from quant_common import FLANK, read_fasta, run_minimap2_bam, run_pblat, parse_bam_hits, parse_psl_hits, classify_exclusive


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', required=True, help='Locus catalog TSV (bsj_id/chrom/start/end/strand)')
    p.add_argument('--tier1_counts', required=True, help='This sample\'s tier1 counts TSV (quant_chunk_remap.py)')
    p.add_argument('--genome_bam', required=True, help='This sample\'s genome-wide alignment BAM')
    p.add_argument('--refs_fasta', required=True, help='circle_refs.fa (quant_build_refs.py)')
    p.add_argument('--refs_lengths', required=True, help='circle_refs_lengths.tsv (quant_build_refs.py)')
    p.add_argument('--sample', required=True)
    p.add_argument('--minimap2', default='minimap2')
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--pblat', default='pblat')
    p.add_argument('--threads', type=int, default=16)
    return p.parse_args()


def build_overlap_clusters(meta: pd.DataFrame) -> dict:
    parent = {bid: bid for bid in meta["bsj_id"]}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for (_chrom, _strand), grp in meta.groupby(["chrom", "strand"]):
        grp = grp.sort_values("start")
        active_end = None
        active_id = None
        for row in grp.itertuples(index=False):
            if active_end is not None and row.start < active_end:
                union(row.bsj_id, active_id)
                active_end = max(active_end, row.end)
            else:
                active_end = row.end
                active_id = row.bsj_id

    clusters = defaultdict(list)
    for bid in parent:
        clusters[find(bid)].append(bid)
    return clusters


def main():
    args = parse_args()
    meta = pd.read_csv(args.catalog, sep='\t')
    tier1 = pd.read_csv(args.tier1_counts, sep='\t')
    tier1_ids = set(tier1["bsj_id"])
    tier1_by_id = tier1.set_index("bsj_id")["remap_count"]
    meta_idx = meta.set_index("bsj_id")

    clusters = build_overlap_clusters(meta)
    target_clusters = {rep: members for rep, members in clusters.items() if len(members) > 1}
    all_target_members = set(m for members in target_clusters.values() for m in members)
    print(f"[{args.sample}] catalog: {len(meta)} loci, multi-member overlap clusters: {len(target_clusters)}, "
          f"covering {len(all_target_members)} loci")

    final = {}
    for row in meta.itertuples(index=False):
        final[row.bsj_id] = int(tier1_by_id.get(row.bsj_id, 0)) if row.bsj_id in tier1_ids else 0

    if not target_clusters:
        out_df = pd.DataFrame({"bsj_id": list(final.keys()), "remap_count": list(final.values())})
        out_df.to_csv(f'{args.sample}_tier1p5_counts.tsv', sep='\t', index=False)
        print(f"[{args.sample}] no multi-member overlap clusters, tier1 counts unchanged")
        return

    bam = pysam.AlignmentFile(str(args.genome_bam), "rb")
    reads = {}
    for _rep, members in target_clusters.items():
        sub = meta_idx.loc[members]
        chrom = sub["chrom"].iloc[0]
        start = int(sub["start"].min()) - FLANK
        end = int(sub["end"].max()) + FLANK
        for read in bam.fetch(chrom, max(0, start), end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.query_name not in reads and read.query_sequence:
                reads[read.query_name] = read.query_sequence
    print(f"[{args.sample}] region-queried {len(target_clusters)} cluster spans, "
          f"{len(reads):,} unique candidate reads fetched")

    workdir = Path(f'{args.sample}_overlap_rescue_work')
    workdir.mkdir(parents=True, exist_ok=True)

    reads_fa = workdir / 'candidate_reads.fa'
    with open(reads_fa, 'w') as out:
        for name, seq in reads.items():
            out.write(f">{name}\n{seq}\n")

    ref_seqs = read_fasta(args.refs_fasta)
    lengths_full = pd.read_csv(args.refs_lengths, sep='\t')
    cand_lengths = lengths_full[lengths_full["bsj_id"].isin(all_target_members)]

    refs_fa = workdir / 'cluster_refs.fa'
    rows = []
    with open(refs_fa, 'w') as out:
        for i, row in enumerate(cand_lengths.itertuples(index=False)):
            seq = ref_seqs.get(row.safe_id, "")
            if not seq:
                continue
            safe_id = f"clu{i}"
            out.write(f">{safe_id}\n")
            for j in range(0, len(seq), 80):
                out.write(seq[j:j + 80] + "\n")
            rows.append({"safe_id": safe_id, "bsj_id": row.bsj_id, "join_pos": row.join_pos})
    lens_df = pd.DataFrame(rows)
    print(f"[{args.sample}] cluster reference: {len(rows)} loci -> {refs_fa}")

    if not reads or lens_df.empty:
        out_df = pd.DataFrame({"bsj_id": list(final.keys()), "remap_count": list(final.values())})
        out_df.to_csv(f'{args.sample}_tier1p5_counts.tsv', sep='\t', index=False)
        return

    bam_out = run_minimap2_bam(args.minimap2, args.samtools, refs_fa, reads_fa,
                                workdir / 'cluster_remap.bam', args.threads,
                                preset_args=("-ax", "map-ont", "--eqx", "-N", "10"))
    psl_path = run_pblat(args.pblat, refs_fa, reads_fa, workdir / 'cluster_remap.psl', args.threads,
                          min_score=50, min_identity=90)

    join_pos_by_ref = dict(zip(lens_df["safe_id"], lens_df["join_pos"]))
    safe_to_bsj = dict(zip(lens_df["safe_id"], lens_df["bsj_id"]))
    support = classify_exclusive(parse_bam_hits(bam_out), parse_psl_hits(psl_path), join_pos_by_ref)

    cluster_counts = {safe_to_bsj[s]: len(qs) for s, qs in support.items()}
    print(f"[{args.sample}] overlap rescue: {sum(cluster_counts.values())} reads assigned across "
          f"{sum(1 for v in cluster_counts.values() if v > 0)} loci (of {len(all_target_members)} cluster members)")

    for bid in all_target_members:
        final[bid] = cluster_counts.get(bid, 0)

    out_df = pd.DataFrame({"bsj_id": list(final.keys()), "remap_count": list(final.values())})
    out_df.to_csv(f'{args.sample}_tier1p5_counts.tsv', sep='\t', index=False)
    print(f"-> {args.sample}_tier1p5_counts.tsv")


if __name__ == '__main__':
    main()

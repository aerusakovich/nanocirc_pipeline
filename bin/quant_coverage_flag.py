#!/usr/bin/env python3
"""
quant_coverage_flag.py

Flags loci for targeted (tier2) rescue by comparing the current remap
count (tier1, post overlap-cluster-rescue) against an independent, older
read-support estimate (quant_tool_readset_union.py) for the same locus.

Two earlier flagging approaches were tried and rejected on real data:
  1. BLAT-vs-minimap2 ratio within tier1's own chunk. Rejected because
     strict PSL filtering also suppresses genuine fragmented BLAT
     evidence for real hard cases, making the ratio unreliable.
  2. A flat "tier1 count < 200" cutoff. Rejected because it flagged 68%
     of loci on real data. Most real circRNAs are naturally
     low-expression, so this was not a meaningful signal.

The rule that worked: flag a locus only if an independent older count is
itself substantial (at or above --min_old_tool_count) and the current
remap count is lower than that. On real data this narrowed 34,762
candidate loci down to 44 flagged.

Usage:
    quant_coverage_flag.py \\
        --deduped_metadata sample_deduped_metadata.tsv \\
        --cluster_map      sample_cluster_map.tsv \\
        --remap_counts     sample1_tier1p5_counts.tsv \\
        --old_tool_counts  sample1_old_tool_counts.tsv \\
        --sample           sample1 \\
        --min_old_tool_count 50
"""
import argparse

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--deduped_metadata', required=True, help='quant_locus_dedup.py deduped catalog (needs bsj_id)')
    p.add_argument('--cluster_map', required=True, help='quant_locus_dedup.py cluster map (bsj_id, representative_bsj_id)')
    p.add_argument('--remap_counts', required=True, help='Current remap counts TSV (bsj_id, remap_count)')
    p.add_argument('--old_tool_counts', required=True, help='quant_tool_readset_union.py output (bsj_id, old_tool_count)')
    p.add_argument('--sample', required=True)
    p.add_argument('--min_old_tool_count', type=float, default=50)
    return p.parse_args()


def flag_low_coverage_loci(deduped_metadata: pd.DataFrame, cluster_map: pd.DataFrame,
                            remap_counts: pd.DataFrame, old_tool_counts: pd.DataFrame,
                            min_old_tool_count: float) -> pd.DataFrame:
    merged = deduped_metadata[["bsj_id"]].merge(remap_counts, on="bsj_id", how="left")
    merged["remap_count"] = merged["remap_count"].fillna(0)

    corrob = cluster_map.merge(old_tool_counts[["bsj_id", "old_tool_count"]], on="bsj_id", how="left")
    best = (corrob.groupby("representative_bsj_id")["old_tool_count"].max()
            .reset_index().rename(columns={"representative_bsj_id": "bsj_id"}))

    full = merged.merge(best, on="bsj_id", how="left")
    full["old_tool_count"] = full["old_tool_count"].fillna(0)
    full["gap"] = full["old_tool_count"] - full["remap_count"]
    full["flagged"] = (full["old_tool_count"] >= min_old_tool_count) & (full["remap_count"] < full["old_tool_count"])

    out = full.rename(columns={"remap_count": "total_count"})
    out["mm2_count"] = 0  # placeholder for backward-compat schema; not used meaningfully downstream
    return out


def main():
    args = parse_args()
    deduped_metadata = pd.read_csv(args.deduped_metadata, sep='\t')
    cluster_map = pd.read_csv(args.cluster_map, sep='\t')
    remap_counts = pd.read_csv(args.remap_counts, sep='\t')
    old_tool_counts = pd.read_csv(args.old_tool_counts, sep='\t')

    out = flag_low_coverage_loci(deduped_metadata, cluster_map, remap_counts, old_tool_counts,
                                  args.min_old_tool_count)

    out_path = f'{args.sample}_flagged_loci.tsv'
    out.to_csv(out_path, sep='\t', index=False)

    n_flagged = out["flagged"].sum()
    print(f"[{args.sample}] {n_flagged} / {len(out)} loci flagged for tier2 rescue "
          f"(old_tool_count >= {args.min_old_tool_count} AND remap_count < old_tool_count)")
    if n_flagged:
        print(out.loc[out["flagged"]].sort_values("gap", ascending=False)
              [["bsj_id", "total_count", "old_tool_count", "gap"]].head(20).to_string(index=False))
    print(f"-> {out_path}")


if __name__ == '__main__':
    main()

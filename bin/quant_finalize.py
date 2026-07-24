#!/usr/bin/env python3
"""
quant_finalize.py

Combines tier1.5 (post overlap-cluster-rescue), tier2 and tier3 into one
final per-locus count, and records which tier each locus's final count
came from.

Rule: tier1.5 stands unless tier2 has a count (tier2 overrides it). All
gene-family loci default to 'tier1_gene_family_unresolved' with
low_confidence=True. Tier3 is applied last and always wins over both.

Usage:
    quant_finalize.py \\
        --deduped_metadata   deduped_metadata.tsv \\
        --tier1_counts       sample1_tier1p5_counts.tsv \\
        --flagged_similarity sample1_flagged_loci_similarity.tsv \\
        --tier2_counts       sample1_tier2_counts.tsv \\
        --tier3_counts       sample1_tier3_counts.tsv \\
        --sample             sample1
"""
import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--deduped_metadata', required=True)
    p.add_argument('--tier1_counts', required=True, help='Post overlap-cluster-rescue (tier1.5) counts TSV')
    p.add_argument('--flagged_similarity', required=True, help='needs is_gene_family')
    p.add_argument('--tier2_counts', default=None)
    p.add_argument('--tier3_counts', default=None)
    p.add_argument('--sample', required=True)
    return p.parse_args()


def finalize_counts(deduped_metadata: pd.DataFrame, tier1: pd.DataFrame,
                     flagged_similarity: pd.DataFrame, tier2: pd.DataFrame, tier3: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in ["bsj_id", "chrom", "start", "end", "strand", "gene_id"] if c in deduped_metadata.columns]
    final = deduped_metadata[cols].merge(tier1, on="bsj_id", how="left")
    final["remap_count"] = final["remap_count"].fillna(0)
    final["quant_tier"] = "tier1"
    final["low_confidence"] = False

    tier2_map = dict(zip(tier2["bsj_id"], tier2["tier2_count"]))
    for bsj_id, count in tier2_map.items():
        mask = final["bsj_id"] == bsj_id
        final.loc[mask, "remap_count"] = count
        final.loc[mask, "quant_tier"] = "tier2"

    gene_family_ids = set(flagged_similarity.loc[flagged_similarity["is_gene_family"], "bsj_id"])
    mask = final["bsj_id"].isin(gene_family_ids)
    final.loc[mask, "low_confidence"] = True
    final.loc[mask, "quant_tier"] = "tier1_gene_family_unresolved"

    tier3_map = dict(zip(tier3["bsj_id"], tier3["tier3_count"]))
    for bsj_id, count in tier3_map.items():
        mask = final["bsj_id"] == bsj_id
        final.loc[mask, "remap_count"] = count
        final.loc[mask, "quant_tier"] = "tier3"
        final.loc[mask, "low_confidence"] = False

    return final


def main():
    args = parse_args()
    deduped_metadata = pd.read_csv(args.deduped_metadata, sep='\t')
    tier1 = pd.read_csv(args.tier1_counts, sep='\t')
    flagged_similarity = pd.read_csv(args.flagged_similarity, sep='\t')
    tier2 = pd.read_csv(args.tier2_counts, sep='\t') if args.tier2_counts and Path(args.tier2_counts).exists() \
        else pd.DataFrame(columns=["bsj_id", "tier2_count"])
    tier3 = pd.read_csv(args.tier3_counts, sep='\t') if args.tier3_counts and Path(args.tier3_counts).exists() \
        else pd.DataFrame(columns=["bsj_id", "tier3_count"])

    final = finalize_counts(deduped_metadata, tier1, flagged_similarity, tier2, tier3)

    out_path = f'{args.sample}_quant_final_counts.tsv'
    final.to_csv(out_path, sep='\t', index=False)
    print(f"[{args.sample}] final counts by tier:\n{final['quant_tier'].value_counts().to_string()}")
    print(f"-> {out_path}")


if __name__ == '__main__':
    main()

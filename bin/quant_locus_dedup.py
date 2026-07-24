#!/usr/bin/env python3
"""
quant_locus_dedup.py

Clusters near-duplicate loci in a discovery-set catalog before building
remap references.

Requires both start and end to be within `--tolerance` bp, not either
alone. An OR-based merge was found to wrongly chain together biologically
distinct isoforms that only share one splice boundary.

Runs once per catalog, not per sample. The discovery catalog is shared
across every sample being quantified.

Usage:
    quant_locus_dedup.py \\
        --catalog  catalog.tsv \\
        --tolerance 10 \\
        --prefix   discovery
"""
import argparse

import pandas as pd


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', required=True,
                    help='Locus catalog TSV (bsj_id/chrom/start/end/strand + rank columns)')
    p.add_argument('--tolerance', type=int, default=10,
                    help='bp tolerance for both start AND end to be considered the same locus')
    p.add_argument('--rank-cols', nargs='+', default=['bsj_confidence', 'n_samples_max'],
                    help='Columns (descending) used to pick each cluster representative; '
                         'columns not present in --catalog are ignored')
    p.add_argument('--prefix', required=True, help='Output file prefix')
    return p.parse_args()


def _pick_representative(g: pd.DataFrame, rank_cols) -> pd.Series:
    """Highest-ranked row wins ties in order: rank_cols (each descending),
    then shortest span (most precisely-called boundary)."""
    ascending = [False] * len(rank_cols) + [True]
    rep = g.sort_values(list(rank_cols) + ["span"], ascending=ascending).iloc[0].copy()
    rep["cluster_root"] = g.name
    return rep


def cluster_duplicate_loci(loci_catalog: pd.DataFrame, tolerance: int, rank_cols):
    rank_cols = [c for c in rank_cols if c in loci_catalog.columns]
    metadata = loci_catalog.reset_index(drop=True).copy()
    n = len(metadata)
    uf = UnionFind(n)

    for _, g in metadata.groupby(["chrom", "strand"]):
        idxs = g.sort_values("start").index.tolist()
        for i in range(len(idxs)):
            for j in range(i + 1, len(idxs)):
                if metadata.at[idxs[j], "start"] - metadata.at[idxs[i], "start"] > tolerance:
                    break
                if abs(metadata.at[idxs[j], "end"] - metadata.at[idxs[i], "end"]) <= tolerance:
                    uf.union(idxs[i], idxs[j])

    metadata["cluster_root"] = [uf.find(i) for i in range(n)]
    metadata["span"] = metadata["end"] - metadata["start"]

    reps = metadata.groupby("cluster_root", group_keys=False).apply(
        lambda g: _pick_representative(g, rank_cols))
    root_to_rep = dict(zip(reps["cluster_root"], reps["bsj_id"]))

    cluster_map = metadata[["bsj_id", "cluster_root"]].copy()
    cluster_map["representative_bsj_id"] = cluster_map["cluster_root"].map(root_to_rep)
    cluster_map = cluster_map.drop(columns=["cluster_root"])

    deduped = reps.drop(columns=["cluster_root", "span"]).reset_index(drop=True)

    n_clusters = metadata["cluster_root"].nunique()
    n_multi = (metadata.groupby("cluster_root").size() > 1).sum()
    print(f"{n_clusters} clusters from {n} loci ({n_multi} clusters have >1 member, "
          f"{n - n_clusters} loci collapsed)")
    return deduped, cluster_map


def main():
    args = parse_args()
    catalog = pd.read_csv(args.catalog, sep='\t')
    deduped, cluster_map = cluster_duplicate_loci(catalog, args.tolerance, args.rank_cols)
    deduped.to_csv(f'{args.prefix}_deduped_metadata.tsv', sep='\t', index=False)
    cluster_map.to_csv(f'{args.prefix}_cluster_map.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
quant_tier3_rescue.py

Tier 3: resolves gene-family and repeat-cluster loci (quant_cross_locus_
similarity.py's is_gene_family==True) instead of excluding them. It
builds one combined reference from every flagged gene-family locus plus
all its transitively-linked genuine_family_siblings (union-find over the
sibling graph, since sibling lists can overlap), aligns the full raw
fastq against that combined reference, and uses classify_boundary_aware
(quant_common.py) to split reads among the competing family members: a
read is handed whole to one paralog only if it crosses that paralog's
diagnostic (paralog-unique) sequence, otherwise it contributes a
fractional 1/N split across every family member it qualifies for. This
resolves boundary-variant loci that otherwise show only 15-58% pairwise
read overlap despite being the same junction, without letting reads that
never leave the shared sequence decide the split by total-span noise
alone.

Usage:
    quant_tier3_rescue.py \\
        --flagged_similarity sample1_flagged_loci_similarity.tsv \\
        --deduped_metadata   deduped_metadata.tsv \\
        --reads              sample1.fastq.gz \\
        --sample             sample1 \\
        --minimap2 minimap2 --samtools samtools --pblat pblat --threads 16
"""
import argparse
import gzip
from pathlib import Path
from collections import defaultdict

import pandas as pd

from quant_common import (MIN_SEGMENT_MATCH, read_fasta, run_minimap2_bam, run_pblat,
                           parse_bam_hits, parse_psl_hits, build_diagnostic_regions, classify_boundary_aware)
from quant_build_refs import build_circle_references


class UnionFind:
    """Dict-keyed union-find with lazy insertion and path halving."""
    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
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


def _diagnostic_regions_by_cluster(ref_lengths: pd.DataFrame, ref_seqs: dict, clusters: dict) -> dict:
    """build_diagnostic_regions, applied separately within each gene-family
    union-find cluster, so a paralog's diagnostic sequence is only defined
    relative to the other members of its own cluster, not to unrelated
    families sharing this rescue call."""
    safe_id_by_bsj = dict(zip(ref_lengths["bsj_id"], ref_lengths["safe_id"]))
    diagnostic = {}
    for members in clusters.values():
        safe_ids = [safe_id_by_bsj[m] for m in members if m in safe_id_by_bsj]
        seqs = {sid: ref_seqs[sid] for sid in safe_ids if sid in ref_seqs}
        if len(seqs) < 2:
            diagnostic.update({sid: [] for sid in seqs})
            continue
        diagnostic.update(build_diagnostic_regions(seqs))
    return diagnostic


def tier3_rescue(flagged_similarity: pd.DataFrame, deduped_metadata: pd.DataFrame, genome_fasta,
                  reads_fq: Path, sample: str, minimap2_bin, samtools_bin, pblat_bin, threads=16) -> pd.DataFrame:
    gf = flagged_similarity[flagged_similarity["is_gene_family"]]
    if gf.empty:
        out = pd.DataFrame(columns=["bsj_id", "tier3_count"])
        out.to_csv(f'{sample}_tier3_counts.tsv', sep="\t", index=False)
        print(f"[{sample}] tier3: no gene-family loci")
        return out

    all_loci = set()
    uf_items = set()
    edges = []
    for row in gf.itertuples(index=False):
        uf_items.add(row.bsj_id)
        all_loci.add(row.bsj_id)
        sibs = str(row.genuine_family_siblings).split(",") if row.genuine_family_siblings else []
        for sib in sibs:
            sib = sib.strip()
            if not sib:
                continue
            uf_items.add(sib)
            all_loci.add(sib)
            edges.append((row.bsj_id, sib))

    uf = UnionFind(uf_items)
    for a, b in edges:
        uf.union(a, b)
    clusters = defaultdict(set)
    for locus in all_loci:
        clusters[uf.find(locus)].add(locus)
    n_clusters = len(clusters)
    print(f"[{sample}] tier3: {len(all_loci)} loci across {n_clusters} gene-family clusters")

    workdir = Path(f'{sample}_tier3_work')
    workdir.mkdir(parents=True, exist_ok=True)
    candidates = deduped_metadata[deduped_metadata["bsj_id"].isin(all_loci)]
    ref_lengths = build_circle_references(genome_fasta, candidates, workdir)
    ref_fa = workdir / "circle_refs.fa"

    bam_path = run_minimap2_bam(minimap2_bin, samtools_bin, ref_fa, reads_fq,
                                 workdir / f"{sample}.tier3.minimap2.bam", threads)

    reads_fa = workdir / f"{sample}.reads.fa"
    _fastq_to_fasta(reads_fq, reads_fa)
    psl_path = run_pblat(pblat_bin, ref_fa, reads_fa, workdir / f"{sample}.tier3.psl", threads,
                          min_score=MIN_SEGMENT_MATCH, min_identity=90)

    join_pos_by_ref = dict(zip(ref_lengths["safe_id"], ref_lengths["join_pos"]))
    safe_to_bsj = dict(zip(ref_lengths["safe_id"], ref_lengths["bsj_id"]))
    ref_seqs = read_fasta(ref_fa)
    diagnostic_regions = _diagnostic_regions_by_cluster(ref_lengths, ref_seqs, clusters)
    support = classify_boundary_aware(parse_bam_hits(bam_path), parse_psl_hits(psl_path),
                                       join_pos_by_ref, diagnostic_regions)

    rows = [{"bsj_id": safe_to_bsj[ref], "tier3_count": support.get(ref, 0.0)}
            for ref in join_pos_by_ref]
    out = pd.DataFrame(rows).sort_values("tier3_count", ascending=False)
    out_path = f'{sample}_tier3_counts.tsv'
    out.to_csv(out_path, sep="\t", index=False)
    print(f"[{sample}] tier3: {len(out)} candidates rescued -> {out_path}")
    return out


def main():
    args = parse_args()
    flagged_similarity = pd.read_csv(args.flagged_similarity, sep='\t')
    deduped_metadata = pd.read_csv(args.deduped_metadata, sep='\t')
    tier3_rescue(flagged_similarity, deduped_metadata, args.genome_fasta, Path(args.reads), args.sample,
                 args.minimap2, args.samtools, args.pblat, args.threads)


if __name__ == '__main__':
    main()

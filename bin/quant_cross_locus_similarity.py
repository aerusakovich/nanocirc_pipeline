#!/usr/bin/env python3
"""
quant_cross_locus_similarity.py

Cross-locus similarity check for this sample's flagged (low-coverage)
loci, against the shared catalog-wide similarity DB built once by
quant_build_similarity_db.py.

Some loci that are individually unique against the whole genome still
show spurious read inflation. This happens when they sit in gene-family
or segmental-duplication regions, where several different, individually
real circRNA loci share a lot of sequence with each other. A
locus-vs-genome uniqueness check cannot see this.

Keeps secondary alignments deliberately: since each flagged locus is
itself a member of the full similarity database, its own trivial
100%-identity self-hit is always reported as PRIMARY by minimap2, so any
genuine match to a *different* locus is necessarily reported as SECONDARY.

A same-chromosome hit within SAME_LOCUS_WINDOW bp is treated as a boundary
variant of the same molecule (not a genuine family match).

Usage:
    quant_cross_locus_similarity.py \\
        --flagged_loci    sample1_flagged_loci.tsv \\
        --refs_fasta      circle_refs.fa \\
        --refs_lengths    circle_refs_lengths.tsv \\
        --similarity_db   locus_similarity_db.mmi \\
        --deduped_metadata deduped_metadata.tsv \\
        --sample          sample1
"""
import argparse
import subprocess
from pathlib import Path
from collections import defaultdict

import pysam
import pandas as pd

from quant_common import read_fasta

MIN_OTHER_HIT_FRAC = 0.1   # fraction of query length matched elsewhere in the catalog -> family member
SAME_LOCUS_WINDOW = 20000  # bp; same-chrom hit within this window is a boundary variant, not a repeat family


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--flagged_loci', required=True, help='quant_coverage_flag.py output for this sample')
    p.add_argument('--refs_fasta', required=True, help='circle_refs.fa (quant_build_refs.py)')
    p.add_argument('--refs_lengths', required=True, help='circle_refs_lengths.tsv (quant_build_refs.py)')
    p.add_argument('--similarity_db', required=True, help='locus_similarity_db.mmi (quant_build_similarity_db.py)')
    p.add_argument('--deduped_metadata', required=True, help='needs bsj_id, chrom, start, end')
    p.add_argument('--sample', required=True)
    p.add_argument('--minimap2', default='minimap2')
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--threads', type=int, default=8)
    return p.parse_args()


def is_genuine_family_hit(target_bsj, query_coords, coords: dict) -> bool:
    """coords: bsj_id -> (chrom, start, end). Cross-chromosome or a
    same-chrom hit far outside SAME_LOCUS_WINDOW is genuine; a same-chrom
    hit within the window is a boundary variant of the same molecule."""
    if target_bsj not in coords:
        return True  # conservative: count it
    q_chrom, q_start, q_end = query_coords
    t_chrom, t_start, t_end = coords[target_bsj]
    if q_chrom != t_chrom:
        return True
    return not (t_start < q_end + SAME_LOCUS_WINDOW and t_end > q_start - SAME_LOCUS_WINDOW)


def main():
    args = parse_args()
    flagged_loci = pd.read_csv(args.flagged_loci, sep='\t')
    circle_refs_lengths = pd.read_csv(args.refs_lengths, sep='\t')
    deduped_metadata = pd.read_csv(args.deduped_metadata, sep='\t')
    seqs = read_fasta(args.refs_fasta)

    flagged = flagged_loci[flagged_loci["flagged"]].copy()

    query_fa = Path(f'{args.sample}_similarity_query.fa')
    flagged_set = set(flagged["bsj_id"])
    query_lengths = {}
    with open(query_fa, "w") as out:
        for row in circle_refs_lengths.itertuples(index=False):
            if row.bsj_id not in flagged_set:
                continue
            seq = seqs[row.safe_id][:row.join_pos]
            if len(seq) < 30:
                continue
            query_lengths[row.bsj_id] = len(seq)
            out.write(f">{row.bsj_id}\n")
            for j in range(0, len(seq), 80):
                out.write(seq[j:j + 80] + "\n")

    if not query_lengths:
        flagged["sibling_loci"] = ""
        flagged["genuine_family_siblings"] = ""
        flagged["n_sibling_loci"] = 0
        flagged["is_gene_family"] = False
        flagged["tier2_candidate"] = True
        out_path = f'{args.sample}_flagged_loci_similarity.tsv'
        flagged.to_csv(out_path, sep='\t', index=False)
        print(f"[{args.sample}] no flagged loci, skipping similarity check")
        return

    bam_path = Path(f'{args.sample}_similarity.bam')
    minimap2 = subprocess.Popen(
        [str(args.minimap2), "-ax", "map-ont", "-N", "20", "-t", str(args.threads),
         str(args.similarity_db), str(query_fa)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    with open(bam_path, "wb") as out:
        samtools = subprocess.Popen([str(args.samtools), "view", "-b", "-@", str(args.threads)],
                                     stdin=minimap2.stdout, stdout=out, stderr=subprocess.DEVNULL)
        minimap2.stdout.close()
        samtools.communicate()
    minimap2.wait()

    coords = {row.bsj_id: (row.chrom, row.start, row.end) for row in deduped_metadata.itertuples(index=False)}

    siblings = defaultdict(set)
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped:
            continue
        qname = r.query_name
        target_bsj = r.reference_name
        if qname not in query_lengths or target_bsj == qname:
            continue
        aligned_len = sum(l for op, l in r.cigartuples if op in (0, 7, 8))
        if aligned_len / query_lengths[qname] >= MIN_OTHER_HIT_FRAC:
            siblings[qname].add(target_bsj)

    flagged["sibling_loci"] = flagged["bsj_id"].map(lambda b: ",".join(sorted(siblings.get(b, set()))))

    def genuine_siblings(bsj_id):
        if bsj_id not in coords:
            return []
        q_coords = coords[bsj_id]
        return sorted(t for t in siblings.get(bsj_id, set()) if is_genuine_family_hit(t, q_coords, coords))

    flagged["genuine_family_siblings"] = flagged["bsj_id"].map(lambda b: ",".join(genuine_siblings(b)))
    flagged["n_sibling_loci"] = flagged["genuine_family_siblings"].map(lambda s: len(s.split(",")) if s else 0)
    flagged["is_gene_family"] = flagged["n_sibling_loci"] > 0
    flagged["tier2_candidate"] = ~flagged["is_gene_family"]

    out_path = f'{args.sample}_flagged_loci_similarity.tsv'
    flagged.to_csv(out_path, sep='\t', index=False)
    print(f"[{args.sample}] {flagged['is_gene_family'].sum()} / {len(flagged)} flagged loci are gene-family/repeat "
          f"members (routed to tier3); {flagged['tier2_candidate'].sum()} remain tier2 candidates")
    print(f"-> {out_path}")


if __name__ == '__main__':
    main()

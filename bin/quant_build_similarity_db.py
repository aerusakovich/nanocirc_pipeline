#!/usr/bin/env python3
"""
quant_build_similarity_db.py

Builds the catalog-wide minimap2 similarity database used to flag
gene-family and repeat-cluster loci. This is the catalog-level,
sample-independent step: it runs once and its output is shared with
every sample's quant_cross_locus_similarity.py, instead of being rebuilt
per sample. Uses a single, undoubled copy of each catalog locus's own
circle sequence.

Usage:
    quant_build_similarity_db.py \\
        --refs_fasta   circle_refs.fa \\
        --refs_lengths circle_refs_lengths.tsv \\
        --minimap2 minimap2 --threads 8
"""
import argparse
import subprocess
from pathlib import Path

import pandas as pd

from quant_common import read_fasta, run_checked


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--refs_fasta', required=True, help='circle_refs.fa (quant_build_refs.py)')
    p.add_argument('--refs_lengths', required=True, help='circle_refs_lengths.tsv (quant_build_refs.py)')
    p.add_argument('--minimap2', default='minimap2')
    p.add_argument('--threads', type=int, default=8)
    return p.parse_args()


def main():
    args = parse_args()
    seqs = read_fasta(args.refs_fasta)
    circle_refs_lengths = pd.read_csv(args.refs_lengths, sep='\t')

    db_fa = Path('locus_similarity_db.fa')
    with open(db_fa, "w") as out:
        for row in circle_refs_lengths.itertuples(index=False):
            seq = seqs[row.safe_id][:row.join_pos]  # single (undoubled) copy
            if len(seq) < 30:
                continue
            out.write(f">{row.bsj_id}\n")
            for j in range(0, len(seq), 80):
                out.write(seq[j:j + 80] + "\n")

    db_mmi = Path('locus_similarity_db.mmi')
    run_checked([str(args.minimap2), "-x", "map-ont", "-t", str(args.threads),
                 "-d", str(db_mmi), str(db_fa)],
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print(f"built similarity DB from {len(circle_refs_lengths)} catalog loci -> {db_fa}, {db_mmi}")


if __name__ == '__main__':
    main()

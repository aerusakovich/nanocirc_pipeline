#!/usr/bin/env python3
"""
quant_chunk_prep.py

Splits circle_refs.fa/_lengths.tsv into fixed-size reference chunks for
chunked remap-and-classify at scale. Runs once per catalog, not per
sample, and its output is shared with every sample's
quant_chunk_remap.py.

Why chunk at all: BLAT/pblat runtime scales badly with reference sequence
count, so a catalog with tens of thousands of loci cannot run as one BLAT
job. Chunks are built from a fixed-seed shuffle of the locus catalog, not
genomic order, so a gene's overlapping isoforms do not land in the same
chunk and compete with each other. This fix took one locus's read count
from 1,319 to 5,066 once it was separated from 18 overlapping neighbors.

Usage:
    quant_chunk_prep.py \\
        --catalog       deduped_metadata.tsv \\
        --refs_fasta    circle_refs.fa \\
        --refs_lengths  circle_refs_lengths.tsv \\
        --chunk_size    100 \\
        --seed          42 \\
        --outdir        chunks/
"""
import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', required=True, help='Locus catalog TSV (needs bsj_id/chrom/start/end)')
    p.add_argument('--refs_fasta', required=True, help='circle_refs.fa (quant_build_refs.py)')
    p.add_argument('--refs_lengths', required=True, help='circle_refs_lengths.tsv (quant_build_refs.py)')
    p.add_argument('--chunk_size', type=int, default=100)
    p.add_argument('--seed', type=int, default=42)
    p.add_argument('--outdir', default='chunks')
    return p.parse_args()


def _read_fasta(path) -> dict:
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def shuffle_and_split_references(loci_catalog: pd.DataFrame, refs_fasta: Path, refs_lengths: Path,
                                  chunks_dir: Path, chunk_size=100, seed=42) -> int:
    lengths = pd.read_csv(refs_lengths, sep="\t")
    lengths = lengths.merge(loci_catalog[["bsj_id", "chrom", "start", "end"]], on="bsj_id", how="left")
    seqs = _read_fasta(refs_fasta)

    lengths = lengths.sample(frac=1, random_state=seed).reset_index(drop=True)

    chunks_dir.mkdir(parents=True, exist_ok=True)
    for f in chunks_dir.glob("chunk_*"):
        f.unlink()

    n_chunks = 0
    for start in range(0, len(lengths), chunk_size):
        idx = start // chunk_size
        chunk_lengths = lengths.iloc[start:start + chunk_size]
        with open(chunks_dir / f"chunk_{idx:04d}.fa", "w") as out_fa:
            for safe_id in chunk_lengths["safe_id"]:
                seq = seqs[safe_id]
                out_fa.write(f">{safe_id}\n")
                for j in range(0, len(seq), 80):
                    out_fa.write(seq[j:j + 80] + "\n")
        chunk_lengths.to_csv(chunks_dir / f"chunk_{idx:04d}_lengths.tsv", sep="\t", index=False)
        n_chunks += 1

    print(f"{n_chunks} reference chunks of up to {chunk_size} loci -> {chunks_dir}")
    return n_chunks


def main():
    args = parse_args()
    catalog = pd.read_csv(args.catalog, sep='\t')
    shuffle_and_split_references(catalog, Path(args.refs_fasta), Path(args.refs_lengths),
                                  Path(args.outdir), args.chunk_size, args.seed)


if __name__ == '__main__':
    main()

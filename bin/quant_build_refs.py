#!/usr/bin/env python3
"""
quant_build_refs.py

Builds synthetic tandem-duplicated circle references for every locus in a
catalog.

For each locus:
  1. Take the spliced exon structure (BED12 blockSizes/blockStarts) and
     extract and join the exonic sequence from the reference genome,
     reverse complemented if strand == '-'. A locus with a single block
     becomes a flat genomic span automatically.
  2. Tandem duplicate the resulting spliced sequence (circle+circle) so a
     linear aligner can find reads crossing the back splice junction, no
     matter where in the circle they start. Always use the full duplicate,
     never a short junction-flank probe: testing showed a probe drops
     recall to 36.5% vs 59.8% for the full reference, with more false
     positives, not fewer.

Usage:
    quant_build_refs.py \\
        --catalog     deduped_metadata.tsv \\
        --genome_fasta genome.fa \\
        --outdir      .
"""
import argparse
from pathlib import Path

import pysam
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', required=True,
                    help='Locus catalog TSV (bsj_id/chrom/start/strand/blockSizes/blockStarts)')
    p.add_argument('--genome_fasta', required=True, help='Reference genome FASTA (indexed or not)')
    p.add_argument('--outdir', default='.', help='Output directory for circle_refs.fa/_lengths.tsv')
    return p.parse_args()


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def build_circle_seq(fasta: pysam.FastaFile, chrom: str, start: int, sizes: list, starts: list, strand: str) -> str:
    """Concatenate exon blocks (BED12 blockSizes/blockStarts, both relative
    to `start`) into the circle's own spliced (mature) sequence."""
    exon_seqs = [fasta.fetch(chrom, start + s, start + s + sz) for s, sz in zip(starts, sizes)]
    seq = "".join(exon_seqs).upper()
    if strand == "-":
        seq = revcomp(seq)
    return seq


def build_circle_references(genome_fasta, gt_loci: pd.DataFrame, out_dir: Path) -> pd.DataFrame:
    loci = gt_loci.drop_duplicates(subset="bsj_id")[
        ["bsj_id", "chrom", "start", "strand", "blockSizes", "blockStarts"]
    ].reset_index(drop=True)
    print(f"Building circle references for {len(loci)} loci ...")

    fasta = pysam.FastaFile(str(genome_fasta))
    fasta_chroms = set(fasta.references)

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ref_fa_path = out_dir / "circle_refs.fa"
    lengths_rows = []
    n_ok = n_skip = 0
    with open(ref_fa_path, "w") as out_fa:
        for i, row in enumerate(loci.itertuples(index=False)):
            if row.chrom not in fasta_chroms:
                n_skip += 1
                continue
            sizes = [int(x) for x in str(row.blockSizes).rstrip(",").split(",") if x]
            starts = [int(x) for x in str(row.blockStarts).rstrip(",").split(",") if x]
            try:
                seq = build_circle_seq(fasta, row.chrom, row.start, sizes, starts, row.strand)
            except Exception:
                n_skip += 1
                continue
            if len(seq) < 30:
                n_skip += 1
                continue
            circle_length = len(seq)
            safe_id = f"circ{i}"
            ref_seq = seq + seq
            join_pos = circle_length  # BSJ = crossing this point in the doubled seq

            out_fa.write(f">{safe_id}\n")
            for j in range(0, len(ref_seq), 80):
                out_fa.write(ref_seq[j:j + 80] + "\n")
            lengths_rows.append({
                "safe_id": safe_id, "bsj_id": row.bsj_id, "circle_length": circle_length,
                "ref_length": len(ref_seq), "join_pos": join_pos,
            })
            n_ok += 1

    lengths_df = pd.DataFrame(lengths_rows)
    lengths_df.to_csv(out_dir / "circle_refs_lengths.tsv", sep="\t", index=False)
    print(f"Built {n_ok} full tandem-duplicated circle references ({n_skip} skipped) -> {ref_fa_path}")
    if n_ok:
        print(f"  total reference size: {lengths_df['ref_length'].sum():,} bp "
              f"(longest circle: {lengths_df['circle_length'].max():,} bp)")
    return lengths_df


def main():
    args = parse_args()
    catalog = pd.read_csv(args.catalog, sep='\t')
    build_circle_references(args.genome_fasta, catalog, Path(args.outdir))


if __name__ == '__main__':
    main()

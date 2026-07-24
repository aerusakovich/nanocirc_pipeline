#!/usr/bin/env python3
"""
quant_genome_align.py

One genome-wide minimap2 splice alignment per sample against its full raw
fastq (not an enriched or subsetted read set), producing a
coordinate-sorted, indexed BAM. This is a fixed per-sample cost, so BLAT
rescue scaling does not depend on reference size: reads overlapping any
locus's genomic span become the rescue candidate pool
(quant_rescue_pool.py), regardless of tool or locus count.

Usage:
    quant_genome_align.py \\
        --genome_index genome_splice.mmi \\
        --reads sample.fastq.gz \\
        --sample sample1 \\
        --minimap2 minimap2 --samtools samtools \\
        --threads 16
"""
import argparse
import subprocess
from pathlib import Path

from quant_common import run_checked


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--genome_index', required=True, help='Pre-built minimap2 splice-mode .mmi (quant_genome_index.py)')
    p.add_argument('--reads', required=True, help='Sample FASTQ (full raw reads, not a subset)')
    p.add_argument('--sample', required=True, help='Sample name, used for output file naming')
    p.add_argument('--minimap2', default='minimap2')
    p.add_argument('--samtools', default='samtools')
    p.add_argument('--threads', type=int, default=16)
    return p.parse_args()


def align_run_to_genome(minimap2_bin, samtools_bin, mmi, reads_fq, out_bam: Path, threads=16) -> Path:
    minimap2 = subprocess.Popen(
        [str(minimap2_bin), "-ax", "splice", "-t", str(threads),
         str(mmi), str(reads_fq)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    samtools_sort = subprocess.Popen(
        [str(samtools_bin), "sort", "-@", str(threads), "-m", "2G",
         "-o", str(out_bam), "-"],
        stdin=minimap2.stdout, stderr=subprocess.DEVNULL,
    )
    minimap2.stdout.close()
    samtools_sort.communicate()
    minimap2.wait()
    run_checked([str(samtools_bin), "index", str(out_bam)])
    return out_bam


def main():
    args = parse_args()
    out_bam = Path(f'{args.sample}.genome_align.sorted.bam')
    align_run_to_genome(args.minimap2, args.samtools, args.genome_index, args.reads, out_bam, args.threads)
    print(f"[{args.sample}] genome-wide alignment -> {out_bam}")


if __name__ == '__main__':
    main()

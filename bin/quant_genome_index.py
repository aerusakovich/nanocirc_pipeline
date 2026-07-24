#!/usr/bin/env python3
"""
quant_genome_index.py

One-time minimap2 splice-mode index of the reference genome, reused by
every sample's genome-wide alignment pass. Runs once per pipeline run
and its output is shared with every sample, not rebuilt per sample.

Builds to a per-process temp file first, then atomically renames it into
place (os.replace, a same-directory rename, atomic on POSIX). This
guards against concurrent builds of the shared index file.

Usage:
    quant_genome_index.py --genome_fasta genome.fa --minimap2 minimap2 --output genome_splice.mmi
"""
import argparse
import os
from pathlib import Path

from quant_common import run_checked


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--genome_fasta', required=True)
    p.add_argument('--minimap2', default='minimap2')
    p.add_argument('--output', required=True, help='Output .mmi index path')
    p.add_argument('--threads', type=int, default=8)
    return p.parse_args()


def build_genome_index(genome_fasta, minimap2_bin, out: Path, threads=8) -> Path:
    if out.exists():
        return out
    tmp = out.with_suffix(f".mmi.tmp{os.getpid()}")
    run_checked(
        [str(minimap2_bin), "-x", "splice", "-t", str(threads),
         "-d", str(tmp), str(genome_fasta)],
    )
    if out.exists():
        # another process finished first while we were building, use theirs
        tmp.unlink()
    else:
        os.replace(tmp, out)
    return out


def main():
    args = parse_args()
    build_genome_index(args.genome_fasta, args.minimap2, Path(args.output), args.threads)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
quant_catalog_from_bed12.py

Converts the pipeline's discovery BED12 and confidence TSV (see
docs/output.md's Confidence TSV format) into the flat locus-catalog TSV
schema the quantification modules expect: bsj_id, chrom, start, end,
strand, blockSizes, blockStarts, bsj_confidence.

Usage:
    quant_catalog_from_bed12.py \\
        --bed  sample_discovery.bed12 \\
        --tsv  sample_discovery_confidence.tsv \\
        --output sample_catalog.tsv
"""
import argparse
import csv


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bed',    required=True, help='Discovery BED12')
    p.add_argument('--tsv',    required=True, help='Matching *_confidence.tsv')
    p.add_argument('--output', required=True, help='Output catalog TSV path')
    return p.parse_args()


def load_bsj_confidence(tsv_path):
    conf = {}
    with open(tsv_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if header is None:
                header = line.lstrip('#').split('\t')
                idx = {c: i for i, c in enumerate(header)}
                continue
            cols = line.split('\t')
            conf[cols[idx['bsj_id']]] = cols[idx['bsj_confidence']] if 'bsj_confidence' in idx else ''
    return conf


def main():
    args = parse_args()
    bsj_confidence = load_bsj_confidence(args.tsv)

    with open(args.bed) as fh, open(args.output, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['bsj_id', 'chrom', 'start', 'end', 'strand', 'blockSizes', 'blockStarts', 'bsj_confidence'])
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith(('#', 'track', 'browser')):
                continue
            cols = line.split('\t')
            chrom, start, end, name, _score, strand = cols[0:6]
            block_sizes, block_starts = cols[10], cols[11]
            writer.writerow([name, chrom, start, end, strand, block_sizes, block_starts,
                              bsj_confidence.get(name, '')])


if __name__ == '__main__':
    main()

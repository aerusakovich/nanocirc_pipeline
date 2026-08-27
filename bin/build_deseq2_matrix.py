#!/usr/bin/env python3
"""
build_deseq2_matrix.py

Wide sample x isoform read-count matrix for DESeq2 (or any other
count-based DE tool), built from every run's own <tier>_clean_with_counts.tsv
(the pipeline's final, post-filter per-run quantification output for that
tier). A feature row is one isoform: same BSJ but a different exon
structure gets its own row. Runs quantified against the same crossrun
catalog already share row identity; runs from different catalogs (or no
crossrun merge at all) are unioned here, 0-filled where a run's own catalog
never called that isoform.

Usage:
    build_deseq2_matrix.py --tier high_confidence \
        --counts_tsv run1=run1_high_confidence_clean_with_counts.tsv \
                     run2=run2_high_confidence_clean_with_counts.tsv \
        --groups run1=groupA run2=groupB \
        --out_counts   deseq2_counts_high_confidence.tsv \
        --out_coldata  deseq2_coldata_high_confidence.tsv \
        --out_features deseq2_features_high_confidence.tsv
"""
import argparse

import pandas as pd

FEATURE_COLS = ['chrom', 'start', 'end', 'strand', 'sel_block_count',
                'sel_block_sizes', 'sel_block_starts', 'bsj_id', 'type']


def parse_kv_list(pairs):
    out = {}
    for pair in pairs:
        k, v = pair.split('=', 1)
        out[k] = v
    return out


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tier', required=True)
    p.add_argument('--counts_tsv', nargs='+', required=True, help='sample=path pairs')
    p.add_argument('--groups', nargs='+', default=[], help='sample=group pairs')
    p.add_argument('--out_counts', required=True)
    p.add_argument('--out_coldata', required=True)
    p.add_argument('--out_features', required=True)
    return p.parse_args()


def isoform_id(df):
    return (df['bsj_id'].astype(str) + '|' +
            df['sel_block_sizes'].astype(str) + '|' +
            df['sel_block_starts'].astype(str))


def main():
    args = parse_args()
    sample_paths = parse_kv_list(args.counts_tsv)
    sample_groups = parse_kv_list(args.groups)

    counts = {}
    features = {}
    for sample, path in sample_paths.items():
        df = pd.read_csv(path, sep='\t')
        df = df.rename(columns={'#chrom': 'chrom'})
        df['isoform_id'] = isoform_id(df)
        df = df.drop_duplicates('isoform_id')
        counts[sample] = df.set_index('isoform_id')['nanocirc_quant_reads']
        for _, row in df[['isoform_id'] + FEATURE_COLS].iterrows():
            features.setdefault(row['isoform_id'], row[FEATURE_COLS].to_dict())

    counts_df = pd.DataFrame(counts).fillna(0).astype(int)
    counts_df.index.name = 'isoform_id'
    counts_df = counts_df.sort_index()
    counts_df.to_csv(args.out_counts, sep='\t')

    features_df = pd.DataFrame.from_dict(features, orient='index')
    features_df.index.name = 'isoform_id'
    features_df = features_df.loc[counts_df.index]
    features_df.to_csv(args.out_features, sep='\t')

    coldata_rows = [{'sample': s, 'group': sample_groups.get(s, s)} for s in counts_df.columns]
    pd.DataFrame(coldata_rows).to_csv(args.out_coldata, sep='\t', index=False)

    print(f"[{args.tier}] {counts_df.shape[0]} isoforms x {counts_df.shape[1]} samples")
    print(f"-> {args.out_counts}")
    print(f"-> {args.out_coldata}")
    print(f"-> {args.out_features}")


if __name__ == '__main__':
    main()

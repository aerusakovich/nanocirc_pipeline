#!/usr/bin/env python3
"""
append_run_counts_matrix.py

Fills out clean_with_counts.tsv with every run's quantified
read count per isoform, reusing the already-built DESeq2 count matrix
(build_deseq2_matrix.py's deseq2_counts_<tier>.tsv: isoform_id x run,
values from each run's own nanocirc_quant_reads). Runs after BUILD_DESEQ2_MATRIX
in the pipeline and replaces that same clean_with_counts.tsv file with a
fuller version, rather than publishing a second file alongside it.

Adds one new column per run in the matrix, named <run>_reads. The calling
run's own column duplicates its existing nanocirc_quant_reads value, kept
so every enriched clean.tsv has the same run-column set regardless of
which run it belongs to.

Usage:
    append_run_counts_matrix.py \\
        --clean_tsv     run1_balanced_precision_clean_with_counts.tsv \\
        --counts_matrix deseq2_counts_balanced_precision.tsv \\
        --output        run1_balanced_precision_clean_with_counts.tsv
"""
import argparse

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--clean_tsv', required=True)
    p.add_argument('--counts_matrix', required=True)
    p.add_argument('--output', required=True)
    return p.parse_args()


def isoform_id(df):
    return (df['bsj_id'].astype(str) + '|' +
            df['sel_block_sizes'].astype(str) + '|' +
            df['sel_block_starts'].astype(str))


def main():
    args = parse_args()
    clean = pd.read_csv(args.clean_tsv, sep='\t')
    matrix = pd.read_csv(args.counts_matrix, sep='\t', index_col='isoform_id')

    original_cols = list(clean.columns)
    clean['isoform_id'] = isoform_id(clean)

    assert matrix.index.is_unique, 'counts_matrix has duplicate isoform_id rows, would multiply clean_tsv rows on merge'

    run_cols = {c: f'{c}_reads' for c in matrix.columns}
    joined = clean.merge(matrix.rename(columns=run_cols), on='isoform_id', how='left')

    assert len(joined) == len(clean), (
        f'row count changed: {len(clean)} filtered rows in, {len(joined)} out. '
        'This step must only append columns, never rows.'
    )
    assert set(original_cols) <= set(joined.columns), 'a column from the filtered input was lost'

    joined[list(run_cols.values())] = joined[list(run_cols.values())].fillna(0).astype(int)
    joined = joined.drop(columns=['isoform_id'])
    joined.to_csv(args.output, sep='\t', index=False)

    print(f'[append_run_counts_matrix] {len(joined)} rows, appended {len(run_cols)} '
          f'per-run count columns ({", ".join(run_cols.values())}) -> {args.output}')


if __name__ == '__main__':
    main()

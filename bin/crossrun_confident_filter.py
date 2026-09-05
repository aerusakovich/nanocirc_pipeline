#!/usr/bin/env python3
"""
crossrun_confident_filter.py

Cross-run merging (crossrun_annotate.py) draws its candidate loci from
each run's confidence-label-filtered tier output, before the
per-run, quantification-based read-count guard
(FILTER_CONFIDENT_DISCOVERY / filter_confident_discovery.py) has run --
that guard needs quantified read counts, which are not available until
after quantification, and drops circNICK-LRS-only calls (discovery,
balanced) or weak isoCirc-only calls (high_confidence) with read support
at or below --circrna_confident_min_reads. Left alone, a locus the guard
would drop can still enter the published crossrun catalog, and its
crossrun-level read count then silently reads 0 (no per-sample count
survives to feed the crossrun count matrix for it), rather than the
locus being absent as it would be in the per-sample tier output it is
supposed to represent.

This script re-applies that guard at the crossrun level: a crossrun
locus is kept only if it also survived FILTER_CONFIDENT_DISCOVERY in at
least one of the group's own runs (the same any-one-run-supports-it
logic the rest of crossrun merging already uses). Filters both the
group's clean TSV and its matching BED12 to the same kept set.

Usage:
    crossrun_confident_filter.py \\
        --group_clean_tsv  GROUP_TIER_crossrun_clean.tsv \\
        --group_bed12      GROUP_TIER_crossrun.bed12 \\
        --run_filtered_tsvs run1_TIER_clean_with_counts.tsv run2_TIER_clean_with_counts.tsv \\
        --out_tsv   GROUP_TIER_crossrun_clean.tsv \\
        --out_bed12 GROUP_TIER_crossrun.bed12
"""
import argparse
import csv


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--group_clean_tsv', required=True, help='Crossrun-merged group clean TSV (crossrun_annotate.py output)')
    p.add_argument('--group_bed12', required=True, help='Matching crossrun-merged group BED12')
    p.add_argument('--run_filtered_tsvs', nargs='+', required=True,
                    help='Per-run FILTER_CONFIDENT_DISCOVERY clean_with_counts.tsv paths, '
                         'each over this same group+tier locus set')
    p.add_argument('--out_tsv', required=True)
    p.add_argument('--out_bed12', required=True)
    return p.parse_args()


def load_bsj_ids(path):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        return {row['bsj_id'] for row in reader if row.get('bsj_id')}


def main():
    args = parse_args()

    confident_ids = set()
    for path in args.run_filtered_tsvs:
        confident_ids |= load_bsj_ids(path)

    # Read both inputs fully before writing anything: --out_tsv/--out_bed12
    # share a filename with --group_clean_tsv/--group_bed12 (the crossrun
    # merge output gets republished under the same name), so opening an
    # output path for writing while still reading the same path would
    # truncate it out from under the read.
    with open(args.group_clean_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fieldnames = reader.fieldnames
        all_rows = list(reader)
    rows = [row for row in all_rows if row.get('bsj_id') in confident_ids]

    with open(args.group_bed12) as fh:
        bed_lines = fh.readlines()
    bed_out = [line for line in bed_lines
               if len(line.rstrip('\n').split('\t')) >= 4 and line.rstrip('\n').split('\t')[3] in confident_ids]

    with open(args.out_tsv, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    with open(args.out_bed12, 'w') as fh:
        fh.writelines(bed_out)

    print(f'[crossrun_confident_filter] {len(all_rows)} -> {len(rows)} clean TSV rows kept '
          f'({len(all_rows) - len(rows)} dropped, no confident read support in any run), '
          f'{len(bed_out)} BED12 rows kept -> {args.out_tsv}, {args.out_bed12}')


if __name__ == '__main__':
    main()

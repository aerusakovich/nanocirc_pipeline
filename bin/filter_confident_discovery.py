#!/usr/bin/env python3
"""
filter_confident_discovery.py
------------------------------
Drops low-confidence patterns from a quantified clean TSV
(circrna_clean.py / crossrun_annotate.py's _clean_with_counts.tsv, after
quant_append_counts.py) and its matching BED12. Found by benchmarking
against simulated ground truth:

  1. CircNick-LRS-only loci with weak read support. CircNick-LRS' own
     boundary calls are often unsupported once checked again
     (see nanocirc_quant_reads), even inside a real gene. This is the
     single biggest false-positive source in the unfiltered catalog.

Needs read support at or below --min_reads. A missing nanocirc_quant_reads
value also fails this check. Applied to the discovery and balanced
tiers. high_confidence already excludes this pattern through its own
stricter multi-tool rules (checked: 0 rows affected there).

  2. IsoCirc-only calls with weak read support, guarding the isocirc_only
     exception in filter_confidence.py (which lets IsoCirc's own
     Low-confidence calls into the high_confidence tier). Only applied
     when --category high_confidence: IsoCirc-only calls are otherwise a
     strong signal (~0.96-0.97 precision unfiltered, see the balanced
     per-tool breakdown) and are intentionally left unfiltered on
     discovery/balanced to avoid re-litigating their already-validated
     numbers -- this rule exists only to keep high_confidence's near-zero
     false-positive rate intact now that a Low-confidence door has been
     opened there for IsoCirc.

Usage:
    filter_confident_discovery.py \\
        --tsv         sample_discovery_clean_with_counts.tsv \\
        --bed         sample_discovery.bed12 \\
        --min_reads   2 \\
        --category    discovery \\
        --out_tsv     sample_discovery_clean_with_counts.tsv \\
        --out_bed     sample_discovery_clean.bed12
"""
import argparse


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--tsv',       required=True, help='clean_with_counts.tsv (has type, supporting_tools, nanocirc_quant_reads)')
    p.add_argument('--bed',       required=True, help='matching BED12')
    p.add_argument('--min_reads', type=int, required=True, help='drop candidates at or below this read count (params.circrna_confident_min_reads)')
    p.add_argument('--category',  default='', help='tier name (meta.category); enables the isocirc-only guard when set to high_confidence')
    p.add_argument('--out_tsv',   required=True, help='output TSV path')
    p.add_argument('--out_bed',   required=True, help='output BED12 path')
    return p.parse_args()


def weak_read_support(reads_str, min_reads):
    try:
        return int(reads_str) <= min_reads
    except (TypeError, ValueError):
        return True  # missing or blank count, treat as unsupported


def should_drop(row, min_reads, category=''):
    tools = row.get('supporting_tools', '')
    tool_list = [t for t in tools.split(',') if t]
    low_read = weak_read_support(row.get('nanocirc_quant_reads', ''), min_reads)

    circnick_only = category in ('discovery', 'balanced') and tool_list == ['circnick']
    isocirc_only_weak = category == 'high_confidence' and tool_list == ['isocirc']

    return low_read and (circnick_only or isocirc_only_weak)


def main():
    args = parse_args()

    import csv
    with open(args.tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fieldnames = reader.fieldnames
        rows = list(reader)

    kept_ids = set()
    kept_rows = []
    for row in rows:
        if should_drop(row, args.min_reads, args.category):
            continue
        kept_ids.add(row.get('bsj_id', ''))
        kept_rows.append(row)

    with open(args.out_tsv, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(kept_rows)

    n_bed = 0
    with open(args.bed) as fi, open(args.out_bed, 'w') as fo:
        for line in fi:
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 4 and cols[3] in kept_ids:
                fo.write(line)
                n_bed += 1

    print(f'[filter_confident_discovery] {len(rows)} -> {len(kept_rows)} rows kept '
          f'({len(rows) - len(kept_rows)} dropped), {n_bed} BED12 rows kept '
          f'-> {args.out_tsv}, {args.out_bed}')


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
quant_append_counts.py

Appends quantification read counts onto a per-sample, per-tier clean TSV
(circrna_clean.py's output), for downstream analysis such as differential
expression. This is the join step connecting quant_finalize.py's output
back to the rest of the pipeline's annotated circRNA tables. Getting this
join wrong causes about 19% of loci to silently get 0 counts.

Why the join cannot go directly from clean.tsv's bsj_id to quant_final_
counts.tsv's bsj_id: quant_final_counts.tsv is keyed by the
representative locus id from the (possibly per-group) deduped catalog
(quant_locus_dedup.py collapses near-duplicate loci to one
representative). clean.tsv's bsj_id is in the full pre-dedup catalog's
key space. Joining clean.tsv directly against quant_final_counts.tsv
would silently return 0 or missing for every non-representative locus.
The join must hop through cluster_map.tsv (bsj_id maps to
representative_bsj_id, covering every locus in the full catalog, not
just representatives): look up clean.tsv's bsj_id there first, falling
back to the row's own bsj_id if absent (it is already its own
representative, a singleton locus with no near-duplicates), then look up
that representative id in quant_final_counts.tsv.

bsj_id is used exactly as it appears everywhere in this chain. Isoform
rows (bsj_id suffixed, e.g. '|iso1') are never stripped to a base id.
Each isoform is its own distinct catalog locus with its own circle
reference and its own representative in cluster_map.tsv, since that is
how the catalog was built (quant_catalog_from_bed12.py copies the BED12
name field verbatim). Stripping the suffix here would look up the wrong
key against cluster_map.tsv and quant_final_counts.tsv, which both still
have the full, un-stripped id.

Usage:
    quant_append_counts.py \\
        --clean_tsv    sample1_discovery_clean.tsv \\
        --cluster_map  unit_cluster_map.tsv \\
        --quant_counts sample1_quant_final_counts.tsv \\
        --output       sample1_discovery_clean_with_counts.tsv
"""
import argparse
import csv


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--clean_tsv', required=True, help="circrna_clean.py's <sample>_<tier>_clean.tsv")
    p.add_argument('--cluster_map', required=True, help='quant_locus_dedup.py cluster_map (bsj_id, representative_bsj_id)')
    p.add_argument('--quant_counts', required=True, help='quant_finalize.py final counts (bsj_id, remap_count, quant_tier, low_confidence)')
    p.add_argument('--output', required=True)
    return p.parse_args()


def load_cluster_map(path):
    rep_of = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            rep_of[row['bsj_id']] = row['representative_bsj_id']
    return rep_of


def load_quant_counts(path):
    counts = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            counts[row['bsj_id']] = {
                'remap_count': row.get('remap_count', '0'),
                'quant_tier': row.get('quant_tier', ''),
                'low_confidence': row.get('low_confidence', ''),
            }
    return counts


def main():
    args = parse_args()
    rep_of = load_cluster_map(args.cluster_map)
    counts = load_quant_counts(args.quant_counts)

    with open(args.clean_tsv) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        fieldnames = list(reader.fieldnames) + [
            'nanocirc_quant_reads', 'nanocirc_quant_tier', 'nanocirc_quant_low_confidence',
        ]
        rows = list(reader)

    n_via_map = n_fallback = n_missing = 0
    with open(args.output, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            bsj_id = row.get('bsj_id', '')

            if bsj_id in rep_of:
                rep_id = rep_of[bsj_id]
                n_via_map += 1
            else:
                rep_id = bsj_id  # singleton locus, already its own representative
                n_fallback += 1

            hit = counts.get(rep_id)
            if hit is None:
                n_missing += 1
                row['nanocirc_quant_reads'] = 0
                row['nanocirc_quant_tier'] = ''
                row['nanocirc_quant_low_confidence'] = ''
            else:
                row['nanocirc_quant_reads'] = hit['remap_count']
                row['nanocirc_quant_tier'] = hit['quant_tier']
                row['nanocirc_quant_low_confidence'] = hit['low_confidence']

            writer.writerow(row)

    print(f'[quant_append_counts] {len(rows)} rows: {n_via_map} via cluster_map, '
          f'{n_fallback} singleton fallback, {n_missing} with no quant count found '
          f'(defaulted to 0) -> {args.output}')


if __name__ == '__main__':
    main()

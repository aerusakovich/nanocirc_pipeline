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

Isoform identity is genomic structure, not the raw blockSizes/blockStarts
TEXT: two tools/samples can write the identical single-block structure as
"486,"/"0," or "486"/"0" (isocirc/circFL omit the trailing comma, our own
CIRILONG_TO_BED12/CIRCNICK_TO_BED12 conversions include it -- see
smart_merge.py's normalize_block_field). Concatenating that raw text into
the join key used to split one real isoform into two rows whenever
different samples' winning source tool happened to differ -- on real
data, a locus with real counts only in "resistant" samples and only in
"sensitive" samples under the other text form, which looks exactly like a
100%-switching isoform when it is really one isoform counted once each
way. Fixed here in two steps: (1) normalize blockSizes/blockStarts before
building the key, so pure text variants collapse immediately; (2) beyond
that, genuinely near-duplicate structures within the same BSJ (a few bp
of real boundary-calling noise, same fix as smart_merge.py's
abs_struct_similar) are also clustered and their counts summed, keeping
only structures that differ by more than that tolerance as separate rows.

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
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from smart_merge import normalize_block_field, abs_struct_similar

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
    p.add_argument('--merge_tolerance', type=int, default=10,
                   help='bp tolerance for clustering near-duplicate exon structures '
                        'within the same BSJ (default matches quant_locus_dedup.py)')
    p.add_argument('--out_counts', required=True)
    p.add_argument('--out_coldata', required=True)
    p.add_argument('--out_features', required=True)
    return p.parse_args()


def isoform_id(df):
    return (df['bsj_id'].astype(str) + '|' +
            df['sel_block_sizes'].astype(str) + '|' +
            df['sel_block_starts'].astype(str))


def base_bsj(isoform_ids: pd.Series) -> pd.Series:
    """bsj_id without an "|isoN" suffix: chrom:start-end:strand."""
    return isoform_ids.str.split('|').str[0]


def base_locus(isoform_ids: pd.Series) -> pd.Series:
    """chrom:start-end, strand stripped -- the broader grouping key used to
    catch a "." (unresolved-strand) call that is really a duplicate of a
    resolved-strand ("+"/"-") call at the identical coordinates (same root
    cause as classify_types()'s own -s/-S blind spot, see
    sanity_check_isoform_types.py): a strand-unresolved catalog entry can
    never be compared against a resolved one under base_bsj alone, since
    the strand character itself makes the two base_bsj strings differ."""
    return base_bsj(isoform_ids).str.rsplit(':', n=1).str[0]


def _abs_coords(sizes_str, starts_str):
    sizes = [int(x) for x in str(sizes_str).rstrip(',').split(',') if x]
    starts = [int(x) for x in str(starts_str).rstrip(',').split(',') if x]
    return tuple((s, sz) for s, sz in zip(starts, sizes))


def cluster_isoforms(features_df: pd.DataFrame, tolerance: int) -> dict:
    """{isoform_id: representative_isoform_id}, clustering by exon-structure
    similarity within each locus (chrom:start-end, strand-agnostic). A "."
    (unresolved-strand) member can merge into a real-strand ("+"/"-")
    cluster; two members with different REAL strands never merge (a
    genuine antisense pair is a real biological conflict, not a
    duplicate) -- same rule smart_merge.py's group_relaxed() already uses
    at discovery time. The representative is whichever member has the
    most exon blocks with a resolved real strand if any exist (so the
    kept row's own strand/type is the confirmed one, not "."), else
    lexicographic order."""
    remap = {}
    for _locus, grp in features_df.groupby(base_locus(features_df.index.to_series())):
        structs = [(iid, _abs_coords(row.sel_block_sizes, row.sel_block_starts), row.strand)
                   for iid, row in grp.iterrows()]
        clusters = []  # [ [struct, {strands seen}, [isoform_id, ...]] ]
        for iid, s, strand in structs:
            placed = False
            for cl in clusters:
                if not abs_struct_similar(s, cl[0], tolerance):
                    continue
                real_strands = cl[1] - {'.'}
                if strand != '.' and real_strands and strand not in real_strands:
                    continue  # genuine antisense conflict, not a duplicate
                cl[1].add(strand)
                cl[2].append(iid)
                placed = True
                break
            if not placed:
                clusters.append([s, {strand}, [iid]])
        for _struct, _cl_strands, members in clusters:
            resolved = [m for m in members if not m.split('|')[0].endswith(':.')]
            rep = sorted(resolved)[0] if resolved else sorted(members)[0]
            for m in members:
                remap[m] = rep
    return remap


def main():
    args = parse_args()
    sample_paths = parse_kv_list(args.counts_tsv)
    sample_groups = parse_kv_list(args.groups)

    counts = {}
    features = {}
    for sample, path in sample_paths.items():
        df = pd.read_csv(path, sep='\t')
        df = df.rename(columns={'#chrom': 'chrom'})
        df['sel_block_sizes'] = df['sel_block_sizes'].map(normalize_block_field)
        df['sel_block_starts'] = df['sel_block_starts'].map(normalize_block_field)
        df['isoform_id'] = isoform_id(df)
        # Sum, not drop_duplicates: a sample whose own catalog somehow carries
        # two rows for the same normalized isoform_id should have its real
        # read support combined, not silently truncated to the first row.
        counts[sample] = df.groupby('isoform_id')['nanocirc_quant_reads'].sum()
        for _, row in df.drop_duplicates('isoform_id')[['isoform_id'] + FEATURE_COLS].iterrows():
            features.setdefault(row['isoform_id'], row[FEATURE_COLS].to_dict())

    counts_df = pd.DataFrame(counts).fillna(0).astype(int)
    counts_df.index.name = 'isoform_id'
    features_df = pd.DataFrame.from_dict(features, orient='index')
    features_df.index.name = 'isoform_id'
    features_df = features_df.loc[counts_df.index]

    remap = cluster_isoforms(features_df, args.merge_tolerance)
    n_merged = sum(1 for iid, rep in remap.items() if iid != rep)
    if n_merged:
        print(f"[{args.tier}] merging {n_merged} near-duplicate isoform(s) into "
              f"{len(set(remap.values()))} representative structures (tolerance {args.merge_tolerance}bp)")
        # Representative's own row (its own block sizes/starts), not an
        # arbitrary merged sibling's, is what the merged row's feature
        # columns should report.
        features_df = features_df.loc[sorted(set(remap.values()))]
        counts_df.index = counts_df.index.map(remap)
        counts_df = counts_df.groupby(counts_df.index).sum()

    counts_df = counts_df.sort_index()
    counts_df.to_csv(args.out_counts, sep='\t')

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

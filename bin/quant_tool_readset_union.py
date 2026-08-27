#!/usr/bin/env python3
"""
quant_tool_readset_union.py

Per-sample, per-tool read-level BSJ support, and the uniform read-ID-union
count derived from it. Each native tool file is an optional CLI arg: a
tool that was not run for this sample simply has no file to parse.

Why: each tool's own "expression" or readCount column is not comparable
across tools, since each has different junction-calling stringency and
counting logic. isocirc, circfl and ciri-long all keep the actual
supporting long-read IDs per BSJ, so one uniform metric, the number of
unique reads supporting a BSJ locus, can be recomputed by taking the
union of read-ID sets across whichever of these three tools called that
locus. circnick does not expose per-read IDs anywhere in its outputs
(only an aggregate BSJ_reads count), so it cannot be deduplicated against
the other tools' reads. It is tracked separately and used as a fallback
only for loci no read-based tool detected at all.

This module's tier1/tier1.5/tier2/tier3 remap-based quantification is the
primary, trusted count for downstream analysis. This union-of-tool-reads
count exists only as an independent older baseline that
quant_coverage_flag.py compares tier1 against.

Usage:
    quant_tool_readset_union.py \\
        --catalog   catalog.tsv \\
        --tsv       sample_discovery_confidence.tsv \\
        --sample    sample1 \\
        --isocirc_out    isocirc.out \\
        --circfl_pass    circFL_Normal_pass.txt \\
        --cirilong_reads sample1.reads \\
        --circnick_annotated sample1.circRNA_candidates.annotated.txt
"""
import argparse

import pandas as pd

STRAND_VARIANTS = ["+", "-", ".", "U"]  # tools encode "unstranded" differently ("." vs "U")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog', required=True, help='Locus catalog TSV (bsj_id/chrom/start/end/strand)')
    p.add_argument('--tsv', required=True, help='Discovery confidence TSV (per-tool 0/1 presence flags)')
    p.add_argument('--sample', required=True)
    p.add_argument('--isocirc_out', default=None, help='Native isocirc.out (omit if isocirc was not run)')
    p.add_argument('--circfl_pass', default=None, help='Native circFL_Normal_pass.txt (omit if circfl was not run)')
    p.add_argument('--cirilong_reads', default=None, help='Native CIRI-long .reads file (omit if cirilong was not run)')
    p.add_argument('--circnick_annotated', default=None, help='Native circnick annotated.txt (omit if circnick was not run)')
    return p.parse_args()


def parse_isocirc(f) -> dict:
    df = pd.read_csv(f, sep="\t", skiprows=1, dtype=str)
    df.columns = [c.lstrip("#") for c in df.columns]
    df["start"] = df["startCoor0based"].astype(int)
    df["end"] = df["endCoor"].astype(int)
    out = {}
    for (chrom, start, end, strand), g in df.groupby(["chrom", "start", "end", "geneStrand"]):
        reads = set()
        for ids in g["readIDs"].dropna():
            reads.update(x for x in ids.split(",") if x)
        out[(chrom, int(start), int(end), strand)] = reads
    return out


def parse_circfl(f) -> dict:
    df = pd.read_csv(f, sep="\t", dtype=str)
    df["start0"] = df["start"].astype(int) - 1
    df["end"] = df["end"].astype(int)
    out = {}
    for (chrom, start, end, strand), g in df.groupby(["chr", "start0", "end", "strand"]):
        reads = set()
        for ids in g["readID"].dropna():
            reads.update(x for x in ids.split(",") if x)
        out[(chrom, int(start), int(end), strand)] = reads
    return out


def parse_cirilong(f) -> dict:
    df = pd.read_csv(f, sep="\t", dtype=str)
    loc = df["circ_id"].str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
    df = df.join(loc)
    df["start0"] = df["start"].astype(int) - 1
    df["end"] = df["end"].astype(int)
    out = {}
    for (chrom, start, end, strand), g in df.groupby(["chrom", "start0", "end", "strand"]):
        out[(chrom, int(start), int(end), strand)] = set(g["read_id"])
    return out


def parse_circnick_native_counts(f) -> dict:
    df = pd.read_csv(f, sep="\t", dtype=str)
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["BSJ_reads"] = df["BSJ_reads"].astype(int)
    out = {}
    for _, row in df.iterrows():
        key = (row["chr"], row["start"], row["end"], row["strand"])
        out[key] = max(out.get(key, 0), row["BSJ_reads"])
    return out


def _lookup(readset: dict, chrom, start, end, strand):
    key = (chrom, start, end, strand)
    if key in readset:
        return readset[key]
    for s in STRAND_VARIANTS:
        if s == strand:
            continue
        alt = (chrom, start, end, s)
        if alt in readset:
            return readset[alt]
    return None


def load_tool_presence(tsv_path) -> pd.DataFrame:
    rows = []
    header_idx = {}
    with open(tsv_path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            if not header_idx:
                header = line.lstrip('#').split('\t')
                header_idx = {c: i for i, c in enumerate(header)}
                continue
            cols = line.split('\t')
            row = {'bsj_id': cols[header_idx['bsj_id']]}
            for tool in ('isocirc', 'circfl', 'cirilong', 'circnick'):
                row[tool] = int(cols[header_idx[tool]]) if tool in header_idx and cols[header_idx[tool]] not in ('', 'NA') else 0
            rows.append(row)
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    catalog = pd.read_csv(args.catalog, sep='\t')
    tool_presence = load_tool_presence(args.tsv)

    readsets = {"isocirc": {}, "circfl": {}, "cirilong": {}, "circnick": {}}
    if args.isocirc_out:
        readsets["isocirc"] = parse_isocirc(args.isocirc_out)
    if args.circfl_pass:
        readsets["circfl"] = parse_circfl(args.circfl_pass)
    if args.cirilong_reads:
        readsets["cirilong"] = parse_cirilong(args.cirilong_reads)
    if args.circnick_annotated:
        readsets["circnick"] = parse_circnick_native_counts(args.circnick_annotated)

    for tool in ("isocirc", "circfl", "cirilong"):
        n_reads = sum(len(v) for v in readsets[tool].values())
        print(f"  [{args.sample}] {tool:9s}: {len(readsets[tool]):6d} loci, {n_reads:7d} read-locus assignments")
    print(f"  [{args.sample}] circnick : {len(readsets['circnick']):6d} loci (native count only, no read IDs)")

    READ_TOOLS = ["isocirc", "circfl", "cirilong"]
    # A bsj_id can repeat across multiple isoform rows (multi-isoform recovery);
    # chrom/start/end/strand are identical across those rows, so any one suffices.
    coords = catalog.drop_duplicates(subset="bsj_id").set_index("bsj_id")[["chrom", "start", "end", "strand"]]

    union_counts, circnick_counts = {}, {}
    for row in tool_presence.itertuples(index=False):
        if row.bsj_id not in coords.index:
            continue
        chrom, start, end, strand = coords.loc[row.bsj_id]
        union_reads = set()
        any_read_tool_flagged = False
        for tool in READ_TOOLS:
            if getattr(row, tool, 0) == 1:
                any_read_tool_flagged = True
                rs = _lookup(readsets[tool], chrom, int(start), int(end), strand)
                if rs is not None:
                    union_reads |= rs
        if union_reads or any_read_tool_flagged:
            union_counts[row.bsj_id] = len(union_reads)
        if getattr(row, "circnick", 0) == 1:
            cn = _lookup(readsets["circnick"], chrom, int(start), int(end), strand)
            if cn is not None:
                circnick_counts[row.bsj_id] = cn

    # old_tool_count: row-wise max across the read-ID union and circnick's own
    # aggregate count. Either source can be the best available independent
    # estimate, depending on which tools called a given locus.
    all_ids = set(union_counts) | set(circnick_counts)
    rows = [{"bsj_id": bid, "old_tool_count": max(union_counts.get(bid, 0), circnick_counts.get(bid, 0))}
            for bid in all_ids]
    out_df = pd.DataFrame(rows)
    out_path = f'{args.sample}_old_tool_counts.tsv'
    out_df.to_csv(out_path, sep='\t', index=False)
    print(f"[{args.sample}] old-tool-count baseline: {len(out_df)} loci -> {out_path}")


if __name__ == '__main__':
    main()

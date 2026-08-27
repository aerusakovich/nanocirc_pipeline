#!/usr/bin/env python3
"""
filter_confidence.py
--------------------
Filters a scored confidence TSV + matching BED12 using the two independent
confidence axes produced by add_isoform_confidence.py:

    bsj_consensus     : quality of BSJ detection across tools (Low/Medium/High)
    isoform_consensus : quality of exon-structure agreement across tools

Filter modes:
  trusted_only  : Balanced mode.
                  Removes Low entries on each axis unless the call comes from
                  a trusted tool (--trusted_tools, default cirilong,isocirc,circfl),
                  kept even when confidence is Low.
                  - bsj_consensus Low: keep if bsj_source in trusted_tools
                  - isoform_consensus Low: keep if any tool in isoform_tools
                    is in trusted_tools
                  Other tools' Low calls are removed on both axes.

  no_low        : Removes all entries with Low on either axis, no exceptions.
                  Used only by the --run_benchmark_modes diagnostic filters.

  isocirc_only  : High-confidence mode (this pipeline's default precision tier).
                  Same as no_low, except IsoCirc is given a narrow exception
                  (IsoCirc alone is far more precise than any other single
                  tool -- see the balanced per-tool precision breakdown that
                  motivated this mode: isocirc-only calls run ~0.96-0.97
                  precision, vs. ~0.56 for cirilong-only, so a general
                  trusted_only exception here would just re-create balanced
                  and give back the precision gain).
                  - bsj_consensus Low: keep only if bsj_source == isocirc
                  - isoform_consensus Low: keep only if isocirc is present
                    in isoform_tools (isocirc is the sole priority tool on
                    this axis, not merely one of several trusted ones)
                  CIRI-long, circFL, CircNick-LRS Low calls are all removed
                  (CIRI-long is trusted elsewhere but not given a pass here,
                  since this mode's whole point is IsoCirc's own precision,
                  not a general trusted-tool exception).
                  Guard: the pipeline additionally routes this tier's output
                  through FILTER_CONFIDENT_DISCOVERY so any IsoCirc-only Low
                  call admitted here still needs read support above
                  params.circrna_confident_min_reads to survive.

  high_only     : Strictest filter, used only by the --run_benchmark_modes
                  diagnostic filters. Keeps only entries where both axes are
                  'High'. Medium is removed.

Usage:
    filter_confidence.py \\
        --bed    sample_smart_consensus.bed12 \\
        --tsv    sample_smart_consensus_confidence.tsv \\
        --mode   trusted_only \\
        --prefix sample_smart_consensus_filtered
"""

import argparse
import sys
import os

ISOCIRC_ONLY_TOOLS = {'isocirc'}
ISOCIRC_EXCEPTION_MODES = {'isocirc_only'}
NO_EXCEPTION_MODES = {'no_low', 'high_only'}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bed',    required=True, help='BED12 input file')
    p.add_argument('--tsv',    required=True, help='Scored confidence TSV input file')
    p.add_argument('--mode',   required=True,
                   choices=['trusted_only', 'no_low', 'isocirc_only', 'high_only'],
                   help='Filter mode')
    p.add_argument('--trusted_tools', default='cirilong,isocirc,circfl',
                   help='Comma-separated tool names trusted_only rescues Low calls '
                        'from (params.circrna_trusted_tools). Only affects '
                        'trusted_only; isocirc_only always trusts IsoCirc alone, '
                        'regardless of this flag.')
    p.add_argument('--prefix', required=True,
                   help='Output file prefix (outputs: prefix.bed12, prefix_confidence.tsv)')
    return p.parse_args()


def _bsj_low_passes(cols, idx, mode, trusted_tools):
    """
    Decide whether a Low bsj_consensus record is kept.
    no_low / high_only → always False.
    trusted_only → True only if the BSJ source is in trusted_tools.
    isocirc_only → True only if the BSJ source is IsoCirc.
    """
    if mode in NO_EXCEPTION_MODES:
        return False
    trusted = ISOCIRC_ONLY_TOOLS if mode in ISOCIRC_EXCEPTION_MODES else trusted_tools
    # check bsj_source (smart TSVs) or the relevant tool presence flag(s)
    if 'bsj_source' in idx:
        return cols[idx['bsj_source']] in trusted
    return any(t in idx and cols[idx[t]] == '1' for t in trusted)


def _isoform_low_passes(cols, idx, mode, trusted_tools):
    """
    Decide whether a Low isoform_consensus record is kept.
    no_low / high_only → always False.
    trusted_only → True only if isoform_tools contains a trusted tool.
    isocirc_only → True only if isoform_tools contains IsoCirc.
    """
    if mode in NO_EXCEPTION_MODES:
        return False
    trusted = ISOCIRC_ONLY_TOOLS if mode in ISOCIRC_EXCEPTION_MODES else trusted_tools
    # check isoform_tools column (smart TSVs)
    if 'isoform_tools' in idx:
        tools = set(cols[idx['isoform_tools']].split(','))
        return bool(tools & trusted)
    # fallback for non-smart TSVs: relevant tool presence flag(s)
    return any(t in idx and cols[idx[t]] == '1' for t in trusted)


def passes_filter(cols, idx, mode, trusted_tools):
    """
    Return True if this TSV row should be kept.
    idx: dict mapping column name → column index.

    Both axes must independently pass:
      bsj_consensus      != Low  (or Low + trusted/isocirc exception)
      isoform_consensus  != Low  (or Low + trusted/isocirc exception)

    high_only: Medium is always removed on both axes; only High passes
    (no exception modes apply, so a Low axis never survives here).
    """
    bsj_cons = cols[idx['bsj_consensus']]    if 'bsj_consensus'     in idx else 'NA'
    iso_cons = cols[idx['isoform_consensus']] if 'isoform_consensus' in idx else 'NA'

    if mode == 'high_only':
        bsj_ok = (bsj_cons == 'High') or (bsj_cons == 'Low' and _bsj_low_passes(cols, idx, mode, trusted_tools))
        iso_ok = (iso_cons == 'High') or (iso_cons == 'Low' and _isoform_low_passes(cols, idx, mode, trusted_tools))
        return bsj_ok and iso_ok

    bsj_ok = (bsj_cons != 'Low') or _bsj_low_passes(cols, idx, mode, trusted_tools)
    iso_ok = (iso_cons != 'Low') or _isoform_low_passes(cols, idx, mode, trusted_tools)

    return bsj_ok and iso_ok


def main():
    args = parse_args()
    trusted_tools = {t.strip() for t in args.trusted_tools.split(',') if t.strip()}

    kept_ids     = set()
    tsv_out      = []
    header_idx   = {}

    # ── Pass 1: filter TSV ────────────────────────────────────────────────────
    with open(args.tsv) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue

            # Header detection: starts with '#' or is the first non-empty line
            if line.startswith('#') or not header_idx:
                header = line.lstrip('#').split('\t')
                # Re-attach '#' to first column for output fidelity
                raw_header = line
                header_idx = {col: i for i, col in enumerate(header)}
                # Also index the raw '#chrom' as 'chrom'
                if line.startswith('#chrom'):
                    header_idx['chrom'] = 0
                tsv_out.append(raw_header)
                continue

            if 'bsj_consensus' not in header_idx and 'isoform_consensus' not in header_idx:
                print(
                    '[filter_confidence] WARNING: bsj_consensus / isoform_consensus '
                    'columns not found, TSV may not have been scored by '
                    'add_isoform_confidence.py. Passing all rows through.',
                    file=sys.stderr
                )
                tsv_out.append(line)
                continue

            cols   = line.split('\t')
            bsj_id = cols[4]   # bsj_id is always column 4 in all confidence TSVs

            if passes_filter(cols, header_idx, args.mode, trusted_tools):
                kept_ids.add(bsj_id)
                tsv_out.append(line)

    # ── Pass 2: filter BED12 by kept bsj_ids ─────────────────────────────────
    bed_out = []
    with open(args.bed) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line or line.startswith(('#', 'track', 'browser')):
                continue
            cols = line.split('\t')
            if cols[3] in kept_ids:   # BED12 name field (col 4, 0-indexed 3) = bsj_id
                bed_out.append(line)

    # ── Write outputs ─────────────────────────────────────────────────────────
    with open(f'{args.prefix}.bed12', 'w') as fh:
        if bed_out:
            fh.write('\n'.join(bed_out) + '\n')

    with open(f'{args.prefix}_confidence.tsv', 'w') as fh:
        if tsv_out:
            fh.write('\n'.join(tsv_out) + '\n')

    print(
        f'[filter_confidence] {args.mode}: '
        f'{len(kept_ids)} entries kept → {args.prefix}.bed12',
        file=sys.stderr
    )


if __name__ == '__main__':
    main()

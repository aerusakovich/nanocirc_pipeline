#!/usr/bin/env python
"""
smart_merge.py
--------------
Merges circRNA calls from several tools into one consensus BED12 file.
Calls are grouped by relaxed BSJ (backsplice junction), meaning start and
end can differ by up to tolerance bp, and must be on the same strand.

Four modes:

  smart_consensus         BSJ by majority vote (ties go to BSJ priority).
                          Structure by majority vote among tools that share
                          the winning BSJ (ties go to struct priority).

  smart_consensus_xstruct Same BSJ vote as smart_consensus. Structure vote
                          uses ALL tools in the group, comparing exon
                          positions in genome coordinates (within
                          struct_tolerance bp per boundary). A tool with a
                          slightly different BSJ can still vote for or
                          against the winning structure. If the winning
                          structure has no tool at the exact winning BSJ,
                          it is rebased (block_starts shifted by the small
                          BSJ offset, at most struct_tolerance bp).

  smart_consensus_hybrid  Same BSJ vote. Structure vote uses only tools at
                          the exact winning BSJ (no rebasing). Tools with a
                          different BSJ get their own isoform entry. This
                          mode makes the pipeline's 'discovery' output.

  smart_priority          BSJ always from the BSJ priority tool. Structure
                          always from the struct priority tool.

Within each BSJ group, every distinct (BSJ, structure) combination is kept.
The winning one is labelled 'main'. The rest are labelled 'iso1', 'iso2',
and so on. The BED12 name field shows this: 'chr:start-end:strand' for
main, 'chr:start-end:strand|iso1' for the others.

BSJ priority order:
    IsoCirc > circFL > CircNick-LRS > CIRI-long

Exon structure priority order:
    IsoCirc > circFL > CIRI-long > CircNick-LRS

Usage:
    smart_merge.py \\
        --sample   SAMPLE_ID \\
        --tool_names  cirilong isocirc circfl \\
        --bed_files   a.bed    b.bed   c.bed  \\
        --tolerance   5 \\
        --n_active    3 \\
        --outdir      results/
"""

import os
import sys
import argparse
from collections import defaultdict


# ── Priority tables ────────────────────────────────────────────────────────────
BSJ_PRIORITY    = ['isocirc', 'circfl', 'circnick', 'cirilong']
STRUCT_PRIORITY = ['isocirc', 'circfl', 'cirilong', 'circnick']


def parse_args(args=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample',     required=True)
    p.add_argument('--tool_names', required=True, nargs='+',
                   help='Tool names in the same order as --bed_files')
    p.add_argument('--bed_files',  required=True, nargs='+',
                   help='BED12 paths in the same order as --tool_names')
    p.add_argument('--tolerance',       type=int, default=5,
                   help='BSJ coordinate tolerance in bp (default: 5)')
    p.add_argument('--struct_tolerance', type=int, default=None,
                   help='Per-exon-boundary tolerance for cross-BSJ structure '
                        'voting in smart_consensus_xstruct (bp). '
                        'Defaults to --tolerance if not set.')
    p.add_argument('--n_active',   type=int, default=None,
                   help='Total active tools (used only for confidence info)')
    p.add_argument('--conf_tsvs',  nargs='+', default=None,
                   help='Confidence TSVs in the same order as --bed_files, one per '
                        'entry (each entry is one run\'s own prior smart_merge.py '
                        'output). Turns on cross-run merge mode for the '
                        'smart_consensus_hybrid output. See cross_run_hybrid_entries().')
    p.add_argument('--min_corroboration', type=int, default=2,
                   help='Cross-run merge mode only (--conf_tsvs given). Drop a '
                        'structure seen in only 1 run unless that run\'s own tool '
                        'agreement meets this number. A structure seen in 2+ runs is '
                        'always kept. Default 2.')
    p.add_argument('--outdir',     default='.',
                   help='Output directory (default: current dir)')
    return p.parse_args(args)


def print_error(msg):
    print(f'ERROR: {msg}', file=sys.stderr)
    sys.exit(1)


# ── I/O ────────────────────────────────────────────────────────────────────────

def read_bed12(path, tool_name):
    """Read BED12 and return list of record dicts."""
    records = []
    if not path or not os.path.exists(path):
        return records
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip('\n')
            if not line or line.startswith(('#', 'track', 'browser')):
                continue
            cols = line.split('\t')
            if len(cols) < 12:
                print_error(f'{tool_name} line {lineno}: expected 12 columns, got {len(cols)}')
            records.append({
                'chrom':        cols[0],
                'start':        int(cols[1]),
                'end':          int(cols[2]),
                'name':         cols[3],
                'score':        cols[4],
                'strand':       cols[5] if cols[5] in ('+', '-', '.') else '.',
                'thick_start':  cols[6],
                'thick_end':    cols[7],
                'rgb':          cols[8],
                'block_count':  cols[9],
                'block_sizes':  cols[10],
                'block_starts': cols[11],
                'tool':         tool_name,
            })
    return records


def load_conf_lookup(path):
    """Read one entry's confidence TSV (its 'main' row and every 'isoN' row)
    into {bsj_id: {'struct_agree': int, 'struct_source': str}}.
    The bsj_id key only works inside this same file. It was made together
    with the matching bed12 in the same run (see cross_run_hybrid_entries())."""
    lookup = {}
    if not path or not os.path.exists(path):
        return lookup
    import csv
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            bsj_id = row.get('bsj_id', '')
            if not bsj_id:
                continue
            try:
                struct_agree = int(row.get('struct_agree_count', '') or 0)
            except ValueError:
                struct_agree = 0
            lookup[bsj_id] = {
                'struct_agree': struct_agree,
                'struct_source': row.get('struct_source', '') or None,
            }
    return lookup


# ── Grouping ───────────────────────────────────────────────────────────────────

def group_relaxed(records, tolerance):
    """
    Group records by relaxed BSJ: same chrom and strand, start within
    tolerance, end within tolerance. Uses union-find over a start-sorted
    list, so a chain of near-duplicate calls (each one close to the next,
    even if the first and last are far apart) becomes one group instead of
    splitting into several.
    Returns {(chrom, start, end, strand): {tool: [record, ...]}}
    """
    parent = list(range(len(records)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    by_chrom_strand = defaultdict(list)
    for i, rec in enumerate(records):
        by_chrom_strand[(rec['chrom'], rec['strand'])].append(i)

    for idxs in by_chrom_strand.values():
        idxs_sorted = sorted(idxs, key=lambda i: records[i]['start'])
        for a in range(len(idxs_sorted)):
            for b in range(a + 1, len(idxs_sorted)):
                i, j = idxs_sorted[a], idxs_sorted[b]
                if records[j]['start'] - records[i]['start'] > tolerance:
                    break
                if abs(records[j]['end'] - records[i]['end']) <= tolerance:
                    union(i, j)

    groups = {}
    key_of_root = {}
    for i, rec in enumerate(records):
        root = find(i)
        if root not in key_of_root:
            key_of_root[root] = (rec['chrom'], rec['start'], rec['end'], rec['strand'])
        key = key_of_root[root]
        groups.setdefault(key, {}).setdefault(rec['tool'], []).append(rec)
    return groups


# ── Helpers ────────────────────────────────────────────────────────────────────

def best_record(recs):
    """Highest-score record from a list (read count proxy)."""
    def _score(r):
        try:
            return int(r['score'])
        except ValueError:
            return 0
    return max(recs, key=_score)


def bsj_key(rec):
    return (rec['start'], rec['end'])


def struct_key(rec):
    return (rec['block_count'], rec['block_sizes'], rec['block_starts'])


def make_bsj_id(chrom, start, end, strand, suffix=None):
    base = f'{chrom}:{start}-{end}:{strand}'
    return f'{base}|{suffix}' if suffix else base


def max_score_of(recs_dict):
    """Largest integer score across a {tool: record} dict."""
    best = 0
    for rec in recs_dict.values():
        try:
            best = max(best, int(rec['score']))
        except ValueError:
            pass
    return best


def _priority_tool(tools_present, priority_list):
    """Return highest-priority tool from a collection, or first if none listed."""
    for t in priority_list:
        if t in tools_present:
            return t
    return next(iter(tools_present))


# ── Absolute-coordinate structure helpers ──────────────────────────────────────

def absolute_exon_coords(rec):
    """Convert BED12 block_starts (relative) to (genome_start, size) pairs
    for each exon."""
    base   = rec['start']
    sizes  = [int(x) for x in rec['block_sizes'].split(',')  if x.strip()]
    starts = [int(x) for x in rec['block_starts'].split(',') if x.strip()]
    return tuple((base + s, sz) for s, sz in zip(starts, sizes))


def abs_struct_similar(a, b, tolerance):
    """Two structures are similar if they have the same number of exons,
    each exon start is within tolerance bp, and each exon size is equal."""
    if len(a) != len(b):
        return False
    return all(abs(ea[0] - eb[0]) <= tolerance and ea[1] == eb[1]
               for ea, eb in zip(a, b))


def group_by_abs_struct(tool_best, tolerance):
    """Group tools by exon structure similarity (genome coordinates).
    Returns a list of [representative_structure, [tool, ...]].
    The first tool seen in a group becomes its representative."""
    groups = []  # [ [abs_struct, [tools]] ]
    for tool, rec in tool_best.items():
        abs_struct = absolute_exon_coords(rec)
        matched = False
        for grp in groups:
            if abs_struct_similar(abs_struct, grp[0], tolerance):
                grp[1].append(tool)
                matched = True
                break
        if not matched:
            groups.append([abs_struct, [tool]])
    return groups


def vote_struct_groups(groups):
    """Vote over a list of [structure, [tools]] pairs.
    Most votes wins. Ties go to STRUCT_PRIORITY.
    Returns (winning_index, agree_count)."""
    max_count = max(len(tools) for _, tools in groups)
    winners   = [i for i, (_, tools) in enumerate(groups)
                 if len(tools) == max_count]

    if len(winners) == 1:
        return winners[0], max_count

    for prio_tool in STRUCT_PRIORITY:
        for idx in winners:
            if prio_tool in groups[idx][1]:
                return idx, max_count

    return winners[0], max_count


def rebase_struct(rec, new_start, new_end):
    """Return a copy of rec with block_starts shifted so the record's start
    becomes new_start. Exon genome positions stay the same, only the
    relative offsets change. Used when the winning structure has a
    slightly different BSJ than the winning BSJ (offset at most
    struct_tolerance bp)."""
    offset     = rec['start'] - new_start        # positive if rec is to the right
    old_starts = [int(x) for x in rec['block_starts'].split(',') if x.strip()]
    new_starts = [s + offset for s in old_starts]
    new_rec    = dict(rec)
    new_rec['start']        = new_start
    new_rec['end']          = new_end
    new_rec['block_starts'] = ','.join(str(s) for s in new_starts)
    return new_rec


# ── Voting ─────────────────────────────────────────────────────────────────────

def vote_majority(votes_dict, priority_list):
    """Given {key: [tool, ...]} vote counts, return (winning_key, agree_count).
    A strict majority (more than half) wins. If there is a tie, or no
    majority, fall back to priority_list."""
    n = sum(len(tools) for tools in votes_dict.values())
    max_count = max(len(tools) for tools in votes_dict.values())

    if max_count > n / 2:
        winners = [k for k, tools in votes_dict.items() if len(tools) == max_count]
        if len(winners) == 1:
            return winners[0], max_count
        # tie among top keys, use priority tool
        for tool in priority_list:
            for k in winners:
                if tool in votes_dict[k]:
                    return k, max_count

    # no majority, fall back to priority list
    all_tools_present = {t for tools in votes_dict.values() for t in tools}
    for tool in priority_list:
        if tool in all_tools_present:
            for k, tools in votes_dict.items():
                if tool in tools:
                    return k, len(tools)

    k = next(iter(votes_dict))
    return k, len(votes_dict[k])


# ── Entry collection ───────────────────────────────────────────────────────────

def collect_entries_consensus(tool_best):
    """Collect all isoform entries for smart_consensus mode.

    For each relaxed group:
      1. Majority vote on (start, end) gives the winning BSJ.
      2. Among tools that agree on this BSJ, majority vote on structure
         gives the main entry. Each other structure at this BSJ is its
         own iso entry.
      3. Tools with a different BSJ (within tolerance) get one iso entry
         per distinct BSJ, using the priority structure for each.

    Returns a list of entry dicts:
      start, end, struct_rec, bsj_src, struct_src,
      bsj_agree, struct_agree, isoform_label, isoform_tools
    """
    entries = []

    bsj_votes = defaultdict(list)
    for tool, rec in tool_best.items():
        bsj_votes[bsj_key(rec)].append(tool)

    winning_bsj, bsj_agree = vote_majority(bsj_votes, BSJ_PRIORITY)
    bsj_src = _priority_tool(set(bsj_votes[winning_bsj]), BSJ_PRIORITY)

    agree_map    = {t: r for t, r in tool_best.items() if bsj_key(r) == winning_bsj}
    disagree_bsj = defaultdict(dict)
    for t, r in tool_best.items():
        if bsj_key(r) != winning_bsj:
            disagree_bsj[bsj_key(r)][t] = r

    struct_votes = defaultdict(list)
    for t, r in agree_map.items():
        struct_votes[struct_key(r)].append(t)

    winning_sk, struct_agree = vote_majority(struct_votes, STRUCT_PRIORITY)

    main_tools = struct_votes[winning_sk]
    struct_src = _priority_tool(set(main_tools), STRUCT_PRIORITY)
    entries.append({
        'start':          winning_bsj[0],
        'end':            winning_bsj[1],
        'struct_rec':     agree_map[struct_src],
        'bsj_src':        bsj_src,
        'struct_src':     struct_src,
        'bsj_agree':      bsj_agree,
        'struct_agree':   struct_agree,
        'isoform_label':  'main',
        'isoform_tools':  list(main_tools),
    })

    iso_n = 0

    for sk, tools in struct_votes.items():
        if sk == winning_sk:
            continue
        iso_n += 1
        src = _priority_tool(set(tools), STRUCT_PRIORITY)
        entries.append({
            'start':         winning_bsj[0],
            'end':           winning_bsj[1],
            'struct_rec':    agree_map[src],
            'bsj_src':       bsj_src,
            'struct_src':    src,
            'bsj_agree':     bsj_agree,
            'struct_agree':  len(tools),
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(tools),
        })

    for bk, bk_tools in sorted(disagree_bsj.items()):
        iso_n += 1
        bsj_src_minor    = _priority_tool(set(bk_tools), BSJ_PRIORITY)
        struct_src_minor  = _priority_tool(set(bk_tools), STRUCT_PRIORITY)
        entries.append({
            'start':         bk[0],
            'end':           bk[1],
            'struct_rec':    bk_tools[struct_src_minor],
            'bsj_src':       bsj_src_minor,
            'struct_src':    struct_src_minor,
            'bsj_agree':     len(bk_tools),
            'struct_agree':  1,
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(bk_tools.keys()),
        })

    return entries


def collect_entries_consensus_xstruct(tool_best, struct_tolerance):
    """Like collect_entries_consensus, but the structure vote uses ALL
    tools in the relaxed BSJ group. Compares exon positions in genome
    coordinates, within struct_tolerance bp.

    1. BSJ vote gives the winning BSJ (same as smart_consensus).
    2. Group ALL tools by exon structure similarity.
    3. Most-votes wins on those groups, giving the winning structure.
    4. Main entry uses the winning BSJ coords and winning structure.
       - If a tool at the exact winning BSJ is in the winning group,
         use its record directly (no rebasing).
       - Otherwise, rebase the highest-priority tool from that group.
    5. Each other structure group becomes one isoform entry, placed at the
       winning BSJ coords (rebased if needed). There are no separate
       minority-BSJ isoforms here, since every tool is already placed by
       structure similarity.
    """
    entries = []

    bsj_votes = defaultdict(list)
    for tool, rec in tool_best.items():
        bsj_votes[bsj_key(rec)].append(tool)

    winning_bsj, bsj_agree = vote_majority(bsj_votes, BSJ_PRIORITY)
    bsj_src   = _priority_tool(set(bsj_votes[winning_bsj]), BSJ_PRIORITY)
    agree_map = {t: r for t, r in tool_best.items() if bsj_key(r) == winning_bsj}

    abs_groups = group_by_abs_struct(tool_best, struct_tolerance)

    winning_idx, struct_agree = vote_struct_groups(abs_groups)
    winning_sg_tools = abs_groups[winning_idx][1]

    bsj_tools_in_winner = [t for t in winning_sg_tools if t in agree_map]
    if bsj_tools_in_winner:
        struct_src = _priority_tool(set(bsj_tools_in_winner), STRUCT_PRIORITY)
        struct_rec = agree_map[struct_src]
    else:
        # rebase, winning structure has no tool at the exact winning BSJ
        struct_src = _priority_tool(set(winning_sg_tools), STRUCT_PRIORITY)
        struct_rec = rebase_struct(tool_best[struct_src],
                                   winning_bsj[0], winning_bsj[1])

    entries.append({
        'start':         winning_bsj[0],
        'end':           winning_bsj[1],
        'struct_rec':    struct_rec,
        'bsj_src':       bsj_src,
        'struct_src':    struct_src,
        'bsj_agree':     bsj_agree,
        'struct_agree':  struct_agree,
        'isoform_label': 'main',
        'isoform_tools': list(winning_sg_tools),
    })

    iso_n = 0
    for i, (_, sg_tools) in enumerate(abs_groups):
        if i == winning_idx:
            continue
        iso_n += 1
        bsj_tools_in_sg = [t for t in sg_tools if t in agree_map]
        if bsj_tools_in_sg:
            src  = _priority_tool(set(bsj_tools_in_sg), STRUCT_PRIORITY)
            srec = agree_map[src]
        else:
            src  = _priority_tool(set(sg_tools), STRUCT_PRIORITY)
            srec = rebase_struct(tool_best[src],
                                 winning_bsj[0], winning_bsj[1])
        entries.append({
            'start':         winning_bsj[0],
            'end':           winning_bsj[1],
            'struct_rec':    srec,
            'bsj_src':       bsj_src,
            'struct_src':    src,
            'bsj_agree':     bsj_agree,
            'struct_agree':  len(sg_tools),
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(sg_tools),
        })

    return entries


def collect_entries_consensus_hybrid(tool_best, struct_tolerance):
    """BSJ majority vote, then a structure vote using only the tools at the
    exact winning BSJ (no borrowing from other BSJs, no rebasing).

    1. BSJ vote: majority vote across all tools (same as consensus/xstruct).
    2. Structure vote: group only the winning-BSJ tools by exon structure
       similarity (within struct_tolerance bp), then most-votes wins.
    3. Tools with a different BSJ get their own isoform entry at their own
       BSJ coords (same as smart_consensus, so BSJ diversity is kept).

    Difference from consensus: small exon-boundary differences (within
                                struct_tolerance bp) count as the same
                                structure, not different ones.
    Difference from xstruct:    no rebasing. Tools with a different BSJ
                                never vote on the winning BSJ's structure.
    """
    entries = []

    bsj_votes = defaultdict(list)
    for tool, rec in tool_best.items():
        bsj_votes[bsj_key(rec)].append(tool)

    winning_bsj, bsj_agree = vote_majority(bsj_votes, BSJ_PRIORITY)
    bsj_src = _priority_tool(set(bsj_votes[winning_bsj]), BSJ_PRIORITY)

    agree_map    = {t: r for t, r in tool_best.items() if bsj_key(r) == winning_bsj}
    disagree_bsj = defaultdict(dict)
    for t, r in tool_best.items():
        if bsj_key(r) != winning_bsj:
            disagree_bsj[bsj_key(r)][t] = r

    # structure vote using genome coords, exact-BSJ tools only
    abs_groups = group_by_abs_struct(agree_map, struct_tolerance)
    winning_idx, struct_agree = vote_struct_groups(abs_groups)
    winning_sg_tools = abs_groups[winning_idx][1]

    struct_src = _priority_tool(set(winning_sg_tools), STRUCT_PRIORITY)
    entries.append({
        'start':         winning_bsj[0],
        'end':           winning_bsj[1],
        'struct_rec':    agree_map[struct_src],
        'bsj_src':       bsj_src,
        'struct_src':    struct_src,
        'bsj_agree':     bsj_agree,
        'struct_agree':  struct_agree,
        'isoform_label': 'main',
        'isoform_tools': list(winning_sg_tools),
    })

    iso_n = 0
    for i, (_, sg_tools) in enumerate(abs_groups):
        if i == winning_idx:
            continue
        iso_n += 1
        src = _priority_tool(set(sg_tools), STRUCT_PRIORITY)
        entries.append({
            'start':         winning_bsj[0],
            'end':           winning_bsj[1],
            'struct_rec':    agree_map[src],
            'bsj_src':       bsj_src,
            'struct_src':    src,
            'bsj_agree':     bsj_agree,
            'struct_agree':  len(sg_tools),
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(sg_tools),
        })

    for bk, bk_tools in sorted(disagree_bsj.items()):
        iso_n += 1
        bsj_src_minor    = _priority_tool(set(bk_tools), BSJ_PRIORITY)
        struct_src_minor = _priority_tool(set(bk_tools), STRUCT_PRIORITY)
        entries.append({
            'start':         bk[0],
            'end':           bk[1],
            'struct_rec':    bk_tools[struct_src_minor],
            'bsj_src':       bsj_src_minor,
            'struct_src':    struct_src_minor,
            'bsj_agree':     len(bk_tools),
            'struct_agree':  1,
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(bk_tools.keys()),
        })

    return entries


def collect_entries_priority(tool_best):
    """Collect all isoform entries for smart_priority mode.

    Main entry:
      BSJ from the highest-priority tool present.
      Structure from the highest struct-priority tool present.

    Every other unique (BSJ, structure) pair becomes its own isoform entry.

    Returns the same entry dict format as collect_entries_consensus.
    """
    entries = []

    bsj_src    = _priority_tool(set(tool_best.keys()), BSJ_PRIORITY)
    struct_src = _priority_tool(set(tool_best.keys()), STRUCT_PRIORITY)

    main_bsj = bsj_key(tool_best[bsj_src])
    main_sk  = struct_key(tool_best[struct_src])

    bsj_agree    = sum(1 for r in tool_best.values() if bsj_key(r)    == main_bsj)
    struct_agree = sum(1 for r in tool_best.values() if struct_key(r) == main_sk)

    # isoform_tools = tools that share BOTH winning BSJ AND winning struct
    main_tools = [t for t, r in tool_best.items()
                  if bsj_key(r) == main_bsj and struct_key(r) == main_sk]
    # bsj_src and struct_src must be represented even if they differ
    main_tools_set = set(main_tools) | {bsj_src, struct_src}

    entries.append({
        'start':         main_bsj[0],
        'end':           main_bsj[1],
        'struct_rec':    tool_best[struct_src],
        'bsj_src':       bsj_src,
        'struct_src':    struct_src,
        'bsj_agree':     bsj_agree,
        'struct_agree':  struct_agree,
        'isoform_label': 'main',
        'isoform_tools': list(main_tools_set),
    })

    # Collect remaining unique (bsj_key, struct_key) combos
    combos = defaultdict(list)
    for t, r in tool_best.items():
        combo = (bsj_key(r), struct_key(r))
        if combo != (main_bsj, main_sk):
            combos[combo].append(t)

    iso_n = 0
    for (bk, sk), tools in sorted(combos.items()):
        iso_n += 1
        src          = _priority_tool(set(tools), STRUCT_PRIORITY)
        bsj_src_iso  = _priority_tool(set(tools), BSJ_PRIORITY)
        bsj_agree_iso    = sum(1 for r in tool_best.values() if bsj_key(r) == bk)
        struct_agree_iso = sum(1 for r in tool_best.values() if struct_key(r) == sk)
        entries.append({
            'start':         bk[0],
            'end':           bk[1],
            'struct_rec':    tool_best[src],
            'bsj_src':       bsj_src_iso,
            'struct_src':    src,
            'bsj_agree':     bsj_agree_iso,
            'struct_agree':  struct_agree_iso,
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': list(tools),
        })

    return entries


def cross_run_hybrid_entries(tool_map, struct_tolerance, min_corroboration=2):
    """Cross-run version of collect_entries_consensus_hybrid(). Used when
    --conf_tsvs is given. Here, each entry is a run, and each run's own
    input is already a merged output (main plus isoN rows), not a single
    raw call per tool.

    This function uses two separate votes, like the normal hybrid mode's
    BSJ vote then structure vote:

    1. BSJ vote (position only, ignores structure). The winning BSJ is the
       exact (start, end) position backed by the most different runs. If a
       run gives several of its own records at the same position, that
       counts once, not once per record. Ties are broken by the summed
       struct_agree_count at that position. bsj_agree is the resulting
       run count. This must stay separate from the structure vote:
       crossrun_annotate.py's tier filter reads bsj_confidence as "how many
       runs support this BSJ", not "how many runs support the winning
       structure". Mixing the two would quietly change balanced_recall and
       high_confidence tier membership in a way nobody checked against
       ground truth.

    2. Structure vote, using only records at the winning BSJ position (no
       borrowing from other positions, same rule collect_entries_consensus_hybrid
       already uses). Every record at that position is grouped by exon
       structure (like group_by_abs_struct). A group's own isoform label
       ('iso1', 'iso2') does not matter here: 'iso1' from one run has no
       link to 'iso1' from another. Each group's rank is the SUM of
       struct_agree_count across all its records (total tool agreement for
       that exact structure, across every run that has it there). This is
       not a raw record count and not raw score, so a structure with
       strong agreement in a few runs can beat one with weak (for example
       single-tool) agreement spread over many runs. Ties go to
       STRUCT_PRIORITY of the best contributing source.

       One weak case is dropped: a structure group backed by only one run,
       where that run's own tool agreement is below min_corroboration.
       Tests showed this cuts false positives by about 34% for a small
       recall cost. A group backed by 2 or more runs, or by min_corroboration or
       more tools within one run, is always kept.

    Minority BSJ positions (every position other than the winner) each get
    their own isoform entry at their own coordinates. This matches
    collect_entries_consensus_hybrid's own handling of BSJ disagreement, so
    BSJ diversity stays in the isoform list. Each minority position picks
    its best-ranked structure and is checked against the same weak-case
    rule above.

    tool_map: {entry_name: [record, ...]}. ALL records for this entry in
    this relaxed BSJ group (not reduced to one via best_record()).

    Returns entries in the same dict format as the other collect_entries_*
    functions.
    """
    all_records = [rec for recs in tool_map.values() for rec in recs]

    bsj_positions = defaultdict(list)
    for rec in all_records:
        bsj_positions[bsj_key(rec)].append(rec)

    def position_runs(recs):
        return {r['tool'] for r in recs}

    def position_weight(recs):
        return sum(r.get('struct_agree', 0) or 0 for r in recs)

    def struct_clusters(recs):
        clusters = []  # [ [abs_struct, [record, ...]] ]
        for rec in recs:
            abs_struct = absolute_exon_coords(rec)
            matched = False
            for grp in clusters:
                if abs_struct_similar(abs_struct, grp[0], struct_tolerance):
                    grp[1].append(rec)
                    matched = True
                    break
            if not matched:
                clusters.append([abs_struct, [rec]])
        return clusters

    def cluster_weight(cl):
        return position_weight(cl[1])

    def cluster_priority_rank(cl):
        best = len(STRUCT_PRIORITY)
        for r in cl[1]:
            src = r.get('struct_source')
            if src in STRUCT_PRIORITY:
                best = min(best, STRUCT_PRIORITY.index(src))
        return best

    def cluster_representative(cl):
        def rec_key(r):
            src = r.get('struct_source')
            src_rank = STRUCT_PRIORITY.index(src) if src in STRUCT_PRIORITY else len(STRUCT_PRIORITY)
            run_rank = BSJ_PRIORITY.index(r['tool']) if r['tool'] in BSJ_PRIORITY else len(BSJ_PRIORITY)
            return (src_rank, -(r.get('struct_agree', 0) or 0), run_rank)
        return min(cl[1], key=rec_key)

    def rank_clusters(recs):
        return sorted(struct_clusters(recs), key=lambda cl: (-cluster_weight(cl), cluster_priority_rank(cl)))

    def weak_evidence(cl):
        return len(position_runs(cl[1])) == 1 and cluster_weight(cl) < min_corroboration

    # ── 1. BSJ vote (position only) ──
    winning_bsj = max(bsj_positions, key=lambda k: (len(position_runs(bsj_positions[k])), position_weight(bsj_positions[k])))
    winning_records = bsj_positions[winning_bsj]
    bsj_agree = len(position_runs(winning_records))
    bsj_src = cluster_representative([None, winning_records])['tool']

    # ── 2. Structure vote, restricted to the winning BSJ's own records ──
    kept = [cl for cl in rank_clusters(winning_records) if not weak_evidence(cl)]

    entries = []
    for i, cl in enumerate(kept):
        rep = cluster_representative(cl)
        runs = sorted(position_runs(cl[1]))
        entries.append({
            'start':         winning_bsj[0],
            'end':           winning_bsj[1],
            'struct_rec':    rep,
            'bsj_src':       bsj_src,
            'struct_src':    rep['tool'],
            'bsj_agree':     bsj_agree,
            'struct_agree':  cluster_weight(cl),
            'isoform_label': 'main' if i == 0 else f'iso{i}',
            'isoform_tools': runs,
        })

    # ── 3. Minority-BSJ positions: one entry each, at their own coords ──
    iso_n = len(entries) - 1
    for k in sorted(bsj_positions):
        if k == winning_bsj:
            continue
        recs = bsj_positions[k]
        if len(position_runs(recs)) == 1 and position_weight(recs) < min_corroboration:
            continue
        best_cl = rank_clusters(recs)[0]
        rep = cluster_representative(best_cl)
        runs = sorted(position_runs(recs))
        iso_n += 1
        entries.append({
            'start':         k[0],
            'end':           k[1],
            'struct_rec':    rep,
            'bsj_src':       rep['tool'],
            'struct_src':    rep['tool'],
            'bsj_agree':     len(runs),
            'struct_agree':  cluster_weight(best_cl),
            'isoform_label': f'iso{iso_n}',
            'isoform_tools': runs,
        })

    return entries


# ── Output ─────────────────────────────────────────────────────────────────────

def bed12_line(chrom, start, end, name, score, strand,
               block_count, block_sizes, block_starts):
    return '\t'.join([
        chrom, str(start), str(end), name, str(score), strand,
        str(start), str(start), '0',
        block_count, block_sizes, block_starts,
    ])


def write_outputs(groups, active_tools, sample, outdir, struct_tolerance, cross_run_mode=False, min_corroboration=2):
    os.makedirs(outdir, exist_ok=True)

    tool_flags  = active_tools
    tool_blocks = []
    for t in active_tools:
        tool_blocks += [f'{t}_block_sizes', f'{t}_block_starts']

    header = (
        ['#chrom', 'start', 'end', 'strand', 'bsj_id', 'bsj_confidence']
        + tool_flags
        + tool_blocks
        + ['isoform_confidence']
        + ['sel_block_count', 'sel_block_sizes', 'sel_block_starts']
        + ['bsj_source', 'struct_source', 'struct_agree_count', 'max_score']
        + ['isoform_label', 'isoform_tools']
    )

    mode_lines = {
        'consensus':         {'bed': [], 'tsv': ['\t'.join(header)]},
        'consensus_xstruct': {'bed': [], 'tsv': ['\t'.join(header)]},
        'consensus_hybrid':  {'bed': [], 'tsv': ['\t'.join(header)]},
        'priority':          {'bed': [], 'tsv': ['\t'.join(header)]},
    }

    for group_key, tool_map in sorted(groups.items()):
        chrom, _, _, strand = group_key

        # Best record per tool (within this group)
        tool_best = {tool: best_record(recs) for tool, recs in tool_map.items()}
        group_score = max_score_of(tool_best)

        # Per-tool flags and block columns (group-level, same for all isoforms)
        flags  = ['1' if t in tool_best else '0' for t in active_tools]
        blocks = []
        for t in active_tools:
            if t in tool_best:
                blocks += [tool_best[t]['block_sizes'], tool_best[t]['block_starts']]
            else:
                blocks += ['.', '.']

        xstruct_fn = lambda tb: collect_entries_consensus_xstruct(tb, struct_tolerance)
        hybrid_fn  = (lambda tb: cross_run_hybrid_entries(tool_map, struct_tolerance, min_corroboration)) if cross_run_mode \
                     else (lambda tb: collect_entries_consensus_hybrid(tb, struct_tolerance))
        for mode, collect_fn in [
            ('consensus',         collect_entries_consensus),
            ('consensus_xstruct', xstruct_fn),
            ('consensus_hybrid',  hybrid_fn),
            ('priority',          collect_entries_priority),
        ]:
            entries = collect_fn(tool_best)

            for entry in entries:
                start   = entry['start']
                end     = entry['end']
                srec    = entry['struct_rec']
                label   = entry['isoform_label']
                iso_tools_str = ','.join(sorted(entry['isoform_tools']))

                bsj_id = make_bsj_id(chrom, start, end, strand,
                                      suffix=(label if label != 'main' else None))

                # score: for cross-run hybrid entries, srec is the actual
                # record that won this structure's group, so use its own
                # score directly. tool_best[t] could be a different record
                # from the same entry, with a different (higher) score,
                # which would wrongly give this structure another one's
                # read support.
                if cross_run_mode and mode == 'consensus_hybrid':
                    try:
                        score = int(srec.get('score', group_score))
                    except (TypeError, ValueError):
                        score = group_score
                else:
                    iso_tool_recs = {t: tool_best[t] for t in entry['isoform_tools']
                                     if t in tool_best}
                    score = max_score_of(iso_tool_recs) if iso_tool_recs else group_score

                bed = bed12_line(
                    chrom, start, end, bsj_id, score, strand,
                    srec['block_count'], srec['block_sizes'], srec['block_starts']
                )

                tsv_row = '\t'.join(
                    [chrom, str(start), str(end), strand, bsj_id,
                     str(entry['bsj_agree'])]    # bsj_confidence = tool-agreement count
                    + flags
                    + blocks
                    + ['']                        # isoform_confidence placeholder
                    + [srec['block_count'], srec['block_sizes'], srec['block_starts']]
                    + [entry['bsj_src'], entry['struct_src'],
                       str(entry['struct_agree']), str(score)]  # max_score last
                    + [label, iso_tools_str]
                )

                mode_lines[mode]['bed'].append(bed)
                mode_lines[mode]['tsv'].append(tsv_row)

    outputs = {
        f'{sample}_smart_consensus.bed12':                mode_lines['consensus']['bed'],
        f'{sample}_smart_consensus_confidence.tsv':       mode_lines['consensus']['tsv'],
        f'{sample}_smart_consensus_xstruct.bed12':        mode_lines['consensus_xstruct']['bed'],
        f'{sample}_smart_consensus_xstruct_confidence.tsv': mode_lines['consensus_xstruct']['tsv'],
        f'{sample}_smart_consensus_hybrid.bed12':         mode_lines['consensus_hybrid']['bed'],
        f'{sample}_smart_consensus_hybrid_confidence.tsv': mode_lines['consensus_hybrid']['tsv'],
        f'{sample}_smart_priority.bed12':                 mode_lines['priority']['bed'],
        f'{sample}_smart_priority_confidence.tsv':        mode_lines['priority']['tsv'],
    }
    for fname, lines in outputs.items():
        with open(os.path.join(outdir, fname), 'w') as fh:
            if lines:
                fh.write('\n'.join(lines) + '\n')

    n_groups = len(groups)
    for mode in ('consensus', 'consensus_xstruct', 'consensus_hybrid', 'priority'):
        n_all = len(mode_lines[mode]['bed'])
        print(
            f'[smart_merge] {sample}: {n_groups} groups → smart_{mode}: {n_all} isoform entries',
            file=sys.stderr
        )


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    if len(args.tool_names) != len(args.bed_files):
        print_error('--tool_names and --bed_files must have the same number of entries')

    cross_run_mode = args.conf_tsvs is not None
    if cross_run_mode and len(args.conf_tsvs) != len(args.bed_files):
        print_error('--conf_tsvs must have the same number of entries as --bed_files')

    all_records  = []
    active_tools = []
    conf_paths = args.conf_tsvs if cross_run_mode else [None] * len(args.bed_files)
    for tool, path, conf_path in zip(args.tool_names, args.bed_files, conf_paths):
        recs = read_bed12(path, tool)
        if cross_run_mode and recs:
            conf_lookup = load_conf_lookup(conf_path)
            for rec in recs:
                info = conf_lookup.get(rec['name'], {})
                rec['struct_agree']  = info.get('struct_agree', 0)
                rec['struct_source'] = info.get('struct_source')
        if recs:
            all_records.extend(recs)
            active_tools.append(tool)

    if len(active_tools) < 2:
        print(
            f'[smart_merge] Only {len(active_tools)} active tool(s), need at least 2. Writing empty outputs.',
            file=sys.stderr
        )
        os.makedirs(args.outdir, exist_ok=True)
        for suffix in [
            '_smart_consensus.bed12',
            '_smart_consensus_confidence.tsv',
            '_smart_consensus_xstruct.bed12',
            '_smart_consensus_xstruct_confidence.tsv',
            '_smart_priority.bed12',
            '_smart_priority_confidence.tsv',
        ]:
            open(os.path.join(args.outdir, args.sample + suffix), 'w').close()
        return

    # Keep active_tools in priority order so TSV columns are consistent
    ordered = [t for t in (BSJ_PRIORITY + [t for t in active_tools if t not in BSJ_PRIORITY])
               if t in active_tools]

    struct_tol = args.struct_tolerance if args.struct_tolerance is not None else args.tolerance
    groups = group_relaxed(all_records, args.tolerance)
    write_outputs(groups, ordered, args.sample, args.outdir, struct_tol,
                  cross_run_mode=cross_run_mode, min_corroboration=args.min_corroboration)


if __name__ == '__main__':
    main()

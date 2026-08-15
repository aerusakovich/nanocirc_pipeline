#!/usr/bin/env python3
"""
quant_common.py

Shared read-classification primitives used by the quant_*.py scripts
(quant_chunk_remap.py, quant_overlap_rescue.py, quant_tier2_rescue.py,
quant_tier3_rescue.py). This file exists so those scripts share one copy
of the classification logic instead of three separate copies.

Not a Nextflow entry point itself. It is imported by the other
quant_*.py scripts, which sit alongside it in bin/ so it resolves via
sys.path[0].
"""
import subprocess
import sys
from pathlib import Path
from collections import defaultdict

import pysam

ANCHOR_BP = 20              # window checked on each side of the join
MIN_ANCHOR_MATCH = 10       # of ANCHOR_BP, min aligned (non-clipped) bp required on each side
MAX_EDITS_PER_SIDE = 15     # max combined mismatches+indels within ANCHOR_BP on each side
MIN_ALIGNED_FRACTION = 0.3  # min fraction of query length that must be aligned (not S/H) overall

FLANK = 3000            # bp; matches circle_refs.py's own read-length assumption
JOIN_PROXIMITY = 50     # bp; how close a split segment's edge must be to join_pos
MIN_SEGMENT_MATCH = 50  # bp; PSL per-block floor for a block to count as evidence


def run_checked(cmd, **kwargs):
    """subprocess.run(check=True) wrapper that propagates a signal-killed
    child's real exit code (128+signum, standard POSIX convention), instead
    of letting an uncaught CalledProcessError mask it as Python's own
    generic exit code 1. This matters because Nextflow's automatic OOM
    retry (conf/base.config's `errorStrategy = { task.exitStatus in
    (130..145)+104+175 ? 'retry' : 'finish' }`) keys off the real exit
    code, for example 137 for a SIGKILL'd minimap2. It never fires if a
    subprocess failure here surfaces as a flat, uninformative 1 instead."""
    try:
        return subprocess.run(cmd, check=True, **kwargs)
    except subprocess.CalledProcessError as e:
        if e.returncode < 0:
            sys.exit(128 - e.returncode)
        raise


def read_fasta(path) -> dict:
    seqs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:], []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def run_minimap2_bam(minimap2_bin, samtools_bin, ref_fa, reads_fa_or_fq, out_bam, threads=4,
                      preset_args=("-ax", "map-ont", "--eqx", "-N", "5")):
    """minimap2 (BAM-mode CIGAR output, --eqx so the anchor filter can tell
    match from mismatch without an MD tag) piped straight into samtools view -b."""
    minimap2 = subprocess.Popen(
        [str(minimap2_bin), *preset_args, "-t", str(threads), str(ref_fa), str(reads_fa_or_fq)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    with open(out_bam, "wb") as out:
        samtools = subprocess.Popen(
            [str(samtools_bin), "view", "-b", "-@", str(threads)],
            stdin=minimap2.stdout, stdout=out, stderr=subprocess.DEVNULL,
        )
        minimap2.stdout.close()
        samtools.communicate()
    minimap2.wait()
    return Path(out_bam)


def run_pblat(pblat_bin, ref_fa, reads_fa, out_psl, threads=4, min_score=50, min_identity=90):
    run_checked(
        [str(pblat_bin), str(ref_fa), str(reads_fa), str(out_psl),
         "-threads=" + str(threads), f"-minScore={min_score}", f"-minIdentity={min_identity}", "-noHead"],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    return Path(out_psl)


def parse_bam_hits(bam_path):
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    for r in bam.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary:
            continue
        yield (r.query_name, r.reference_name, r.reference_start, r.reference_end,
               r.cigartuples, r.query_length or r.infer_read_length() or 0)


def parse_psl_hits(psl_path, min_segment_len=MIN_SEGMENT_MATCH):
    """Parses each individual PSL block (not the best hit, not the whole
    tStart/tEnd span), keeping only blocks >= min_segment_len bp. This
    avoids spurious hits from large gapped records where many tiny
    fragments sum past a per-line floor."""
    if not Path(psl_path).exists():
        return
    with open(psl_path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 21:
                continue
            qname, rname = fields[9], fields[13]
            block_sizes = [int(x) for x in fields[18].strip(",").split(",") if x]
            t_starts = [int(x) for x in fields[20].strip(",").split(",") if x]
            for size, t_start in zip(block_sizes, t_starts):
                if size >= min_segment_len:
                    yield qname, rname, t_start, t_start + size


# CIGAR op codes (pysam): 0=M 1=I 2=D 3=N 4=S 5=H 6=P 7='=' 8='X'
def edits_and_match_in_window(cigartuples, ref_start, win_start, win_end):
    """Walk a segment's CIGAR from ref_start; return (matched_bp, edit_bp) that
    fall within the reference-coordinate window [win_start, win_end)."""
    ref_pos = ref_start
    matched = edits = 0
    for op, length in cigartuples:
        if ref_pos >= win_end:
            break
        if op in (0, 7):  # M (unexpected w/ --eqx, handled defensively) or '='
            matched += max(0, min(ref_pos + length, win_end) - max(ref_pos, win_start))
            ref_pos += length
        elif op in (2, 8):  # D or X
            edits += max(0, min(ref_pos + length, win_end) - max(ref_pos, win_start))
            ref_pos += length
        elif op == 3:  # N
            ref_pos += length
        elif op == 1:  # I, does not consume reference; charge at current ref_pos
            if win_start <= ref_pos < win_end:
                edits += length
        # S(4)/H(5)/P(6): no reference consumption, ignore
    return matched, edits


def aligned_fraction(cigartuples, query_length):
    if not query_length:
        return 0.0
    aligned = sum(length for op, length in cigartuples if op in (0, 1, 7, 8))  # M/I/=/X consume query
    return aligned / query_length


def passes_anchor_filter(segs, jp):
    """segs: list of (ref_start, ref_end, cigartuples, query_length) for one read's
    segments (primary + supplementary) against one reference. Require a clean,
    well-anchored match immediately flanking jp on both sides, plus an overall
    not-mostly-clipped guard.

    Per-window quality is taken from the best single segment covering that
    window, not summed across all segments. Summing double-counts a
    genuine read's overlapping primary and supplementary segments on
    short tandem-duplicated references, which rejected 100% of
    true-origin reads at a real 158bp-circle locus under a summing
    scheme."""
    win_before = (jp - ANCHOR_BP, jp)
    win_after = (jp, jp + ANCHOR_BP)
    m_before = e_before = m_after = e_after = 0
    qlen = 0
    aligned_bp_sum = 0
    for ref_start, ref_end, cigartuples, seg_qlen in segs:
        if cigartuples is None:
            continue
        if ref_end > win_before[0] and ref_start < win_before[1]:
            m, e = edits_and_match_in_window(cigartuples, ref_start, *win_before)
            if m > m_before or (m == m_before and e < e_before):
                m_before, e_before = m, e
        if ref_end > win_after[0] and ref_start < win_after[1]:
            m, e = edits_and_match_in_window(cigartuples, ref_start, *win_after)
            if m > m_after or (m == m_after and e < e_after):
                m_after, e_after = m, e
        qlen = max(qlen, seg_qlen)
        aligned_bp_sum += aligned_fraction(cigartuples, seg_qlen) * seg_qlen
    anchor_ok = (m_before >= MIN_ANCHOR_MATCH and e_before <= MAX_EDITS_PER_SIDE
                 and m_after >= MIN_ANCHOR_MATCH and e_after <= MAX_EDITS_PER_SIDE)
    clip_ok = qlen > 0 and (aligned_bp_sum / qlen) >= MIN_ALIGNED_FRACTION
    return anchor_ok and clip_ok


def proximity_gate(segs, jp):
    """segs: list of (start, end, ...). A segment merely ending before /
    starting after the join, however far away, is not evidence of crossing
    it on a tandem-duplicated reference (copy1/copy2 are identical sequence,
    so an ordinary read overlapping one interior exon can otherwise look
    like it has segments on 'both sides' of the join purely by chance)."""
    crosses = any(s < jp < e for s, e, *_ in segs)
    before_near = any(jp - JOIN_PROXIMITY <= e <= jp for s, e, *_ in segs)
    after_near = any(jp <= s <= jp + JOIN_PROXIMITY for s, e, *_ in segs)
    return crosses or (before_near and after_near)


def anchor_coverage_ok(segs, jp):
    """Coordinate-only (no CIGAR needed) approximation of the anchor-quality
    check, for plain (start,end) segments from either BAM or PSL hits:
    >=MIN_ANCHOR_MATCH bp of block/segment coverage within ANCHOR_BP of the
    join on each side."""
    win_before, win_after = (jp - ANCHOR_BP, jp), (jp, jp + ANCHOR_BP)
    cov_before = sum(max(0, min(e, win_before[1]) - max(s, win_before[0])) for s, e in segs)
    cov_after = sum(max(0, min(e, win_after[1]) - max(s, win_after[0])) for s, e in segs)
    return cov_before >= MIN_ANCHOR_MATCH and cov_after >= MIN_ANCHOR_MATCH


def score_chunk_bam(hits_iter, join_pos_by_ref):
    """BAM/minimap2 hits qualifying under the JOIN_PROXIMITY gate and
    passes_anchor_filter (edit-count/clip quality, the strictest check
    available since real CIGARs are present). Returns
    {(qname, rname): score} (score = total matched span) for every
    qualifying pair, rather than collapsing straight to a per-locus
    supporting set. Keeping a score, instead of a plain boolean
    qualify-or-not, lets quant_chunk_remap.py resolve a read's assignment
    globally across every chunk in the catalog, not just within the one
    chunk this locus happened to land in. This closes the double-crediting
    bug where the same physical read independently qualifies for two
    genomically-related loci that the fixed-seed chunk shuffle happened to
    split apart (see quant_chunk_remap.py's module docstring)."""
    segments = defaultdict(lambda: defaultdict(list))
    for qname, rname, ref_start, ref_end, cigartuples, qlen in hits_iter:
        if rname not in join_pos_by_ref:
            continue
        segments[rname][qname].append((ref_start, ref_end, cigartuples, qlen))

    scores = {}
    for rname, reads in segments.items():
        jp = join_pos_by_ref[rname]
        for qname, segs in reads.items():
            if not proximity_gate(segs, jp):
                continue
            if not passes_anchor_filter(segs, jp):
                continue
            scores[(qname, rname)] = sum(e - s for s, e, *_ in segs)
    return scores


def score_chunk_psl(hits_iter, join_pos_by_ref):
    """BLAT/PSL hits (already per-block, >=MIN_SEGMENT_MATCH bp, from
    parse_psl_hits) qualifying under the JOIN_PROXIMITY gate and a
    block-coverage plus minimum-anchor-bp approximation of the
    anchor-quality check (PSL has no CIGAR). Returns
    {(qname, rname): score}, same rationale as score_chunk_bam."""
    segments = defaultdict(lambda: defaultdict(list))
    for qname, rname, t_start, t_end in hits_iter:
        if rname not in join_pos_by_ref:
            continue
        segments[rname][qname].append((t_start, t_end))

    scores = {}
    for rname, reads in segments.items():
        jp = join_pos_by_ref[rname]
        for qname, segs in reads.items():
            if not proximity_gate(segs, jp):
                continue
            if not anchor_coverage_ok(segs, jp):
                continue
            scores[(qname, rname)] = sum(e - s for s, e in segs)
    return scores


DIAGNOSTIC_KMER = 21  # bp; k-mer length used to find sequence unique to one sibling reference


def _kmer_set(seq: str, k: int) -> set:
    return {seq[i:i + k] for i in range(0, len(seq) - k + 1)}


def build_diagnostic_regions(ref_seqs: dict, k: int = DIAGNOSTIC_KMER) -> dict:
    """For a group of candidate reference sequences that are meant to
    compete for the same reads (same-BSJ sibling isoforms, or gene-family
    paralogs), find the reference-local positions in each one whose 
    sequence is absent from every OTHER member of the group.

    Returns {ref_id: sorted list of (start, end) diagnostic intervals} in
    that ref's own reference-local coordinates (matching the coordinate
    space of the BAM/PSL hits aligned against it)."""
    kmers_by_ref = {rid: _kmer_set(seq, k) for rid, seq in ref_seqs.items()}
    diagnostic = {}
    for rid, seq in ref_seqs.items():
        other_kmers = set()
        for other_id, kset in kmers_by_ref.items():
            if other_id != rid:
                other_kmers |= kset
        intervals = []
        run_start = None
        n_windows = len(seq) - k + 1
        for i in range(max(0, n_windows)):
            unique = seq[i:i + k] not in other_kmers
            if unique and run_start is None:
                run_start = i
            elif not unique and run_start is not None:
                intervals.append((run_start, i + k - 1))
                run_start = None
        if run_start is not None:
            intervals.append((run_start, len(seq)))
        diagnostic[rid] = intervals
    return diagnostic


def _overlap_bp(intervals, s, e):
    return sum(max(0, min(b, e) - max(a, s)) for a, b in intervals)


def _matched_blocks_from_cigar(ref_start, cigartuples):
    """Reference-coordinate blocks actually covered by matched/mismatched
    bases (M/=/X), splitting on deletions and introns (D/N) instead of
    treating them as covered span. A read aligned with e.g. a 300bp deletion has NOT crossed
    those 300bp, it skipped them, so a plain (reference_start,
    reference_end) span would wrongly count that gap as aligned coverage."""

    blocks = []
    pos = ref_start
    block_start = None
    for op, length in cigartuples:
        if op in (0, 7, 8):  # M, '=', X: consumes reference and is aligned
            if block_start is None:
                block_start = pos
            pos += length
        elif op in (2, 3):  # D, N: gap, not covered
            if block_start is not None:
                blocks.append((block_start, pos))
                block_start = None
            pos += length
        # I(1)/S(4)/H(5)/P(6): no reference consumption
    if block_start is not None:
        blocks.append((block_start, pos))
    return blocks


def classify_boundary_aware(bam_hits_iter, psl_hits_iter, join_pos_by_ref, diagnostic_regions_by_ref):
    """Like classify_exclusive (same proximity/anchor qualification gates),
    but a read is only handed whole to the single best-scoring reference
    when it has diagnostic-region coverage (build_diagnostic_regions)
    that actually favors that reference. A read whose aligned span falls
    entirely within sequence shared by two or more of its qualifying
    references carries no information to decide between them: rather than
    winner-take-all by total span, it is split with equal fractional
    weight across all of them (RSEM/Salmon-style multi-mapper handling).

    This targets a specific failure mode of plain total-span scoring:
    for a nested pair of isoforms (one a structural subset of the
    other), most of the shorter isoform's reads map almost
    entirely inside the shared region, where total span differs between
    candidates only by noise (alignment softclipping, indel wobble) and
    not by genuine evidence, so it should not be allowed to decide the
    read's assignment.

    Returns {rname: float_count} (fractional for ambiguous reads),
    unlike classify_exclusive's {rname: set(qname)}."""
    segments = defaultdict(lambda: defaultdict(list))  # rname -> qname -> [(s,e), ...]
    for qname, rname, ref_start, ref_end, cigartuples, qlen in bam_hits_iter:
        if rname not in join_pos_by_ref:
            continue
        if cigartuples:
            segments[rname][qname].extend(_matched_blocks_from_cigar(ref_start, cigartuples))
        else:
            segments[rname][qname].append((ref_start, ref_end))
    for qname, rname, t_start, t_end in psl_hits_iter:
        if rname not in join_pos_by_ref:
            continue
        segments[rname][qname].append((t_start, t_end))

    total_by_qname = defaultdict(dict)  # qname -> rname -> total aligned span
    diag_by_qname = defaultdict(dict)   # qname -> rname -> diagnostic-region span
    for rname, reads in segments.items():
        jp = join_pos_by_ref[rname]
        diag_regions = diagnostic_regions_by_ref.get(rname, [])
        for qname, segs in reads.items():
            if not (proximity_gate(segs, jp) and anchor_coverage_ok(segs, jp)):
                continue
            total_by_qname[qname][rname] = sum(e - s for s, e in segs)
            diag_by_qname[qname][rname] = sum(_overlap_bp(diag_regions, s, e) for s, e in segs)

    counts = defaultdict(float)
    for qname, per_ref_total in total_by_qname.items():
        per_ref_diag = diag_by_qname[qname]
        diag_positive = {r: d for r, d in per_ref_diag.items() if d > 0}
        if diag_positive:
            best = max(diag_positive, key=lambda r: (diag_positive[r], per_ref_total[r]))
            counts[best] += 1.0
        else:
            weight = 1.0 / len(per_ref_total)
            for r in per_ref_total:
                counts[r] += weight
    return dict(counts)


def classify_exclusive(bam_hits_iter, psl_hits_iter, join_pos_by_ref):
    """Best-hit-wins (RSEM/Salmon-style) assignment across BOTH minimap2 and
    BLAT hits combined: each read is assigned to the single QUALIFYING locus
    (JOIN_PROXIMITY gate + anchor-coverage check) where its total matched
    span is greatest, instead of crediting it to every qualifying locus
    independently. Used for tier2/tier3 rescue and overlap-cluster rescue,
    where boundary-variant/gene-family/overlapping loci otherwise show only
    15-58% pairwise read overlap despite being 'the same' junction under an
    independent-crediting scheme."""
    segments = defaultdict(lambda: defaultdict(list))  # rname -> qname -> [(s,e), ...]
    for qname, rname, ref_start, ref_end, cigartuples, qlen in bam_hits_iter:
        if rname not in join_pos_by_ref:
            continue
        segments[rname][qname].append((ref_start, ref_end))
    for qname, rname, t_start, t_end in psl_hits_iter:
        if rname not in join_pos_by_ref:
            continue
        segments[rname][qname].append((t_start, t_end))

    scores = defaultdict(dict)  # qname -> rname -> score
    for rname, reads in segments.items():
        jp = join_pos_by_ref[rname]
        for qname, segs in reads.items():
            if not (proximity_gate(segs, jp) and anchor_coverage_ok(segs, jp)):
                continue
            scores[qname][rname] = sum(e - s for s, e in segs)

    supporting = defaultdict(set)
    for qname, per_ref in scores.items():
        if not per_ref:
            continue
        best_rname = max(per_ref, key=per_ref.get)
        supporting[best_rname].add(qname)
    return supporting

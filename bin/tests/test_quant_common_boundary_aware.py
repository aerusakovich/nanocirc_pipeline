from quant_common import build_diagnostic_regions, classify_boundary_aware, _overlap_bp, _matched_blocks_from_cigar


def test_diagnostic_regions_find_the_differing_middle():
    # shared "AAAA" prefix and "GGGG" suffix, differing middle (CCCC vs TTTT).
    # A k-mer straddling the true edit boundary is itself legitimately
    # unique, so the diagnostic run can bleed up to k-1bp into the
    # nominally-shared flanks; only the single outermost bp (entirely
    # inside a k-mer that is wholly shared) is guaranteed excluded.
    seqs = {"A": "AAAACCCCGGGG", "B": "AAAATTTTGGGG"}
    diag = build_diagnostic_regions(seqs, k=4)
    assert _overlap_bp(diag["A"], 0, 1) == 0
    assert _overlap_bp(diag["B"], 0, 1) == 0
    assert _overlap_bp(diag["A"], 11, 12) == 0
    assert _overlap_bp(diag["B"], 11, 12) == 0
    # the differing middle must be marked diagnostic for both
    assert _overlap_bp(diag["A"], 4, 8) > 0
    assert _overlap_bp(diag["B"], 4, 8) > 0


def test_diagnostic_regions_identical_sequences_are_never_diagnostic():
    seqs = {"A": "ACGTACGTACGT", "B": "ACGTACGTACGT"}
    diag = build_diagnostic_regions(seqs, k=4)
    assert diag["A"] == []
    assert diag["B"] == []


def test_diagnostic_regions_singleton_group_is_fully_diagnostic():
    seqs = {"A": "ACGTACGTACGT"}
    diag = build_diagnostic_regions(seqs, k=4)
    assert _overlap_bp(diag["A"], 0, 12) == 12


def _psl_hit(qname, rname, start, end):
    return (qname, rname, start, end)


def test_boundary_aware_diagnostic_evidence_overrides_larger_total_span():
    # read1 hits both A and B, crossing jp=100 on both. B's total aligned
    # span (140) is larger than A's (80), but only A's span reaches into
    # A's diagnostic region [40, 60); B's span (0-140) never touches B's
    # diagnostic region [140, 160). A must win despite the smaller span.
    join_pos_by_ref = {"A": 100, "B": 100}
    diagnostic_regions = {"A": [(40, 60)], "B": [(140, 160)]}
    psl_hits = [
        _psl_hit("read1", "A", 40, 120),
        _psl_hit("read1", "B", 0, 140),
    ]
    support = classify_boundary_aware([], psl_hits, join_pos_by_ref, diagnostic_regions)
    assert support == {"A": 1.0}


def test_boundary_aware_no_diagnostic_evidence_splits_proportionally():
    # read2 hits both A and B entirely within the shared region; neither
    # diagnostic region is reached by either alignment.
    join_pos_by_ref = {"A": 100, "B": 100}
    diagnostic_regions = {"A": [(40, 60)], "B": [(140, 160)]}
    psl_hits = [
        _psl_hit("read2", "A", 85, 115),
        _psl_hit("read2", "B", 85, 115),
    ]
    support = classify_boundary_aware([], psl_hits, join_pos_by_ref, diagnostic_regions)
    assert support == {"A": 0.5, "B": 0.5}


def test_boundary_aware_combined_reads_accumulate_correctly():
    join_pos_by_ref = {"A": 100, "B": 100}
    diagnostic_regions = {"A": [(40, 60)], "B": [(140, 160)]}
    psl_hits = [
        _psl_hit("read1", "A", 40, 120),
        _psl_hit("read1", "B", 0, 140),
        _psl_hit("read2", "A", 85, 115),
        _psl_hit("read2", "B", 85, 115),
    ]
    support = classify_boundary_aware([], psl_hits, join_pos_by_ref, diagnostic_regions)
    assert support == {"A": 1.5, "B": 0.5}


def test_boundary_aware_single_reference_no_split_needed():
    join_pos_by_ref = {"A": 100}
    diagnostic_regions = {"A": []}
    psl_hits = [_psl_hit("read1", "A", 85, 115)]
    support = classify_boundary_aware([], psl_hits, join_pos_by_ref, diagnostic_regions)
    assert support == {"A": 1.0}


def test_matched_blocks_from_cigar_excludes_deletion_span():
    # 40M, 130D, 50M, 20D, 80M: two deletions must not appear as covered span
    cigar = [(0, 40), (2, 130), (0, 50), (2, 20), (0, 80)]
    blocks = _matched_blocks_from_cigar(80, cigar)
    assert blocks == [(80, 120), (250, 300), (320, 400)]


def test_boundary_aware_bam_deletion_spanning_diagnostic_region_is_not_credited():
    # real-world case this reproduces: a base-isoform read aligned against
    # its longer sibling (which contains an extra cassette exon) with a
    # deletion standing in for that exon. The read's raw (ref_start,
    # ref_end) span (80, 400) would fully cover the sibling's diagnostic
    # region (300, 320), but the read never actually matched bases there;
    # it must not be credited as diagnostic evidence for the sibling.
    join_pos_by_ref = {"A": 100, "B": 100}
    diagnostic_regions = {"A": [], "B": [(300, 320)]}
    cigar = [(0, 40), (2, 130), (0, 50), (2, 20), (0, 80)]  # matches test above
    bam_hits = [("read1", "B", 80, 400, cigar, 170)]
    psl_hits = [("read1", "A", 85, 115)]  # competing hit, crosses jp on A too

    support = classify_boundary_aware(bam_hits, psl_hits, join_pos_by_ref, diagnostic_regions)
    # neither ref has diagnostic evidence for this read -> ambiguous, split 50/50
    assert support == {"A": 0.5, "B": 0.5}


def test_boundary_aware_read_failing_gates_is_excluded():
    # segment nowhere near jp=100 fails proximity_gate entirely
    join_pos_by_ref = {"A": 100, "B": 100}
    diagnostic_regions = {"A": [], "B": []}
    psl_hits = [_psl_hit("read1", "A", 500, 520)]
    support = classify_boundary_aware([], psl_hits, join_pos_by_ref, diagnostic_regions)
    assert support == {}

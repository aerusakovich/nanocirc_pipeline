import os

import smart_merge as sm


def rec(chrom='chr1', start=0, end=100, strand='+', tool='t', score='0',
        block_sizes='100,', block_starts='0,', block_count='1', name=None):
    return {
        'chrom': chrom, 'start': start, 'end': end,
        'name': name or f'{chrom}:{start}-{end}',
        'score': score, 'strand': strand,
        'thick_start': str(start), 'thick_end': str(end), 'rgb': '0',
        'block_count': block_count, 'block_sizes': block_sizes,
        'block_starts': block_starts, 'tool': tool,
    }


# ── small pure helpers ──────────────────────────────────────────────────────

def test_bsj_key():
    assert sm.bsj_key(rec(start=10, end=20)) == (10, 20)


def test_struct_key():
    r = rec(block_count='2', block_sizes='10,20,', block_starts='0,50,')
    assert sm.struct_key(r) == ('2', '10,20,', '0,50,')


def test_make_bsj_id_without_suffix():
    assert sm.make_bsj_id('chr1', 10, 20, '+') == 'chr1:10-20:+'


def test_make_bsj_id_with_suffix():
    assert sm.make_bsj_id('chr1', 10, 20, '+', suffix='iso1') == 'chr1:10-20:+|iso1'


def test_best_record_picks_highest_score():
    recs = [rec(tool='a', score='5'), rec(tool='b', score='12'), rec(tool='c', score='3')]
    assert sm.best_record(recs)['tool'] == 'b'


def test_best_record_treats_non_numeric_score_as_zero():
    recs = [rec(tool='a', score='NA'), rec(tool='b', score='1')]
    assert sm.best_record(recs)['tool'] == 'b'


def test_max_score_of():
    d = {'a': rec(score='4'), 'b': rec(score='9'), 'c': rec(score='bad')}
    assert sm.max_score_of(d) == 9


def test_priority_tool_picks_first_in_priority_order():
    assert sm._priority_tool({'cirilong', 'circnick', 'circfl'}, sm.BSJ_PRIORITY) == 'circfl'


def test_priority_tool_falls_back_to_arbitrary_when_none_listed():
    present = {'onlytool'}
    assert sm._priority_tool(present, sm.BSJ_PRIORITY) == 'onlytool'


# ── absolute-coordinate structure helpers ───────────────────────────────────

def test_absolute_exon_coords():
    r = rec(start=100, block_sizes='10,20,', block_starts='0,50,')
    assert sm.absolute_exon_coords(r) == ((100, 10), (150, 20))


def test_abs_struct_similar_within_tolerance():
    a = ((100, 10), (150, 20))
    b = ((102, 10), (151, 20))
    assert sm.abs_struct_similar(a, b, tolerance=5) is True


def test_abs_struct_similar_rejects_size_mismatch():
    a = ((100, 10), (150, 20))
    b = ((100, 11), (150, 20))
    assert sm.abs_struct_similar(a, b, tolerance=5) is False


def test_abs_struct_similar_rejects_different_exon_count():
    a = ((100, 10), (150, 20))
    b = ((100, 10),)
    assert sm.abs_struct_similar(a, b, tolerance=5) is False


def test_abs_struct_similar_rejects_beyond_tolerance():
    a = ((100, 10),)
    b = ((110, 10),)
    assert sm.abs_struct_similar(a, b, tolerance=5) is False


def test_group_by_abs_struct_merges_similar_and_splits_different():
    tool_best = {
        'isocirc': rec(start=100, block_sizes='10,', block_starts='0,'),
        'circfl':  rec(start=102, block_sizes='10,', block_starts='0,'),  # similar to isocirc
        'cirilong': rec(start=500, block_sizes='10,', block_starts='0,'),  # different locus
    }
    groups = sm.group_by_abs_struct(tool_best, tolerance=5)
    assert len(groups) == 2
    tool_sets = sorted((sorted(tools) for _, tools in groups), key=len)
    assert tool_sets[0] == ['cirilong']
    assert tool_sets[1] == ['circfl', 'isocirc']


def test_vote_struct_groups_unique_majority():
    groups = [['A', ['isocirc', 'circfl']], ['B', ['cirilong']]]
    idx, count = sm.vote_struct_groups(groups)
    assert (idx, count) == (0, 2)


def test_vote_struct_groups_tie_uses_struct_priority():
    groups = [['A', ['circnick']], ['B', ['cirilong']]]
    idx, count = sm.vote_struct_groups(groups)
    # STRUCT_PRIORITY = isocirc > circfl > cirilong > circnick, so cirilong (group B) wins the tie
    assert (idx, count) == (1, 1)


def test_rebase_struct_preserves_absolute_exon_positions():
    r = rec(start=110, block_sizes='10,20,', block_starts='0,50,')
    original_coords = sm.absolute_exon_coords(r)
    rebased = sm.rebase_struct(r, new_start=100, new_end=170)
    assert rebased['start'] == 100
    assert rebased['end'] == 170
    assert sm.absolute_exon_coords(rebased) == original_coords


# ── voting ───────────────────────────────────────────────────────────────

def test_vote_majority_unique_majority_wins():
    votes = {'A': ['isocirc', 'circfl'], 'B': ['cirilong']}
    key, count = sm.vote_majority(votes, sm.BSJ_PRIORITY)
    assert (key, count) == ('A', 2)


def test_vote_majority_no_majority_falls_back_to_priority_order():
    # two singleton groups tied at 1/2: no strict majority, priority list breaks it
    votes = {'B': ['circnick'], 'A': ['cirilong']}
    key, count = sm.vote_majority(votes, sm.BSJ_PRIORITY)
    # BSJ_PRIORITY = isocirc > circfl > cirilong > circnick, cirilong wins over circnick
    assert key == 'A'


def test_vote_majority_priority_order_wins_over_dict_insertion_order():
    # circnick inserted first, but isocirc should still win via priority
    votes = {'circnick_key': ['circnick'], 'isocirc_key': ['isocirc']}
    key, count = sm.vote_majority(votes, sm.BSJ_PRIORITY)
    assert key == 'isocirc_key'


# ── grouping (union-find over relaxed BSJ) ──────────────────────────────────

def test_group_relaxed_merges_records_within_tolerance():
    records = [rec(tool='a', start=100, end=200), rec(tool='b', start=103, end=204)]
    groups = sm.group_relaxed(records, tolerance=5)
    assert len(groups) == 1
    tools = next(iter(groups.values()))
    assert set(tools) == {'a', 'b'}


def test_group_relaxed_keeps_records_beyond_tolerance_apart():
    records = [rec(tool='a', start=100, end=200), rec(tool='b', start=110, end=210)]
    groups = sm.group_relaxed(records, tolerance=5)
    assert len(groups) == 2


def test_group_relaxed_chains_through_an_intermediate_record():
    # A and C are 8bp apart (out of tolerance directly), but B bridges them
    records = [
        rec(tool='a', start=100, end=200),
        rec(tool='b', start=104, end=204),
        rec(tool='c', start=108, end=208),
    ]
    groups = sm.group_relaxed(records, tolerance=5)
    assert len(groups) == 1
    tools = next(iter(groups.values()))
    assert set(tools) == {'a', 'b', 'c'}


def test_group_relaxed_requires_same_chrom_and_strand():
    records = [
        rec(tool='a', chrom='chr1', strand='+', start=100, end=200),
        rec(tool='b', chrom='chr2', strand='+', start=100, end=200),
        rec(tool='c', chrom='chr1', strand='-', start=100, end=200),
    ]
    groups = sm.group_relaxed(records, tolerance=5)
    assert len(groups) == 3


def test_group_relaxed_merges_unstranded_with_real_strand():
    records = [
        rec(tool='a', chrom='chr1', strand='+', start=100, end=200),
        rec(tool='b', chrom='chr1', strand='.', start=100, end=200),
    ]
    groups = sm.group_relaxed(records, tolerance=5)
    assert len(groups) == 1
    (key, tool_map), = groups.items()
    assert key == ('chr1', 100, 200, '+')  # real strand wins over '.'
    assert set(tool_map) == {'a', 'b'}


def test_group_relaxed_still_keeps_conflicting_real_strands_apart():
    records = [
        rec(tool='a', chrom='chr1', strand='+', start=100, end=200),
        rec(tool='b', chrom='chr1', strand='-', start=100, end=200),
        rec(tool='c', chrom='chr1', strand='.', start=100, end=200),
    ]
    groups = sm.group_relaxed(records, tolerance=5)
    # 'c' can bridge to either real-strand group, but '+' and '-' must
    # still never end up in the same group as each other.
    assert len(groups) == 2
    strands = {key[3] for key in groups}
    assert strands == {'+', '-'}


# ── consensus_hybrid minority-BSJ folding ───────────────────────────────────

def test_collect_entries_consensus_hybrid_folds_minority_bsj_within_tolerance():
    tool_best = {
        'cirilong': rec(tool='cirilong', start=100, end=200, score='5'),
        'circfl':   rec(tool='circfl',   start=100, end=200, score='3'),
        'isocirc':  rec(tool='isocirc',  start=102, end=200, score='2'),  # 2bp off, within tolerance
    }
    entries = sm.collect_entries_consensus_hybrid(tool_best, struct_tolerance=5)
    assert len(entries) == 1
    main = entries[0]
    assert main['isoform_label'] == 'main'
    assert set(main['isoform_tools']) == {'cirilong', 'circfl', 'isocirc'}
    assert main['bsj_agree'] == 3


def test_collect_entries_consensus_hybrid_keeps_minority_bsj_beyond_tolerance():
    tool_best = {
        'cirilong': rec(tool='cirilong', start=100, end=200, score='5'),
        'circfl':   rec(tool='circfl',   start=100, end=200, score='3'),
        'isocirc':  rec(tool='isocirc',  start=500, end=600, score='2'),  # far away
    }
    entries = sm.collect_entries_consensus_hybrid(tool_best, struct_tolerance=5)
    assert len(entries) == 2
    labels = {e['isoform_label'] for e in entries}
    assert labels == {'main', 'iso1'}
    iso1 = next(e for e in entries if e['isoform_label'] == 'iso1')
    assert iso1['isoform_tools'] == ['isocirc']


def test_collect_entries_consensus_hybrid_folds_chained_minority_bsj():
    # 3-tool strict majority wins the BSJ vote outright at 100-200. A minority
    # at 110-200 is beyond tolerance from the winner, so it gets its own
    # entry (iso1). A second minority at 113-200 is within tolerance of that
    # iso1 position (3bp) but NOT of the winner (13bp): it should fold into
    # iso1, not spawn a separate iso2 and not fold into main.
    tool_best = {
        'a': rec(tool='a', start=100, end=200, score='5'),
        'b': rec(tool='b', start=100, end=200, score='4'),
        'c': rec(tool='c', start=100, end=200, score='3'),
        'd': rec(tool='d', start=110, end=200, score='2'),
        'e': rec(tool='e', start=113, end=200, score='1'),
    }
    entries = sm.collect_entries_consensus_hybrid(tool_best, struct_tolerance=5)
    assert len(entries) == 2
    main = next(e for e in entries if e['isoform_label'] == 'main')
    iso1 = next(e for e in entries if e['isoform_label'] == 'iso1')
    assert set(main['isoform_tools']) == {'a', 'b', 'c'}
    assert set(iso1['isoform_tools']) == {'d', 'e'}


# ── I/O ──────────────────────────────────────────────────────────────────

def test_read_bed12_parses_fields(tmp_path):
    bed = tmp_path / 'a.bed'
    bed.write_text(
        '# comment\n'
        '\n'
        'chr1\t10\t20\tname1\t5\t+\t10\t20\t0\t1\t10,\t0,\n'
    )
    records = sm.read_bed12(str(bed), 'isocirc')
    assert len(records) == 1
    r = records[0]
    assert r['chrom'] == 'chr1'
    assert r['start'] == 10
    assert r['end'] == 20
    assert r['tool'] == 'isocirc'


def test_read_bed12_missing_file_returns_empty_list():
    assert sm.read_bed12('/no/such/file.bed', 'isocirc') == []


def test_read_bed12_rejects_short_lines(tmp_path):
    bed = tmp_path / 'bad.bed'
    bed.write_text('chr1\t10\t20\tname1\n')
    import pytest
    with pytest.raises(SystemExit):
        sm.read_bed12(str(bed), 'isocirc')


def test_load_conf_lookup(tmp_path):
    tsv = tmp_path / 'conf.tsv'
    tsv.write_text(
        'bsj_id\tstruct_agree_count\tstruct_source\n'
        'chr1:10-20:+\t3\tisocirc\n'
    )
    lookup = sm.load_conf_lookup(str(tsv))
    assert lookup['chr1:10-20:+'] == {'struct_agree': 3, 'struct_source': 'isocirc'}


def test_load_conf_lookup_missing_file_returns_empty_dict():
    assert sm.load_conf_lookup('/no/such/file.tsv') == {}


# ── regression: main() must write all 8 output files when no tool has calls ──
# (earlier bug: an incomplete suffix list here silently skipped the
# _smart_consensus_hybrid.* files when fewer than 2 tools were active)

def test_main_writes_all_empty_outputs_when_no_tool_has_calls(tmp_path, monkeypatch):
    outdir = tmp_path / 'out'
    monkeypatch.setattr('sys.argv', [
        'smart_merge.py',
        '--sample', 'testsample',
        '--tool_names', 'cirilong', 'isocirc',
        '--bed_files', str(tmp_path / 'missing_a.bed'), str(tmp_path / 'missing_b.bed'),
        '--outdir', str(outdir),
    ])
    sm.main()

    expected = [
        '_smart_consensus.bed12', '_smart_consensus_confidence.tsv',
        '_smart_consensus_xstruct.bed12', '_smart_consensus_xstruct_confidence.tsv',
        '_smart_consensus_hybrid.bed12', '_smart_consensus_hybrid_confidence.tsv',
        '_smart_priority.bed12', '_smart_priority_confidence.tsv',
    ]
    for suffix in expected:
        path = outdir / f'testsample{suffix}'
        assert path.exists(), f'missing {path}'
        assert path.stat().st_size == 0


# ── splice-motif strand check ────────────────────────────────────────────────

def make_fasta(tmp_path, chroms):
    """chroms: {name: sequence}. Writes name.fa (one line per record) + its
    .fai, returns the .fa path. Mutating positions in `sequence` via slicing
    (a mutable bytearray/list would help, but plain string concatenation of
    short fixture sequences is simpler to read in each test)."""
    fasta_path = tmp_path / 'genome.fa'
    fai_lines = []
    offset = 0
    with open(fasta_path, 'w') as fh:
        for name, seq in chroms.items():
            header = f'>{name}\n'
            fh.write(header)
            offset += len(header)
            fai_lines.append(f'{name}\t{len(seq)}\t{offset}\t{len(seq)}\t{len(seq) + 1}')
            fh.write(seq + '\n')
            offset += len(seq) + 1
    with open(str(fasta_path) + '.fai', 'w') as fh:
        fh.write('\n'.join(fai_lines) + '\n')
    return str(fasta_path)


def seq_with(length, positions):
    """A length-N sequence of 'N' with positions {offset: 2-char string}
    overlaid at those 0-based offsets."""
    chars = ['N'] * length
    for pos, s in positions.items():
        chars[pos] = s[0]
        chars[pos + 1] = s[1]
    return ''.join(chars)


def test_resolve_strand_by_motif_unresolved_by_tool(tmp_path):
    fasta = make_fasta(tmp_path, {'chr1': seq_with(70, {8: 'AG', 60: 'GT'})})
    fa = sm.FastaRandomAccess(fasta)
    strand, status = sm.resolve_strand_by_motif(fa, 'chr1', 10, 60, '.', '50,', '0,')
    assert strand == '.'
    assert status == 'unresolved_by_tool'


def test_resolve_strand_by_motif_single_exon_confirmed(tmp_path):
    # outer BSJ only (single exon): acceptor just before start, donor just after end
    fasta = make_fasta(tmp_path, {'chr1': seq_with(70, {8: 'AG', 60: 'GT'})})
    fa = sm.FastaRandomAccess(fasta)
    strand, status = sm.resolve_strand_by_motif(fa, 'chr1', 10, 60, '+', '50,', '0,')
    assert (strand, status) == ('+', 'confirmed')


def test_resolve_strand_by_motif_single_exon_never_hard_flips(tmp_path):
    # motif is a clean, unambiguous MINUS signal, but the call is single-exon
    # ('+' called): must downgrade to ambiguous, never flip on 1 junction alone.
    fasta = make_fasta(tmp_path, {'chr1': seq_with(70, {8: 'AC', 60: 'CT'})})
    fa = sm.FastaRandomAccess(fasta)
    strand, status = sm.resolve_strand_by_motif(fa, 'chr1', 10, 60, '+', '50,', '0,')
    assert (strand, status) == ('.', 'ambiguous')


def test_resolve_strand_by_motif_single_exon_ambiguous_no_signal(tmp_path):
    fasta = make_fasta(tmp_path, {'chr1': seq_with(70, {})})  # all 'N', no motif either way
    fa = sm.FastaRandomAccess(fasta)
    strand, status = sm.resolve_strand_by_motif(fa, 'chr1', 10, 60, '+', '50,', '0,')
    assert (strand, status) == ('.', 'ambiguous')


def test_resolve_strand_by_motif_multi_exon_confirmed(tmp_path):
    # start=10, end=60, exons [10,20) and [50,60), internal intron (20,50).
    # outer: acceptor@8, donor@60. internal: donor@20, acceptor@48.
    fasta = make_fasta(tmp_path, {'chr1': seq_with(70, {
        8: 'AG', 60: 'GT',    # outer BSJ, '+'
        20: 'GT', 48: 'AG',   # internal junction, '+'
    })})
    fa = sm.FastaRandomAccess(fasta)
    strand, status = sm.resolve_strand_by_motif(
        fa, 'chr1', 10, 60, '+', '10,10,', '0,40,')
    assert (strand, status) == ('+', 'confirmed')


def test_resolve_strand_by_motif_multi_exon_flips_on_unanimous_disagreement(tmp_path):
    fasta = make_fasta(tmp_path, {'chr1': seq_with(70, {
        8: 'AC', 60: 'CT',    # outer BSJ, '-'
        20: 'CT', 48: 'AC',   # internal junction, '-'
    })})
    fa = sm.FastaRandomAccess(fasta)
    strand, status = sm.resolve_strand_by_motif(
        fa, 'chr1', 10, 60, '+', '10,10,', '0,40,')
    assert (strand, status) == ('-', 'flipped')


def test_write_outputs_never_emits_duplicate_bsj_id_from_cross_group_flip(tmp_path):
    # write_outputs()'s own defensive fallback for groups that arrive
    # already split by conflicting strand (main() now pre-resolves each
    # record's strand by motif before grouping, so this specific circfl/
    # circnick shape no longer reaches write_outputs split in two via the
    # real main() path -- see test_main_merges_opposite_strand_tools_calling_
    # identical_structure below for that combined outcome). This test calls
    # write_outputs directly with pre-built groups, covering any other
    # caller (or a genuine antisense pair, where the two records' OWN
    # motif evidence disagrees and group_relaxed is right to keep them
    # apart) that can still reach write_outputs in this shape: if both
    # groups converge on the same final strand, they'd collide on the exact
    # same bsj_id -- which broke a downstream script's
    # catalog.set_index('bsj_id') uniqueness assumption on real human data.
    # circfl: single-exon '+' call, motif confirms '+'.
    # circnick: multi-exon '-' call at the SAME coordinates (kept as its own
    # group since '+' vs '-' never merge), whose own junctions unanimously
    # support '+' instead -- eligible to flip since it has 2+ junctions.
    fasta = make_fasta(tmp_path, {'chr1': seq_with(70, {
        8: 'AG', 60: 'GT',    # outer BSJ (shared coords), '+'
        20: 'GT', 48: 'AG',   # circnick's internal junction, also '+'
    })})

    circfl_rec = rec(tool='circfl', chrom='chr1', strand='+', start=10, end=60,
                      block_count='1', block_sizes='50,', block_starts='0,')
    circnick_rec = rec(tool='circnick', chrom='chr1', strand='-', start=10, end=60,
                        block_count='2', block_sizes='10,10,', block_starts='0,40,')

    groups = sm.group_relaxed([circfl_rec, circnick_rec], tolerance=5)
    assert len(groups) == 2  # confirms they're still kept apart by group_relaxed itself

    outdir = tmp_path / 'out'
    sm.write_outputs(groups, ['circfl', 'circnick'], 'testsample', str(outdir),
                      struct_tolerance=5, fasta=fasta)

    bed_lines = (outdir / 'testsample_smart_consensus_hybrid.bed12').read_text().strip().split('\n')
    names = [line.split('\t')[3] for line in bed_lines]
    assert len(names) == len(set(names)), f'duplicate bsj_id in output: {names}'

    strands = sorted(line.split('\t')[5] for line in bed_lines)
    assert strands == ['+', '.'], (
        "circfl keeps its confirmed '+'; circnick's flip to the same '+' "
        "must be refused and downgraded to '.' instead of colliding"
    )


def test_write_outputs_disambiguates_two_groups_both_landing_on_ambiguous(tmp_path):
    # Real crash on mouse species benchmark data: two separate groups
    # (again kept apart by group_relaxed since they conflict on strand)
    # each independently resolve to '.' (no clean single-exon signal either
    # way, not a flip). Downgrading a flip to '.' isn't enough here since
    # NEITHER side ever had a confident strand to begin with -- both want
    # '.' from the start. Falls back to a disambiguating suffix.
    fasta = make_fasta(tmp_path, {'chr1': seq_with(70, {})})  # all 'N': no motif signal at all

    circfl_rec = rec(tool='circfl', chrom='chr1', strand='+', start=10, end=60,
                      block_count='1', block_sizes='50,', block_starts='0,')
    circnick_rec = rec(tool='circnick', chrom='chr1', strand='-', start=10, end=60,
                        block_count='1', block_sizes='50,', block_starts='0,')

    groups = sm.group_relaxed([circfl_rec, circnick_rec], tolerance=5)
    assert len(groups) == 2

    outdir = tmp_path / 'out'
    sm.write_outputs(groups, ['circfl', 'circnick'], 'testsample', str(outdir),
                      struct_tolerance=5, fasta=fasta)

    bed_lines = (outdir / 'testsample_smart_consensus_hybrid.bed12').read_text().strip().split('\n')
    names = [line.split('\t')[3] for line in bed_lines]
    assert len(names) == len(set(names)), f'duplicate bsj_id in output: {names}'
    assert sorted(line.split('\t')[5] for line in bed_lines) == ['.', '.']


def test_resolve_strand_by_motif_multi_exon_mixed_signal_stays_ambiguous(tmp_path):
    # outer junction supports '+', internal junction supports '-': mixed,
    # never a hard flip on a non-unanimous signal.
    fasta = make_fasta(tmp_path, {'chr1': seq_with(70, {
        8: 'AG', 60: 'GT',    # outer BSJ, '+'
        20: 'CT', 48: 'AC',   # internal junction, '-'
    })})
    fa = sm.FastaRandomAccess(fasta)
    strand, status = sm.resolve_strand_by_motif(
        fa, 'chr1', 10, 60, '+', '10,10,', '0,40,')
    assert (strand, status) == ('.', 'ambiguous')


def write_bed12(path, records):
    with open(path, 'w') as fh:
        for r in records:
            fh.write('\t'.join([
                r['chrom'], str(r['start']), str(r['end']), r.get('name', 'x'),
                r.get('score', '0'), r['strand'], str(r['start']), str(r['start']),
                '0', r['block_count'], r['block_sizes'], r['block_starts'],
            ]) + '\n')


def test_main_merges_opposite_strand_tools_calling_identical_structure(tmp_path, monkeypatch):
    # End-to-end reproduction of a real gap found on ciri_long_run1
    # (chr10:46580893-46581036: cirilong called '+', circnick called the
    # identical structure '.', each ending up with bsj_confidence=1 instead
    # of a combined 2). Two tools call the exact same single-exon structure
    # on opposite strands; the motif unambiguously supports '+'. Pre-
    # resolving each raw record's strand before grouping (main()) should
    # let group_relaxed see them as strand-compatible and combine them into
    # one entry, instead of splitting into two weaker ones.
    fasta = make_fasta(tmp_path, {'chr1': seq_with(70, {8: 'AG', 60: 'GT'})})

    write_bed12(tmp_path / 'toolA.bed12', [dict(
        chrom='chr1', start=10, end=60, strand='+',
        block_count='1', block_sizes='50,', block_starts='0,', score='10')])
    write_bed12(tmp_path / 'toolB.bed12', [dict(
        chrom='chr1', start=10, end=60, strand='-',
        block_count='1', block_sizes='50,', block_starts='0,', score='8')])

    outdir = tmp_path / 'out'
    monkeypatch.setattr('sys.argv', [
        'smart_merge.py',
        '--sample', 'testsample',
        '--tool_names', 'toolA', 'toolB',
        '--bed_files', str(tmp_path / 'toolA.bed12'), str(tmp_path / 'toolB.bed12'),
        '--tolerance', '5', '--struct_tolerance', '5', '--n_active', '2',
        '--fasta', fasta,
        '--outdir', str(outdir),
    ])
    sm.main()

    bed_lines = (outdir / 'testsample_smart_consensus_hybrid.bed12').read_text().strip().split('\n')
    assert len(bed_lines) == 1, f'expected one merged entry, got: {bed_lines}'
    assert bed_lines[0].split('\t')[5] == '+'

    tsv_path = outdir / 'testsample_smart_consensus_hybrid_confidence.tsv'
    header, row = tsv_path.read_text().strip().split('\n')
    cols = dict(zip(header.lstrip('#').split('\t'), row.split('\t')))
    assert cols['bsj_confidence'] == '2'
    assert cols['toolA'] == '1' and cols['toolB'] == '1'
    assert cols['strand_status'] == 'confirmed'


def test_write_outputs_preserves_isoform_label_over_generic_dup_suffix(tmp_path):
    # Real case (ciri_long_run1 chr1:92102472-92108001): a BSJ with 3+
    # distinct isoform structures from multi-isoform recovery. Without a
    # motif check to distinguish them, every isoform's resolved_strand stays
    # identical to the group's raw strand, so the 2nd+ isoform always
    # collides with an earlier one at the exact same (coords, strand). The
    # collision must not erase iso2's real isoform_label by falling back to
    # a generic 'dup2' -- group-level bsj_confidence/tool flags are already
    # shared correctly across every isoform row; only the bsj_id suffix was
    # ever wrong.
    isocirc_recs = [
        rec(tool='isocirc', chrom='chr1', strand='+', start=10, end=60, score='30',
            block_count='2', block_sizes='20,20,', block_starts='0,30,'),
        rec(tool='isocirc', chrom='chr1', strand='+', start=10, end=60, score='20',
            block_count='2', block_sizes='15,15,', block_starts='0,35,'),
        rec(tool='isocirc', chrom='chr1', strand='+', start=10, end=60, score='10',
            block_count='2', block_sizes='10,10,', block_starts='0,40,'),
    ]

    groups = sm.group_relaxed(isocirc_recs, tolerance=5)
    assert len(groups) == 1

    outdir = tmp_path / 'out'
    sm.write_outputs(groups, ['isocirc'], 'testsample', str(outdir), struct_tolerance=5)

    tsv_lines = (outdir / 'testsample_smart_consensus_hybrid_confidence.tsv').read_text().strip().split('\n')
    header = tsv_lines[0].lstrip('#').split('\t')
    rows = [dict(zip(header, line.split('\t'))) for line in tsv_lines[1:]]

    labels = sorted(r['isoform_label'] for r in rows)
    assert labels == ['iso1', 'iso2', 'main']

    suffixes = sorted(r['bsj_id'].split('|')[-1] if '|' in r['bsj_id'] else 'main' for r in rows)
    assert suffixes == ['iso1', 'iso2', 'main'], f'dup{{N}} fallback used instead of real isoform label: {suffixes}'

    # Group-level support (1 supporting tool, isocirc) is identical across
    # all 3 isoform rows either way -- the collision never fragmented it.
    assert {r['bsj_confidence'] for r in rows} == {'1'}
    assert {r['isocirc'] for r in rows} == {'1'}

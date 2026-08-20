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

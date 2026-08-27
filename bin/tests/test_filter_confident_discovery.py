import filter_confident_discovery as fcd


def test_weak_read_support_below_threshold():
    assert fcd.weak_read_support('1', min_reads=2) is True


def test_weak_read_support_above_threshold():
    assert fcd.weak_read_support('5', min_reads=2) is False


def test_weak_read_support_at_threshold_counts_as_weak():
    assert fcd.weak_read_support('2', min_reads=2) is True


def test_weak_read_support_missing_or_invalid_value_is_weak():
    assert fcd.weak_read_support('', min_reads=2) is True
    assert fcd.weak_read_support(None, min_reads=2) is True
    assert fcd.weak_read_support('NA', min_reads=2) is True


def test_drop_circnick_only_low_read_in_discovery():
    row = {'supporting_tools': 'circnick', 'nanocirc_quant_reads': '0'}
    assert fcd.should_drop(row, min_reads=2, category='discovery') is True


def test_drop_circnick_only_low_read_in_balanced():
    row = {'supporting_tools': 'circnick', 'nanocirc_quant_reads': '1'}
    assert fcd.should_drop(row, min_reads=2, category='balanced') is True


def test_keep_circnick_only_with_enough_reads():
    row = {'supporting_tools': 'circnick', 'nanocirc_quant_reads': '5'}
    assert fcd.should_drop(row, min_reads=2, category='discovery') is False


def test_keep_circnick_corroborated_by_another_tool_even_with_zero_reads():
    row = {'supporting_tools': 'circnick,isocirc', 'nanocirc_quant_reads': '0'}
    assert fcd.should_drop(row, min_reads=2, category='discovery') is False


def test_circnick_only_rule_does_not_apply_outside_discovery_and_balanced():
    row = {'supporting_tools': 'circnick', 'nanocirc_quant_reads': '0'}
    assert fcd.should_drop(row, min_reads=2, category='high_confidence') is False
    assert fcd.should_drop(row, min_reads=2, category='') is False


def test_drop_isocirc_only_low_read_in_high_confidence():
    row = {'supporting_tools': 'isocirc', 'nanocirc_quant_reads': '0'}
    assert fcd.should_drop(row, min_reads=2, category='high_confidence') is True


def test_isocirc_only_rule_does_not_apply_outside_high_confidence():
    row = {'supporting_tools': 'isocirc', 'nanocirc_quant_reads': '0'}
    assert fcd.should_drop(row, min_reads=2, category='discovery') is False
    assert fcd.should_drop(row, min_reads=2, category='balanced') is False


def test_never_drop_when_read_support_is_strong():
    row = {'supporting_tools': 'circnick', 'nanocirc_quant_reads': '100'}
    assert fcd.should_drop(row, min_reads=2, category='discovery') is False


def test_should_drop_missing_supporting_tools_column_is_safe():
    row = {'nanocirc_quant_reads': '0'}
    assert fcd.should_drop(row, min_reads=2, category='discovery') is False


def test_main_end_to_end(tmp_path, monkeypatch):
    tsv_in = tmp_path / 'in.tsv'
    tsv_in.write_text(
        'bsj_id\tsupporting_tools\tnanocirc_quant_reads\n'
        'chr1:1-100:+\tcircnick\t0\n'
        'chr1:200-300:+\tcircnick,isocirc\t0\n'
        'chr1:400-500:+\tcirilong\t10\n'
    )
    bed_in = tmp_path / 'in.bed12'
    bed_in.write_text(
        'chr1\t0\t100\tchr1:1-100:+\t0\t+\t0\t100\t0\t1\t100,\t0,\n'
        'chr1\t199\t300\tchr1:200-300:+\t0\t+\t199\t300\t0\t1\t101,\t0,\n'
        'chr1\t399\t500\tchr1:400-500:+\t0\t+\t399\t500\t0\t1\t101,\t0,\n'
    )
    out_tsv = tmp_path / 'out.tsv'
    out_bed = tmp_path / 'out.bed12'

    monkeypatch.setattr('sys.argv', [
        'filter_confident_discovery.py',
        '--tsv', str(tsv_in), '--bed', str(bed_in),
        '--min_reads', '2', '--category', 'discovery',
        '--out_tsv', str(out_tsv), '--out_bed', str(out_bed),
    ])
    fcd.main()

    tsv_text = out_tsv.read_text()
    kept_lines = tsv_text.strip().splitlines()
    assert len(kept_lines) == 3  # header + 2 surviving rows
    assert 'chr1:1-100:+' not in tsv_text
    assert 'chr1:200-300:+' in tsv_text
    assert 'chr1:400-500:+' in tsv_text

    bed_lines = out_bed.read_text().strip().splitlines()
    assert len(bed_lines) == 2

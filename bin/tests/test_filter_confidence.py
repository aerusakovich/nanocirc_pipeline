import filter_confidence as fc

TRUSTED = {'cirilong', 'isocirc', 'circfl'}


def make_row(bsj_consensus, isoform_consensus, bsj_source=None, isoform_tools=None):
    """Build (cols, idx) the way passes_filter expects: idx maps column
    name -> position in cols, matching a real smart-merge confidence TSV."""
    cols = []
    idx = {}

    def add(name, value):
        idx[name] = len(cols)
        cols.append(value)

    add('bsj_consensus', bsj_consensus)
    add('isoform_consensus', isoform_consensus)
    if bsj_source is not None:
        add('bsj_source', bsj_source)
    if isoform_tools is not None:
        add('isoform_tools', isoform_tools)
    return cols, idx


# ── no_low / high_only: no exceptions ───────────────────────────────────────

def test_no_low_drops_low_bsj_regardless_of_source():
    cols, idx = make_row('Low', 'High', bsj_source='isocirc')
    assert fc.passes_filter(cols, idx, 'no_low', TRUSTED) is False


def test_no_low_keeps_high_high():
    cols, idx = make_row('High', 'High')
    assert fc.passes_filter(cols, idx, 'no_low', TRUSTED) is True


def test_no_low_keeps_medium_medium():
    cols, idx = make_row('Medium', 'Medium')
    assert fc.passes_filter(cols, idx, 'no_low', TRUSTED) is True


def test_high_only_drops_medium():
    cols, idx = make_row('Medium', 'High')
    assert fc.passes_filter(cols, idx, 'high_only', TRUSTED) is False


def test_high_only_keeps_high_high():
    cols, idx = make_row('High', 'High')
    assert fc.passes_filter(cols, idx, 'high_only', TRUSTED) is True


# ── trusted_only ─────────────────────────────────────────────────────────

def test_trusted_only_rescues_low_bsj_from_trusted_source():
    cols, idx = make_row('Low', 'High', bsj_source='cirilong')
    assert fc.passes_filter(cols, idx, 'trusted_only', TRUSTED) is True


def test_trusted_only_drops_low_bsj_from_untrusted_source():
    cols, idx = make_row('Low', 'High', bsj_source='circnick')
    assert fc.passes_filter(cols, idx, 'trusted_only', TRUSTED) is False


def test_trusted_only_rescues_low_isoform_if_any_tool_trusted():
    cols, idx = make_row('High', 'Low', isoform_tools='circnick,circfl')
    assert fc.passes_filter(cols, idx, 'trusted_only', TRUSTED) is True


def test_trusted_only_drops_low_isoform_if_no_tool_trusted():
    cols, idx = make_row('High', 'Low', isoform_tools='circnick')
    assert fc.passes_filter(cols, idx, 'trusted_only', TRUSTED) is False


# ── isocirc_only: narrow isocirc-only exception (this pipeline's high_confidence tier) ─

def test_isocirc_only_rescues_low_bsj_from_isocirc():
    cols, idx = make_row('Low', 'High', bsj_source='isocirc')
    assert fc.passes_filter(cols, idx, 'isocirc_only', TRUSTED) is True


def test_isocirc_only_does_not_rescue_low_bsj_from_other_trusted_tool():
    # cirilong is trusted_only-trusted, but isocirc_only only ever trusts isocirc
    cols, idx = make_row('Low', 'High', bsj_source='cirilong')
    assert fc.passes_filter(cols, idx, 'isocirc_only', TRUSTED) is False


def test_isocirc_only_rescues_low_isoform_when_isocirc_present():
    cols, idx = make_row('High', 'Low', isoform_tools='isocirc,circnick')
    assert fc.passes_filter(cols, idx, 'isocirc_only', TRUSTED) is True


# ── fallback for non-smart TSVs (tool presence flags instead of source cols) ─

def test_trusted_only_fallback_uses_tool_presence_flags():
    cols = ['Low', 'High', '1', '0']
    idx = {'bsj_consensus': 0, 'isoform_consensus': 1, 'cirilong': 2, 'circnick': 3}
    assert fc.passes_filter(cols, idx, 'trusted_only', TRUSTED) is True

    cols2 = ['Low', 'High', '0', '1']
    assert fc.passes_filter(cols2, idx, 'trusted_only', TRUSTED) is False

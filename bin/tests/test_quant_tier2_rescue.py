import pandas as pd

from quant_tier2_rescue import _expand_to_siblings


def test_flagged_isoform_pulls_in_unflagged_sibling():
    deduped = pd.DataFrame({
        "bsj_id": ["chr19:35329175-35336184:+", "chr19:35329175-35336184:+|iso1",
                   "chr1:100-200:+"],
    })
    out = _expand_to_siblings(["chr19:35329175-35336184:+|iso1"], deduped)
    assert set(out) == {"chr19:35329175-35336184:+", "chr19:35329175-35336184:+|iso1"}


def test_no_siblings_returns_only_flagged():
    deduped = pd.DataFrame({"bsj_id": ["chr1:100-200:+", "chr2:500-600:-"]})
    out = _expand_to_siblings(["chr1:100-200:+"], deduped)
    assert out == ["chr1:100-200:+"]


def test_multiple_flagged_loci_each_pull_their_own_siblings():
    deduped = pd.DataFrame({
        "bsj_id": ["chr1:100-200:+", "chr1:100-200:+|iso1", "chr1:100-200:+|iso2",
                   "chr2:500-600:-", "chr2:500-600:-|iso1"],
    })
    out = _expand_to_siblings(["chr1:100-200:+|iso1", "chr2:500-600:-"], deduped)
    assert set(out) == set(deduped["bsj_id"])


def test_empty_flagged_ids_returns_empty():
    deduped = pd.DataFrame({"bsj_id": ["chr1:100-200:+", "chr1:100-200:+|iso1"]})
    assert _expand_to_siblings([], deduped) == []

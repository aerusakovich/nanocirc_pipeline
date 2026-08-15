import pandas as pd

from quant_finalize import finalize_counts


def _deduped():
    return pd.DataFrame({
        "bsj_id": ["chr1:100-200:+", "chr1:100-200:+|iso1", "chr2:500-600:-"],
        "chrom": ["chr1", "chr1", "chr2"],
        "start": [100, 100, 500],
        "end": [200, 200, 600],
        "strand": ["+", "+", "-"],
    })


def test_finalize_accepts_fractional_tier2_override():
    # tier1 counts are plain ints (pandas infers an int64 remap_count column);
    # classify_boundary_aware's proportional split can hand tier2/tier3 a
    # fractional count, which must not raise a pandas dtype error on assignment.
    deduped = _deduped()
    tier1 = pd.DataFrame({"bsj_id": deduped["bsj_id"], "remap_count": [10, 5, 20]})
    flagged = pd.DataFrame({"bsj_id": [], "is_gene_family": []})
    tier2 = pd.DataFrame({"bsj_id": ["chr1:100-200:+|iso1"], "tier2_count": [238.5]})
    tier3 = pd.DataFrame(columns=["bsj_id", "tier3_count"])

    final = finalize_counts(deduped, tier1, flagged, tier2, tier3)

    row = final.loc[final["bsj_id"] == "chr1:100-200:+|iso1"].iloc[0]
    assert row["remap_count"] == 238.5
    assert row["quant_tier"] == "tier2"


def test_finalize_accepts_fractional_tier3_override():
    deduped = _deduped()
    tier1 = pd.DataFrame({"bsj_id": deduped["bsj_id"], "remap_count": [10, 5, 20]})
    flagged = pd.DataFrame({"bsj_id": ["chr2:500-600:-"], "is_gene_family": [True]})
    tier2 = pd.DataFrame(columns=["bsj_id", "tier2_count"])
    tier3 = pd.DataFrame({"bsj_id": ["chr2:500-600:-"], "tier3_count": [12.25]})

    final = finalize_counts(deduped, tier1, flagged, tier2, tier3)

    row = final.loc[final["bsj_id"] == "chr2:500-600:-"].iloc[0]
    assert row["remap_count"] == 12.25
    assert row["quant_tier"] == "tier3"
    assert row["low_confidence"] == False

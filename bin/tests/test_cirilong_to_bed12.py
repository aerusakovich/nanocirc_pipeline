from cirilong_to_bed12 import convert, load_isoform_fractions


def _info_line(circ_id, start, end, score, isoform_field, strand="+"):
    return (f'chr1\tCIRI-long\tcircRNA\t{start}\t{end}\t{score}\t{strand}\t.\t'
            f'circ_id "{circ_id}"; isoform "{isoform_field}"; gene_id "G1";\n')


def _read_bed12(path):
    rows = []
    with open(path) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            rows.append({"name": cols[3], "score": int(cols[4])})
    return rows


def test_load_isoform_fractions_parses_key_and_value(tmp_path):
    iso_path = tmp_path / "sample.isoforms"
    iso_path.write_text(
        "isoform_ID\tsample\n"
        "chr1:100-200|100-150,180-200\t0.6666666666666666\n"
        "chr1:100-200|100-200\t0.3333333333333333\n"
    )
    fractions = load_isoform_fractions(str(iso_path))
    assert fractions[("chr1:100-200", "100-150,180-200")] == 0.6666666666666666
    assert fractions[("chr1:100-200", "100-200")] == 0.3333333333333333


def test_load_isoform_fractions_missing_file_returns_empty(tmp_path):
    assert load_isoform_fractions(str(tmp_path / "does_not_exist.isoforms")) == {}


def test_load_isoform_fractions_no_file_sentinel_returns_empty():
    assert load_isoform_fractions("NO_FILE") == {}
    assert load_isoform_fractions("") == {}


def test_convert_splits_ambiguous_isoform_score_by_fraction(tmp_path):
    # Real observed case: chr1:1040881-1041172, .info score=3 shared across
    # both interpretations, .isoforms resolves it to 2 reads / 1 read.
    info = tmp_path / "sample.info"
    info.write_text(_info_line(
        "chr1:1040881-1041172", 1040881, 1041172, 3,
        "1040881-1041017,1041093-1041172|1040881-1041172",
    ))
    isoforms = tmp_path / "sample.isoforms"
    isoforms.write_text(
        "isoform_ID\tsample\n"
        "chr1:1040881-1041172|1040881-1041017,1041093-1041172\t0.6666666666666666\n"
        "chr1:1040881-1041172|1040881-1041172\t0.3333333333333333\n"
    )
    out = tmp_path / "out.bed12"
    convert(str(info), str(isoforms), str(out))

    rows = {r["name"]: r["score"] for r in _read_bed12(out)}
    assert rows["chr1:1040881-1041172_iso1"] == 2
    assert rows["chr1:1040881-1041172_iso2"] == 1


def test_convert_single_isoform_locus_keeps_shared_score(tmp_path):
    info = tmp_path / "sample.info"
    info.write_text(_info_line("chr1:100-500", 100, 500, 7, "100-500"))
    isoforms = tmp_path / "sample.isoforms"
    isoforms.write_text("isoform_ID\tsample\nchr1:100-500|100-500\t1.0\n")
    out = tmp_path / "out.bed12"
    convert(str(info), str(isoforms), str(out))

    rows = _read_bed12(out)
    assert len(rows) == 1
    assert rows[0]["name"] == "chr1:100-500"
    assert rows[0]["score"] == 7


def test_convert_falls_back_to_shared_score_when_isoforms_missing(tmp_path):
    info = tmp_path / "sample.info"
    info.write_text(_info_line(
        "chr1:1040881-1041172", 1040881, 1041172, 3,
        "1040881-1041017,1041093-1041172|1040881-1041172",
    ))
    out = tmp_path / "out.bed12"
    convert(str(info), "NO_FILE", str(out))

    rows = {r["name"]: r["score"] for r in _read_bed12(out)}
    assert rows["chr1:1040881-1041172_iso1"] == 3
    assert rows["chr1:1040881-1041172_iso2"] == 3


def test_convert_falls_back_when_fraction_missing_for_this_locus(tmp_path):
    # .isoforms file exists but doesn't cover this specific locus (e.g.
    # stale/mismatched inputs); must not crash, must fall back per-row.
    info = tmp_path / "sample.info"
    info.write_text(_info_line(
        "chr1:1040881-1041172", 1040881, 1041172, 3,
        "1040881-1041017,1041093-1041172|1040881-1041172",
    ))
    isoforms = tmp_path / "sample.isoforms"
    isoforms.write_text("isoform_ID\tsample\nchr2:9999-9999|9999-9999\t1.0\n")
    out = tmp_path / "out.bed12"
    convert(str(info), str(isoforms), str(out))

    rows = {r["name"]: r["score"] for r in _read_bed12(out)}
    assert rows["chr1:1040881-1041172_iso1"] == 3
    assert rows["chr1:1040881-1041172_iso2"] == 3

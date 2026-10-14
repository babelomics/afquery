"""Genotype tallies at multi-allelic sites.

Each ALT allele at a position is its own row. A sample carrying only another
allele there is neither a carrier nor hom-ref for this one, so it is left out of
every tally and N_HET + N_HOM_ALT + N_HOM_REF + N_FAIL + N_NO_COVERAGE falls
short of n_eligible by exactly that number of samples.
"""
import csv
import io
import warnings

import cyvcf2
import pytest

from afquery.models import AfqueryWarning

from test_variant_info import _build_db

SITE = ("chr1", 5000)

COHORT = [
    ("A_HET", "female", "WGS", [("chr1", 5000, "G", "A", "0/1")]),
    ("A_HOM", "male",   "WGS", [("chr1", 5000, "G", "A", "1/1")]),
    ("T_HET", "female", "WGS", [("chr1", 5000, "G", "T", "0/1")]),
    ("AT",    "male",   "WGS", [("chr1", 5000, "G", "A,T", "1/2")]),
    ("REF1",  "female", "WGS", [("chr1", 100, "C", "G", "0/1")]),
    ("REF2",  "male",   "WGS", [("chr1", 100, "C", "G", "0/1")]),
]

# (N_HET, N_HOM_ALT, N_HOM_REF, samples carrying only the other allele)
EXPECTED = {
    "A": (2, 1, 2, 1),   # het A_HET, AT; hom A_HOM; hom-ref REF1, REF2; other T_HET
    "T": (2, 0, 2, 2),   # het T_HET, AT; hom-ref REF1, REF2; other A_HET, A_HOM
}


@pytest.fixture(scope="module")
def db(tmp_path_factory):
    return _build_db(tmp_path_factory.mktemp("multiallelic"), COHORT)


def _tallies(r):
    return (r.N_HET, r.N_HOM_ALT, r.N_HOM_REF)


def _check(results):
    by_alt = {r.variant.alt: r for r in results if (r.variant.chrom, r.variant.pos) == SITE}
    assert set(by_alt) == set(EXPECTED)
    for alt, (het, hom, hom_ref, other) in EXPECTED.items():
        r = by_alt[alt]
        assert _tallies(r) == (het, hom, hom_ref), alt
        assert r.n_samples_eligible == 6
        total = r.N_HET + r.N_HOM_ALT + r.N_HOM_REF + r.N_FAIL + r.N_NO_COVERAGE
        assert total == r.n_samples_eligible - other, alt


def test_point_query(db):
    _check(db.query(*SITE))


def test_region_query(db):
    _check(db.query_region("chr1", 1, 10000))


def test_batch_query_one_allele_requested(db):
    """The allele not asked for still keeps its carriers out of hom-ref."""
    [r] = db.query_batch("chr1", [(5000, "G", "A")])
    assert _tallies(r) == EXPECTED["A"][:3]


def test_batch_multi_query(db):
    results = db.query_batch_multi([("chr1", 5000, "G", "T"), ("chr1", 5000, "G", "A")])
    _check(results)


def test_dump_matches_query(db):
    buf = io.StringIO()
    db.dump(output=buf, by_sex=True)
    rows =[r for r in csv.DictReader(io.StringIO(buf.getvalue())) if r["pos"] == "5000"]
    assert {r["alt"] for r in rows} == set(EXPECTED)
    for row in rows:
        [base] = db.query(*SITE, alt=row["alt"])
        assert int(row["N_HOM_REF"]) == base.N_HOM_REF == EXPECTED[row["alt"]][2]
        for sex in ("male", "female"):
            [g] = db.query(*SITE, alt=row["alt"], sex=sex)
            assert int(row[f"N_HOM_REF_{sex}"]) == g.N_HOM_REF
            assert int(row[f"N_HET_{sex}"]) == g.N_HET


def test_annotate_matches_query(db, tmp_path):
    """C is not stored at the site: its hom-ref excludes carriers of A and T."""
    vcf_in = tmp_path / "in.vcf"
    vcf_in.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t5000\t.\tG\tA,T,C\t.\tPASS\t.\n"
    )
    out = tmp_path / "out.vcf"
    db.annotate_vcf(str(vcf_in), str(out))
    [v] = list(cyvcf2.VCF(str(out)))
    assert tuple(v.INFO.get("AFQUERY_N_HOM_REF")) == (2, 2, 2)
    assert tuple(v.INFO.get("AFQUERY_N_HET")) == (2, 2, 0)


def test_variant_info_agrees(db):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AfqueryWarning)
        carriers = db.variant_info(*SITE)
    names = sorted(c.sample_name for c in carriers)
    assert names == ["AT", "AT", "A_HET", "A_HOM", "T_HET"]


# ---------------------------------------------------------------------------
# Coverage evidence is judged per position, not per allele
# ---------------------------------------------------------------------------

PANEL_COHORT = COHORT + [
    ("P_T",  "female", "PANEL", [("chr1", 5000, "G", "T", "0/1")]),
    ("P_R1", "female", "PANEL", [("chr1", 100, "C", "G", "0/1")]),
    ("P_R2", "male",   "PANEL", [("chr1", 100, "C", "G", "0/1")]),
]


@pytest.fixture(scope="module")
def panel_db(tmp_path_factory):
    return _build_db(
        tmp_path_factory.mktemp("multiallelic_panel"), PANEL_COHORT,
        beds={"PANEL": "chr1\t0\t6000\n"},
    )


def test_min_pass_counts_any_allele_as_evidence(panel_db):
    """The panel has no A carrier but one T carrier: it was sequenced here, so
    with --min-pass 1 its non-carriers stay hom-ref for A as well."""
    by_alt = {r.variant.alt: r for r in panel_db.query(*SITE, min_pass=1)}
    assert by_alt["A"].N_NO_COVERAGE == 0
    assert by_alt["A"].N_HOM_REF == 4     # REF1, REF2, P_R1, P_R2
    assert by_alt["T"].N_NO_COVERAGE == 0


def test_min_pass_failure_never_marks_carriers_uncovered(panel_db):
    """With --min-pass 2 the panel fails the gate; P_T has a call, so only the
    panel's true non-carriers move to N_NO_COVERAGE, for every allele."""
    by_alt = {r.variant.alt: r for r in panel_db.query(*SITE, min_pass=2)}
    for alt in ("A", "T"):
        assert by_alt[alt].N_NO_COVERAGE == 2, alt   # P_R1, P_R2
    assert by_alt["A"].N_HOM_REF == 2                 # REF1, REF2
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AfqueryWarning)
        carriers = panel_db.variant_info(*SITE, alt="A", min_pass=2)
    no_cov = sorted(c.sample_name for c in carriers if c.genotype == "no_coverage")
    assert no_cov == ["P_R1", "P_R2"]


def test_min_covered_judged_per_position(tmp_path):
    """--min-covered counts quality-passing carriers of every allele at the
    position: the panel carries A, so it was sequenced here for T as well."""
    db = _build_db(tmp_path, [
        ("W_T",  "female", "WGS",   [("chr1", 5000, "G", "T", "0/1")]),
        ("W_R",  "female", "WGS",   [("chr1", 100, "C", "G", "0/1")]),
        ("P_A1", "female", "PANEL", [("chr1", 5000, "G", "A", "0/1")]),
        ("P_A2", "male",   "PANEL", [("chr1", 5000, "G", "A", "0/1")]),
        ("P_R1", "female", "PANEL", [("chr1", 100, "C", "G", "0/1")]),
        ("P_R2", "male",   "PANEL", [("chr1", 100, "C", "G", "0/1")]),
    ], beds={"PANEL": "chr1\t0\t6000\n"}, min_covered=1)
    by_alt = {r.variant.alt: r for r in db.query(*SITE)}
    for alt in ("A", "T"):
        assert by_alt[alt].N_NO_COVERAGE == 0, alt
    assert by_alt["A"].N_HOM_REF == 3    # W_R, P_R1, P_R2
    assert by_alt["T"].N_HOM_REF == 3


def test_min_covered_still_flags_unevidenced_panel(tmp_path):
    """With no quality-passing carrier of any allele in the panel, its
    non-carriers stay uncovered for every allele at the position."""
    db = _build_db(tmp_path, [
        ("W_A",  "female", "WGS",   [("chr1", 5000, "G", "A", "0/1")]),
        ("W_T",  "female", "WGS",   [("chr1", 5000, "G", "T", "0/1")]),
        ("P_R1", "female", "PANEL", [("chr1", 100, "C", "G", "0/1")]),
        ("P_R2", "male",   "PANEL", [("chr1", 100, "C", "G", "0/1")]),
    ], beds={"PANEL": "chr1\t0\t6000\n"}, min_covered=1)
    by_alt = {r.variant.alt: r for r in db.query(*SITE)}
    for alt in ("A", "T"):
        assert by_alt[alt].N_NO_COVERAGE == 2, alt


def test_annotate_absent_allele_shares_site_coverage(panel_db, tmp_path):
    """C is not stored: under --min-pass 2 the panel's non-carriers are as
    uncovered for C as they are for the stored A and T."""
    vcf_in = tmp_path / "in.vcf"
    vcf_in.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t5000\t.\tG\tA,T,C\t.\tPASS\t.\n"
    )
    out = tmp_path / "out.vcf"
    panel_db.annotate_vcf(str(vcf_in), str(out), min_pass=2)
    [v] = list(cyvcf2.VCF(str(out)))
    assert tuple(v.INFO.get("AFQUERY_N_NO_COVERAGE")) == (2, 2, 2)
    assert tuple(v.INFO.get("AFQUERY_N_HOM_REF")) == (2, 2, 2)


# ---------------------------------------------------------------------------
# Failed calls are judged per position
# ---------------------------------------------------------------------------

FAIL_COHORT = [
    ("AT",      "female", "WGS", [("chr1", 5000, "G", "A,T", "1/2")]),
    ("AHET",    "female", "WGS", [("chr1", 5000, "G", "A", "0/1")]),
    ("THOM",    "male",   "WGS", [("chr1", 5000, "G", "T", "1/1")]),
    # As written by normalize_vcf.sh: split records, GT masked at a failed site
    ("F_SPLIT", "male",   "WGS", [("chr1", 5000, "G", "A", "./.", "LowQual"),
                                  ("chr1", 5000, "G", "T", "./.", "LowQual")]),
    ("F_MULTI", "male",   "WGS", [("chr1", 5000, "G", "A,T", "1/2", "LowQual")]),
    ("F_A",     "female", "WGS", [("chr1", 5000, "G", "A", "0/1", "LowQual")]),
    ("REF1",    "female", "WGS", [("chr1", 100, "C", "G", "0/1")]),
    ("REF2",    "male",   "WGS", [("chr1", 100, "C", "G", "0/1")]),
]


@pytest.fixture(scope="module")
def fail_db(tmp_path_factory):
    return _build_db(tmp_path_factory.mktemp("multiallelic_fail"), FAIL_COHORT)


def test_failed_calls_leave_an_for_every_allele(fail_db):
    """A failed call for any allele is out of AN: AN = AC_ref + sum of AC."""
    by_alt = {r.variant.alt: r for r in fail_db.query(*SITE)}
    assert (by_alt["A"].AC, by_alt["T"].AC) == (2, 3)
    for r in by_alt.values():
        assert r.AN == 10          # 5 called diploid samples
        assert r.N_FAIL == 3       # F_SPLIT, F_MULTI, F_A: the same for both alleles
        assert r.n_samples_eligible == 8
    assert (by_alt["A"].N_HOM_REF, by_alt["T"].N_HOM_REF) == (2, 2)
    ac_ref = 1 + 2 + 2             # AHET, REF1, REF2
    assert by_alt["A"].AN == ac_ref + by_alt["A"].AC + by_alt["T"].AC


def test_annotate_counts_each_failed_sample_once(fail_db, tmp_path):
    vcf_in = tmp_path / "in.vcf"
    vcf_in.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t5000\t.\tG\tA,T\t.\tPASS\t.\n"
    )
    out = tmp_path / "out.vcf"
    fail_db.annotate_vcf(str(vcf_in), str(out))
    [v] = list(cyvcf2.VCF(str(out)))
    assert v.INFO.get("AFQUERY_N_FAIL") == 3
    assert v.INFO.get("AFQUERY_AN") == 10


def test_variant_info_lists_site_failures_for_each_allele(fail_db):
    """F_A failed at the position, so it is listed as alt/FAIL for T as well,
    matching N_FAIL."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AfqueryWarning)
        carriers = fail_db.variant_info(*SITE, alt="T")
    failed = sorted(c.sample_name for c in carriers if c.genotype == "alt")
    assert failed == ["F_A", "F_MULTI", "F_SPLIT"]


def test_all_eligible_failed_gives_an_zero(tmp_path):
    db = _build_db(tmp_path, [
        ("F1", "female", "WGS", [("chr1", 5000, "G", "A", "0/1", "LowQual")]),
        ("F2", "male",   "WGS", [("chr1", 5000, "G", "A", "./.", "LowQual")]),
    ])
    [r] = db.query(*SITE)
    assert (r.AC, r.AN, r.AF, r.N_FAIL) == (0, 0, None, 2)
    vcf_in = tmp_path / "in.vcf"
    vcf_in.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t5000\t.\tG\tA\t.\tPASS\t.\n"
    )
    out = tmp_path / "out.vcf"
    db.annotate_vcf(str(vcf_in), str(out))
    [v] = list(cyvcf2.VCF(str(out)))
    assert v.INFO.get("AFQUERY_AN") == 0
    assert v.INFO.get("AFQUERY_N_FAIL") == 2
    assert "AFQUERY_AF=." in str(v)


def test_fail_as_ref_counts_failed_calls_in_an(fail_db):
    """--fail-as-ref only changes AN: the 3 failed samples count as hom-ref."""
    by_alt = {r.variant.alt: r for r in fail_db.query(*SITE, fail_as_ref=True)}
    for r in by_alt.values():
        assert r.AN == 16
        assert r.N_FAIL == 3
    assert (by_alt["A"].AC, by_alt["T"].AC) == (2, 3)
    assert (by_alt["A"].N_HOM_REF, by_alt["T"].N_HOM_REF) == (2, 2)


def test_fail_as_ref_in_annotate_and_dump(fail_db, tmp_path):
    vcf_in = tmp_path / "in.vcf"
    vcf_in.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t5000\t.\tG\tA,T\t.\tPASS\t.\n"
    )
    out = tmp_path / "out.vcf"
    fail_db.annotate_vcf(str(vcf_in), str(out), fail_as_ref=True)
    [v] = list(cyvcf2.VCF(str(out)))
    assert v.INFO.get("AFQUERY_AN") == 16
    assert v.INFO.get("AFQUERY_N_FAIL") == 3

    buf = io.StringIO()
    fail_db.dump(output=buf, by_sex=True, fail_as_ref=True)
    rows = [r for r in csv.DictReader(io.StringIO(buf.getvalue())) if r["pos"] == "5000"]
    for row in rows:
        [want] = fail_db.query(*SITE, alt=row["alt"], fail_as_ref=True)
        assert int(row["AN"]) == want.AN == 16
        for sex in ("male", "female"):
            [g] = fail_db.query(*SITE, alt=row["alt"], sex=sex, fail_as_ref=True)
            assert int(row[f"AN_{sex}"]) == g.AN

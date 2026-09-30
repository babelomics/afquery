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

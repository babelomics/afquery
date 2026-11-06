import json
from pathlib import Path

import pytest

from afquery import Database

DATA_DIR = Path(__file__).parent / "data"


def load_expected():
    return json.loads((DATA_DIR / "expected_results.json").read_text())


@pytest.mark.parametrize("case", load_expected())
def test_query_expected(test_db, case):
    db = Database(test_db)
    results = db.query(
        chrom=case["chrom"],
        pos=case["pos"],
        phenotype=case["phenotype"],
        sex=case["sex"],
    )

    # Find the result matching the expected alt
    match = next((r for r in results if r.variant.alt == case["alt"]), None)
    assert match is not None, (
        f"No result for alt={case['alt']} at {case['chrom']}:{case['pos']}"
    )
    assert match.AC == case["AC"], f"{case['description']}: AC mismatch"
    assert match.AN == case["AN"], f"{case['description']}: AN mismatch"
    assert match.n_samples_eligible == case["n_eligible"], \
        f"{case['description']}: n_eligible mismatch"
    for field in ("N_HET", "N_HOM_ALT", "N_HEMI"):
        if field in case:
            assert getattr(match, field) == case[field], f"{case['description']}: {field}"

    if case["AF"] is None:
        assert match.AF is None
    else:
        assert abs(match.AF - case["AF"]) < 1e-9, f"{case['description']}: AF mismatch"


def test_query_an_zero_returns_empty(test_db):
    # C50={2,4,6}, male filter → only sample 4 (male, WES_A tech_id=1)
    # pos=3500 is in WES_B only (3000-4000) → sample 4 NOT covered → eligible={}
    db = Database(test_db)
    results = db.query(chrom="chr1", pos=3500, phenotype=["C50"], sex="male")
    assert results == []


def test_query_no_variant_at_position(test_db):
    # Position that is covered but has no variant in the Parquet
    db = Database(test_db)
    results = db.query(chrom="chr1", pos=1000, phenotype=["E11.9"], sex="both")
    assert results == []


def test_query_unknown_phenotype(test_db):
    db = Database(test_db)
    results = db.query(chrom="chr1", pos=1500, phenotype=["UNKNOWN"], sex="both")
    assert results == []


def test_query_chrom_normalization(test_db):
    # "1" should work the same as "chr1"
    db = Database(test_db)
    r1 = db.query(chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both")
    r2 = db.query(chrom="1",    pos=1500, phenotype=["E11.9"], sex="both")
    assert len(r1) == len(r2)
    assert r1[0].AC == r2[0].AC


def test_query_multiple_phenotype_union(test_db):
    # E11.9 + I10 should give more eligible than E11.9 alone (larger union)
    db = Database(test_db)
    r_e119 = db.query(chrom="chr1", pos=1500, phenotype=["E11.9"],        sex="both")
    r_both = db.query(chrom="chr1", pos=1500, phenotype=["E11.9", "I10"], sex="both")
    assert r_both[0].n_samples_eligible >= r_e119[0].n_samples_eligible


def test_query_nonexistent_chrom(test_db):
    db = Database(test_db)
    results = db.query(chrom="chr22", pos=1000, phenotype=["E11.9"], sex="both")
    assert results == []


def test_info(test_db):
    db = Database(test_db)
    data = db.info()
    assert data["genome_build"] == "GRCh37"
    assert data["sample_count"] == 10


def test_query_ref_alt_filter_exact(test_db):
    # Query at position 1500 which has multiple alleles, then filter to specific ref/alt
    db = Database(test_db)
    # Get all results first
    all_results = db.query(chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both")
    assert len(all_results) > 0, "Test position should have at least one variant"

    # Get the first result's ref/alt
    first_ref = all_results[0].variant.ref
    first_alt = all_results[0].variant.alt

    # Query with that specific ref/alt
    filtered_results = db.query(
        chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both",
        ref=first_ref, alt=first_alt
    )

    # Should get exactly one result matching that allele
    assert len(filtered_results) == 1
    assert filtered_results[0].variant.ref == first_ref
    assert filtered_results[0].variant.alt == first_alt
    assert filtered_results[0].AC == all_results[0].AC
    assert filtered_results[0].AN == all_results[0].AN


def test_query_ref_alt_filter_no_match(test_db):
    # Query with ref/alt that doesn't exist at the position
    db = Database(test_db)
    results = db.query(
        chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both",
        ref="NONEXIST", alt="ALLELE"
    )
    assert results == []


def test_query_empty_phenotype_all_samples(test_db):
    # Query with no phenotype filter should return results with all samples eligible
    db = Database(test_db)
    results = db.query(chrom="chr1", pos=1500, phenotype=[], sex="both")
    assert len(results) > 0, "Query with empty phenotype should return results"
    # With no phenotype filter, all samples should be eligible (or a larger subset)
    assert results[0].n_samples_eligible > 0


def test_query_empty_phenotype_vs_filtered(test_db):
    # Empty phenotype should include more eligible samples than a single phenotype filter
    db = Database(test_db)
    results_no_pheno = db.query(chrom="chr1", pos=1500, phenotype=[], sex="both")
    results_e119 = db.query(chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both")

    assert len(results_no_pheno) > 0, "Empty phenotype query should return results"
    assert len(results_e119) > 0, "E11.9 phenotype query should return results"

    # The unfiltered query should have more or equal eligible samples
    assert results_no_pheno[0].n_samples_eligible >= results_e119[0].n_samples_eligible

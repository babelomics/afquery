"""Tests for SampleFilter, tech filtering, phenotype exclusion, and CLI helpers."""
import json
from pathlib import Path

import pytest

from afquery import Database
from afquery.models import SampleFilter
from afquery.cli import _expand_tokens

# Test DB setup (from conftest.py):
# SAMPLES: {0,1,2,3}=WGS, {4,5,6}=WES_kit_A, {7,8,9}=WES_kit_B
# E11.9={0,1,2,5,6,7}, I10={1,3,5,8,9}, C50={2,4,6}, J45={0,3,7,8,9}, G20={1,2,3,5}
# chr1 pos=1500: het=[5], hom=[2], fail=[0]
#   BED coverage: WES_kit_A covers chr1 1000-1999, WES_kit_B covers chr1 3000-3999
#   So covered at 1500: WGS={0,1,2,3} + WES_kit_A={4,5,6} = {0,1,2,3,4,5,6}


# ---------------------------------------------------------------------------
# SampleFilter.parse unit tests
# ---------------------------------------------------------------------------

def test_sample_filter_parse_include_only():
    sf = SampleFilter.parse(["E11.9", "I10"], [], sex="both")
    assert sf.phenotype_include == ["E11.9", "I10"]
    assert sf.phenotype_exclude == []
    assert sf.tech_include == []
    assert sf.tech_exclude == []
    assert sf.sex == "both"


def test_sample_filter_parse_exclude_prefix():
    sf = SampleFilter.parse(["E11.9", "^I10"], ["WGS", "^WES_kit_A"], sex="male")
    assert sf.phenotype_include == ["E11.9"]
    assert sf.phenotype_exclude == ["I10"]
    assert sf.tech_include == ["WGS"]
    assert sf.tech_exclude == ["WES_kit_A"]
    assert sf.sex == "male"


def test_sample_filter_parse_all_exclude():
    sf = SampleFilter.parse(["^E11.9", "^I10"], ["^WES_kit_A"], sex="female")
    assert sf.phenotype_include == []
    assert sf.phenotype_exclude == ["E11.9", "I10"]
    assert sf.tech_include == []
    assert sf.tech_exclude == ["WES_kit_A"]
    assert sf.sex == "female"


def test_sample_filter_parse_empty():
    sf = SampleFilter.parse([], [], sex="both")
    assert sf.phenotype_include == []
    assert sf.phenotype_exclude == []
    assert sf.tech_include == []
    assert sf.tech_exclude == []


# ---------------------------------------------------------------------------
# Phenotype exclusion via ^ prefix
# ---------------------------------------------------------------------------

def test_filter_phenotype_exclusion(test_db):
    # Exclude E11.9={0,1,2,5,6,7} → remaining {3,4,8,9}
    # covered at 1500: {0..6} → eligible = {3,4}
    # het=[5], hom=[2], fail=[0] → none in {3,4} → AC=0, AN=4 (2 diploid samples)
    db = Database(test_db)
    results = db.query(chrom="chr1", pos=1500, phenotype=["^E11.9"], sex="both")
    assert len(results) > 0, "Should return results when AN > 0"
    r = results[0]
    assert r.n_samples_eligible == 2, "Only {3,4} covered and not in E11.9"
    assert r.AC == 0, "None of {3,4} carry the variant"
    assert r.AN == 4  # 2 diploid autosomes


def test_filter_phenotype_include_and_exclude(test_db):
    # E11.9={0,1,2,5,6,7}, I10={1,3,5,8,9}
    # phenotype_include=["E11.9"], exclude=["I10"] → {0,2,6,7}
    # covered at 1500: {0..6} → eligible = {0,2,6}
    db = Database(test_db)
    r_combined = db.query(chrom="chr1", pos=1500, phenotype=["E11.9", "^I10"], sex="both")
    r_e119_only = db.query(chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both")

    assert len(r_combined) > 0
    # Excluding I10 from E11.9 should reduce eligible samples
    assert r_combined[0].n_samples_eligible < r_e119_only[0].n_samples_eligible


def test_filter_phenotype_include_exclude_same(test_db):
    # Include and exclude same code → empty bitmap → AN=0 → no results
    db = Database(test_db)
    results = db.query(chrom="chr1", pos=1500, phenotype=["E11.9", "^E11.9"], sex="both")
    assert results == []


# ---------------------------------------------------------------------------
# Tech filter: inclusion
# ---------------------------------------------------------------------------

def test_filter_tech_inclusion(test_db):
    # tech=["WGS"] → sample_bm = {0,1,2,3}
    # covered at 1500 = {0..6} → eligible = {0,1,2,3}
    db = Database(test_db)
    r_wgs = db.query(chrom="chr1", pos=1500, sex="both", tech=["WGS"])
    r_all = db.query(chrom="chr1", pos=1500, sex="both")

    assert len(r_wgs) > 0
    # WGS-only should have fewer eligible than no-filter (which includes WES_A too)
    assert r_wgs[0].n_samples_eligible < r_all[0].n_samples_eligible
    assert r_wgs[0].n_samples_eligible == 4  # {0,1,2,3}


def test_filter_tech_inclusion_wes_a(test_db):
    # tech=["WES_kit_A"] → sample_bm = {4,5,6}
    # covered at 1500: WES_kit_A covers 1000-1999 → {4,5,6} covered
    # eligible = {4,5,6}
    db = Database(test_db)
    results = db.query(chrom="chr1", pos=1500, sex="both", tech=["WES_kit_A"])
    assert len(results) > 0
    assert results[0].n_samples_eligible == 3  # {4,5,6}


def test_filter_tech_inclusion_wes_b_no_coverage(test_db):
    # tech=["WES_kit_B"] → sample_bm = {7,8,9}
    # WES_kit_B covers pos 3000-3999, NOT pos 1500 → eligible = {}
    db = Database(test_db)
    results = db.query(chrom="chr1", pos=1500, sex="both", tech=["WES_kit_B"])
    assert results == [], "WES_kit_B does not cover pos=1500"


# ---------------------------------------------------------------------------
# Tech filter: exclusion
# ---------------------------------------------------------------------------

def test_filter_tech_exclusion(test_db):
    # tech=["^WES_kit_A"] → exclude {4,5,6} → sample_bm = {0,1,2,3,7,8,9}
    # covered at 1500: {0..6} → eligible = {0,1,2,3}
    db = Database(test_db)
    r_excl = db.query(chrom="chr1", pos=1500, sex="both", tech=["^WES_kit_A"])
    r_all  = db.query(chrom="chr1", pos=1500, sex="both")

    assert len(r_excl) > 0
    assert r_excl[0].n_samples_eligible < r_all[0].n_samples_eligible
    assert r_excl[0].n_samples_eligible == 4  # {0,1,2,3} — WES_B not covered at 1500


# ---------------------------------------------------------------------------
# Combined phenotype + tech filter
# ---------------------------------------------------------------------------

def test_filter_combined(test_db):
    # phenotype=["E11.9", "^I10"], tech=["WGS"]
    # E11.9={0,1,2,5,6,7} minus I10∩E11.9={1,5} → ph_bm={0,2,6,7}
    # tech WGS → {0,1,2,3}
    # sample_bm = {0,2,6,7} & {0,1,2,3} = {0,2}
    # covered at 1500: {0..6} → eligible = {0,2}
    # fail=[0]→{0}: S00's call failed → out of AC and AN
    # hom=[2]→{2} in eligible → AC = 2, AN = 2
    db = Database(test_db)
    results = db.query(
        chrom="chr1", pos=1500,
        phenotype=["E11.9", "^I10"], sex="both", tech=["WGS"],
    )
    assert len(results) > 0
    r = results[0]
    assert r.n_samples_eligible == 2  # {0,2}
    assert r.AC == 2  # sample 2 hom(+2)
    assert r.AN == 2  # sample 2 only; sample 0 has a failed call
    assert r.N_FAIL == 1


# ---------------------------------------------------------------------------
# _expand_tokens unit tests
# ---------------------------------------------------------------------------

def test_expand_tokens_empty():
    assert _expand_tokens(()) == []


def test_expand_tokens_single():
    assert _expand_tokens(("E11.9",)) == ["E11.9"]


def test_expand_tokens_comma():
    assert _expand_tokens(("E11.9,I10",)) == ["E11.9", "I10"]


def test_expand_tokens_multiple_args():
    assert _expand_tokens(("E11.9", "I10")) == ["E11.9", "I10"]


def test_expand_tokens_mixed():
    assert _expand_tokens(("E11.9,I10", "^C50")) == ["E11.9", "I10", "^C50"]


def test_expand_tokens_with_spaces():
    assert _expand_tokens(("E11.9, I10",)) == ["E11.9", "I10"]


def test_expand_tokens_exclusion():
    assert _expand_tokens(("^E11.9",)) == ["^E11.9"]


# ---------------------------------------------------------------------------
# annotate without --phenotype uses all samples
# ---------------------------------------------------------------------------

def test_annotate_no_phenotype_uses_all_samples(test_db, tmp_path):
    """annotate_vcf with phenotype=None should use all samples (no restriction)."""
    db = Database(test_db)
    # Query all samples explicitly
    r_all = db.query(chrom="chr1", pos=1500, sex="both")
    # Query with phenotype=None explicitly
    r_none = db.query(chrom="chr1", pos=1500, phenotype=None, sex="both")
    assert len(r_all) > 0
    assert r_all[0].n_samples_eligible == r_none[0].n_samples_eligible
    assert r_all[0].AC == r_none[0].AC
    assert r_all[0].AN == r_none[0].AN


# ---------------------------------------------------------------------------
# CLI integration tests for --tech and --phenotype ^
# ---------------------------------------------------------------------------

def test_cli_tech_filter(test_db, tmp_path):
    """CLI --tech WGS should restrict to WGS samples only."""
    from click.testing import CliRunner
    from afquery.cli import cli

    runner = CliRunner()
    result_wgs = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500",
        "--tech", "WGS", "--format", "json",
    ])
    result_all = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500",
        "--format", "json",
    ])
    assert result_wgs.exit_code == 0, result_wgs.output
    assert result_all.exit_code == 0, result_all.output

    data_wgs = json.loads(result_wgs.output)
    data_all = json.loads(result_all.output)
    assert len(data_wgs) > 0
    assert data_wgs[0]["n_eligible"] < data_all[0]["n_eligible"]


def test_cli_phenotype_exclusion(test_db, tmp_path):
    """CLI --phenotype ^E11.9 should exclude E11.9 samples."""
    from click.testing import CliRunner
    from afquery.cli import cli

    runner = CliRunner()
    result_excl = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500",
        "--phenotype", "^E11.9", "--format", "json",
    ])
    result_all = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500",
        "--format", "json",
    ])
    assert result_excl.exit_code == 0, result_excl.output
    assert result_all.exit_code == 0, result_all.output

    data_excl = json.loads(result_excl.output)
    data_all = json.loads(result_all.output)
    assert len(data_excl) > 0
    # Excluding E11.9 should reduce eligible samples
    assert data_excl[0]["n_eligible"] < data_all[0]["n_eligible"]


def test_cli_comma_phenotype(test_db):
    """CLI --phenotype E11.9,I10 should work the same as --phenotype E11.9 --phenotype I10."""
    from click.testing import CliRunner
    from afquery.cli import cli

    runner = CliRunner()
    result_comma = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500",
        "--phenotype", "E11.9,I10", "--format", "json",
    ])
    result_repeat = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500",
        "--phenotype", "E11.9", "--phenotype", "I10", "--format", "json",
    ])
    assert result_comma.exit_code == 0
    assert result_repeat.exit_code == 0

    data_comma = json.loads(result_comma.output)
    data_repeat = json.loads(result_repeat.output)
    assert len(data_comma) > 0
    assert data_comma[0]["n_eligible"] == data_repeat[0]["n_eligible"]
    assert data_comma[0]["AC"] == data_repeat[0]["AC"]

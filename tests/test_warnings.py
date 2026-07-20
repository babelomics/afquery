import shutil
import warnings

import pytest
from click.testing import CliRunner

from afquery import Database, AfqueryWarning
from afquery.capture import CaptureIndex
from afquery.cli import query as query_cmd


def test_warn_unknown_phenotype_include(test_db):
    db = Database(test_db)
    with pytest.warns(AfqueryWarning, match="FAKE_CODE.*include will match 0 samples"):
        results = db.query(chrom="chr1", pos=1500, phenotype=["FAKE_CODE"], sex="both")
    assert results == []


def test_warn_unknown_phenotype_exclude(test_db):
    db = Database(test_db)
    with pytest.warns(AfqueryWarning, match="FAKE_CODE.*exclude has no effect"):
        results = db.query(chrom="chr1", pos=1500, phenotype=["^FAKE_CODE"], sex="both")
    # Exclude of unknown phenotype should not affect results
    assert len(results) > 0


def test_warn_unknown_tech_include(test_db):
    db = Database(test_db)
    with pytest.warns(AfqueryWarning, match="FAKE_TECH.*include will match 0 samples"):
        results = db.query(chrom="chr1", pos=1500, phenotype=[], sex="both", tech=["FAKE_TECH"])
    assert results == []


def test_warn_unknown_tech_exclude(test_db):
    db = Database(test_db)
    with pytest.warns(AfqueryWarning, match="FAKE_TECH.*exclude has no effect"):
        results = db.query(chrom="chr1", pos=1500, phenotype=[], sex="both", tech=["^FAKE_TECH"])
    assert len(results) > 0


def test_warn_contradictory_phenotype_filter(test_db):
    db = Database(test_db)
    # Include and exclude the same code — produces empty eligible set
    with pytest.warns(AfqueryWarning, match="empty eligible set"):
        results = db.query(chrom="chr1", pos=1500, phenotype=["E11.9", "^E11.9"], sex="both")
    assert results == []


def test_warn_chrom_not_in_db_query(test_db):
    db = Database(test_db)
    with pytest.warns(AfqueryWarning, match="chr22.*has no data"):
        results = db.query(chrom="chr22", pos=1000, phenotype=[], sex="both")
    assert results == []


def test_warn_chrom_not_in_db_query_batch(test_db):
    db = Database(test_db)
    with pytest.warns(AfqueryWarning, match="chr22.*has no data"):
        results = db.query_batch(chrom="chr22", variants=[(1000, "A", "T")], phenotype=[], sex="both")
    assert results == []


def test_warn_chrom_not_in_db_query_region(test_db):
    db = Database(test_db)
    with pytest.warns(AfqueryWarning, match="chr22.*has no data"):
        results = db.query_region(chrom="chr22", start=1000, end=2000, phenotype=[], sex="both")
    assert results == []


def test_invalid_sex_raises_valueerror(test_db):
    db = Database(test_db)
    with pytest.raises(ValueError, match="Invalid sex"):
        db.query(chrom="chr1", pos=1500, phenotype=[], sex="xyz")


def test_no_warn_suppresses(test_db):
    runner = CliRunner()
    result = runner.invoke(query_cmd, [
        "--db", test_db,
        "--locus", "chr22:1000",
        "--no-warn",
    ])
    assert result.exit_code == 0
    assert "Warning" not in (result.output or "")


def test_existing_tests_not_broken_by_warnings(test_db):
    # Ensure previously passing test cases don't emit unexpected warnings
    db = Database(test_db)
    with warnings.catch_warnings():
        warnings.simplefilter("error", AfqueryWarning)
        results = db.query(chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both")
    assert len(results) > 0


def test_warn_capture_bed_matches_no_known_chrom(test_db, tmp_path):
    # A capture BED whose chrom names match nothing would silently drop that
    # technology's samples from AN at every position — it must be loud.
    db_copy = tmp_path / "db_bad_bed"
    shutil.copytree(test_db, db_copy)
    bad_bed = tmp_path / "bad.bed"
    bad_bed.write_text("contigZ\t100\t200\n")
    CaptureIndex.from_bed(str(bad_bed)).save(str(db_copy / "capture" / "tech_1.pickle"))

    with pytest.warns(AfqueryWarning, match="WES_kit_A.*match no known chromosome"):
        Database(str(db_copy))


def test_no_capture_warning_for_valid_db(test_db):
    with warnings.catch_warnings():
        warnings.simplefilter("error", AfqueryWarning)
        Database(test_db)


def test_warn_chrom_message_includes_available(test_db):
    db = Database(test_db)
    with pytest.warns(AfqueryWarning) as record:
        db.query(chrom="chr99", pos=1000, phenotype=[], sex="both")
    msg = str(record[0].message)
    assert "chr99" in msg
    assert "Available" in msg

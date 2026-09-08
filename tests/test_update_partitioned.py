"""add-samples against a bucketed database.

Every other add_samples test runs on the hand-built flat fixture in conftest.py,
which is the one layout create-db never produces. These build a real database
through run_preprocess, so the merge is exercised against the layout it will
actually meet.

Cohort (tests/data): S00-S03 wgs, S04-S06 wes_kit_a (chr1:1000-2000,
chrX:1000-2000), S07-S09 wes_kit_b (chr1:3000-4000). chr1:1500 A>T is carried
by S00 het, S02 hom and S05 het.
"""
import shutil
import sqlite3

import pyarrow.parquet as pq
import pytest

from afquery import storage
from afquery.bitmaps import deserialize
from afquery.database import Database
from afquery.preprocess import run_preprocess
from afquery.preprocess.update import add_samples, check_database

from test_update import write_manifest, write_vcf


# ---------------------------------------------------------------------------
# Fixtures — build once per module, copy per test
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def partitioned_src(tmp_path_factory, data_dir):
    db = tmp_path_factory.mktemp("partitioned_src")
    run_preprocess(
        manifest_path=str(data_dir / "manifest.tsv"), output_dir=str(db),
        genome_build="GRCh37", bed_dir=str(data_dir / "beds"), threads=2,
    )
    return str(db)


@pytest.fixture
def partitioned_db(partitioned_src, tmp_path):
    dest = tmp_path / "db"
    shutil.copytree(partitioned_src, dest)
    return str(dest)


@pytest.fixture(scope="module")
def covered_src(tmp_path_factory, data_dir):
    """Same cohort with a coverage-evidence threshold, so filtered_bitmap is live."""
    db = tmp_path_factory.mktemp("covered_src")
    run_preprocess(
        manifest_path=str(data_dir / "manifest.tsv"), output_dir=str(db),
        genome_build="GRCh37", bed_dir=str(data_dir / "beds"), threads=2,
        min_covered=2,
    )
    return str(db)


@pytest.fixture
def covered_db(covered_src, tmp_path):
    dest = tmp_path / "covered"
    shutil.copytree(covered_src, dest)
    return str(dest)


def _add_one(db, tmp_path, name, variants, tech="wgs", sex="male",
             phenotype="ZZNEW", bed_dir=None):
    vcf = str(tmp_path / f"{name}.vcf")
    write_vcf(vcf, name, variants)
    manifest = str(tmp_path / f"{name}.tsv")
    write_manifest(manifest, [(name, sex, tech, vcf, phenotype)])
    return add_samples(db, manifest, threads=1, bed_dir=bed_dir)


def _row(db, chrom, bucket, pos, ref, alt):
    """One row of a bucket file as deserialized bitmaps, or None."""
    path = storage.bucket_path(f"{db}/variants", chrom, bucket)
    if not path.exists():
        return None
    t = pq.read_table(str(path))
    for i in range(len(t)):
        if (t["pos"][i].as_py(), t["ref"][i].as_py(), t["alt"][i].as_py()) == (pos, ref, alt):
            return {
                name: deserialize(t[name][i].as_py())
                for name in ("het_bitmap", "hom_bitmap", "fail_bitmap",
                             "filtered_bitmap", "quality_pass_bitmap")
            }
    return None


# ---------------------------------------------------------------------------
# The regression: added samples must be visible to queries
# ---------------------------------------------------------------------------

def test_add_samples_partitioned_new_variant_is_queryable(partitioned_db, tmp_path):
    _add_one(partitioned_db, tmp_path, "S10", [("chr1", 7000, "A", "T", "0/1")])

    db = Database(partitioned_db)
    results = db.query(chrom="chr1", pos=7000, sex="both")
    assert results, "variant added by add-samples is not visible to the query engine"
    assert any(r.AC >= 1 for r in results)


def test_add_samples_partitioned_existing_variant_gains_carrier(partitioned_db, tmp_path):
    before = Database(partitioned_db).query(chrom="chr1", pos=1500, sex="both")
    ac_before = before[0].AC
    an_before = before[0].AN

    # S10 is hom-alt, so it contributes 2 to AC and 2 to AN
    _add_one(partitioned_db, tmp_path, "S10", [("chr1", 1500, "A", "T", "1/1")])

    after = Database(partitioned_db).query(chrom="chr1", pos=1500, sex="both")
    assert after[0].AC == ac_before + 2
    assert after[0].AN == an_before + 2


def test_add_samples_partitioned_sets_carrier_bit(partitioned_db, tmp_path):
    """The new sample id must be a set bit in the bucket, not merely counted."""
    _add_one(partitioned_db, tmp_path, "S10", [("chr1", 1500, "A", "T", "1/1")])

    con = sqlite3.connect(f"{partitioned_db}/metadata.sqlite")
    sid = con.execute("SELECT sample_id FROM samples WHERE sample_name='S10'").fetchone()[0]
    con.close()

    row = _row(partitioned_db, "chr1", 0, 1500, "A", "T")
    assert row is not None
    assert sid in row["hom_bitmap"]
    assert sid not in row["het_bitmap"]


def test_add_samples_partitioned_phenotype_isolated(partitioned_db, tmp_path):
    """Filtered to the new sample alone the counts must describe that sample."""
    _add_one(partitioned_db, tmp_path, "S10", [("chr1", 1500, "A", "T", "1/1")])

    res = Database(partitioned_db).query(
        chrom="chr1", pos=1500, phenotype=["ZZNEW"], sex="both"
    )
    assert len(res) == 1
    r = res[0]
    assert (r.AC, r.AN, r.N_HET, r.N_HOM_ALT, r.N_HOM_REF) == (2, 2, 0, 1, 0)


# ---------------------------------------------------------------------------
# Layout is preserved, and never split
# ---------------------------------------------------------------------------

def test_add_samples_preserves_chrom_layouts(partitioned_db, tmp_path):
    variants = f"{partitioned_db}/variants"
    before = {c: storage.chrom_layout(variants, c)
              for c in storage.partitioned_chroms(variants) | storage.flat_chroms(variants)}

    _add_one(partitioned_db, tmp_path, "S10", [("chr1", 7000, "A", "T", "0/1")])

    after = {c: storage.chrom_layout(variants, c) for c in before}
    assert after == before
    assert storage.mixed_layout_chroms(variants) == []


def test_add_samples_partitioned_creates_new_bucket(partitioned_db, tmp_path):
    assert 2 not in storage.existing_bucket_ids(f"{partitioned_db}/variants", "chr1")

    _add_one(partitioned_db, tmp_path, "S10", [("chr1", 2_500_000, "C", "G", "0/1")])

    assert 2 in storage.existing_bucket_ids(f"{partitioned_db}/variants", "chr1")
    assert _row(partitioned_db, "chr1", 2, 2_500_000, "C", "G") is not None
    res = Database(partitioned_db).query(chrom="chr1", pos=2_500_000, sex="both")
    assert res and res[0].AC == 1


def test_add_samples_partitioned_new_chromosome_uses_buckets(partitioned_db, tmp_path):
    variants = f"{partitioned_db}/variants"
    assert "chr2" not in storage.partitioned_chroms(variants)

    _add_one(partitioned_db, tmp_path, "S10", [("chr2", 1000, "G", "A", "0/1")])

    assert storage.chrom_layout(variants, "chr2") == storage.PARTITIONED
    assert storage.existing_bucket_ids(variants, "chr2") == [0]
    assert not storage.flat_path(variants, "chr2").exists()
    res = Database(partitioned_db).query(chrom="chr2", pos=1000, sex="both")
    assert res and res[0].AC == 1


def test_add_samples_partitioned_check_database_clean(partitioned_db, tmp_path):
    _add_one(partitioned_db, tmp_path, "S10",
             [("chr1", 1500, "A", "T", "1/1"), ("chr1", 2_500_000, "C", "G", "0/1")])

    errors = [r for r in check_database(partitioned_db) if r.severity == "error"]
    assert errors == []


# ---------------------------------------------------------------------------
# Phase 2: filtered_bitmap depends on the cohort, not on which bucket changed
# ---------------------------------------------------------------------------

def test_add_samples_partitioned_recomputes_filtered_bitmap_in_untouched_buckets(
    covered_db, tmp_path, data_dir
):
    """A sample whose only variant is in bucket 2 still changes bucket 0.

    filtered_bitmap is (tech_bm - carriers) for every WES tech short of
    min_covered quality-passing carriers. Enlarging wes_kit_a therefore moves it
    at chr1:1500, which lives in a bucket the new sample contributed nothing to.
    """
    before = _row(covered_db, "chr1", 0, 1500, "A", "T")
    assert before is not None
    assert before["filtered_bitmap"], "fixture should have a live filtered_bitmap"

    _add_one(covered_db, tmp_path, "S10", [("chr1", 2_500_000, "C", "G", "0/1")],
             tech="wes_kit_a", bed_dir=str(data_dir / "beds"))

    con = sqlite3.connect(f"{covered_db}/metadata.sqlite")
    sid = con.execute("SELECT sample_id FROM samples WHERE sample_name='S10'").fetchone()[0]
    con.close()

    after = _row(covered_db, "chr1", 0, 1500, "A", "T")
    assert sid in after["filtered_bitmap"], (
        "bucket 0 kept a stale filtered_bitmap: the new wes_kit_a sample is a "
        "non-carrier of an under-covered tech and must be marked"
    )
    assert before["filtered_bitmap"] < after["filtered_bitmap"]


def test_add_samples_partitioned_untouched_chrom_is_not_rewritten(covered_db, tmp_path, data_dir):
    """Chromosomes the batch never touches keep their files as they were."""
    chrX = storage.bucket_path(f"{covered_db}/variants", "chrX", 5)
    before = chrX.read_bytes()

    _add_one(covered_db, tmp_path, "S10", [("chr1", 2_500_000, "C", "G", "0/1")],
             tech="wes_kit_a", bed_dir=str(data_dir / "beds"))

    assert chrX.read_bytes() == before

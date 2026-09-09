"""Phase 4 tests: add_samples, remove_samples, check_database."""
import json
import os
import shutil
import sqlite3

import pyarrow.parquet as pq
import pytest

from afquery.bitmaps import deserialize
from afquery.database import Database
from afquery.preprocess.update import (
    CheckResult,
    UpdateError,
    add_samples,
    check_database,
    remove_samples,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def write_vcf(path: str, sample_name: str, variants: list) -> None:
    """Write a minimal single-sample VCF. variants: [(chrom, pos, ref, alt, gt), ...]"""
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
        for chrom, pos, ref, alt, gt in variants:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")


def write_manifest(path: str, entries: list) -> None:
    """Write a TSV manifest. entries: [(name, sex, tech, vcf_path, phenotype_csv), ...]"""
    with open(path, "w") as f:
        f.write("sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n")
        for name, sex, tech, vcf, codes in entries:
            f.write(f"{name}\t{sex}\t{tech}\t{vcf}\t{codes}\n")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def fresh_db(test_db, tmp_path):
    """A mutable per-test copy of the session test_db."""
    dest = tmp_path / "db"
    shutil.copytree(test_db, dest)
    return str(dest)


# ---------------------------------------------------------------------------
# add_samples tests
# ---------------------------------------------------------------------------

def test_add_samples_basic(fresh_db, tmp_path):
    vcf_s10 = str(tmp_path / "S10.vcf")
    vcf_s11 = str(tmp_path / "S11.vcf")
    write_vcf(vcf_s10, "S10", [("chr1", 7000, "A", "T", "0/1")])
    write_vcf(vcf_s11, "S11", [("chr1", 9000, "G", "C", "0/1")])

    manifest = str(tmp_path / "manifest.tsv")
    write_manifest(manifest, [
        ("S10", "male",   "wgs", vcf_s10, "E11.9"),
        ("S11", "female", "wgs", vcf_s11, "I10"),
    ])

    result = add_samples(fresh_db, manifest, threads=1)

    assert result["new_samples"] == 2

    con = sqlite3.connect(os.path.join(fresh_db, "metadata.sqlite"))
    count = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    s10_id = con.execute(
        "SELECT sample_id FROM samples WHERE sample_name='S10'"
    ).fetchone()[0]
    s11_id = con.execute(
        "SELECT sample_id FROM samples WHERE sample_name='S11'"
    ).fetchone()[0]
    con.close()

    assert count == 12
    assert s10_id == 10
    assert s11_id == 11

    manifest_data = json.loads(open(os.path.join(fresh_db, "manifest.json")).read())
    assert manifest_data["sample_count"] == 12


def test_add_samples_new_tech(fresh_db, tmp_path):
    bed_dir = str(tmp_path / "beds")
    os.makedirs(bed_dir)
    with open(os.path.join(bed_dir, "wes_kit_c.bed"), "w") as f:
        f.write("chr1\t1000\t10000\n")

    vcf_s10 = str(tmp_path / "S10.vcf")
    write_vcf(vcf_s10, "S10", [("chr1", 5000, "T", "G", "0/1")])

    manifest = str(tmp_path / "manifest.tsv")
    write_manifest(manifest, [("S10", "male", "wes_kit_c", vcf_s10, "E11.9")])

    add_samples(fresh_db, manifest, threads=1, bed_dir=bed_dir)

    con = sqlite3.connect(os.path.join(fresh_db, "metadata.sqlite"))
    row = con.execute(
        "SELECT tech_id FROM technologies WHERE tech_name='wes_kit_c'"
    ).fetchone()
    con.close()

    assert row is not None
    assert row[0] == 3  # Fourth tech, 0-indexed


def test_add_samples_existing_tech(fresh_db, tmp_path):
    vcf_s10 = str(tmp_path / "S10.vcf")
    write_vcf(vcf_s10, "S10", [("chr1", 7000, "A", "T", "0/1")])

    manifest = str(tmp_path / "manifest.tsv")
    write_manifest(manifest, [("S10", "male", "WGS", vcf_s10, "E11.9")])

    add_samples(fresh_db, manifest, threads=1)

    con = sqlite3.connect(os.path.join(fresh_db, "metadata.sqlite"))
    tech_count = con.execute("SELECT COUNT(*) FROM technologies").fetchone()[0]
    con.close()

    assert tech_count == 3  # Still 3 technologies; WGS reused


def test_add_samples_id_monotonicity(fresh_db, tmp_path):
    vcf_s10 = str(tmp_path / "S10.vcf")
    write_vcf(vcf_s10, "S10", [])
    manifest1 = str(tmp_path / "manifest1.tsv")
    write_manifest(manifest1, [("S10", "male", "wgs", vcf_s10, "E11.9")])
    add_samples(fresh_db, manifest1, threads=1)

    vcf_s11 = str(tmp_path / "S11.vcf")
    write_vcf(vcf_s11, "S11", [])
    manifest2 = str(tmp_path / "manifest2.tsv")
    write_manifest(manifest2, [("S11", "female", "wgs", vcf_s11, "I10")])
    add_samples(fresh_db, manifest2, threads=1)

    con = sqlite3.connect(os.path.join(fresh_db, "metadata.sqlite"))
    s10_id = con.execute(
        "SELECT sample_id FROM samples WHERE sample_name='S10'"
    ).fetchone()[0]
    s11_id = con.execute(
        "SELECT sample_id FROM samples WHERE sample_name='S11'"
    ).fetchone()[0]
    con.close()

    assert s10_id == 10
    assert s11_id == 11


def test_add_samples_new_variant(fresh_db, tmp_path):
    vcf_s10 = str(tmp_path / "S10.vcf")
    write_vcf(vcf_s10, "S10", [("chr1", 7000, "A", "T", "0/1")])

    manifest = str(tmp_path / "manifest.tsv")
    write_manifest(manifest, [("S10", "male", "wgs", vcf_s10, "E11.9")])

    add_samples(fresh_db, manifest, threads=1)

    db = Database(fresh_db)
    results = db.query(chrom="chr1", pos=7000, phenotype=["E11.9"], sex="both")
    assert results
    assert any(r.AC >= 1 for r in results)


def test_add_samples_existing_variant(fresh_db, tmp_path):
    # chr1:1500 A>T has het=[0,5], hom=[2] in the base DB
    # S10 (male, wgs, E11.9) also carries chr1:1500 het → AC should increase
    db_before = Database(fresh_db)
    results_before = db_before.query(chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both")
    ac_before = results_before[0].AC if results_before else 0

    vcf_s10 = str(tmp_path / "S10.vcf")
    write_vcf(vcf_s10, "S10", [("chr1", 1500, "A", "T", "0/1")])

    manifest = str(tmp_path / "manifest.tsv")
    write_manifest(manifest, [("S10", "male", "wgs", vcf_s10, "E11.9")])

    add_samples(fresh_db, manifest, threads=1)

    db_after = Database(fresh_db)
    results_after = db_after.query(chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both")
    ac_after = results_after[0].AC if results_after else 0

    assert ac_after > ac_before


def test_add_samples_wrong_genome_build(fresh_db, tmp_path):
    vcf_s10 = str(tmp_path / "S10.vcf")
    write_vcf(vcf_s10, "S10", [("chr1", 7000, "A", "T", "0/1")])

    manifest = str(tmp_path / "manifest.tsv")
    write_manifest(manifest, [("S10", "male", "wgs", vcf_s10, "E11.9")])

    with pytest.raises(UpdateError, match="genome_build mismatch"):
        add_samples(fresh_db, manifest, threads=1, genome_build="GRCh38")


def test_add_samples_duplicate_name(fresh_db, tmp_path):
    vcf = str(tmp_path / "S00.vcf")
    write_vcf(vcf, "S00", [])  # S00 already exists in the DB

    manifest = str(tmp_path / "manifest.tsv")
    write_manifest(manifest, [("S00", "male", "wgs", vcf, "E11.9")])

    with pytest.raises(UpdateError, match="already in database"):
        add_samples(fresh_db, manifest, threads=1)


# ---------------------------------------------------------------------------
# remove_samples tests
# ---------------------------------------------------------------------------

def test_remove_samples_basic(fresh_db):
    result = remove_samples(fresh_db, ["S09"])

    assert result["removed"] == ["S09"]

    con = sqlite3.connect(os.path.join(fresh_db, "metadata.sqlite"))
    count = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    s09 = con.execute(
        "SELECT sample_id FROM samples WHERE sample_name='S09'"
    ).fetchone()
    con.close()

    assert count == 9
    assert s09 is None

    manifest_data = json.loads(open(os.path.join(fresh_db, "manifest.json")).read())
    assert manifest_data["sample_count"] == 9


def test_remove_samples_id_not_reused(fresh_db, tmp_path):
    remove_samples(fresh_db, ["S09"])

    vcf_s10 = str(tmp_path / "S10.vcf")
    write_vcf(vcf_s10, "S10", [])
    manifest = str(tmp_path / "manifest.tsv")
    write_manifest(manifest, [("S10", "male", "wgs", vcf_s10, "E11.9")])
    add_samples(fresh_db, manifest, threads=1)

    con = sqlite3.connect(os.path.join(fresh_db, "metadata.sqlite"))
    s10_id = con.execute(
        "SELECT sample_id FROM samples WHERE sample_name='S10'"
    ).fetchone()[0]
    con.close()

    assert s10_id == 10  # Not 9 (removed S09's ID must not be reused)


def test_remove_samples_query_impact(fresh_db):
    # S00 (id=0, male, wgs, E11.9+J45) is het-carrier of chr1:1500 A>T
    db_before = Database(fresh_db)
    results_before = db_before.query(chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both")
    ac_before = results_before[0].AC if results_before else 0

    remove_samples(fresh_db, ["S00"])

    db_after = Database(fresh_db)
    results_after = db_after.query(chrom="chr1", pos=1500, phenotype=["E11.9"], sex="both")
    ac_after = results_after[0].AC if results_after else 0

    assert ac_after < ac_before


def test_remove_nonexistent_raises(fresh_db):
    with pytest.raises(UpdateError, match="not found"):
        remove_samples(fresh_db, ["NONEXISTENT_SAMPLE"])


def test_remove_updates_precomputed_bitmaps(fresh_db):
    # S09 (id=9) is female; after removal, female bitmap must not contain 9
    remove_samples(fresh_db, ["S09"])

    con = sqlite3.connect(os.path.join(fresh_db, "metadata.sqlite"))
    row = con.execute(
        "SELECT bitmap_data FROM precomputed_bitmaps "
        "WHERE bitmap_type='sex' AND bitmap_key='female'"
    ).fetchone()
    con.close()

    bm = deserialize(bytes(row[0]))
    assert 9 not in bm


def test_remove_leaves_empty_bitmap_row(fresh_db):
    # chr1:3500 G>C has het=[7] — S07 is the only carrier
    # After removing S07, the row must still exist (with empty bitmaps)
    remove_samples(fresh_db, ["S07"])

    table = pq.read_table(os.path.join(fresh_db, "variants", "chr1.parquet"))
    pos_col = table["pos"].to_pylist()
    het_col = table["het_bitmap"].to_pylist()
    hom_col = table["hom_bitmap"].to_pylist()

    found = False
    for i, pos in enumerate(pos_col):
        if pos == 3500:
            found = True
            het_bm = deserialize(het_col[i])
            hom_bm = deserialize(hom_col[i])
            assert 7 not in het_bm
            assert len(het_bm) == 0
            assert len(hom_bm) == 0

    assert found, "Row at pos=3500 should still exist after removing its only carrier"


# ---------------------------------------------------------------------------
# check_database tests
# ---------------------------------------------------------------------------

def test_check_valid_db(fresh_db):
    results = check_database(fresh_db)
    errors = [r for r in results if r.severity == "error"]
    assert not errors


def test_check_missing_parquet(fresh_db):
    # Corrupt a Parquet file to make it unreadable
    pq_path = os.path.join(fresh_db, "variants", "chr1.parquet")
    with open(pq_path, "wb") as f:
        f.write(b"not a valid parquet file")

    results = check_database(fresh_db)
    errors = [r for r in results if r.severity == "error"]
    assert any("chr1" in r.message for r in errors)


def test_check_missing_capture(fresh_db):
    os.remove(os.path.join(fresh_db, "capture", "tech_0.pickle"))

    results = check_database(fresh_db)
    errors = [r for r in results if r.severity == "error"]
    assert any("tech_0" in r.message for r in errors)


def test_check_invalid_manifest_json(fresh_db):
    with open(os.path.join(fresh_db, "manifest.json"), "w") as f:
        f.write("{ invalid json !!!")

    results = check_database(fresh_db)
    errors = [r for r in results if r.severity == "error"]
    assert errors


def test_check_missing_manifest_key(fresh_db):
    manifest = json.loads(open(os.path.join(fresh_db, "manifest.json")).read())
    del manifest["genome_build"]
    with open(os.path.join(fresh_db, "manifest.json"), "w") as f:
        json.dump(manifest, f)

    results = check_database(fresh_db)
    errors = [r for r in results if r.severity == "error"]
    assert any("genome_build" in r.message for r in errors)


def test_check_after_add(fresh_db, tmp_path):
    vcf_s10 = str(tmp_path / "S10.vcf")
    write_vcf(vcf_s10, "S10", [("chr1", 7000, "A", "T", "0/1")])
    manifest = str(tmp_path / "manifest.tsv")
    write_manifest(manifest, [("S10", "male", "wgs", vcf_s10, "E11.9")])

    add_samples(fresh_db, manifest, threads=1)

    results = check_database(fresh_db)
    errors = [r for r in results if r.severity == "error"]
    assert not errors


def test_check_after_remove(fresh_db):
    remove_samples(fresh_db, ["S09"])

    results = check_database(fresh_db)
    errors = [r for r in results if r.severity == "error"]
    assert not errors


def test_check_sample_count_mismatch(fresh_db):
    # Manually set wrong sample_count in manifest
    manifest = json.loads(open(os.path.join(fresh_db, "manifest.json")).read())
    manifest["sample_count"] = 999
    with open(os.path.join(fresh_db, "manifest.json"), "w") as f:
        json.dump(manifest, f)

    results = check_database(fresh_db)
    non_info = [r for r in results if r.severity in ("warning", "error")]
    assert any("sample_count" in r.message or "mismatch" in r.message for r in non_info)


# ---------------------------------------------------------------------------
# Mixed layout: a chromosome stored in both layouts at once
# ---------------------------------------------------------------------------

def _split_chr1(db: str) -> None:
    """Put chr1 in both layouts, the state an older add-samples used to leave."""
    variants = os.path.join(db, "variants")
    os.makedirs(os.path.join(variants, "chr1"), exist_ok=True)
    shutil.copy(
        os.path.join(variants, "chr1.parquet"),
        os.path.join(variants, "chr1", "bucket_0.parquet"),
    )


def test_check_database_reports_mixed_layout(fresh_db):
    _split_chr1(fresh_db)

    errors = [r for r in check_database(fresh_db) if r.severity == "error"]
    assert any("chr1.parquet" in r.message for r in errors)


def test_add_samples_refuses_mixed_layout(fresh_db, tmp_path):
    _split_chr1(fresh_db)

    vcf = str(tmp_path / "S10.vcf")
    write_vcf(vcf, "S10", [("chr1", 7000, "A", "T", "0/1")])
    manifest = str(tmp_path / "manifest.tsv")
    write_manifest(manifest, [("S10", "male", "wgs", vcf, "E11.9")])

    with pytest.raises(UpdateError, match="Mixed variant layout"):
        add_samples(fresh_db, manifest, threads=1)


def test_add_samples_manifest_sample_count_after_removal(fresh_db, tmp_path):
    """sample_count is a count, not the next free id — they diverge after a removal."""
    remove_samples(fresh_db, ["S00"])

    vcf = str(tmp_path / "S10.vcf")
    write_vcf(vcf, "S10", [("chr1", 7000, "A", "T", "0/1")])
    manifest = str(tmp_path / "manifest.tsv")
    write_manifest(manifest, [("S10", "male", "wgs", vcf, "E11.9")])
    add_samples(fresh_db, manifest, threads=1)

    con = sqlite3.connect(os.path.join(fresh_db, "metadata.sqlite"))
    real_count = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
    con.close()

    data = json.loads(open(os.path.join(fresh_db, "manifest.json")).read())
    assert data["sample_count"] == real_count == 10
    assert data["next_sample_id"] == 11

    assert [r for r in check_database(fresh_db) if r.severity == "error"] == []

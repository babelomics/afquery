"""Query results checked against an independent reimplementation.

The oracle (tests/oracle.py) derives the expected counts from the VCF, BED and
manifest text without importing afquery. Anything the database and the oracle
disagree about is a real counting error, not a fixture that drifted.
"""
import shutil

import pytest

from afquery.database import Database
from afquery.preprocess import run_preprocess
from afquery.preprocess.build import build_all_parquets
from afquery.preprocess.update import add_samples

import oracle
from test_update import write_manifest, write_vcf


def _assert_matches(db_dir, cohort, samples=None, phenotype=None, label=""):
    db = Database(str(db_dir))
    checked = 0
    for chrom, pos, ref, alt in cohort.variants():
        want = cohort.expect(chrom, pos, ref, alt, samples)
        kwargs = {"chrom": chrom, "pos": pos, "sex": "both"}
        if phenotype is not None:
            kwargs["phenotype"] = phenotype
        got = [r for r in db.query(**kwargs)
               if (r.variant.ref, r.variant.alt) == (ref, alt)]

        if want["AN"] == 0:
            assert not got or got[0].AC == 0
            continue

        assert got, f"{label} {chrom}:{pos} {ref}>{alt} returned nothing, expected {want}"
        r = got[0]
        actual = {
            "AC": r.AC, "AN": r.AN, "N_HET": r.N_HET, "N_HOM_ALT": r.N_HOM_ALT,
            "N_FAIL": r.N_FAIL, "N_HOM_REF": r.N_HOM_REF,
        }
        assert actual == want, f"{label} {chrom}:{pos} {ref}>{alt}"
        checked += 1
    assert checked > 0, f"{label}: oracle checked nothing"


# ---------------------------------------------------------------------------
# Fixture cohort, both layouts
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def fixture_cohort(data_dir):
    return oracle.Cohort(data_dir / "manifest.tsv", bed_dir=data_dir / "beds")


@pytest.fixture(scope="module")
def built_db(tmp_path_factory, data_dir):
    db = tmp_path_factory.mktemp("oracle_db")
    run_preprocess(
        manifest_path=str(data_dir / "manifest.tsv"), output_dir=str(db),
        genome_build="GRCh37", bed_dir=str(data_dir / "beds"), threads=2,
    )
    return db


def test_partitioned_build_matches_oracle(built_db, fixture_cohort):
    _assert_matches(built_db, fixture_cohort, label="partitioned")


def test_flat_build_matches_oracle(tmp_path, built_db, fixture_cohort, data_dir):
    """The legacy layout must agree with the oracle too, and so with the buckets."""
    import sqlite3

    from afquery.models import Sample
    from afquery.preprocess.ingest import ingest_all

    flat = tmp_path / "flat"
    shutil.copytree(built_db, flat)
    shutil.rmtree(flat / "variants")
    (flat / "variants").mkdir()

    con = sqlite3.connect(flat / "metadata.sqlite")
    rows = con.execute("SELECT sample_id, sample_name, sex, tech_id FROM samples").fetchall()
    con.close()
    samples = [Sample(*r) for r in rows]
    vcfs = [str(data_dir / "vcfs" / f"{r[1]}.vcf") for r in rows]

    tmp = tmp_path / "ingest"
    tmp.mkdir()
    ingest_all(samples, vcfs, str(tmp), n_workers=1)
    build_all_parquets(str(tmp), str(flat / "variants"), n_workers=1, partitioned=False)

    _assert_matches(flat, fixture_cohort, label="flat")


def test_phenotype_subset_matches_oracle(built_db, fixture_cohort):
    subset = fixture_cohort.with_phenotype("E11.9")
    assert subset
    _assert_matches(built_db, fixture_cohort, samples=subset,
                    phenotype=["E11.9"], label="phenotype E11.9")


# ---------------------------------------------------------------------------
# After add-samples — the case the bug lived in
# ---------------------------------------------------------------------------

@pytest.fixture
def grown_db(built_db, tmp_path, data_dir):
    """Fixture cohort plus one wgs and one wes_kit_a sample, added via update."""
    db = tmp_path / "grown"
    shutil.copytree(built_db, db)

    entries = []
    for name, tech, sex, variants in [
        ("S10", "wgs", "male", [("chr1", 1500, "A", "T", "1/1"),
                                ("chr1", 7000, "A", "T", "0/1"),
                                ("chrX", 5000000, "A", "G", "0/1")]),
        ("S11", "wes_kit_a", "female", [("chr1", 1500, "A", "T", "0/1"),
                                        ("chr1", 1800, "C", "G", "1/1")]),
    ]:
        vcf = tmp_path / f"{name}.vcf"
        write_vcf(str(vcf), name, variants)
        entries.append((name, sex, tech, str(vcf), "ZZNEW"))

    manifest = tmp_path / "grow.tsv"
    write_manifest(str(manifest), entries)
    add_samples(str(db), str(manifest), threads=1, bed_dir=str(data_dir / "beds"))
    return db, manifest


def test_after_add_samples_matches_oracle(grown_db, data_dir, tmp_path):
    db, grow_manifest = grown_db

    # Oracle over the union of the original manifest and the added one
    combined = tmp_path / "combined.tsv"
    original = (data_dir / "manifest.tsv").read_text().rstrip("\n").splitlines()
    added = grow_manifest.read_text().rstrip("\n").splitlines()[1:]
    header = original[0].split("\t")

    out = [original[0]]
    for line in original[1:]:
        row = dict(zip(header, line.split("\t")))
        # the fixture manifest's paths are relative to tests/data
        row["vcf_path"] = str(data_dir / row["vcf_path"])
        out.append("\t".join(row[c] for c in header))
    for line in added:
        # write_manifest emits a different column order
        row = dict(zip(["sample_name", "sex", "tech_name", "vcf_path", "phenotype_codes"],
                       line.split("\t")))
        out.append("\t".join(row[c] for c in header))
    combined.write_text("\n".join(out) + "\n")

    cohort = oracle.Cohort(combined, bed_dir=data_dir / "beds")
    for name in ("S10", "S11"):
        assert name in cohort.samples

    _assert_matches(db, cohort, label="after add")


def test_after_add_samples_new_phenotype_matches_oracle(grown_db, data_dir, tmp_path):
    """Restricted to the added samples, the numbers must describe only them."""
    db, grow_manifest = grown_db
    cohort = oracle.Cohort(grow_manifest, bed_dir=data_dir / "beds")
    _assert_matches(db, cohort, samples=["S10", "S11"],
                    phenotype=["ZZNEW"], label="added-only")

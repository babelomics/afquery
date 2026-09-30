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


# ---------------------------------------------------------------------------
# Multi-allelic cohort
# ---------------------------------------------------------------------------

def test_multiallelic_cohort_matches_oracle(tmp_path):
    """Two ALT alleles at one position, a 1/2 carrier, partial capture and a
    FILTER failure: every allele's tallies agree with the oracle."""
    calls = {
        # name: (sex, tech, [(chrom, pos, ref, alt, gt, filter)])
        "A_HET":  ("female", "WGS",   [("chr1", 5000, "G", "A", "0/1", "PASS")]),
        "A_HOM":  ("male",   "WGS",   [("chr1", 5000, "G", "A", "1/1", "PASS")]),
        "T_HET":  ("female", "WGS",   [("chr1", 5000, "G", "T", "0/1", "PASS")]),
        "AT_HET": ("male",   "WGS",   [("chr1", 5000, "G", "A,T", "1/2", "PASS")]),
        "T_FAIL": ("female", "WGS",   [("chr1", 5000, "G", "T", "0/1", "LowQual")]),
        "REF":    ("female", "WGS",   [("chr1", 100, "C", "G", "0/1", "PASS")]),
        "P_T":    ("female", "PANEL", [("chr1", 5000, "G", "T", "1/1", "PASS")]),
        "P_REF":  ("male",   "PANEL", [("chr1", 100, "C", "G", "0/1", "PASS")]),
        "P_OFF":  ("female", "PANEL", [("chr1", 9000, "C", "G", "0/1", "PASS"),
                                       ("chr1", 9000, "C", "T", "0/1", "PASS")]),
        "X_A":    ("male",   "WGS",   [("chrX", 5000000, "A", "G,C", "1", "PASS")]),
        "X_C":    ("female", "WGS",   [("chrX", 5000000, "A", "C", "0/1", "PASS")]),
    }
    beds = tmp_path / "beds"
    beds.mkdir()
    (beds / "PANEL.bed").write_text("chr1\t0\t6000\n")
    manifest = ["sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes"]
    for name, (sex, tech, records) in calls.items():
        vcf = tmp_path / f"{name}.vcf"
        with open(vcf, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write('##FILTER=<ID=LowQual,Description="Low quality">\n')
            for contig in sorted({r[0] for r in records}):
                f.write(f"##contig=<ID={contig}>\n")
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}\n")
            for chrom, pos, ref, alt, gt, flt in records:
                f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t.\tGT\t{gt}\n")
        manifest.append(f"{name}\t{sex}\t{tech}\t{vcf}\tE11.9")
    manifest_path = tmp_path / "manifest.tsv"
    manifest_path.write_text("\n".join(manifest) + "\n")

    db = tmp_path / "db"
    run_preprocess(
        manifest_path=str(manifest_path), output_dir=str(db),
        genome_build="GRCh37", bed_dir=str(beds), threads=1,
    )
    cohort = oracle.Cohort(manifest_path, bed_dir=beds)
    _assert_matches(db, cohort, label="multiallelic")

    # Spot-check chr1:5000, where all 11 samples are eligible (the panel BED
    # covers it). Hom-ref for both alleles: REF, P_REF, P_OFF, X_A, X_C.
    got = {r.variant.alt: r for r in Database(str(db)).query("chr1", 5000)}
    assert got["A"].n_samples_eligible == 11
    # A: A_HET, AT_HET het; A_HOM hom; T_HET, T_FAIL, P_T carry only T.
    assert (got["A"].N_HET, got["A"].N_HOM_ALT, got["A"].N_HOM_REF) == (2, 1, 5)
    # T: T_HET, AT_HET het; P_T hom; T_FAIL failed; A_HET, A_HOM carry only A.
    assert (got["T"].N_HET, got["T"].N_HOM_ALT, got["T"].N_FAIL, got["T"].N_HOM_REF) == (2, 1, 1, 5)

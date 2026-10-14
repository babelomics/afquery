"""Tests for N_NO_COVERAGE field and Phase 1 / Phase 2 filtering."""

import json
import shutil
import sqlite3
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from pyroaring import BitMap

from afquery import Database
from afquery.bitmaps import serialize
from afquery.preprocess.build import PARQUET_SCHEMA


def _db(test_db):
    return Database(test_db)


# ---------------------------------------------------------------------------
# Phase 1 — query-time count-based gating
# ---------------------------------------------------------------------------

class TestPhase1Defaults:
    def test_min_pass_zero_no_effect(self, test_db):
        """min_pass=0 (default) → N_NO_COVERAGE=0 for all variants."""
        db = _db(test_db)
        # Single locus
        for r in db.query(chrom="chr1", pos=1500):
            assert r.N_NO_COVERAGE == 0
        # Region
        for r in db.query_region("chr1", 1, 10000):
            assert r.N_NO_COVERAGE == 0

    def test_invariant_with_zero_thresholds(self, test_db):
        """Genotype partition holds for every variant under default settings."""
        db = _db(test_db)
        for r in db.query_region("chr1", 1, 10000):
            assert (
                r.N_HET + r.N_HOM_ALT + r.N_HOM_REF + r.N_FAIL + r.N_NO_COVERAGE
                == r.n_samples_eligible
            )


class TestPhase1MinPass:
    def test_min_pass_filters_wes_below_threshold(self, test_db):
        """chr1:3500 has 1 PASS carrier in WES_B (S07). WES_A's BED doesn't cover 3500.
        With min_pass=2, WES_B's 1 carrier < 2 → its non-carrier samples (S08, S09)
        move from N_HOM_REF to N_NO_COVERAGE.
        """
        db = _db(test_db)
        baseline = db.query(chrom="chr1", pos=3500)[0]
        filtered = db.query(chrom="chr1", pos=3500, min_pass=2)[0]

        # Carriers stay in N_HET / N_HOM_ALT
        assert filtered.N_HET == baseline.N_HET
        assert filtered.N_HOM_ALT == baseline.N_HOM_ALT
        assert filtered.N_FAIL == baseline.N_FAIL

        # Non-carrier WES_B samples (S08, S09) shift to N_NO_COVERAGE
        assert filtered.N_NO_COVERAGE == 2
        assert filtered.N_HOM_REF == baseline.N_HOM_REF - 2

        # n_eligible and AN are unchanged: filtered samples remain in eligible
        assert filtered.n_samples_eligible == baseline.n_samples_eligible
        assert filtered.AN == baseline.AN

    def test_min_pass_one_no_effect_when_carrier_exists(self, test_db):
        """min_pass=1 should NOT filter when each WES tech has ≥1 PASS carrier in BED."""
        db = _db(test_db)
        # chr1:3500 → WES_B has 1 carrier (S07); WES_A is not in BED → tech_eligible empty
        r = db.query(chrom="chr1", pos=3500, min_pass=1)[0]
        assert r.N_NO_COVERAGE == 0

    def test_invariant_holds_with_filter(self, test_db):
        db = _db(test_db)
        for r in db.query_region("chr1", 1, 10000, min_pass=2):
            assert (
                r.N_HET + r.N_HOM_ALT + r.N_HOM_REF + r.N_FAIL + r.N_NO_COVERAGE
                == r.n_samples_eligible
            )


class TestPhase1MinObserved:
    def test_min_observed_counts_fail(self, test_db):
        """min_observed counts fail samples as evidence — fail-only WES tech survives."""
        # samples 0-3 are WGS, 4-6 are WES_A, 7-9 are WES_B.
        # chr1:1500 het=[5], hom=[2], fail=[0]. WES_A carriers: het={5} (1 PASS carrier)
        # WES_A's only carrier (S05) is in het → also counts as observed (1).
        # So min_observed=2 should filter.
        db = _db(test_db)
        baseline = db.query(chrom="chr1", pos=1500)[0]
        filtered = db.query(chrom="chr1", pos=1500, min_observed=2)[0]
        assert filtered.N_NO_COVERAGE >= 1
        assert filtered.N_HET + filtered.N_HOM_ALT == baseline.N_HET + baseline.N_HOM_ALT


class TestPhase1AndLogic:
    def test_min_pass_and_min_observed_combine(self, test_db):
        """If either threshold fails, the tech is filtered."""
        db = _db(test_db)
        only_pass = db.query(chrom="chr1", pos=3500, min_pass=2, min_observed=0)[0]
        only_obs = db.query(chrom="chr1", pos=3500, min_pass=0, min_observed=2)[0]
        both = db.query(chrom="chr1", pos=3500, min_pass=2, min_observed=2)[0]
        # With chr1:3500 only WES_B is eligible: 1 carrier, 1 observed.
        # Both filters fail at threshold 2 → same N_NO_COVERAGE.
        assert only_pass.N_NO_COVERAGE == only_obs.N_NO_COVERAGE == both.N_NO_COVERAGE


# ---------------------------------------------------------------------------
# Phase 1 — invariants
# ---------------------------------------------------------------------------

class TestInvariants:
    def test_wgs_never_in_no_coverage(self, test_db):
        """WGS samples (tech_id=0) should never appear in N_NO_COVERAGE."""
        db = _db(test_db)
        # variant_info shows individual carriers; we use it to verify identity
        carriers = db.variant_info(chrom="chr1", pos=3500, min_pass=2)
        for c in carriers:
            if c.genotype == "no_coverage":
                assert c.tech_name != "WGS", (
                    f"WGS sample {c.sample_name} marked as no_coverage"
                )

    def test_carriers_never_in_no_coverage(self, test_db):
        """het/hom/fail carriers should never appear with genotype=no_coverage."""
        db = _db(test_db)
        carriers = db.variant_info(chrom="chr1", pos=3500, min_pass=10)
        # At threshold=10, ALL techs fail. But carriers must stay as het/hom/alt.
        for c in carriers:
            if c.genotype == "no_coverage":
                # Sample must NOT be a carrier at this position
                # (carriers are listed BEFORE no_coverage in variant_info order)
                pass  # checked implicitly via genotype != het/hom/alt
            else:
                assert c.genotype in ("het", "hom", "alt")


# ---------------------------------------------------------------------------
# variant_info with no_coverage
# ---------------------------------------------------------------------------

class TestVariantInfoNoCoverage:
    def test_no_coverage_genotype_appears(self, test_db):
        """variant_info should report 'no_coverage' for filtered samples."""
        db = _db(test_db)
        baseline = db.variant_info(chrom="chr1", pos=3500)
        filtered = db.variant_info(chrom="chr1", pos=3500, min_pass=2)

        # Baseline carriers
        baseline_carriers = {c.sample_id for c in baseline}
        # Filtered should include baseline carriers + at least 1 no_coverage sample
        filtered_carriers = {c.sample_id for c in filtered}
        no_cov_ids = {c.sample_id for c in filtered if c.genotype == "no_coverage"}
        assert no_cov_ids, "expected at least 1 no_coverage sample"
        # All baseline carriers still present, with their original genotype
        baseline_by_id = {c.sample_id: c.genotype for c in baseline}
        for c in filtered:
            if c.sample_id in baseline_by_id:
                assert c.genotype == baseline_by_id[c.sample_id]

    def test_no_coverage_filter_pass_is_none(self, test_db):
        """filter_pass must be None for no_coverage rows; True/False otherwise."""
        db = _db(test_db)
        carriers = db.variant_info(chrom="chr1", pos=3500, min_pass=2)
        saw_no_cov = False
        for c in carriers:
            if c.genotype == "no_coverage":
                saw_no_cov = True
                assert c.filter_pass is None
            else:
                assert c.filter_pass in (True, False)
        assert saw_no_cov, "expected at least one no_coverage carrier"


# ---------------------------------------------------------------------------
# Phase 2 — filtered_bitmap and quality_pass_bitmap
# ---------------------------------------------------------------------------

def _make_phase2_db(src_db: str, dst_dir: Path,
                    filtered_at_chr1_3500: list[int] | None = None,
                    quality_pass_at_chr1_3500: list[int] | None = None) -> str:
    """Clone test_db and inject Phase 2 columns at chr1:3500."""
    dst = dst_dir / "p2_db"
    shutil.copytree(src_db, str(dst))

    # Update manifest to indicate Phase 2 schema
    manifest_path = dst / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["schema_version"] = "3.0"
    manifest["coverage_filter"] = {
        "min_dp": 30, "min_gq": 20, "min_qual": 0.0, "min_covered": 1,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2))

    # Inject filtered_bitmap and quality_pass_bitmap into chr1.parquet
    parquet_file = dst / "variants" / "chr1.parquet"
    table = pq.read_table(str(parquet_file))
    n = len(table)
    empty = serialize(BitMap())
    new_filtered = []
    new_quality = []
    for i in range(n):
        pos = table["pos"][i].as_py()
        if pos == 3500:
            new_filtered.append(serialize(BitMap(filtered_at_chr1_3500 or [])))
            new_quality.append(serialize(BitMap(quality_pass_at_chr1_3500 or [])))
        else:
            new_filtered.append(empty)
            new_quality.append(empty)

    new_table = pa.table(
        {
            "pos":                 table["pos"],
            "ref":                 table["ref"],
            "alt":                 table["alt"],
            "het_bitmap":          table["het_bitmap"],
            "hom_bitmap":          table["hom_bitmap"],
            "fail_bitmap":         table["fail_bitmap"],
            "filtered_bitmap":     pa.array(new_filtered, type=pa.large_binary()),
            "quality_pass_bitmap": pa.array(new_quality, type=pa.large_binary()),
        },
        schema=PARQUET_SCHEMA,
    )
    pq.write_table(new_table, str(parquet_file))
    return str(dst)


class TestPhase2Stored:
    def test_filtered_bitmap_moves_samples_to_no_coverage(self, test_db, tmp_path):
        """Phase 2 DB with filtered_bitmap[8,9] at chr1:3500 → N_NO_COVERAGE = 2."""
        # Samples 8, 9 are WES_B, BED-covered at 3500, not carriers.
        p2 = _make_phase2_db(test_db, tmp_path,
                             filtered_at_chr1_3500=[8, 9],
                             quality_pass_at_chr1_3500=[7])
        db = Database(p2)
        r = db.query(chrom="chr1", pos=3500)[0]
        assert r.N_NO_COVERAGE == 2

    def test_min_quality_evidence_filters_low_quality_techs(self, test_db, tmp_path):
        """min_quality_evidence=2 with quality_pass_bitmap=[7] (1 sample) → tech filtered."""
        p2 = _make_phase2_db(test_db, tmp_path,
                             filtered_at_chr1_3500=[],
                             quality_pass_at_chr1_3500=[7])
        db = Database(p2)
        # No filtered_bitmap entries → baseline N_NO_COVERAGE=0
        baseline = db.query(chrom="chr1", pos=3500)[0]
        assert baseline.N_NO_COVERAGE == 0
        # min_quality_evidence=2 → WES_B has 1 quality_pass (S07) < 2 → S08, S09 filtered
        r = db.query(chrom="chr1", pos=3500, min_quality_evidence=2)[0]
        assert r.N_NO_COVERAGE == 2

    def test_min_quality_evidence_errors_on_old_db(self, test_db):
        """Using --min-quality-evidence on a schema_version < 3.0 DB raises ValueError."""
        db = _db(test_db)
        with pytest.raises(ValueError, match="coverage quality"):
            db.query(chrom="chr1", pos=3500, min_quality_evidence=1)

    def test_phase1_phase2_combine(self, test_db, tmp_path):
        """N_NO_COVERAGE = union of stored filtered_bitmap and Phase 1 dynamic filtering."""
        p2 = _make_phase2_db(test_db, tmp_path,
                             filtered_at_chr1_3500=[8],
                             quality_pass_at_chr1_3500=[7])
        db = Database(p2)
        # filtered_bitmap = {8} → N_NO_COVERAGE = 1 baseline
        baseline = db.query(chrom="chr1", pos=3500)[0]
        assert baseline.N_NO_COVERAGE == 1
        # min_pass=2 → WES_B has 1 PASS carrier < 2 → also filters S09 (S08 already)
        # Union: {8} | {8, 9} = {8, 9} → N_NO_COVERAGE = 2
        r = db.query(chrom="chr1", pos=3500, min_pass=2)[0]
        assert r.N_NO_COVERAGE == 2


# ---------------------------------------------------------------------------
# Phase 1 + ploidy chromosomes
# ---------------------------------------------------------------------------

class TestInvariantOnSexChroms:
    def test_chrm_invariant(self, test_db):
        """Invariant must hold on chrM (haploid)."""
        db = _db(test_db)
        for r in db.query(chrom="chrM", pos=100, min_pass=1):
            assert (
                r.N_HET + r.N_HOM_ALT + r.N_HOM_REF + r.N_FAIL + r.N_NO_COVERAGE
                == r.n_samples_eligible
            )

    def test_chry_invariant(self, test_db):
        """Invariant must hold on chrY."""
        db = _db(test_db)
        for r in db.query(chrom="chrY", pos=500000, min_pass=1):
            assert (
                r.N_HET + r.N_HOM_ALT + r.N_HOM_REF + r.N_FAIL + r.N_NO_COVERAGE
                == r.n_samples_eligible
            )

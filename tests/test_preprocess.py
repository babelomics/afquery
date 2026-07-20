import json
import logging
import os
import pickle
import sqlite3
import time as _time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from pyroaring import BitMap

from afquery.bitmaps import deserialize
from afquery.capture import CaptureIndex
from afquery.database import Database
from afquery.models import Sample, Technology
from afquery.preprocess import run_preprocess, _write_sqlite, _write_manifest
from afquery.preprocess.build import (
    _build_chrom_worker,
    _get_chroms_from_consolidated,
    _get_chroms_from_file,
    build_all_parquets,
    build_chromosome_parquet,
    consolidate_temp_files,
    get_chroms_in_temp_files,
)
from afquery.preprocess.ingest import INGEST_SCHEMA, IngestError, ingest_sample
from afquery.preprocess.manifest import ManifestError, parse_manifest


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_vcf(path: Path, sample_name: str, records: list[tuple]) -> None:
    """Write a minimal VCF. records = [(chrom, pos, ref, alt, gt_str), ...]"""
    contigs = sorted({r[0] for r in records}) if records else ["chr1"]
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        for contig in contigs:
            f.write(f"##contig=<ID={contig}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        )
        for chrom, pos, ref, alt, gt in records:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")


def _write_vcf_with_filter(path: Path, sample_name: str, records: list[tuple]) -> None:
    """Write a minimal VCF. records = [(chrom, pos, ref, alt, gt_str, filter_str), ...]"""
    contigs = sorted({r[0] for r in records}) if records else ["chr1"]
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##FILTER=<ID=LowQual,Description="Low quality">\n')
        for contig in contigs:
            f.write(f"##contig=<ID={contig}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        )
        for chrom, pos, ref, alt, gt, flt in records:
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t.\tGT\t{gt}\n")


def _write_parquet(path: Path, rows: list[tuple]) -> None:
    """Write ingest-schema Parquet. rows = [(chrom, pos, ref, alt, gt_ac, sample_id[, filter_pass]), ...]
    filter_pass defaults to True if not provided. Quality columns (dp/gq/qual) default to None."""
    n = len(rows)
    table = pa.table(
        {
            "chrom":       pa.array([r[0] for r in rows], type=pa.utf8()),
            "pos":         pa.array([r[1] for r in rows], type=pa.uint32()),
            "ref":         pa.array([r[2] for r in rows], type=pa.utf8()),
            "alt":         pa.array([r[3] for r in rows], type=pa.utf8()),
            "gt_ac":       pa.array([r[4] for r in rows], type=pa.uint8()),
            "sample_id":   pa.array([r[5] for r in rows], type=pa.uint32()),
            "filter_pass": pa.array([r[6] if len(r) > 6 else True for r in rows], type=pa.bool_()),
            "dp":          pa.array([None] * n, type=pa.int32()),
            "gq":          pa.array([None] * n, type=pa.int32()),
            "qual":        pa.array([None] * n, type=pa.float32()),
        },
        schema=INGEST_SCHEMA,
    )
    pq.write_table(table, path)


# ---------------------------------------------------------------------------
# Manifest unit tests
# ---------------------------------------------------------------------------

class TestManifest:
    def _make_manifest(self, tmp_path, content: str, vcfs: list[str] | None = None) -> Path:
        if vcfs is None:
            vcfs = []
        for v in vcfs:
            vcf_path = tmp_path / v
            vcf_path.parent.mkdir(parents=True, exist_ok=True)
            _write_vcf(vcf_path, "SAMPLE", [])
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text(content)
        return manifest

    def test_valid_wgs_manifest(self, tmp_path):
        manifest = self._make_manifest(
            tmp_path,
            "sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n"
            "S00\tmale\twgs\tS00.vcf\tE11.9,J45\n"
            "S01\tfemale\twgs\tS01.vcf\tI10\n",
            vcfs=["S00.vcf", "S01.vcf"],
        )
        samples, techs = parse_manifest(str(manifest))
        assert len(samples) == 2
        assert samples[0].sample_name == "S00"
        assert samples[0].sex == "male"
        assert samples[0].tech_name == "wgs"
        assert samples[0].phenotype_codes == ["E11.9", "J45"]
        assert len(techs) == 1
        assert techs[0].tech_name == "wgs"
        assert techs[0].bed_path is None

    def test_wes_tech_with_bed(self, tmp_path):
        bed = tmp_path / "wes_kit_a.bed"
        bed.write_text("chr1\t999\t2000\n")
        manifest = self._make_manifest(
            tmp_path,
            "sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n"
            "S00\tmale\twes_kit_a\tS00.vcf\tE11.9\n",
            vcfs=["S00.vcf"],
        )
        samples, techs = parse_manifest(str(manifest), bed_dir=str(tmp_path))
        assert techs[0].bed_path == str(tmp_path / "wes_kit_a.bed")

    def test_missing_required_column_raises(self, tmp_path):
        manifest = self._make_manifest(
            tmp_path,
            "sample_name\tsex\tvcf_path\tphenotype_codes\n"  # missing tech_name
            "S00\tmale\tS00.vcf\tE11.9\n",
            vcfs=["S00.vcf"],
        )
        with pytest.raises(ManifestError, match="Missing required columns"):
            parse_manifest(str(manifest))

    def test_duplicate_sample_name_raises(self, tmp_path):
        manifest = self._make_manifest(
            tmp_path,
            "sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n"
            "S00\tmale\twgs\tS00.vcf\tE11.9\n"
            "S00\tfemale\twgs\tS01.vcf\tI10\n",
            vcfs=["S00.vcf", "S01.vcf"],
        )
        with pytest.raises(ManifestError, match="Duplicate sample_name"):
            parse_manifest(str(manifest))

    def test_invalid_sex_raises(self, tmp_path):
        manifest = self._make_manifest(
            tmp_path,
            "sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n"
            "S00\tunknown\twgs\tS00.vcf\tE11.9\n",
            vcfs=["S00.vcf"],
        )
        with pytest.raises(ManifestError, match="invalid sex"):
            parse_manifest(str(manifest))

    def test_nonexistent_vcf_raises(self, tmp_path):
        manifest = self._make_manifest(
            tmp_path,
            "sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n"
            "S00\tmale\twgs\tmissing.vcf\tE11.9\n",
        )
        with pytest.raises(ManifestError, match="vcf_path does not exist"):
            parse_manifest(str(manifest))

    def test_empty_phenotype_raises(self, tmp_path):
        manifest = self._make_manifest(
            tmp_path,
            "sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n"
            "S00\tmale\twgs\tS00.vcf\t\n",
            vcfs=["S00.vcf"],
        )
        with pytest.raises(ManifestError, match="phenotype_codes is empty"):
            parse_manifest(str(manifest))

    def test_whitespace_stripping(self, tmp_path):
        manifest = self._make_manifest(
            tmp_path,
            "sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n"
            "S00 \t male \t wgs \t S00.vcf \t E11.9 , J45 \n",
            vcfs=["S00.vcf"],
        )
        samples, _ = parse_manifest(str(manifest))
        assert samples[0].sample_name == "S00"
        assert samples[0].sex == "male"
        assert samples[0].phenotype_codes == ["E11.9", "J45"]

    def test_wes_missing_bed_raises(self, tmp_path):
        manifest = self._make_manifest(
            tmp_path,
            "sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n"
            "S00\tmale\twes_kit_x\tS00.vcf\tE11.9\n",
            vcfs=["S00.vcf"],
        )
        with pytest.raises(ManifestError, match="BED file not found"):
            parse_manifest(str(manifest), bed_dir=str(tmp_path))

    def test_relative_vcf_path_resolved(self, tmp_path):
        vcf = tmp_path / "sub" / "sample.vcf"
        vcf.parent.mkdir()
        _write_vcf(vcf, "S00", [])
        manifest = tmp_path / "manifest.tsv"
        manifest.write_text(
            "sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n"
            "S00\tmale\twgs\tsub/sample.vcf\tE11.9\n"
        )
        samples, _ = parse_manifest(str(manifest))
        assert os.path.isabs(samples[0].vcf_path)
        assert os.path.exists(samples[0].vcf_path)

    def test_technology_deduplication(self, tmp_path):
        manifest = self._make_manifest(
            tmp_path,
            "sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes\n"
            "S00\tmale\twgs\tS00.vcf\tE11.9\n"
            "S01\tfemale\twgs\tS01.vcf\tI10\n"
            "S02\tmale\twgs\tS02.vcf\tJ45\n",
            vcfs=["S00.vcf", "S01.vcf", "S02.vcf"],
        )
        _, techs = parse_manifest(str(manifest))
        assert len(techs) == 1


# ---------------------------------------------------------------------------
# Regions unit tests
# ---------------------------------------------------------------------------

class TestRegions:
    def test_wgs_tech_covers_all(self, tmp_path, data_dir):
        from afquery.models import Technology
        from afquery.preprocess.regions import build_capture_indices

        tech = Technology(tech_id=0, tech_name="wgs", bed_path=None)
        build_capture_indices([tech], str(tmp_path))

        idx = CaptureIndex.load(str(tmp_path / "tech_0.pickle"))
        assert idx.covers("chr1", 1)
        assert idx.covers("chrX", 5_000_000)
        assert idx.covers("chrY", 999_999_999)

    def test_wes_tech_covers_bed_only(self, tmp_path, data_dir):
        from afquery.models import Technology
        from afquery.preprocess.regions import build_capture_indices

        bed_path = str(data_dir / "beds" / "wes_kit_a.bed")
        tech = Technology(tech_id=1, tech_name="wes_kit_a", bed_path=bed_path)
        build_capture_indices([tech], str(tmp_path))

        idx = CaptureIndex.load(str(tmp_path / "tech_1.pickle"))
        assert idx.covers("chr1", 1500)       # chr1:999-2000 covers pos 1500
        assert not idx.covers("chr1", 5000)   # outside bed
        assert not idx.covers("chrY", 500_000)


# ---------------------------------------------------------------------------
# Ingest unit tests
# ---------------------------------------------------------------------------

class TestIngest:
    def test_het_variant(self, tmp_path):
        vcf = tmp_path / "s0.vcf"
        _write_vcf(vcf, "S0", [("chr1", 1000, "A", "T", "0/1")])
        out, _ = ingest_sample(0, str(vcf), str(tmp_path))
        tbl = pq.read_table(out)
        assert len(tbl) == 1
        assert tbl["gt_ac"][0].as_py() == 1
        assert tbl["chrom"][0].as_py() == "chr1"

    def test_hom_alt_variant(self, tmp_path):
        vcf = tmp_path / "s0.vcf"
        _write_vcf(vcf, "S0", [("chr1", 1000, "A", "T", "1/1")])
        out, _ = ingest_sample(0, str(vcf), str(tmp_path))
        tbl = pq.read_table(out)
        assert len(tbl) == 1
        assert tbl["gt_ac"][0].as_py() == 2

    def test_hom_ref_produces_no_rows(self, tmp_path):
        vcf = tmp_path / "s0.vcf"
        _write_vcf(vcf, "S0", [("chr1", 1000, "A", "T", "0/0")])
        out, _ = ingest_sample(0, str(vcf), str(tmp_path))
        tbl = pq.read_table(out)
        assert len(tbl) == 0

    def test_missing_genotype_produces_no_rows(self, tmp_path):
        vcf = tmp_path / "s0.vcf"
        _write_vcf(vcf, "S0", [("chr1", 1000, "A", "T", "./.")])
        out, _ = ingest_sample(0, str(vcf), str(tmp_path))
        tbl = pq.read_table(out)
        assert len(tbl) == 0

    def test_multiallelic_produces_two_rows(self, tmp_path):
        vcf = tmp_path / "s0.vcf"
        _write_vcf(vcf, "S0", [("chr1", 1000, "A", "T,G", "1/2")])
        out, _ = ingest_sample(0, str(vcf), str(tmp_path))
        tbl = pq.read_table(out)
        assert len(tbl) == 2
        alts = sorted(tbl["alt"].to_pylist())
        assert alts == ["G", "T"]

    def test_chrM_het(self, tmp_path):
        vcf = tmp_path / "s0.vcf"
        _write_vcf(vcf, "S0", [("chrM", 100, "C", "A", "0/1")])
        out, _ = ingest_sample(0, str(vcf), str(tmp_path))
        tbl = pq.read_table(out)
        assert len(tbl) == 1
        assert tbl["chrom"][0].as_py() == "chrM"
        assert tbl["gt_ac"][0].as_py() == 1

    def test_missing_gt_pass_produces_no_rows(self, tmp_path):
        """GT=./. with FILTER=PASS → still no row (no allele information)."""
        vcf = tmp_path / "s0.vcf"
        _write_vcf_with_filter(vcf, "S0", [("chr1", 1000, "A", "T", "./.", "PASS")])
        out, _ = ingest_sample(0, str(vcf), str(tmp_path))
        tbl = pq.read_table(out)
        assert len(tbl) == 0

    def test_missing_gt_failed_filter_produces_fail_row(self, tmp_path):
        """GT=./. with FILTER=LowQual → one row with gt_ac=0, filter_pass=False."""
        vcf = tmp_path / "s0.vcf"
        _write_vcf_with_filter(vcf, "S0", [("chr1", 1000, "A", "T", "./.", "LowQual")])
        out, _ = ingest_sample(0, str(vcf), str(tmp_path))
        tbl = pq.read_table(out)
        assert len(tbl) == 1
        assert tbl["gt_ac"][0].as_py() == 0
        assert tbl["filter_pass"][0].as_py() == False
        assert tbl["chrom"][0].as_py() == "chr1"
        assert tbl["pos"][0].as_py() == 1000
        assert tbl["alt"][0].as_py() == "T"

    def test_missing_gt_failed_multiallelic_produces_rows_per_alt(self, tmp_path):
        """GT=./. with FILTER=LowQual at multiallelic site → one row per ALT."""
        vcf = tmp_path / "s0.vcf"
        _write_vcf_with_filter(vcf, "S0", [("chr1", 1000, "A", "T,G", "./.", "LowQual")])
        out, _ = ingest_sample(0, str(vcf), str(tmp_path))
        tbl = pq.read_table(out)
        assert len(tbl) == 2
        assert all(tbl["gt_ac"].to_pylist()[i] == 0 for i in range(2))
        assert all(not tbl["filter_pass"].to_pylist()[i] for i in range(2))
        assert sorted(tbl["alt"].to_pylist()) == ["G", "T"]

    def test_chrom_normalization(self, tmp_path):
        vcf = tmp_path / "s0.vcf"
        _write_vcf(vcf, "S0", [("1", 1000, "A", "T", "0/1")])
        out, _ = ingest_sample(0, str(vcf), str(tmp_path))
        tbl = pq.read_table(out)
        assert tbl["chrom"][0].as_py() == "chr1"


# ---------------------------------------------------------------------------
# Build unit tests
# ---------------------------------------------------------------------------

class TestBuild:
    def test_two_samples_correct_bitmaps(self, tmp_path):
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),  # het, pass
        ])
        _write_parquet(tmp_path / "sample_1.parquet", [
            ("chr1", 1000, "A", "T", 2, 1),  # hom, pass
        ])
        count = build_chromosome_parquet("chr1", str(tmp_path), str(variants_dir))
        assert count == 1

        tbl = pq.read_table(variants_dir / "chr1.parquet")
        het_bm = deserialize(bytes(tbl["het_bitmap"][0].as_py()))
        hom_bm = deserialize(bytes(tbl["hom_bitmap"][0].as_py()))
        fail_bm = deserialize(bytes(tbl["fail_bitmap"][0].as_py()))
        assert list(het_bm) == [0]
        assert list(hom_bm) == [1]
        assert list(fail_bm) == []

    def test_output_sorted_by_pos_alt(self, tmp_path):
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 2000, "G", "C", 1, 0),
            ("chr1", 1000, "A", "T", 1, 0),
            ("chr1", 1000, "A", "G", 1, 0),
        ])
        build_chromosome_parquet("chr1", str(tmp_path), str(variants_dir))
        tbl = pq.read_table(variants_dir / "chr1.parquet")
        positions = tbl["pos"].to_pylist()
        alts = tbl["alt"].to_pylist()
        assert positions == [1000, 1000, 2000]
        assert alts[:2] == ["G", "T"]  # sorted by alt at same pos

    def test_multi_chromosome_separate_files(self, tmp_path):
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1",  1000, "A", "T", 1, 0),
            ("chrX",  5000, "G", "C", 1, 0),
        ])
        result = build_all_parquets(str(tmp_path), str(variants_dir))
        assert "chr1" in result
        assert "chrX" in result
        # Default is partitioned=True: bucket files in chrom subdirectories
        assert (variants_dir / "chr1").is_dir()
        assert (variants_dir / "chrX").is_dir()
        assert any((variants_dir / "chr1").glob("bucket_*.parquet"))
        assert any((variants_dir / "chrX").glob("bucket_*.parquet"))

    def test_no_variants_no_parquet(self, tmp_path):
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        _write_parquet(tmp_path / "sample_0.parquet", [])
        result = build_all_parquets(str(tmp_path), str(variants_dir))
        assert result == {}
        assert not list(variants_dir.glob("*.parquet"))

    def test_partitioned_flat_disabled(self, tmp_path):
        """build_chromosome_parquet with partitioned=False writes a flat .parquet file."""
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
        ])
        count = build_chromosome_parquet("chr1", str(tmp_path), str(variants_dir), partitioned=False)
        assert count == 1
        assert (variants_dir / "chr1.parquet").exists()
        assert not (variants_dir / "chr1").is_dir()

    def test_partitioned_creates_bucket_files(self, tmp_path):
        """build_chromosome_parquet with partitioned=True writes bucket files."""
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),           # bucket 0
            ("chr1", 1_500_000, "G", "C", 1, 0),      # bucket 1
        ])
        count = build_chromosome_parquet("chr1", str(tmp_path), str(variants_dir), partitioned=True)
        assert count == 2
        assert (variants_dir / "chr1" / "bucket_0.parquet").exists()
        assert (variants_dir / "chr1" / "bucket_1.parquet").exists()
        assert not (variants_dir / "chr1.parquet").exists()

    def test_parallel_build_matches_serial(self, tmp_path):
        """Parallel build (n_workers > 1) produces same results as serial (n_workers=1)."""
        variants_dir_s = tmp_path / "serial"
        variants_dir_p = tmp_path / "parallel"
        variants_dir_s.mkdir()
        variants_dir_p.mkdir()

        # Write multi-chrom temp parquet
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
            ("chr1", 2000, "G", "C", 2, 0),
            ("chr2", 500, "T", "A", 1, 0),
        ])

        # Serial build (partitioned=False to compare flat files easily)
        build_all_parquets(str(tmp_path), str(variants_dir_s), n_workers=1, partitioned=False)
        # Parallel build
        build_all_parquets(str(tmp_path), str(variants_dir_p), n_workers=2, partitioned=False)

        for chrom in ("chr1", "chr2"):
            s_file = variants_dir_s / f"{chrom}.parquet"
            p_file = variants_dir_p / f"{chrom}.parquet"
            assert s_file.exists(), f"serial {chrom} missing"
            assert p_file.exists(), f"parallel {chrom} missing"
            s_tbl = pq.read_table(s_file)
            p_tbl = pq.read_table(p_file)
            assert s_tbl["pos"].to_pylist() == p_tbl["pos"].to_pylist()
            assert s_tbl["alt"].to_pylist() == p_tbl["alt"].to_pylist()

    def test_parallel_build_default_workers(self, tmp_path):
        """build_all_parquets with n_workers=None uses all CPUs without error."""
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
        ])
        result = build_all_parquets(str(tmp_path), str(variants_dir), n_workers=None)
        assert "chr1" in result

    def test_consolidate_creates_partitioned_dir(self, tmp_path):
        """consolidate_temp_files() creates chrom=X/ subdirectories (Hive partitioning)."""
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
            ("chr2", 2000, "G", "C", 2, 0),
        ])
        _write_parquet(tmp_path / "sample_1.parquet", [
            ("chr1", 1500, "C", "G", 1, 1),
        ])
        consolidated = consolidate_temp_files(str(tmp_path))
        assert os.path.isdir(consolidated), "consolidated path must be a directory"
        subdirs = {e.name for e in os.scandir(consolidated) if e.is_dir()}
        assert "chrom=chr1" in subdirs
        assert "chrom=chr2" in subdirs
        # Each subdir should contain at least one parquet file
        assert any((tmp_path / "consolidated" / "chrom=chr1").glob("*.parquet"))
        assert any((tmp_path / "consolidated" / "chrom=chr2").glob("*.parquet"))

    def test_build_from_partitioned_dir(self, tmp_path):
        """build_chromosome_parquet() reads correctly from a Hive-partitioned directory."""
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),  # het sample 0, pass
        ])
        _write_parquet(tmp_path / "sample_1.parquet", [
            ("chr1", 1000, "A", "T", 2, 1),  # hom sample 1, pass
        ])
        consolidated = consolidate_temp_files(str(tmp_path))
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        count = build_chromosome_parquet(
            "chr1", str(tmp_path), str(variants_dir),
            partitioned=False, consolidated_path=consolidated,
        )
        assert count == 1
        tbl = pq.read_table(variants_dir / "chr1.parquet")
        het_bm = deserialize(bytes(tbl["het_bitmap"][0].as_py()))
        hom_bm = deserialize(bytes(tbl["hom_bitmap"][0].as_py()))
        fail_bm = deserialize(bytes(tbl["fail_bitmap"][0].as_py()))
        assert list(het_bm) == [0]
        assert list(hom_bm) == [1]
        assert list(fail_bm) == []

    def test_build_all_parquets_partitioned_workers(self, tmp_path):
        """With Hive-partitioned consolidated_path, build_all_parquets uses n_chroms workers (no cap of 4)."""
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
            ("chr2",  500, "G", "C", 1, 0),
        ])
        consolidated = consolidate_temp_files(str(tmp_path))
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        # n_workers=8 > 4 (old cap) but only 2 chroms present → should use 2 workers, not fail
        result = build_all_parquets(
            str(tmp_path), str(variants_dir),
            n_workers=8, partitioned=False, consolidated_path=consolidated,
        )
        assert "chr1" in result
        assert "chr2" in result


# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def preprocessed_db(tmp_path_factory):
    db_path = tmp_path_factory.mktemp("preprocess_db")
    data_dir = Path(__file__).parent / "data"
    run_preprocess(
        manifest_path=str(data_dir / "manifest.tsv"),
        output_dir=str(db_path),
        genome_build="GRCh37",
        bed_dir=str(data_dir / "beds"),
        threads=2,
    )
    return str(db_path)


def test_preprocess_creates_files(preprocessed_db):
    db = Path(preprocessed_db)
    assert (db / "manifest.json").exists()
    assert (db / "metadata.sqlite").exists()
    assert (db / "capture").is_dir()
    assert (db / "variants").is_dir()
    # Parquets may be in subdirectories (partitioned format) — use rglob
    assert any((db / "variants").rglob("*.parquet"))
    assert any((db / "capture").glob("*.pickle"))


def test_preprocess_manifest_content(preprocessed_db):
    manifest = json.loads((Path(preprocessed_db) / "manifest.json").read_text())
    assert manifest["genome_build"] == "GRCh37"
    assert manifest["sample_count"] == 10


def test_preprocess_query_matches_expected(preprocessed_db, data_dir):
    expected = json.loads((data_dir / "expected_results.json").read_text())
    db = Database(preprocessed_db)

    for case in expected:
        results = db.query(
            chrom=case["chrom"],
            pos=case["pos"],
            phenotype=case["phenotype"],
            sex=case["sex"],
        )
        matching = [
            r for r in results
            if r.variant.ref == case["ref"] and r.variant.alt == case["alt"]
        ]
        assert matching, f"No result for {case['description']}"
        r = matching[0]
        assert r.AC == case["AC"], f"{case['description']}: AC {r.AC} != {case['AC']}"
        assert r.AN == case["AN"], f"{case['description']}: AN {r.AN} != {case['AN']}"
        assert abs(r.AF - case["AF"]) < 1e-9, f"{case['description']}: AF mismatch"
        assert r.n_samples_eligible == case["n_eligible"], (
            f"{case['description']}: n_eligible {r.n_samples_eligible} != {case['n_eligible']}"
        )


# ---------------------------------------------------------------------------
# _get_chroms_from_file unit tests
# ---------------------------------------------------------------------------

class TestGetChromsFromFile:
    def test_single_chrom(self, tmp_path):
        _write_parquet(tmp_path / "s.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
        ])
        chroms = _get_chroms_from_file(str(tmp_path / "s.parquet"))
        assert chroms == ["chr1"]

    def test_multiple_chroms(self, tmp_path):
        _write_parquet(tmp_path / "s.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
            ("chrX", 5000, "G", "C", 1, 0),
        ])
        chroms = _get_chroms_from_file(str(tmp_path / "s.parquet"))
        assert "chr1" in chroms
        assert "chrX" in chroms

    def test_invalid_path_returns_empty(self, tmp_path):
        chroms = _get_chroms_from_file(str(tmp_path / "nonexistent.parquet"))
        assert chroms == []

    def test_unknown_contig_filtered(self, tmp_path):
        _write_parquet(tmp_path / "s.parquet", [
            ("chrFAKE", 1000, "A", "T", 1, 0),
        ])
        chroms = _get_chroms_from_file(str(tmp_path / "s.parquet"))
        assert "chrFAKE" not in chroms

    def test_unknown_contigs_reported_once(self, tmp_path, caplog):
        # Regression: the "Skipping unknown contig" warning sat downstream of this
        # filter, so it could never fire and contigs vanished with no signal.
        _write_parquet(tmp_path / "s.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
            ("chrFAKE", 1000, "A", "T", 1, 0),
            ("chrFAKE", 2000, "A", "T", 1, 0),
            ("chrOTHER", 3000, "A", "T", 1, 0),
        ])
        with caplog.at_level(logging.WARNING, logger="afquery.preprocess.build"):
            chroms = _get_chroms_from_file(str(tmp_path / "s.parquet"))
        assert chroms == ["chr1"]
        records = [r for r in caplog.records if "unknown contig" in r.message]
        assert len(records) == 1
        message = records[0].getMessage()
        assert "chrFAKE" in message
        assert "chrOTHER" in message


# ---------------------------------------------------------------------------
# _get_chroms_from_consolidated unit tests
# ---------------------------------------------------------------------------

class TestGetChromsFromConsolidated:
    def test_hive_directory_returns_chroms(self, tmp_path):
        (tmp_path / "chrom=chr1").mkdir()
        (tmp_path / "chrom=chr1" / "data.parquet").touch()
        (tmp_path / "chrom=chrX").mkdir()
        (tmp_path / "chrom=chrX" / "data.parquet").touch()
        chroms = _get_chroms_from_consolidated(str(tmp_path))
        assert "chr1" in chroms
        assert "chrX" in chroms

    def test_hive_no_chrom_dirs_returns_empty(self, tmp_path):
        (tmp_path / "other_dir").mkdir()
        (tmp_path / "another").mkdir()
        chroms = _get_chroms_from_consolidated(str(tmp_path))
        assert chroms == []

    def test_hive_unknown_chroms_filtered(self, tmp_path):
        (tmp_path / "chrom=chrFAKE").mkdir()
        (tmp_path / "chrom=chrFAKE" / "data.parquet").touch()
        chroms = _get_chroms_from_consolidated(str(tmp_path))
        assert "chrFAKE" not in chroms

    def test_hive_unknown_chroms_reported(self, tmp_path, caplog):
        (tmp_path / "chrom=chr1").mkdir()
        (tmp_path / "chrom=chrFAKE").mkdir()
        with caplog.at_level(logging.WARNING, logger="afquery.preprocess.build"):
            chroms = _get_chroms_from_consolidated(str(tmp_path))
        assert chroms == ["chr1"]
        assert any("chrFAKE" in r.getMessage() for r in caplog.records)

    def test_single_file_delegates(self, tmp_path):
        _write_parquet(tmp_path / "consolidated.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
        ])
        chroms = _get_chroms_from_consolidated(str(tmp_path / "consolidated.parquet"))
        assert "chr1" in chroms


# ---------------------------------------------------------------------------
# _build_chrom_worker unit tests
# ---------------------------------------------------------------------------

class TestBuildChromWorker:
    def test_returns_tuple(self, tmp_path):
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
        ])
        result = _build_chrom_worker(
            "chr1", str(tmp_path), str(variants_dir), 100_000, False
        )
        chrom, count, elapsed = result
        assert chrom == "chr1"
        assert count == 1
        assert elapsed >= 0.0

    def test_empty_source_returns_zero(self, tmp_path):
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        _write_parquet(tmp_path / "sample_0.parquet", [])
        result = _build_chrom_worker(
            "chr1", str(tmp_path), str(variants_dir), 100_000, False
        )
        chrom, count, elapsed = result
        assert chrom == "chr1"
        assert count == 0


# ---------------------------------------------------------------------------
# build_all_parquets: resume and force behaviors
# ---------------------------------------------------------------------------

class TestBuildAllParquetsResume:
    def test_resume_skips_existing_file(self, tmp_path):
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
        ])
        build_all_parquets(str(tmp_path), str(variants_dir), partitioned=False)
        out_file = variants_dir / "chr1.parquet"
        assert out_file.exists()
        mtime_before = out_file.stat().st_mtime

        # Second call with resume=True should skip
        result2 = build_all_parquets(
            str(tmp_path), str(variants_dir), partitioned=False, resume=True
        )
        assert result2 == {}
        assert out_file.stat().st_mtime == mtime_before

    def test_resume_false_rebuilds(self, tmp_path):
        variants_dir = tmp_path / "variants"
        variants_dir.mkdir()
        _write_parquet(tmp_path / "sample_0.parquet", [
            ("chr1", 1000, "A", "T", 1, 0),
        ])
        build_all_parquets(str(tmp_path), str(variants_dir), partitioned=False)
        out_file = variants_dir / "chr1.parquet"
        mtime_before = out_file.stat().st_mtime

        _time.sleep(0.05)
        result2 = build_all_parquets(
            str(tmp_path), str(variants_dir), partitioned=False, resume=False
        )
        assert "chr1" in result2
        assert out_file.stat().st_mtime > mtime_before


# ---------------------------------------------------------------------------
# _write_sqlite and _write_manifest unit tests
# ---------------------------------------------------------------------------

class TestWriteSqliteAndManifest:
    def test_write_sqlite_creates_all_tables(self, tmp_path):
        samples = [Sample(0, "S00", "male", 0), Sample(1, "S01", "female", 0)]
        technologies = [Technology(0, "WGS", None)]
        pairs = [(0, "E11.9"), (1, "I10")]
        ingested_at = "2026-01-01T00:00:00+00:00"

        _write_sqlite(str(tmp_path), samples, technologies, pairs,
                      vcf_paths=None, ingested_at=ingested_at)

        con = sqlite3.connect(str(tmp_path / "metadata.sqlite"))
        tables = {r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'"
        ).fetchall()}
        con.close()

        assert "samples" in tables
        assert "technologies" in tables
        assert "sample_phenotype" in tables
        assert "precomputed_bitmaps" in tables
        assert "changelog" in tables

    def test_write_sqlite_inserts_samples(self, tmp_path):
        samples = [Sample(0, "S00", "male", 0), Sample(1, "S01", "female", 0)]
        technologies = [Technology(0, "WGS", None)]
        pairs = [(0, "E11.9")]

        _write_sqlite(str(tmp_path), samples, technologies, pairs)

        con = sqlite3.connect(str(tmp_path / "metadata.sqlite"))
        count = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        con.close()
        assert count == 2

    def test_write_sqlite_adds_changelog_entry(self, tmp_path):
        samples = [Sample(0, "S00", "male", 0)]
        technologies = [Technology(0, "WGS", None)]
        pairs = [(0, "E11.9")]

        _write_sqlite(str(tmp_path), samples, technologies, pairs)

        con = sqlite3.connect(str(tmp_path / "metadata.sqlite"))
        row = con.execute(
            "SELECT event_type, notes FROM changelog WHERE event_type='preprocess'"
        ).fetchone()
        con.close()
        assert row is not None
        assert "1 samples ingested" in row[1]

    def test_write_manifest_fields(self, tmp_path):
        _write_manifest(str(tmp_path), "GRCh37", 5, db_version="1.0")
        manifest = json.loads((tmp_path / "manifest.json").read_text())
        assert manifest["genome_build"] == "GRCh37"
        assert manifest["schema_version"] == "2.0"
        assert manifest["pass_only_filter"] is True
        assert manifest["sample_count"] == 5
        assert "created_at" in manifest

    def test_write_manifest_custom_version(self, tmp_path):
        _write_manifest(str(tmp_path), "GRCh38", 3, db_version="2.5")
        manifest = json.loads((tmp_path / "manifest.json").read_text())
        assert manifest["db_version"] == "2.5"
        assert manifest["genome_build"] == "GRCh38"

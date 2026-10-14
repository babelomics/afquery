"""Tests for the dump command (bulk allele frequency export)."""
import csv
import io
import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from afquery.cli import cli
from afquery.database import Database
from afquery.dump import _build_groups, _dump_bucket_worker, _discover_flat_buckets
from afquery.models import SampleFilter
from afquery.query import QueryEngine


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_csv(text: str) -> list[dict]:
    """Parse CSV string into list of dicts."""
    reader = csv.DictReader(io.StringIO(text))
    return list(reader)


def _check_ac_zero_af(row: dict, suffix: str) -> None:
    """An AC=0 column set has AF=0, or an empty AF when AN=0 (only failed calls)."""
    if row[f"AN{suffix}"] == "0":
        assert row[f"AF{suffix}"] == ""
    else:
        assert float(row[f"AF{suffix}"]) == 0.0


def _dump_to_string(db_path, **kwargs) -> tuple[str, dict]:
    """Run Database.dump() into a StringIO and return (csv_text, stats)."""
    buf = io.StringIO()
    database = Database(db_path)
    stats = database.dump(output=buf, **kwargs)
    return buf.getvalue(), stats


# ---------------------------------------------------------------------------
# Unit tests: _build_groups
# ---------------------------------------------------------------------------

class TestBuildGroups:
    def test_no_flags_returns_empty(self, test_db):
        engine = QueryEngine(test_db)
        base_sf = SampleFilter()
        groups = _build_groups(engine, base_sf, False, False, [], False)
        assert groups == []

    def test_by_sex(self, test_db):
        engine = QueryEngine(test_db)
        base_sf = SampleFilter()
        groups = _build_groups(engine, base_sf, True, False, [], False)
        labels = [g[0] for g in groups]
        assert "male" in labels
        assert "female" in labels
        assert len(groups) == 2

    def test_by_sex_sfs(self, test_db):
        engine = QueryEngine(test_db)
        base_sf = SampleFilter()
        groups = _build_groups(engine, base_sf, True, False, [], False)
        male_sf = next(sf for lbl, sf in groups if lbl == "male")
        female_sf = next(sf for lbl, sf in groups if lbl == "female")
        assert male_sf.sex == "male"
        assert female_sf.sex == "female"

    def test_by_phenotype(self, test_db):
        engine = QueryEngine(test_db)
        base_sf = SampleFilter()
        groups = _build_groups(engine, base_sf, False, False, ["E11.9", "I10"], False)
        labels = [g[0] for g in groups]
        assert "E11.9" in labels
        assert "I10" in labels
        assert len(groups) == 2

    def test_by_sex_and_phenotype_cartesian(self, test_db):
        engine = QueryEngine(test_db)
        base_sf = SampleFilter()
        groups = _build_groups(engine, base_sf, True, False, ["E11.9"], False)
        labels = [g[0] for g in groups]
        assert "male_E11.9" in labels
        assert "female_E11.9" in labels
        assert len(groups) == 2

    def test_by_tech(self, test_db):
        engine = QueryEngine(test_db)
        base_sf = SampleFilter()
        groups = _build_groups(engine, base_sf, False, True, [], False)
        labels = [g[0] for g in groups]
        # Test DB has WGS, WES_kit_A, WES_kit_B
        assert len(groups) == 3

    def test_all_groups_includes_phenotypes_from_db(self, test_db):
        engine = QueryEngine(test_db)
        base_sf = SampleFilter()
        # all_groups includes sex × tech × phenotype Cartesian product
        groups = _build_groups(engine, base_sf, False, False, [], True)
        labels = [g[0] for g in groups]
        # All phenotype codes from DB must appear in at least one label
        for code in ["C50", "E11.9", "G20", "I10", "J45"]:
            assert any(code in lbl for lbl in labels), f"{code} not found in any label"
        # Should have many groups (sex=2 × tech=3 × pheno=5 = 30)
        assert len(groups) == 30

    def test_base_sf_excludes_preserved(self, test_db):
        engine = QueryEngine(test_db)
        base_sf = SampleFilter(phenotype_exclude=["G20"], tech_exclude=["WES_kit_B"])
        groups = _build_groups(engine, base_sf, True, False, [], False)
        for _, sf in groups:
            assert "G20" in sf.phenotype_exclude
            assert "WES_kit_B" in sf.tech_exclude

    def test_base_sf_sex_override(self, test_db):
        """When base_sf has sex=male but by_sex is True, female group should still be created."""
        engine = QueryEngine(test_db)
        base_sf = SampleFilter(sex="male")
        groups = _build_groups(engine, base_sf, True, False, [], False)
        labels = [g[0] for g in groups]
        assert "male" in labels
        assert "female" in labels


# ---------------------------------------------------------------------------
# Functional tests: Database.dump()
# ---------------------------------------------------------------------------

class TestDumpBasic:
    def test_returns_stats(self, test_db):
        _, stats = _dump_to_string(test_db)
        assert "n_rows" in stats
        assert "n_buckets" in stats
        assert "n_chroms" in stats

    def test_csv_header_base_columns(self, test_db):
        text, _ = _dump_to_string(test_db)
        rows = _parse_csv(text)
        assert len(rows) > 0
        assert "chrom" in rows[0]
        assert "pos" in rows[0]
        assert "ref" in rows[0]
        assert "alt" in rows[0]
        assert "AC" in rows[0]
        assert "AN" in rows[0]
        assert "AF" in rows[0]
        assert "N_HET" in rows[0]
        assert "N_HOM_ALT" in rows[0]
        assert "N_HOM_REF" in rows[0]
        assert "N_FAIL" in rows[0]

    def test_ac_positive_filter(self, test_db):
        """By default all exported rows must have AC > 0."""
        text, _ = _dump_to_string(test_db)
        rows = _parse_csv(text)
        for row in rows:
            assert int(row["AC"]) > 0, f"AC=0 in row: {row}"

    def test_include_ac_zero_flag(self, test_db):
        """include_ac_zero=True produces a superset; AC=0 rows have AF=0, or an
        empty AF when AN=0 because every eligible sample has a failed call.

        chr1:3500 G>C has only S07 (phenotype J45) as carrier. When filtered by
        phenotype I10, no I10 sample carries this variant but several are covered
        at that position → AC=0, AN>0.
        """
        # Default: chr1:3500 G>C is omitted for phenotype I10 (no carriers)
        text_default, _ = _dump_to_string(test_db, chrom="chr1", phenotype=["I10"])
        rows_default = _parse_csv(text_default)
        default_positions = {(r["pos"], r["alt"]) for r in rows_default}
        assert ("3500", "C") not in default_positions, "chr1:3500 G>C should be absent by default"

        # With include_ac_zero=True it should appear
        text_all, _ = _dump_to_string(
            test_db, chrom="chr1", phenotype=["I10"], include_ac_zero=True
        )
        rows_all = _parse_csv(text_all)
        assert len(rows_all) >= len(rows_default)

        ac_zero_rows = [r for r in rows_all if int(r["AC"]) == 0]
        assert len(ac_zero_rows) > 0, "Expected at least one AC=0 row"
        for row in ac_zero_rows:
            _check_ac_zero_af(row, "")

        # chr1:5000: both eligible I10 samples (S01, S03) have a failed call there
        [failed_only] = [r for r in rows_all if r["pos"] == "5000"]
        assert failed_only["AN"] == "0"
        assert failed_only["AF"] == ""
        assert failed_only["N_FAIL"] == "2"

    def test_include_ac_zero_cli(self, test_db, tmp_path):
        """--all-variants CLI flag includes AC=0 rows for a phenotype-filtered dump."""
        runner = CliRunner()
        out_default = str(tmp_path / "default.csv")
        out_all = str(tmp_path / "all.csv")

        result_default = runner.invoke(cli, [
            "dump", "--db", test_db, "--chrom", "chr1",
            "--phenotype", "I10", "--output", out_default,
        ])
        assert result_default.exit_code == 0
        rows_default = _parse_csv(Path(out_default).read_text())

        result_all = runner.invoke(cli, [
            "dump", "--db", test_db, "--chrom", "chr1",
            "--phenotype", "I10", "--all-variants", "--output", out_all,
        ])
        assert result_all.exit_code == 0
        rows_all = _parse_csv(Path(out_all).read_text())

        assert len(rows_all) >= len(rows_default)
        ac_zero_rows = [r for r in rows_all if r.get("AC") == "0"]
        assert len(ac_zero_rows) > 0, "Expected AC=0 rows with --all-variants"

    def test_include_ac_zero_with_disaggregation(self, test_db):
        """AC=0 rows appear correctly when combined with by_sex disaggregation.

        All sex-disaggregated AC_ columns must also be 0 for AC=0 base rows.
        """
        text, _ = _dump_to_string(
            test_db, chrom="chr1", phenotype=["I10"],
            by_sex=True, include_ac_zero=True,
        )
        rows = _parse_csv(text)
        ac_zero_rows = [r for r in rows if int(r["AC"]) == 0]
        assert len(ac_zero_rows) > 0, "Expected AC=0 rows with disaggregation"
        for row in ac_zero_rows:
            _check_ac_zero_af(row, "")
            for sex in ("male", "female"):
                _check_ac_zero_af(row, f"_{sex}")
            sex_ac_cols = [k for k in row if k.startswith("AC_")]
            assert sex_ac_cols, "Expected disaggregated AC_ columns"
            for col in sex_ac_cols:
                assert int(row[col]) == 0, f"Expected {col}=0, got {row[col]}"

    def test_values_match_query(self, test_db):
        """AC/AN/AF/N_HET match direct query() results."""
        database = Database(test_db)
        text, _ = _dump_to_string(test_db, chrom="chr1")
        rows = _parse_csv(text)
        for row in rows:
            results = database.query(chrom=row["chrom"], pos=int(row["pos"]),
                                     ref=row["ref"], alt=row["alt"])
            assert len(results) == 1
            r = results[0]
            assert int(row["AC"]) == r.AC
            assert int(row["AN"]) == r.AN
            assert abs(float(row["AF"]) - r.AF) < 1e-9
            assert int(row["N_HET"]) == r.N_HET
            assert int(row["N_HOM_ALT"]) == r.N_HOM_ALT
            assert int(row["N_FAIL"]) == r.N_FAIL

    def test_all_chroms_present(self, test_db):
        """Dump without --chrom includes all chromosomes that have data with AC>0."""
        text, _ = _dump_to_string(test_db)
        rows = _parse_csv(text)
        chroms = {r["chrom"] for r in rows}
        assert "chr1" in chroms

    def test_genomic_order(self, test_db):
        """Rows for chr1 are sorted by (pos, alt)."""
        text, _ = _dump_to_string(test_db, chrom="chr1")
        rows = _parse_csv(text)
        positions = [(int(r["pos"]), r["alt"]) for r in rows]
        assert positions == sorted(positions)


class TestDumpFilters:
    def test_chrom_filter(self, test_db):
        text, _ = _dump_to_string(test_db, chrom="chr1")
        rows = _parse_csv(text)
        assert all(r["chrom"] == "chr1" for r in rows)
        assert len(rows) > 0

    def test_region_filter(self, test_db):
        text, _ = _dump_to_string(test_db, chrom="chr1", start=1000, end=2000)
        rows = _parse_csv(text)
        assert all(1000 <= int(r["pos"]) <= 2000 for r in rows)
        assert len(rows) == 1  # only chr1:1500

    def test_region_filter_excludes_others(self, test_db):
        text, _ = _dump_to_string(test_db, chrom="chr1", start=3000, end=4000)
        rows = _parse_csv(text)
        assert all(3000 <= int(r["pos"]) <= 4000 for r in rows)
        assert len(rows) == 1  # only chr1:3500

    def test_phenotype_filter_reduces_rows(self, test_db):
        """Filtering by phenotype may change AC values or reduce rows."""
        text_all, _ = _dump_to_string(test_db, chrom="chr1")
        text_pheno, _ = _dump_to_string(test_db, chrom="chr1", phenotype=["E11.9"])
        rows_all = _parse_csv(text_all)
        rows_pheno = _parse_csv(text_pheno)
        # Phenotype filter should produce <= rows and <= AC
        assert len(rows_pheno) <= len(rows_all)

    def test_sex_filter_male_only(self, test_db):
        """Male-only dump has different AN (only male samples counted)."""
        database = Database(test_db)
        buf = io.StringIO()
        database.dump(output=buf, chrom="chr1", sex="male")
        rows = _parse_csv(buf.getvalue())
        # Cross-check with query
        for row in rows:
            results = database.query(
                chrom=row["chrom"], pos=int(row["pos"]),
                ref=row["ref"], alt=row["alt"], sex="male",
            )
            assert len(results) == 1
            assert int(row["AC"]) == results[0].AC
            assert int(row["AN"]) == results[0].AN


class TestDumpDisaggregation:
    def test_by_sex_columns_present(self, test_db):
        database = Database(test_db)
        buf = io.StringIO()
        database.dump(output=buf, chrom="chr1", by_sex=True)
        text = buf.getvalue()
        rows = _parse_csv(text)
        assert len(rows) > 0
        assert "AC_male" in rows[0]
        assert "AN_male" in rows[0]
        assert "AF_male" in rows[0]
        assert "AC_female" in rows[0]
        assert "AN_female" in rows[0]

    def test_by_sex_values_match_query(self, test_db):
        """AC_male / AC_female match direct query(sex=male/female)."""
        database = Database(test_db)
        buf = io.StringIO()
        database.dump(output=buf, chrom="chr1", by_sex=True)
        rows = _parse_csv(buf.getvalue())
        for row in rows:
            for sx in ("male", "female"):
                results = database.query(
                    chrom=row["chrom"], pos=int(row["pos"]),
                    ref=row["ref"], alt=row["alt"], sex=sx,
                )
                expected_ac = results[0].AC if results else 0
                assert int(row[f"AC_{sx}"]) == expected_ac

    def test_by_phenotype_columns_present(self, test_db):
        database = Database(test_db)
        buf = io.StringIO()
        database.dump(output=buf, chrom="chr1", by_phenotype=["E11.9", "I10"])
        text = buf.getvalue()
        rows = _parse_csv(text)
        assert len(rows) > 0
        assert "AC_E11.9" in rows[0]
        assert "AC_I10" in rows[0]

    def test_by_phenotype_values_match_query(self, test_db):
        database = Database(test_db)
        buf = io.StringIO()
        database.dump(output=buf, chrom="chr1", by_phenotype=["E11.9"])
        rows = _parse_csv(buf.getvalue())
        for row in rows:
            results = database.query(
                chrom=row["chrom"], pos=int(row["pos"]),
                ref=row["ref"], alt=row["alt"], phenotype=["E11.9"],
            )
            expected_ac = results[0].AC if results else 0
            assert int(row["AC_E11.9"]) == expected_ac

    def test_by_sex_and_phenotype_cartesian_labels(self, test_db):
        database = Database(test_db)
        buf = io.StringIO()
        database.dump(output=buf, chrom="chr1", by_sex=True, by_phenotype=["E11.9"])
        rows = _parse_csv(buf.getvalue())
        assert len(rows) > 0
        assert "AC_male_E11.9" in rows[0]
        assert "AC_female_E11.9" in rows[0]

    def test_n_fail_in_group_columns(self, test_db):
        """N_FAIL_<label> columns present in disaggregated output."""
        database = Database(test_db)
        buf = io.StringIO()
        database.dump(output=buf, chrom="chr1", by_sex=True)
        rows = _parse_csv(buf.getvalue())
        assert len(rows) > 0
        assert "N_FAIL_male" in rows[0]
        assert "N_FAIL_female" in rows[0]

    def test_by_tech_column_count(self, test_db):
        """by_tech produces columns for each technology in the DB."""
        database = Database(test_db)
        buf = io.StringIO()
        database.dump(output=buf, chrom="chr1", by_tech=True)
        rows = _parse_csv(buf.getvalue())
        assert len(rows) > 0
        # Test DB has WGS, WES_kit_A, WES_kit_B
        for tech in ("WGS", "WES_kit_A", "WES_kit_B"):
            assert f"AC_{tech}" in rows[0], f"Missing AC_{tech}"


# ---------------------------------------------------------------------------
# CLI tests
# ---------------------------------------------------------------------------

class TestDumpCLI:
    def test_basic_stdout(self, test_db):
        runner = CliRunner()
        result = runner.invoke(cli, ["dump", "--db", test_db, "--chrom", "chr1"])
        assert result.exit_code == 0
        assert "chrom,pos,ref,alt" in result.output

    def test_output_file(self, test_db, tmp_path):
        runner = CliRunner()
        out_path = str(tmp_path / "out.csv")
        result = runner.invoke(cli, [
            "dump", "--db", test_db, "--chrom", "chr1", "--output", out_path,
        ])
        assert result.exit_code == 0
        with open(out_path) as f:
            content = f.read()
        assert "chrom,pos,ref,alt" in content
        assert "row(s) exported" in result.output  # stats go to stderr (mixed in CliRunner)

    def test_start_without_chrom_fails(self, test_db):
        runner = CliRunner()
        result = runner.invoke(cli, ["dump", "--db", test_db, "--start", "1000"])
        assert result.exit_code != 0

    def test_end_without_chrom_fails(self, test_db):
        runner = CliRunner()
        result = runner.invoke(cli, ["dump", "--db", test_db, "--end", "2000"])
        assert result.exit_code != 0

    def test_start_greater_than_end_fails(self, test_db):
        runner = CliRunner()
        result = runner.invoke(cli, [
            "dump", "--db", test_db, "--chrom", "chr1",
            "--start", "5000", "--end", "1000",
        ])
        assert result.exit_code != 0

    def test_by_sex_header(self, test_db):
        runner = CliRunner()
        result = runner.invoke(cli, [
            "dump", "--db", test_db, "--chrom", "chr1", "--by-sex",
        ])
        assert result.exit_code == 0
        assert "AC_male" in result.output
        assert "AC_female" in result.output

    def test_all_groups_warning_in_help(self, test_db):
        runner = CliRunner()
        result = runner.invoke(cli, ["dump", "--help"])
        assert result.exit_code == 0
        assert "all-groups" in result.output.lower() or "all_groups" in result.output.lower()

    def test_by_phenotype_cli(self, test_db):
        runner = CliRunner()
        result = runner.invoke(cli, [
            "dump", "--db", test_db, "--chrom", "chr1",
            "--by-phenotype", "E11.9",
        ])
        assert result.exit_code == 0
        assert "AC_E11.9" in result.output

    def test_stats_on_stderr(self, test_db, tmp_path):
        """Stats (row count) appear in CLI output (stderr via err=True in click.echo)."""
        runner = CliRunner()
        out_path = str(tmp_path / "out2.csv")
        result = runner.invoke(cli, [
            "dump", "--db", test_db, "--chrom", "chr1", "--output", out_path,
        ])
        assert result.exit_code == 0
        # stats line should be emitted (CliRunner merges stderr into output by default)
        assert "row(s) exported" in result.output


# ---------------------------------------------------------------------------
# Unit tests: _discover_flat_buckets
# ---------------------------------------------------------------------------

class TestDiscoverFlatBuckets:
    def test_basic_returns_bucket_ids(self, test_db):
        flat_path = Path(test_db) / "variants" / "chr1.parquet"
        buckets = _discover_flat_buckets(str(flat_path), None, None)
        assert isinstance(buckets, list)
        assert len(buckets) > 0
        # All chr1 variants (pos 1500, 3500, 5000) are in bucket 0 (< 1M)
        assert 0 in buckets

    def test_range_excludes_variants(self, test_db):
        flat_path = Path(test_db) / "variants" / "chr1.parquet"
        # Start > all chr1 variant positions → no buckets
        buckets = _discover_flat_buckets(str(flat_path), 9_000_000, None)
        assert buckets == []

    def test_range_includes_variants(self, test_db):
        flat_path = Path(test_db) / "variants" / "chr1.parquet"
        # pos_start=1000, pos_end=6000 covers all chr1 variants
        buckets = _discover_flat_buckets(str(flat_path), 1000, 6000)
        assert 0 in buckets


# ---------------------------------------------------------------------------
# Worker path coverage
# ---------------------------------------------------------------------------

class TestDumpWorkerPaths:
    def test_single_worker_explicit(self, test_db):
        """n_workers=1 uses serial path; results match default."""
        buf1 = io.StringIO()
        database = Database(test_db)
        database.dump(output=buf1, chrom="chr1", n_workers=1)
        rows1 = _parse_csv(buf1.getvalue())

        buf2 = io.StringIO()
        database.dump(output=buf2, chrom="chr1")
        rows2 = _parse_csv(buf2.getvalue())

        assert len(rows1) == len(rows2)
        for r1, r2 in zip(rows1, rows2):
            assert r1["pos"] == r2["pos"]
            assert r1["AC"] == r2["AC"]

    def test_multiworker_all_chroms(self, test_db):
        """n_workers=2 across all chroms exercises ProcessPoolExecutor path."""
        buf = io.StringIO()
        database = Database(test_db)
        stats = database.dump(output=buf, n_workers=2)
        rows = _parse_csv(buf.getvalue())
        assert stats["n_rows"] == len(rows)
        assert len(rows) > 0


# ---------------------------------------------------------------------------
# Output and edge-case coverage
# ---------------------------------------------------------------------------

class TestDumpEdgeCases:
    def test_empty_region_returns_zero_rows(self, test_db):
        """A region with no variants yields 0 rows without error."""
        buf = io.StringIO()
        database = Database(test_db)
        stats = database.dump(output=buf, chrom="chr1", start=9_000_000, end=9_100_000)
        assert stats["n_rows"] == 0
        assert stats["n_buckets"] == 0

    def test_dump_to_filelike_object(self, test_db):
        """Passing a StringIO as output writes CSV content correctly."""
        buf = io.StringIO()
        database = Database(test_db)
        stats = database.dump(output=buf, chrom="chr1")
        content = buf.getvalue()
        assert content.startswith("chrom,pos,ref,alt")
        assert stats["n_rows"] > 0

    def test_dump_chrom_filter_cli(self, test_db):
        runner = CliRunner()
        result = runner.invoke(cli, [
            "dump", "--db", test_db, "--chrom", "chr1",
        ])
        assert result.exit_code == 0
        assert "chrom,pos,ref,alt" in result.output
        assert "chr1" in result.output

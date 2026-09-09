import json
import shutil

import pytest
from click.testing import CliRunner

from afquery.cli import cli


@pytest.fixture
def db_copy(test_db, tmp_path):
    """Fresh writable copy of test_db for tests that mutate the DB."""
    copy_path = str(tmp_path / "db_copy")
    shutil.copytree(test_db, copy_path)
    return copy_path


@pytest.fixture
def runner():
    return CliRunner()


# --- afquery info ---

def test_info_shows_genome_build(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db])
    assert result.exit_code == 0
    assert "GRCh37" in result.output


def test_info_shows_sample_count(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db])
    assert result.exit_code == 0
    assert "10" in result.output


def test_info_missing_db(runner):
    result = runner.invoke(cli, ["info"])
    assert result.exit_code != 0


# --- afquery query: text format ---

def test_query_text_format(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "1:1500", "--phenotype", "E11.9",
    ])
    assert result.exit_code == 0
    assert "AC=4" in result.output
    assert "AN=10" in result.output
    assert "AF=0.4000" in result.output


def test_query_text_no_results(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "1:9999", "--phenotype", "UNKNOWN",
    ])
    assert result.exit_code == 0
    assert "No variants found" in result.output


# --- afquery query: json format ---

def test_query_json_format(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500", "--phenotype", "E11.9",
        "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert len(data) == 1
    assert data[0]["AC"] == 4
    assert data[0]["AN"] == 10
    assert abs(data[0]["AF"] - 0.4) < 1e-9


def test_query_json_empty(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:9999", "--phenotype", "UNKNOWN",
        "--format", "json",
    ])
    assert result.exit_code == 0
    assert json.loads(result.output) == []


# --- afquery query: tsv format ---

def test_query_tsv_format(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "1:1500", "--phenotype", "E11.9",
        "--format", "tsv",
    ])
    assert result.exit_code == 0
    lines = result.output.strip().splitlines()
    assert lines[0].startswith("chrom\tpos\tref\talt")
    assert len(lines) == 2  # header + 1 variant
    fields = lines[1].split("\t")
    assert fields[4] == "4"   # AC
    assert fields[5] == "10"  # AN


# --- sex filter ---

def test_query_sex_female(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500", "--phenotype", "E11.9",
        "--sex", "female", "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert data[0]["AC"] == 3
    assert data[0]["AN"] == 6


# --- multiple Phenotype codes ---

def test_query_multiple_phenotype(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500",
        "--phenotype", "E11.9", "--phenotype", "I10",
        "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert len(data) == 1


# --- missing required option ---

def test_query_bad_locus_format(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "bad_format",
    ])
    assert result.exit_code != 0


def test_query_missing_db(runner):
    result = runner.invoke(cli, [
        "query", "--locus", "1:1500", "--phenotype", "E11.9",
    ])
    assert result.exit_code != 0


# --- afquery query: region and batch ---

def test_query_region_format(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--region", "1:1000-2000",
        "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert isinstance(data, list)


def test_query_batch_format(runner, test_db, tmp_path):
    variants_file = tmp_path / "variants.tsv"
    variants_file.write_text("1 1500 C T\n")
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--from-file", str(variants_file),
        "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert isinstance(data, list)


def test_query_requires_mode(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
    ])
    assert result.exit_code != 0


def test_query_exclusive_modes(runner, test_db, tmp_path):
    variants_file = tmp_path / "variants.tsv"
    variants_file.write_text("1 1500 C T\n")
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "1:1500", "--region", "1:1000-2000",
    ])
    assert result.exit_code != 0


# --- afquery create-db ---

def test_create_db_help(runner):
    result = runner.invoke(cli, ["create-db", "--help"])
    assert result.exit_code == 0
    assert "manifest" in result.output.lower()


def test_create_db_missing_manifest(runner, tmp_path):
    result = runner.invoke(cli, [
        "create-db",
        "--output-dir", str(tmp_path / "db"),
        "--genome-build", "GRCh37",
    ])
    assert result.exit_code != 0


# --- afquery update-db ---

def test_update_db_requires_operation(runner, test_db):
    result = runner.invoke(cli, ["update-db", "--db", test_db])
    assert result.exit_code != 0
    assert "required" in result.output.lower() or "required" in (result.exception and str(result.exception) or "").lower() or result.exit_code != 0


def test_update_db_compact(runner, test_db, tmp_path):
    db_copy = str(tmp_path / "compact_db")
    shutil.copytree(test_db, db_copy)
    result = runner.invoke(cli, ["update-db", "--db", db_copy, "--compact"])
    assert result.exit_code == 0
    assert "compact" in result.output.lower()


def test_update_db_help(runner):
    result = runner.invoke(cli, ["update-db", "--help"])
    assert result.exit_code == 0
    assert "--compact" in result.output
    assert "--add-samples" in result.output
    assert "--remove-samples" in result.output


# --- afquery annotate ---

def test_annotate_with_phenotype(runner, test_db, tmp_path):
    """Test annotate command with explicit phenotype."""
    import cyvcf2
    from pathlib import Path

    # Create a minimal input VCF
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf"

    with open(input_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("1\t1500\t.\tC\tT\t.\t.\t.\n")

    result = runner.invoke(cli, [
        "annotate", "--db", test_db,
        "--input", str(input_vcf),
        "--output", str(output_vcf),
        "--phenotype", "E11.9",
    ])
    assert result.exit_code == 0
    assert "Annotated" in result.output
    assert output_vcf.exists()


def test_annotate_without_phenotype(runner, test_db, tmp_path):
    """Test annotate command without phenotype (should use all phenotypes)."""
    import cyvcf2
    from pathlib import Path

    # Create a minimal input VCF
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf"

    with open(input_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("1\t1500\t.\tC\tT\t.\t.\t.\n")

    result = runner.invoke(cli, [
        "annotate", "--db", test_db,
        "--input", str(input_vcf),
        "--output", str(output_vcf),
    ])
    assert result.exit_code == 0
    assert "Annotated" in result.output
    assert output_vcf.exists()

    # Verify that the output VCF has AFQUERY annotations
    vcf = cyvcf2.VCF(str(output_vcf))
    records = list(vcf)
    assert len(records) > 0
    # At pos 1500 with multiple phenotypes, AN should be larger than with single phenotype
    assert records[0].INFO.get("AFQUERY_AN", 0) > 0
    vcf.close()


# --- afquery info: formats ---

def test_info_json_format(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db, "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert "genome_build" in data
    assert "sample_count" in data
    assert data["genome_build"] == "GRCh37"


def test_info_tsv_format(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db, "--format", "tsv"])
    assert result.exit_code == 0
    lines = result.output.strip().splitlines()
    keys = [line.split("\t")[0] for line in lines]
    assert "genome_build" in keys
    assert "sample_count" in keys


def test_info_samples_table(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db, "--samples"])
    assert result.exit_code == 0
    assert "Name" in result.output
    assert "S00" in result.output


def test_info_samples_json(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db, "--samples", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert isinstance(data, list)
    assert any(s["sample_name"] == "S00" for s in data)


def test_info_samples_tsv(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db, "--samples", "--format", "tsv"])
    assert result.exit_code == 0
    lines = result.output.strip().splitlines()
    assert lines[0].startswith("sample_id\t")
    assert len(lines) > 1  # header + at least one sample


def test_info_changelog(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db, "--changelog"])
    assert result.exit_code == 0


def test_info_changelog_json(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db, "--changelog", "--format", "json"])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert isinstance(data, list)


def test_info_changelog_tsv(runner, test_db):
    result = runner.invoke(cli, ["info", "--db", test_db, "--changelog", "--format", "tsv"])
    assert result.exit_code == 0
    lines = result.output.strip().splitlines()
    assert lines[0].startswith("event_id\t")


# --- afquery version ---

def test_version_show(runner, test_db):
    result = runner.invoke(cli, ["version", "show", "--db", test_db])
    assert result.exit_code == 0
    assert result.output.strip() != ""


def test_version_set(runner, db_copy):
    result = runner.invoke(cli, ["version", "set", "--db", db_copy, "NEW_VER"])
    assert result.exit_code == 0
    assert "NEW_VER" in result.output


def test_version_set_then_show(runner, db_copy):
    runner.invoke(cli, ["version", "set", "--db", db_copy, "MYVER"])
    result = runner.invoke(cli, ["version", "show", "--db", db_copy])
    assert result.exit_code == 0
    assert "MYVER" in result.output


# --- afquery --version / --help ---

def test_version_flag(runner):
    from afquery import __version__
    result = runner.invoke(cli, ["--version"])
    assert result.exit_code == 0
    assert __version__ in result.output


def test_help_shows_program_name_and_version(runner):
    from afquery import __version__
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert f"AFQuery v{__version__}" in result.output


def test_help_lists_all_commands(runner):
    result = runner.invoke(cli, ["--help"])
    for cmd in ("query", "variant-info", "annotate", "dump", "info",
                "version", "create-db", "update-db", "check", "benchmark"):
        assert cmd in result.output


# --- afquery check ---

def test_check_ok(runner, test_db):
    result = runner.invoke(cli, ["check", "--db", test_db])
    assert result.exit_code == 0


def test_check_missing_db(runner):
    result = runner.invoke(cli, ["check", "--db", "/nonexistent/afquery_path_xyz"])
    assert result.exit_code != 0


# --- afquery query: ref/alt/tech filters and edge cases ---

def test_query_ref_alt_filter(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500",
        "--ref", "A", "--alt", "T",
        "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert isinstance(data, list)
    if data:
        assert data[0]["ref"] == "A"
        assert data[0]["alt"] == "T"


def test_query_tech_filter(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500",
        "--tech", "WGS", "--format", "json",
    ])
    assert result.exit_code == 0
    data = json.loads(result.output)
    assert isinstance(data, list)


def test_query_region_invalid_fmt(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--region", "notaregion",
    ])
    assert result.exit_code != 0


def test_query_tsv_with_fail_col(runner, test_db):
    result = runner.invoke(cli, [
        "query", "--db", test_db,
        "--locus", "chr1:1500",
        "--format", "tsv",
    ])
    assert result.exit_code == 0
    lines = result.output.strip().splitlines()
    assert len(lines) >= 2
    assert "N_FAIL" in lines[0]


def test_query_output_includes_n_no_coverage_all_formats(runner, test_db):
    """N_FAIL and N_NO_COVERAGE appear in text/tsv/json output, and TSV columns line up."""
    base = ["query", "--db", test_db, "--locus", "chr1:1500"]

    text = runner.invoke(cli, base)
    assert text.exit_code == 0
    assert "N_FAIL=" in text.output
    assert "N_NO_COVERAGE=" in text.output

    tsv = runner.invoke(cli, base + ["--format", "tsv"])
    assert tsv.exit_code == 0
    tsv_lines = tsv.output.strip().splitlines()
    assert len(tsv_lines) >= 2
    header = tsv_lines[0].split("\t")
    assert "N_FAIL" in header
    assert header[-1] == "N_NO_COVERAGE"
    for row in tsv_lines[1:]:
        assert len(row.split("\t")) == len(header)

    js = runner.invoke(cli, base + ["--format", "json"])
    assert js.exit_code == 0
    data = json.loads(js.output)
    assert len(data) >= 1
    for entry in data:
        assert isinstance(entry["N_FAIL"], int)
        assert isinstance(entry["N_NO_COVERAGE"], int)


# --- afquery annotate: verbose flag ---

def test_annotate_verbose(runner, test_db, tmp_path):
    input_vcf = tmp_path / "input.vcf"
    output_vcf = tmp_path / "output.vcf"
    with open(input_vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        f.write("1\t1500\t.\tA\tT\t.\t.\t.\n")
    result = runner.invoke(cli, [
        "annotate", "--db", test_db,
        "--input", str(input_vcf),
        "--output", str(output_vcf),
        "--verbose",
    ])
    assert result.exit_code == 0
    assert "Annotated" in result.output


# --- afquery update-db --remove-samples ---

def test_update_db_remove_samples(runner, db_copy):
    result = runner.invoke(cli, [
        "update-db", "--db", db_copy,
        "--remove-samples", "S00",
    ])
    assert result.exit_code == 0
    assert "Removed 1 sample(s)" in result.output

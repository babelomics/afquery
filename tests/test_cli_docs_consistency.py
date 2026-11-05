"""
Regression tests that enforce consistency between the CLI and documentation.
These tests run in-process using Click's test runner (no subprocess overhead).
"""
import dataclasses
import pytest
from click.testing import CliRunner
from afquery.cli import cli
from afquery.models import QueryResult


runner = CliRunner()


def _help(args: list[str]) -> str:
    result = runner.invoke(cli, args + ["--help"])
    return result.output


# --- CLI option presence / absence ---

def test_create_db_has_no_include_all_filters():
    """--include-all-filters is not yet exposed; docs must not reference it."""
    assert "--include-all-filters" not in _help(["create-db"])


def test_create_db_has_build_memory():
    assert "--build-memory" in _help(["create-db"])


def test_create_db_has_build_threads():
    assert "--build-threads" in _help(["create-db"])


def test_query_has_no_verbose():
    """query command intentionally has no --verbose flag."""
    assert "--verbose" not in _help(["query"])


def test_query_has_region():
    assert "--region" in _help(["query"])


def test_query_has_no_start_end():
    """--start/--end belong to dump, not query."""
    help_text = _help(["query"])
    assert "--start" not in help_text
    assert "--end" not in help_text


@pytest.mark.parametrize("command", ["query", "annotate", "dump"])
def test_fail_as_ref_flag(command):
    assert "--fail-as-ref" in _help([command])


def test_dump_has_start_end():
    help_text = _help(["dump"])
    assert "--start" in help_text
    assert "--end" in help_text


# --- Manifest column names ---

def test_manifest_required_cols_constant():
    """required_cols in manifest.py must contain tech_name and phenotype_codes."""
    import pathlib
    manifest_src = pathlib.Path(
        "src/afquery/preprocess/manifest.py"
    ).read_text()
    assert "tech_name" in manifest_src
    assert "phenotype_codes" in manifest_src


# --- Python API return types ---

def test_query_result_n_fail_is_int_not_optional():
    """N_FAIL must be int with default 0, never Optional."""
    fields = {f.name: f for f in dataclasses.fields(QueryResult)}
    assert "N_FAIL" in fields
    assert fields["N_FAIL"].default == 0


def test_query_result_n_fail_default_zero():
    """Default QueryResult must have N_FAIL=0."""
    from afquery.models import VariantKey
    r = QueryResult(
        variant=VariantKey(chrom="chr1", pos=1, ref="A", alt="T"),
        AC=0, AN=0, AF=0.0, n_samples_eligible=0,
        N_HET=0, N_HOM_ALT=0, N_HOM_REF=0,
    )
    assert r.N_FAIL == 0
    assert r.N_FAIL is not None
    assert r.N_HEMI == 0


# --- Genotype classes in every output, and the docs that describe them ---

import json
import re
from pathlib import Path

DOCS = Path(__file__).resolve().parent.parent / "docs"


def test_query_outputs_report_n_hemi(test_db):
    args = ["query", "--db", test_db, "--locus", "chrX:5000000"]
    tsv = runner.invoke(cli, args + ["--format", "tsv"]).output.splitlines()
    header = tsv[0].split("\t")
    assert header[header.index("N_HOM_ALT") + 1] == "N_HEMI"
    row = dict(zip(header, tsv[1].split("\t")))
    assert (row["N_HET"], row["N_HOM_ALT"], row["N_HEMI"]) == ("1", "0", "1")

    [entry] = json.loads(runner.invoke(cli, args + ["--format", "json"]).output)
    assert (entry["N_HOM_ALT"], entry["N_HEMI"]) == (0, 1)

    text = runner.invoke(cli, args).output
    assert "N_HOM_ALT=0  N_HEMI=1  N_HOM_REF=2" in text


def test_annotate_help_lists_n_hemi():
    assert "AFQUERY_N_HEMI" in _help(["annotate"])


def test_dump_guide_base_columns_match_dump(test_db):
    """The base column line in the dump guide is the header dump writes."""
    import io
    from afquery import Database
    buf = io.StringIO()
    Database(test_db).dump(output=buf, chrom="chr1")
    header = buf.getvalue().splitlines()[0]
    guide = (DOCS / "guides" / "dump-export.md").read_text()
    assert f"```\n{header}\n```" in guide


def test_annotate_guide_lists_every_info_field(test_db, tmp_path):
    """Every AFQUERY_* INFO field annotate adds is in the annotate guide's table."""
    import cyvcf2
    vcf_in = tmp_path / "in.vcf"
    vcf_in.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1500\t.\tA\tT\t.\tPASS\t.\n"
    )
    out = tmp_path / "out.vcf"
    assert runner.invoke(cli, ["annotate", "--db", test_db, "--input", str(vcf_in),
                               "--output", str(out)]).exit_code == 0
    added = set(re.findall(r"ID=(AFQUERY_\w+)", cyvcf2.VCF(str(out)).raw_header))
    documented = set(re.findall(r"^\| `(AFQUERY_\w+)`", (DOCS / "guides" / "annotate-vcf.md").read_text(), re.M))
    assert added == documented

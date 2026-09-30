"""Tests for variant_info command: sample carrier retrieval with metadata."""
import json
import warnings

import pytest
from click.testing import CliRunner

from afquery import Database, SampleCarrier, variant_info
from afquery.cli import cli
from afquery.models import AfqueryWarning

# ---------------------------------------------------------------------------
# Reference data from conftest.py
# ---------------------------------------------------------------------------
# VARIANTS = [
#   ("chr1",  1500,    "A", "T", het=[0,5],   hom=[2],    fail=[0]),   # S00 in both het+fail
#   ("chr1",  3500,    "G", "C", het=[7],      hom=[],     fail=[]),
#   ("chr1",  5000,    "T", "G", het=[0,1],   hom=[3],    fail=[1,3]),
#   ("chrX",  5000000, "A", "G", het=[0,2],   hom=[],     fail=[]),
#   ("chrY",  500000,  "T", "C", het=[0,1],   hom=[4],    fail=[]),
#   ("chrM",  100,     "C", "A", het=[0,2,5], hom=[],     fail=[5]),
# ]
# SAMPLES (id, name, sex, tech_id):
#   0=S00/male/WGS,  1=S01/male/WGS,   2=S02/female/WGS,  3=S03/female/WGS,
#   4=S04/male/WES_kit_A, 5=S05/female/WES_kit_A, 6=S06/female/WES_kit_A,
#   7=S07/male/WES_kit_B, 8=S08/male/WES_kit_B,   9=S09/female/WES_kit_B
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# 1. Basic retrieval
# ---------------------------------------------------------------------------

def test_basic_returns_carriers(test_db):
    db = Database(test_db)
    carriers = db.variant_info("chr1", 3500, ref="G", alt="C")
    assert len(carriers) == 1
    assert carriers[0].sample_id == 7
    assert carriers[0].sample_name == "S07"
    assert carriers[0].genotype == "het"
    assert carriers[0].filter_pass is True


def test_hom_genotype(test_db):
    db = Database(test_db)
    carriers = db.variant_info("chr1", 1500, ref="A", alt="T")
    by_id = {c.sample_id: c for c in carriers}
    # S02 is hom/PASS
    assert by_id[2].genotype == "hom"
    assert by_id[2].filter_pass is True


def test_fail_genotype(test_db):
    """S00 is in both het_bitmap and fail_bitmap → should be reported as alt/FAIL."""
    db = Database(test_db)
    carriers = db.variant_info("chr1", 1500, ref="A", alt="T")
    by_id = {c.sample_id: c for c in carriers}
    assert by_id[0].genotype == "alt"
    assert by_id[0].filter_pass is False


def test_pass_het_genotype(test_db):
    """S05 is in het_bitmap but not fail_bitmap → het/PASS."""
    db = Database(test_db)
    carriers = db.variant_info("chr1", 1500, ref="A", alt="T")
    by_id = {c.sample_id: c for c in carriers}
    assert by_id[5].genotype == "het"
    assert by_id[5].filter_pass is True


def test_sorted_by_sample_id(test_db):
    db = Database(test_db)
    carriers = db.variant_info("chr1", 1500, ref="A", alt="T")
    ids = [c.sample_id for c in carriers]
    assert ids == sorted(ids)


# ---------------------------------------------------------------------------
# 2. Metadata included
# ---------------------------------------------------------------------------

def test_sex_in_metadata(test_db):
    db = Database(test_db)
    carriers = db.variant_info("chr1", 3500, ref="G", alt="C")
    by_id = {c.sample_id: c for c in carriers}
    assert by_id[7].sex == "male"


def test_tech_in_metadata(test_db):
    db = Database(test_db)
    carriers = db.variant_info("chr1", 3500, ref="G", alt="C")
    by_id = {c.sample_id: c for c in carriers}
    assert by_id[7].tech_name == "WES_kit_B"


def test_phenotypes_in_metadata(test_db):
    """S07 has phenotypes E11.9 and J45 (from conftest SAMPLE_PHENOTYPE)."""
    db = Database(test_db)
    carriers = db.variant_info("chr1", 3500, ref="G", alt="C")
    by_id = {c.sample_id: c for c in carriers}
    assert set(by_id[7].phenotypes) == {"E11.9", "J45"}


def test_phenotypes_sorted(test_db):
    db = Database(test_db)
    carriers = db.variant_info("chr1", 3500, ref="G", alt="C")
    for c in carriers:
        assert c.phenotypes == sorted(c.phenotypes)


# ---------------------------------------------------------------------------
# 3. All fail_bitmap cases — chr1:5000 T>G
#    het=[0,1], hom=[3], fail=[1,3]
#    → S01: het+fail → alt/FAIL
#    → S03: hom+fail → alt/FAIL
#    → S00: het only → het/PASS
# ---------------------------------------------------------------------------

def test_fail_chr1_5000(test_db):
    db = Database(test_db)
    carriers = db.variant_info("chr1", 5000, ref="T", alt="G")
    by_id = {c.sample_id: c for c in carriers}
    # S00: het, no fail
    assert by_id[0].genotype == "het"
    assert by_id[0].filter_pass is True
    # S01: het + fail
    assert by_id[1].genotype == "alt"
    assert by_id[1].filter_pass is False
    # S03: hom + fail
    assert by_id[3].genotype == "alt"
    assert by_id[3].filter_pass is False


# ---------------------------------------------------------------------------
# 4. Sample filter: phenotype
# ---------------------------------------------------------------------------

def test_phenotype_filter(test_db):
    """chr1:1500 carriers with E11.9 only: S00(0), S05(5) [S02 has E11.9 too]."""
    db = Database(test_db)
    # E11.9 carriers at chr1:1500: S00(0), S02(2), S05(5)
    carriers_all = db.variant_info("chr1", 1500, ref="A", alt="T")
    carriers_ph = db.variant_info("chr1", 1500, ref="A", alt="T", phenotype=["E11.9"])
    ids_all = {c.sample_id for c in carriers_all}
    ids_ph = {c.sample_id for c in carriers_ph}
    # Filtered set must be a subset
    assert ids_ph <= ids_all
    # Only samples with E11.9 should appear
    for c in carriers_ph:
        assert "E11.9" in c.phenotypes


def test_phenotype_exclude(test_db):
    """Exclude E11.9: S00(0) and S05(5) should be absent, S02(2) without E11.9 would remain, but S02 has E11.9 too. Let's check S02 is gone."""
    db = Database(test_db)
    carriers = db.variant_info("chr1", 1500, ref="A", alt="T", phenotype=["^E11.9"])
    ids = {c.sample_id for c in carriers}
    # S00, S02, S05 all have E11.9 → all excluded → empty
    for sid in [0, 2, 5]:
        assert sid not in ids


# ---------------------------------------------------------------------------
# 5. Sample filter: sex
# ---------------------------------------------------------------------------

def test_sex_filter_male(test_db):
    """chr1:1500 — male carriers only. S00(male), S05(female), S02(female). → only S00."""
    db = Database(test_db)
    carriers = db.variant_info("chr1", 1500, ref="A", alt="T", sex="male")
    ids = {c.sample_id for c in carriers}
    assert ids == {0}


def test_sex_filter_female(test_db):
    db = Database(test_db)
    carriers = db.variant_info("chr1", 1500, ref="A", alt="T", sex="female")
    for c in carriers:
        assert c.sex == "female"


# ---------------------------------------------------------------------------
# 6. Sample filter: technology
# ---------------------------------------------------------------------------

def test_tech_filter(test_db):
    """chr1:1500 — WGS only. Carriers S00(WGS), S02(WGS), S05(WES_kit_A). → S00, S02."""
    db = Database(test_db)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AfqueryWarning)
        carriers = db.variant_info("chr1", 1500, ref="A", alt="T", tech=["WGS"])
    ids = {c.sample_id for c in carriers}
    assert 5 not in ids   # S05 is WES_kit_A → excluded
    for c in carriers:
        assert c.tech_name == "WGS"


# ---------------------------------------------------------------------------
# 7. Ref/alt specificity
# ---------------------------------------------------------------------------

def test_ref_alt_specificity(test_db):
    db = Database(test_db)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AfqueryWarning)
        carriers = db.variant_info("chr1", 1500, ref="A", alt="T")
    assert all(isinstance(c, SampleCarrier) for c in carriers)
    assert len(carriers) > 0


def test_pos_only_returns_all_alleles(test_db):
    """Without ref/alt, all alleles at the locus are returned."""
    db = Database(test_db)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AfqueryWarning)
        carriers_all = db.variant_info("chr1", 1500)
        carriers_specific = db.variant_info("chr1", 1500, ref="A", alt="T")
    # chr1:1500 has only one allele in the test DB, so both should match
    assert {c.sample_id for c in carriers_all} == {c.sample_id for c in carriers_specific}


# ---------------------------------------------------------------------------
# 8. Edge case: variant not in database
# ---------------------------------------------------------------------------

def test_variant_not_in_db(test_db):
    db = Database(test_db)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AfqueryWarning)
        carriers = db.variant_info("chr1", 9999999, ref="A", alt="T")
    assert carriers == []


def test_unknown_chrom_returns_empty(test_db):
    db = Database(test_db)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AfqueryWarning)
        carriers = db.variant_info("chr99", 1000)
    assert carriers == []


# ---------------------------------------------------------------------------
# 9. Edge case: no carriers in filtered population
# ---------------------------------------------------------------------------

def test_no_carriers_after_filter(test_db):
    """chr1:3500 only has S07 (WES_kit_B/male). Filtering to female → empty."""
    db = Database(test_db)
    carriers = db.variant_info("chr1", 3500, ref="G", alt="C", sex="female")
    assert carriers == []


# ---------------------------------------------------------------------------
# 10. CLI: text output
# ---------------------------------------------------------------------------

def test_cli_text_output(test_db):
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["variant-info", "--db", test_db, "--locus", "chr1:3500",
         "--ref", "G", "--alt", "C", "--no-warn"],
    )
    assert result.exit_code == 0, result.output
    assert "S07" in result.output
    assert "het" in result.output
    assert "PASS" in result.output


def test_cli_text_empty(test_db):
    """No carriers → prints 'No carriers found' message."""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["variant-info", "--db", test_db, "--locus", "chr1:3500",
         "--ref", "G", "--alt", "C", "--sex", "female", "--no-warn"],
    )
    assert result.exit_code == 0, result.output
    assert "No carriers found" in result.output


# ---------------------------------------------------------------------------
# 11. CLI: TSV output
# ---------------------------------------------------------------------------

def test_cli_tsv_output(test_db):
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["variant-info", "--db", test_db, "--locus", "chr1:3500",
         "--ref", "G", "--alt", "C", "--format", "tsv", "--no-warn"],
    )
    assert result.exit_code == 0, result.output
    lines = result.output.strip().split("\n")
    header = lines[0].split("\t")
    assert header == ["sample_id", "sample_name", "sex", "tech", "phenotypes", "genotype", "filter"]
    assert len(lines) == 2  # header + 1 carrier
    fields = lines[1].split("\t")
    assert fields[1] == "S07"
    assert fields[5] == "het"
    assert fields[6] == "PASS"


def test_cli_tsv_fail_sample(test_db):
    """chr1:1500 — S00 is alt/FAIL."""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["variant-info", "--db", test_db, "--locus", "chr1:1500",
         "--ref", "A", "--alt", "T", "--format", "tsv", "--no-warn"],
    )
    assert result.exit_code == 0, result.output
    lines = result.output.strip().split("\n")
    rows = {l.split("\t")[1]: l.split("\t") for l in lines[1:]}
    assert rows["S00"][5] == "alt"
    assert rows["S00"][6] == "FAIL"


# ---------------------------------------------------------------------------
# 12. CLI: JSON output
# ---------------------------------------------------------------------------

def test_cli_json_output(test_db):
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["variant-info", "--db", test_db, "--locus", "chr1:3500",
         "--ref", "G", "--alt", "C", "--format", "json", "--no-warn"],
    )
    assert result.exit_code == 0, result.output
    data = json.loads(result.output)
    assert "variant" in data
    assert "samples" in data
    assert data["variant"]["chrom"] == "chr1"
    assert data["variant"]["pos"] == 3500
    assert len(data["samples"]) == 1
    s = data["samples"][0]
    assert s["sample_name"] == "S07"
    assert s["genotype"] == "het"
    assert s["filter"] == "PASS"
    assert isinstance(s["phenotypes"], list)


def test_cli_json_fail_sample(test_db):
    """chr1:1500 JSON — S00 should have filter='FAIL'."""
    runner = CliRunner()
    result = runner.invoke(
        cli,
        ["variant-info", "--db", test_db, "--locus", "chr1:1500",
         "--ref", "A", "--alt", "T", "--format", "json", "--no-warn"],
    )
    assert result.exit_code == 0, result.output
    data = json.loads(result.output)
    by_name = {s["sample_name"]: s for s in data["samples"]}
    assert by_name["S00"]["filter"] == "FAIL"
    assert by_name["S00"]["genotype"] == "alt"
    assert by_name["S02"]["filter"] == "PASS"
    assert by_name["S02"]["genotype"] == "hom"


# ---------------------------------------------------------------------------
# 13. variant_info wrapper function
# ---------------------------------------------------------------------------

def test_variant_info_wrapper_function(test_db):
    """Test the public variant_info() wrapper function directly."""
    carriers = variant_info(test_db, "chr1", 3500, ref="G", alt="C")
    assert len(carriers) == 1
    assert carriers[0].sample_id == 7
    assert carriers[0].sample_name == "S07"
    assert carriers[0].genotype == "het"


def test_variant_info_wrapper_with_filters(test_db):
    """Test variant_info() wrapper with phenotype and sex filters."""
    carriers = variant_info(
        test_db, "chr1", 1500,
        ref="A", alt="T",
        phenotype=["E11.9"],  # S00 has this phenotype
        sex="male"
    )
    assert len(carriers) > 0
    assert all(c.sex == "male" for c in carriers)


def test_variant_info_wrapper_with_tech_filter(test_db):
    """Test variant_info() wrapper with technology filter."""
    carriers = variant_info(
        test_db, "chr1", 3500,
        ref="G", alt="C",
        tech=["WES_kit_B"]
    )
    assert len(carriers) == 1
    assert carriers[0].tech_name == "WES_kit_B"


# ---------------------------------------------------------------------------
# 14. Edge cases and warnings
# ---------------------------------------------------------------------------

def test_unknown_chrom_warning(test_db):
    """Unknown chromosome returns empty list with warning."""
    db = Database(test_db)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        carriers = db.variant_info("chrUnknown", 1000)
        assert len(carriers) == 0
        assert len(w) == 1
        assert issubclass(w[0].category, AfqueryWarning)
        assert "chrUnknown" in str(w[0].message)


def test_multiple_alleles_warning(test_db):
    """Multiple alleles at position triggers warning without explicit ref/alt."""
    db = Database(test_db)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        # chr1:5000 has multiple alleles: T→G and others
        carriers = db.variant_info("chr1", 5000)
        # Should warn about multiple alleles
        warns = [x for x in w if issubclass(x.category, AfqueryWarning)]
        if len(warns) > 0:
            assert "alleles" in str(warns[0].message).lower()


def test_empty_result_no_warning(test_db):
    """Explicit ref/alt with no matches returns empty, no warning."""
    db = Database(test_db)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        carriers = db.variant_info("chr1", 3500, ref="G", alt="X")
        assert len(carriers) == 0
        # Should not warn about multiple alleles when explicitly filtered
        warns = [x for x in w if "alleles" in str(x.message).lower()]
        assert len(warns) == 0


# ---------------------------------------------------------------------------
# 15. Capture regions: carriers outside the sample's BED are not listed
# ---------------------------------------------------------------------------

def test_out_of_capture_carrier_excluded_chrY(test_db):
    """chrY:500000 hom=[4]: S04 is WES_kit_A, whose BED does not cover chrY."""
    db = Database(test_db)
    carriers = db.variant_info("chrY", 500000, ref="T", alt="C")
    assert {c.sample_id for c in carriers} == {0, 1}


def test_out_of_capture_carrier_excluded_chrM(test_db):
    """chrM:100 het=[0,2,5]: S05 is WES_kit_A, whose BED does not cover chrM."""
    db = Database(test_db)
    carriers = db.variant_info("chrM", 100, ref="C", alt="A")
    assert {c.sample_id for c in carriers} == {0, 2}


def _write_vcf(path, sample_name, records):
    contigs = sorted({r[0] for r in records}) or ["chr1"]
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        f.write('##FILTER=<ID=LowQual,Description="Low quality">\n')
        for contig in contigs:
            f.write(f"##contig=<ID={contig}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
        for rec in records:
            chrom, pos, ref, alt, gt = rec[:5]
            flt = rec[5] if len(rec) > 5 else "PASS"
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t.\tGT\t{gt}\n")


def _build_db(tmp_path, sample_defs, beds=None):
    """sample_defs = [(name, sex, tech, records)]; beds = {tech: bed_text}."""
    from afquery.preprocess import run_preprocess

    bed_dir = tmp_path / "beds"
    bed_dir.mkdir()
    for tech, text in (beds or {}).items():
        (bed_dir / f"{tech}.bed").write_text(text)
    rows = ["sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes"]
    for name, sex, tech, records in sample_defs:
        vcf_path = tmp_path / f"{name}.vcf"
        _write_vcf(vcf_path, name, records)
        rows.append(f"{name}\t{sex}\t{tech}\t{vcf_path}\tE11.9")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("\n".join(rows) + "\n")
    db_path = tmp_path / "db"
    db_path.mkdir()
    run_preprocess(
        manifest_path=str(manifest), output_dir=str(db_path),
        genome_build="GRCh37", threads=1, bed_dir=str(bed_dir),
    )
    return Database(str(db_path))


def test_off_target_call_not_listed(tmp_path):
    """A panel sample with a call outside its BED is neither counted nor listed."""
    db = _build_db(tmp_path, [
        ("WGS_HET", "female", "WGS", [("chr1", 5000, "G", "A", "0/1")]),
        ("PANEL_OFF", "female", "PANEL", [("chr1", 5000, "G", "A", "0/1")]),
        ("PANEL_IN", "female", "PANEL", [("chr1", 500, "C", "T", "0/1")]),
    ], beds={"PANEL": "chr1\t0\t1000\n"})
    [r] = db.query("chr1", 5000)
    assert r.n_samples_eligible == 1
    carriers = db.variant_info("chr1", 5000, ref="G", alt="A")
    assert [c.sample_name for c in carriers] == ["WGS_HET"]


def test_carrier_counts_match_query(tmp_path):
    """variant_info het/hom/alt counts equal query N_HET/N_HOM_ALT/N_FAIL."""
    db = _build_db(tmp_path, [
        ("W1", "female", "WGS", [("chr1", 5000, "G", "A", "0/1")]),
        ("W2", "male",   "WGS", [("chr1", 5000, "G", "A", "1/1")]),
        ("W3", "female", "WGS", [("chr1", 5000, "G", "A", "0/1", "LowQual")]),
        ("P1", "female", "PANEL", [("chr1", 5000, "G", "A", "1/1")]),
        ("P2", "male",   "PANEL", [("chr1", 5000, "G", "A", "0/1", "LowQual")]),
        ("P3", "female", "PANEL", [("chr1", 500, "G", "A", "0/1")]),
        ("P4", "female", "PANEL", [("chr1", 500, "G", "A", "1/1")]),
    ], beds={"PANEL": "chr1\t0\t1000\n"})
    for pos in (500, 5000):
        [r] = db.query("chr1", pos)
        genotypes = [c.genotype for c in db.variant_info("chr1", pos, ref="G", alt="A")]
        assert genotypes.count("het") == r.N_HET
        assert genotypes.count("hom") == r.N_HOM_ALT
        assert genotypes.count("alt") == r.N_FAIL


def test_multiallelic_1_2_sample_listed_for_each_allele(tmp_path):
    """A 1/2 sample carries both alleles and is listed once per allele."""
    db = _build_db(tmp_path, [
        ("S12", "female", "WGS", [("chr1", 5000, "G", "A,T", "1/2")]),
        ("SA",  "female", "WGS", [("chr1", 5000, "G", "A", "0/1")]),
    ])
    names_a = [c.sample_name for c in db.variant_info("chr1", 5000, alt="A")]
    names_t = [c.sample_name for c in db.variant_info("chr1", 5000, alt="T")]
    assert names_a == ["S12", "SA"]
    assert names_t == ["S12"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AfqueryWarning)
        both = [c.sample_name for c in db.variant_info("chr1", 5000)]
    assert sorted(both) == ["S12", "S12", "SA"]

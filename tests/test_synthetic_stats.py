"""End-to-end tests: haploid genotype statistics with synthetic VCFs.

Builds real databases via run_preprocess, then queries to verify all 8 metrics.
"""
from pathlib import Path

import pytest

from afquery.database import Database
from afquery.preprocess import run_preprocess


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _write_vcf(path: Path, sample_name: str, records: list[tuple]) -> None:
    """Write a minimal VCF.  records = [(chrom, pos, ref, alt, gt[, filter_str]), ...]"""
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
        for rec in records:
            chrom, pos, ref, alt, gt = rec[:5]
            flt = rec[5] if len(rec) > 5 else "PASS"
            f.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t{flt}\t.\tGT\t{gt}\n")


def _build_db(tmp_path: Path, sample_defs: list[tuple], **kw) -> Database:
    """Build a test database.
    sample_defs = [(name, sex, records), ...]
    """
    rows = ["sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes"]
    for name, sex, records in sample_defs:
        vcf_path = tmp_path / f"{name}.vcf"
        _write_vcf(vcf_path, name, records)
        rows.append(f"{name}\t{sex}\twgs\t{vcf_path}\tE11.9")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text("\n".join(rows) + "\n")

    db_path = tmp_path / "db"
    db_path.mkdir()
    run_preprocess(
        manifest_path=str(manifest),
        output_dir=str(db_path),
        genome_build="GRCh37",
        threads=1,
        **kw,
    )
    return Database(str(db_path))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestSyntheticDiploid:

    def test_basic_diploid_autosome(self, tmp_path):
        """3 samples with GT 0/0, 0/1, 1/1 on chr1."""
        samples = [
            ("S0", "male",   [("chr1", 1000, "A", "T", "0/0")]),
            ("S1", "female", [("chr1", 1000, "A", "T", "0/1")]),
            ("S2", "male",   [("chr1", 1000, "A", "T", "1/1")]),
        ]
        db = _build_db(tmp_path, samples)
        results = db.query(chrom="chr1", pos=1000)
        assert len(results) == 1
        r = results[0]
        assert r.N_HET == 1         # S1
        assert r.N_HOM_ALT == 1     # S2
        assert r.N_HOM_REF == 1     # S0
        assert r.AC == 3            # 1 + 2
        assert r.AN == 6            # 3 * 2
        assert r.n_samples_eligible == 3

    def test_missing_variant_is_hom_ref(self, tmp_path):
        """WGS sample with no variant at position → counted as hom-ref."""
        samples = [
            ("S0", "male",   [("chr1", 1000, "A", "T", "0/1")]),
            ("S1", "female", [("chr1", 1000, "A", "T", "0/0")]),
            ("S2", "male",   []),  # no variant at chr1:1000 → 0/0
        ]
        db = _build_db(tmp_path, samples)
        results = db.query(chrom="chr1", pos=1000)
        assert len(results) == 1
        r = results[0]
        assert r.N_HET == 1         # S0
        assert r.N_HOM_ALT == 0
        assert r.N_HOM_REF == 2     # S1 + S2
        assert r.AC == 1
        assert r.n_samples_eligible == 3


class TestSyntheticHaploid:

    def test_haploid_chrX_males(self, tmp_path):
        """Males with haploid GT=0 and GT=1 on chrX non-PAR."""
        samples = [
            ("S0", "male",   [("chrX", 5000000, "A", "G", "0")]),  # ref
            ("S1", "male",   [("chrX", 5000000, "A", "G", "1")]),  # alt
        ]
        db = _build_db(tmp_path, samples)
        results = db.query(chrom="chrX", pos=5000000)
        assert len(results) == 1
        r = results[0]
        assert r.N_HET == 0         # no diploid samples
        assert r.N_HOM_ALT == 0
        assert r.N_HEMI == 1        # S1 hemizygous ALT
        assert r.N_HOM_REF == 1     # S0 hemizygous REF
        assert r.AC == 1
        assert r.AN == 2            # 2 males * 1

    def test_mixed_ploidy_chrX(self, tmp_path):
        """Male (haploid) + female (diploid) on chrX non-PAR."""
        samples = [
            ("S0", "male",   [("chrX", 5000000, "A", "G", "1")]),    # haploid ALT
            ("S1", "female", [("chrX", 5000000, "A", "G", "0/1")]),  # diploid het
        ]
        db = _build_db(tmp_path, samples)
        results = db.query(chrom="chrX", pos=5000000)
        assert len(results) == 1
        r = results[0]
        assert r.N_HET == 1         # S1 diploid het
        assert r.N_HOM_ALT == 0
        assert r.N_HEMI == 1        # S0 hemizygous ALT
        assert r.AC == 2            # 1 (haploid) + 1 (het)
        assert r.AN == 3            # 1 (male) + 2 (female)
        assert r.n_samples_eligible == 2

    def test_mixed_ploidy_chrX_with_hom(self, tmp_path):
        """Male haploid ALT + female hom-ALT on chrX non-PAR."""
        samples = [
            ("S0", "male",   [("chrX", 5000000, "A", "G", "1")]),    # haploid ALT
            ("S1", "female", [("chrX", 5000000, "A", "G", "1/1")]),  # diploid hom-alt
        ]
        db = _build_db(tmp_path, samples)
        results = db.query(chrom="chrX", pos=5000000)
        assert len(results) == 1
        r = results[0]
        assert r.N_HET == 0
        assert r.N_HOM_ALT == 1     # S1 diploid hom
        assert r.N_HEMI == 1        # S0 hemizygous
        assert r.AC == 3            # 1 (haploid) + 2 (hom)
        assert r.AN == 3

    def test_chrM_all_haploid(self, tmp_path):
        """chrM: all samples haploid regardless of sex."""
        samples = [
            ("S0", "male",   [("chrM", 100, "C", "A", "1")]),
            ("S1", "female", [("chrM", 100, "C", "A", "1")]),
            ("S2", "male",   [("chrM", 100, "C", "A", "0")]),
        ]
        db = _build_db(tmp_path, samples)
        results = db.query(chrom="chrM", pos=100)
        assert len(results) == 1
        r = results[0]
        assert r.N_HET == 0         # all haploid
        assert r.N_HOM_ALT == 2     # S0, S1: chrM carriers are not hemizygous
        assert r.N_HEMI == 0
        assert r.N_HOM_REF == 1     # S2
        assert r.AC == 2
        assert r.AN == 3            # 3 * 1


class TestSyntheticFilter:

    def test_filter_fail_with_haploid(self, tmp_path):
        """FILTER!=PASS on chrX haploid → tracked as N_FAIL, not N_HEMI."""
        samples = [
            ("S0", "male", [("chrX", 5000000, "A", "G", "1", "PASS")]),
            ("S1", "male", [("chrX", 5000000, "A", "G", "1", "LowQual")]),
        ]
        db = _build_db(tmp_path, samples)
        results = db.query(chrom="chrX", pos=5000000)
        assert len(results) == 1
        r = results[0]
        # S0 passes → haploid ALT → N_HEMI
        # S1 fails → N_FAIL (not a carrier for AC)
        assert r.AC == 1
        assert r.N_FAIL == 1
        assert r.N_HEMI == 1        # S0 only
        assert r.N_HOM_ALT == 0


class TestSyntheticHaploidHeterozygous:
    """A male has one allele on chrX non-PAR and chrY: a 0/1 call there cannot be
    read as a haploid genotype and is a failed call. chrM is exempt: a 0/1 there
    is heteroplasmy and the sample is a carrier."""

    @pytest.mark.parametrize("chrom,pos", [("chrX", 5000000), ("chrY", 2700000)])
    def test_male_het_is_failed_call(self, tmp_path, chrom, pos):
        samples = [
            ("M_HET", "male", [(chrom, pos, "A", "G", "0/1")]),
            ("M_ALT", "male", [(chrom, pos, "A", "G", "1")]),
            ("M_REF", "male", []),
        ]
        db = _build_db(tmp_path, samples)
        [r] = db.query(chrom=chrom, pos=pos)
        assert r.N_FAIL == 1          # M_HET
        assert r.N_HEMI == 1          # M_ALT
        assert r.N_HET == 0
        assert r.AC == 1              # M_ALT only
        assert r.AN == 2              # M_ALT + M_REF, one allele each
        assert r.n_samples_eligible == 3
        [r] = db.query(chrom=chrom, pos=pos, fail_as_ref=True)
        assert (r.AC, r.AN) == (1, 3)

    def test_male_hom_call_counts_one_allele(self, tmp_path):
        """1/1 names a single allele: the male carries it, once."""
        samples = [("M", "male", [("chrX", 5000000, "A", "G", "1/1")])]
        db = _build_db(tmp_path, samples)
        [r] = db.query(chrom="chrX", pos=5000000)
        assert (r.AC, r.AN, r.N_FAIL) == (1, 1, 0)
        assert (r.N_HEMI, r.N_HOM_ALT) == (1, 0)

    def test_male_het_in_chrX_par_is_diploid(self, tmp_path):
        samples = [("M", "male", [("chrX", 100000, "A", "G", "0/1")])]
        db = _build_db(tmp_path, samples)
        [r] = db.query(chrom="chrX", pos=100000)
        assert (r.AC, r.AN, r.N_HET, r.N_FAIL) == (1, 2, 1, 0)
        assert r.N_HEMI == 0          # the PAR is diploid: no hemizygotes

    @pytest.mark.parametrize("sex", ["male", "female"])
    def test_chrM_het_is_a_carrier(self, tmp_path, sex):
        samples = [("S", sex, [("chrM", 100, "C", "A", "0/1")])]
        db = _build_db(tmp_path, samples)
        [r] = db.query(chrom="chrM", pos=100)
        assert (r.AC, r.AN, r.N_FAIL) == (1, 1, 0)
        assert (r.N_HOM_ALT, r.N_HEMI, r.N_HET) == (1, 0, 0)


class TestMissingGtFailFilter:

    def test_missing_gt_at_failed_site_counts_as_n_fail(self, tmp_path):
        """GT=./. at FILTER≠PASS → N_FAIL=1, not counted in AC."""
        samples = [
            ("S0", "male",   [("chr1", 1000, "A", "T", "0/1", "PASS")]),
            ("S1", "female", [("chr1", 1000, "A", "T", "./.", "LowQual")]),
        ]
        db = _build_db(tmp_path, samples)
        results = db.query(chrom="chr1", pos=1000)
        assert len(results) == 1
        r = results[0]
        assert r.N_FAIL == 1, f"Expected N_FAIL=1 (missing GT at failed site), got {r.N_FAIL}"
        assert r.AC == 1        # only S0 counts
        assert r.N_HET == 1     # S0
        assert r.N_HOM_REF == 0
        assert r.n_samples_eligible == 2

    def test_missing_gt_pass_does_not_count_as_n_fail(self, tmp_path):
        """GT=./. at FILTER=PASS → still no row, not tracked."""
        samples = [
            ("S0", "male",   [("chr1", 1000, "A", "T", "0/1", "PASS")]),
            ("S1", "female", [("chr1", 1000, "A", "T", "./.", "PASS")]),
        ]
        db = _build_db(tmp_path, samples)
        results = db.query(chrom="chr1", pos=1000)
        assert len(results) == 1
        r = results[0]
        assert r.N_FAIL == 0
        assert r.N_HOM_REF == 1  # S1 counted as hom-ref

    def test_missing_gt_failed_with_haploid(self, tmp_path):
        """GT=. (haploid) at FILTER≠PASS on chrX → N_FAIL=1."""
        samples = [
            ("S0", "male", [("chrX", 5000000, "A", "G", "1", "PASS")]),
            ("S1", "male", [("chrX", 5000000, "A", "G", ".", "LowQual")]),
        ]
        db = _build_db(tmp_path, samples)
        results = db.query(chrom="chrX", pos=5000000)
        assert len(results) == 1
        r = results[0]
        assert r.N_FAIL == 1
        assert r.AC == 1
        assert r.N_HEMI == 1  # S0 hemizygous ALT


class TestSyntheticInvariants:

    def test_all_metrics_consistent(self, tmp_path):
        """Comprehensive check: all 8 metrics consistent across ploidy types."""
        samples = [
            # Autosome
            ("S0", "male",   [("chr1", 1000, "A", "T", "0/1"),
                              ("chrX", 5000000, "C", "G", "1")]),
            ("S1", "female", [("chr1", 1000, "A", "T", "1/1"),
                              ("chrX", 5000000, "C", "G", "0/1")]),
            ("S2", "male",   [("chr1", 1000, "A", "T", "0/0"),
                              ("chrX", 5000000, "C", "G", "0")]),
        ]
        db = _build_db(tmp_path, samples)

        for chrom, pos in [("chr1", 1000), ("chrX", 5000000)]:
            results = db.query(chrom=chrom, pos=pos)
            for r in results:
                # Genotype sum = eligible
                total = r.N_HET + r.N_HOM_ALT + r.N_HEMI + r.N_HOM_REF + (r.N_FAIL or 0)
                assert total == r.n_samples_eligible, \
                    f"Sum mismatch at {chrom}:{pos}"
                # AC <= AN
                assert r.AC <= r.AN, f"AC > AN at {chrom}:{pos}"

    def test_batch_region_consistency(self, tmp_path):
        """query, query_batch, query_region all return same stats."""
        samples = [
            ("S0", "male",   [("chrX", 5000000, "A", "G", "1")]),
            ("S1", "female", [("chrX", 5000000, "A", "G", "0/1")]),
        ]
        db = _build_db(tmp_path, samples)

        point = db.query(chrom="chrX", pos=5000000)
        batch = db.query_batch(chrom="chrX", variants=[(5000000, "A", "G")])
        region = db.query_region(chrom="chrX", start=5000000, end=5000000)

        assert len(point) == len(batch) == len(region) == 1
        for method_result in [batch[0], region[0]]:
            assert point[0].AC == method_result.AC
            assert point[0].AN == method_result.AN
            assert point[0].N_HET == method_result.N_HET
            assert point[0].N_HOM_ALT == method_result.N_HOM_ALT
            assert point[0].N_HEMI == method_result.N_HEMI
            assert point[0].N_HOM_REF == method_result.N_HOM_REF

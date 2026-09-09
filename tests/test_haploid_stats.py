"""Tests for ploidy-aware N_HET / N_HOM_ALT classification.

Haploid ALT samples (chrX non-PAR males, chrY, chrM) should be counted as
N_HOM_ALT, not N_HET, even though they sit in het_bitmap (gt_ac=1).
"""
from afquery import Database


class TestHaploidPointQuery:

    def test_chrX_nonpar_mixed_ploidy(self, test_db):
        """chrX:5000000 het_ids=[0(male),2(female)] → N_HET=1, N_HOM_ALT=1."""
        db = Database(test_db)
        results = db.query(chrom="chrX", pos=5000000)
        assert len(results) == 1
        r = results[0]
        # S00 (male, haploid) → N_HOM_ALT; S02 (female, diploid) → N_HET
        assert r.N_HET == 1, f"Expected N_HET=1 (diploid het S02), got {r.N_HET}"
        assert r.N_HOM_ALT == 1, f"Expected N_HOM_ALT=1 (haploid ALT S00), got {r.N_HOM_ALT}"
        assert r.AC == 2

    def test_chrY_all_haploid(self, test_db):
        """chrY:500000 het_ids=[0,1] — all haploid → N_HET=0."""
        db = Database(test_db)
        results = db.query(chrom="chrY", pos=500000)
        assert len(results) == 1
        r = results[0]
        # S00,S01 in het_bm, haploid → both N_HOM_ALT
        # S04 (male WES_A) not eligible (WES_A has no chrY coverage)
        assert r.N_HET == 0, f"Expected N_HET=0 (chrY all haploid), got {r.N_HET}"
        assert r.N_HOM_ALT == 2, f"Expected N_HOM_ALT=2, got {r.N_HOM_ALT}"
        assert r.AC == 2

    def test_chrM_all_haploid(self, test_db):
        """chrM:100 het_ids=[0,2,5] — all haploid → N_HET=0."""
        db = Database(test_db)
        results = db.query(chrom="chrM", pos=100)
        assert len(results) == 1
        r = results[0]
        # Eligible: WGS{0,1,2,3} (S05 is WES_A, no chrM coverage)
        # het_elig = {0,2}, all haploid on chrM → N_HOM_ALT
        assert r.N_HET == 0, f"Expected N_HET=0 (chrM all haploid), got {r.N_HET}"
        assert r.N_HOM_ALT == 2, f"Expected N_HOM_ALT=2, got {r.N_HOM_ALT}"

    def test_autosome_diploid_regression(self, test_db):
        """chr1:1500 (autosome) — all diploid, unchanged behavior."""
        db = Database(test_db)
        results = db.query(chrom="chr1", pos=1500)
        assert len(results) == 1
        r = results[0]
        # Eligible: WGS{0,1,2,3} ∪ WES_A{4,5,6} (chr1:1500 covered by both)
        # het_ids=[0,5], hom_ids=[2], fail_ids=[0]
        # All diploid → N_HET = 2 (S00,S05), N_HOM_ALT = 1 (S02)
        assert r.N_HET == 2, f"Expected N_HET=2, got {r.N_HET}"
        assert r.N_HOM_ALT == 1, f"Expected N_HOM_ALT=1, got {r.N_HOM_ALT}"
        assert r.AC == 4  # 2 het + 2*1 hom


class TestHaploidInvariants:

    def test_genotype_sum_equals_eligible(self, test_db):
        """N_HET + N_HOM_ALT + N_HOM_REF + N_FAIL + N_NO_COVERAGE == n_eligible for all variants."""
        db = Database(test_db)
        for chrom, pos in [("chr1", 1500), ("chr1", 3500), ("chr1", 5000),
                           ("chrX", 5000000), ("chrY", 500000), ("chrM", 100)]:
            for r in db.query(chrom=chrom, pos=pos):
                total = (
                    r.N_HET + r.N_HOM_ALT + r.N_HOM_REF
                    + (r.N_FAIL or 0) + (r.N_NO_COVERAGE or 0)
                )
                assert total == r.n_samples_eligible, \
                    f"Sum mismatch at {chrom}:{pos}: {total} != {r.n_samples_eligible}"

    def test_ac_le_an(self, test_db):
        """AC <= AN for all variants."""
        db = Database(test_db)
        for chrom, pos in [("chr1", 1500), ("chr1", 3500), ("chr1", 5000),
                           ("chrX", 5000000), ("chrY", 500000), ("chrM", 100)]:
            for r in db.query(chrom=chrom, pos=pos):
                assert r.AC <= r.AN, f"AC ({r.AC}) > AN ({r.AN}) at {chrom}:{pos}"

    def test_autosome_ac_formula(self, test_db):
        """Autosome: AC == N_HET + 2 * N_HOM_ALT (all diploid, ignoring fail)."""
        db = Database(test_db)
        for chrom, pos in [("chr1", 1500), ("chr1", 3500), ("chr1", 5000)]:
            for r in db.query(chrom=chrom, pos=pos):
                expected_ac = r.N_HET + 2 * r.N_HOM_ALT
                assert r.AC == expected_ac, \
                    f"AC mismatch at {chrom}:{pos}: {r.AC} != {expected_ac}"


class TestHaploidSexFiltered:

    def test_chrX_male_only(self, test_db):
        """chrX:5000000 male-only → all haploid → N_HET=0."""
        db = Database(test_db)
        results = db.query(chrom="chrX", pos=5000000, sex="male")
        assert len(results) == 1
        r = results[0]
        assert r.N_HET == 0, f"Expected N_HET=0 (male-only haploid), got {r.N_HET}"
        assert r.N_HOM_ALT == 1  # S00

    def test_chrX_female_only(self, test_db):
        """chrX:5000000 female-only → all diploid → N_HOM_ALT=0."""
        db = Database(test_db)
        results = db.query(chrom="chrX", pos=5000000, sex="female")
        assert len(results) == 1
        r = results[0]
        assert r.N_HET == 1  # S02
        assert r.N_HOM_ALT == 0, f"Expected N_HOM_ALT=0 (female-only diploid), got {r.N_HOM_ALT}"


class TestHaploidBatchAndRegion:

    def test_query_batch_matches_point(self, test_db):
        """query_batch returns same N_HET/N_HOM_ALT as point query."""
        db = Database(test_db)
        point = db.query(chrom="chrX", pos=5000000)
        batch = db.query_batch(chrom="chrX", variants=[(5000000, "A", "G")])
        assert len(point) == len(batch) == 1
        assert point[0].N_HET == batch[0].N_HET
        assert point[0].N_HOM_ALT == batch[0].N_HOM_ALT
        assert point[0].AC == batch[0].AC

    def test_query_region_matches_point(self, test_db):
        """query_region returns same N_HET/N_HOM_ALT as point query."""
        db = Database(test_db)
        point = db.query(chrom="chrX", pos=5000000)
        region = db.query_region(chrom="chrX", start=5000000, end=5000000)
        assert len(point) == len(region) == 1
        assert point[0].N_HET == region[0].N_HET
        assert point[0].N_HOM_ALT == region[0].N_HOM_ALT
        assert point[0].AC == region[0].AC

    def test_batch_chrY_haploid(self, test_db):
        """query_batch on chrY also gives N_HET=0."""
        db = Database(test_db)
        batch = db.query_batch(chrom="chrY", variants=[(500000, "T", "C")])
        assert len(batch) == 1
        assert batch[0].N_HET == 0

    def test_region_chrM_haploid(self, test_db):
        """query_region on chrM also gives N_HET=0."""
        db = Database(test_db)
        region = db.query_region(chrom="chrM", start=100, end=100)
        assert len(region) == 1
        assert region[0].N_HET == 0


class TestChrYEligibility:
    """Females hold no chrY, so they are neither carriers nor homozygous reference there."""

    def test_chrY_excludes_females_from_eligible(self, test_db):
        db = Database(test_db)
        r = db.query(chrom="chrY", pos=500000)[0]

        # The fixture is 5 male and 5 female WGS/WES samples; only the males
        # that the capture covers can be genotyped on chrY.
        assert r.n_samples_eligible == r.AN, (
            "n_samples_eligible counts samples with no chrY: it must agree with AN, "
            "which has always been male-only on this chromosome"
        )

    def test_chrY_counts_partition_the_eligible_set(self, test_db):
        """The per-genotype counts must add up to the eligible set on every chromosome."""
        db = Database(test_db)
        for chrom, pos in (("chr1", 1500), ("chrX", 5000000), ("chrY", 500000), ("chrM", 100)):
            for r in db.query(chrom=chrom, pos=pos):
                total = r.N_HET + r.N_HOM_ALT + r.N_HOM_REF + r.N_FAIL + r.N_NO_COVERAGE
                assert total == r.n_samples_eligible, f"{chrom}:{pos} counts do not partition"

    def test_chrY_no_phantom_hom_ref(self, test_db):
        """S00 and S01 carry chrY:500000; nobody else eligible is homozygous reference."""
        db = Database(test_db)
        r = db.query(chrom="chrY", pos=500000)[0]
        assert r.N_HOM_ALT == 2
        assert r.N_HOM_REF == 0

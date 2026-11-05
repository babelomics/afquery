import pytest
from pyroaring import BitMap
from afquery.ploidy import compute_AN, is_par

# Sample setup: IDs 0-9, males={0,1,4,7,8}, females={2,3,5,6,9}
ALL_10 = BitMap(range(10))
MALES = BitMap([0, 1, 4, 7, 8])
FEMALES = BitMap([2, 3, 5, 6, 9])


def test_is_par_chrX_in_par1_grch37():
    assert is_par("chrX", 100_000, "GRCh37") is True


def test_is_par_chrX_outside_par_grch37():
    assert is_par("chrX", 5_000_000, "GRCh37") is False


def test_is_par_autosome():
    assert is_par("chr1", 1_000_000, "GRCh37") is False


def test_is_par_chrY_in_par1_grch37():
    assert is_par("chrY", 50_000, "GRCh37") is True


def test_is_par_boundary_start():
    # GRCh37 chrX PAR1: 60001-2699520
    assert is_par("chrX", 60_001, "GRCh37") is True


def test_is_par_boundary_end():
    assert is_par("chrX", 2_699_520, "GRCh37") is True


def test_is_par_just_outside_boundary():
    assert is_par("chrX", 60_000, "GRCh37") is False
    assert is_par("chrX", 2_699_521, "GRCh37") is False


# --- compute_AN tests ---

def test_autosome_all_samples():
    AN = compute_AN(ALL_10, MALES, FEMALES, "chr1", 1_000_000, "GRCh37")
    assert AN == 20


def test_autosome_subset():
    subset = BitMap([0, 1, 2, 3, 4])  # 5 samples
    AN = compute_AN(subset, MALES, FEMALES, "chr1", 1_000_000, "GRCh37")
    assert AN == 10


def test_chrX_par_all():
    # PAR1 GRCh37: 60001-2699520
    AN = compute_AN(ALL_10, MALES, FEMALES, "chrX", 100_000, "GRCh37")
    assert AN == 20


def test_chrX_non_par_all():
    # pos 5_000_000 is non-PAR
    # 5 females * 2 + 5 males * 1 = 15
    AN = compute_AN(ALL_10, MALES, FEMALES, "chrX", 5_000_000, "GRCh37")
    assert AN == 15


def test_chrX_non_par_females_only():
    AN = compute_AN(FEMALES, MALES, FEMALES, "chrX", 5_000_000, "GRCh37")
    assert AN == 10  # 5 females * 2


def test_chrX_non_par_males_only():
    AN = compute_AN(MALES, MALES, FEMALES, "chrX", 5_000_000, "GRCh37")
    assert AN == 5  # 5 males * 1


def test_chrY_all():
    # Only males contribute, ploidy=1
    AN = compute_AN(ALL_10, MALES, FEMALES, "chrY", 2_700_000, "GRCh37")
    assert AN == 5  # 5 males


def test_chrY_females_only():
    AN = compute_AN(FEMALES, MALES, FEMALES, "chrY", 2_700_000, "GRCh37")
    assert AN == 0


def test_chrM_all():
    AN = compute_AN(ALL_10, MALES, FEMALES, "chrM", 100, "GRCh37")
    assert AN == 10


def test_chrM_subset():
    subset = BitMap([0, 2, 5])
    AN = compute_AN(subset, MALES, FEMALES, "chrM", 100, "GRCh37")
    assert AN == 3


def test_chrom_normalization():
    # "1" should normalize to "chr1"
    AN = compute_AN(ALL_10, MALES, FEMALES, "1", 1_000_000, "GRCh37")
    assert AN == 20


def test_empty_eligible():
    empty = BitMap()
    assert compute_AN(empty, MALES, FEMALES, "chr1", 1_000_000, "GRCh37") == 0
    assert compute_AN(empty, MALES, FEMALES, "chrX", 5_000_000, "GRCh37") == 0
    assert compute_AN(empty, MALES, FEMALES, "chrY", 2_700_000, "GRCh37") == 0
    assert compute_AN(empty, MALES, FEMALES, "chrM", 100, "GRCh37") == 0


def test_phenotype_subset_eligible():
    # Subset: {0,1,2} — 2 males (0,1), 1 female (2)
    subset = BitMap([0, 1, 2])
    AN = compute_AN(subset, MALES, FEMALES, "chrX", 5_000_000, "GRCh37")
    assert AN == 4  # 1 female*2 + 2 males*1


# --- chrY PARs: counted at their chrX coordinates, nobody genotyped on chrY ---

from afquery.ploidy import sample_ploidy, split_ploidy


@pytest.mark.parametrize("build,pos", [
    ("GRCh37", 10_001), ("GRCh37", 500_000), ("GRCh37", 2_649_520), ("GRCh37", 59_034_050),
    ("GRCh38", 10_001), ("GRCh38", 2_781_479), ("GRCh38", 57_000_000),
])
def test_chrY_par_has_no_alleles(build, pos):
    assert compute_AN(ALL_10, MALES, FEMALES, "chrY", pos, build) == 0
    haploid, diploid = split_ploidy(ALL_10, MALES, FEMALES, "chrY", pos, build)
    assert len(haploid) == len(diploid) == 0
    assert sample_ploidy("chrY", pos, "male", build) == 0
    assert sample_ploidy("chrY", pos, "female", build) == 0


@pytest.mark.parametrize("build,pos", [
    ("GRCh37", 2_649_521), ("GRCh37", 2_700_000), ("GRCh38", 2_781_480),
])
def test_chrY_outside_par_is_haploid_in_males(build, pos):
    assert compute_AN(ALL_10, MALES, FEMALES, "chrY", pos, build) == 5
    haploid, diploid = split_ploidy(ALL_10, MALES, FEMALES, "chrY", pos, build)
    assert haploid == MALES and len(diploid) == 0
    assert sample_ploidy("chrY", pos, "male", build) == 1
    assert sample_ploidy("chrY", pos, "female", build) == 0


def test_chrX_par_stays_diploid_for_everyone():
    """The PAR loci are counted on chrX, where males carry two copies."""
    assert sample_ploidy("chrX", 100_000, "male", "GRCh37") == 2
    haploid, diploid = split_ploidy(ALL_10, MALES, FEMALES, "chrX", 100_000, "GRCh37")
    assert len(haploid) == 0 and diploid == ALL_10


# --- classify_carriers: one source for counts and variant-info labels ---

from afquery.ploidy import classify_carriers


def _classes(chrom, pos, het, hom, build="GRCh37"):
    haploid, diploid = split_ploidy(ALL_10, MALES, FEMALES, chrom, pos, build)
    return classify_carriers(chrom, haploid, diploid, BitMap(het), BitMap(hom))


def test_classify_autosome_by_genotype():
    het, hom, hemi = _classes("chr1", 1000, het=[0, 2], hom=[4])
    assert (list(het), list(hom), list(hemi)) == ([0, 2], [4], [])


def test_classify_chrX_nonpar_males_are_hemizygous():
    # males 0 (GT 1, het bitmap) and 1 (GT 1/1, hom bitmap); female 2 het, 3 hom
    het, hom, hemi = _classes("chrX", 5_000_000, het=[0, 2], hom=[1, 3])
    assert (list(het), list(hom), list(hemi)) == ([2], [3], [0, 1])


def test_classify_chrX_par_is_diploid():
    het, hom, hemi = _classes("chrX", 100_000, het=[0, 2], hom=[1])
    assert (list(het), list(hom), list(hemi)) == ([0, 2], [1], [])


def test_classify_chrY_nonpar_is_hemizygous():
    het, hom, hemi = _classes("chrY", 2_700_000, het=[0], hom=[4])
    assert (list(het), list(hom), list(hemi)) == ([], [], [0, 4])


def test_classify_chrY_par_has_no_carriers():
    het, hom, hemi = _classes("chrY", 500_000, het=[0], hom=[4])
    assert len(het) == len(hom) == len(hemi) == 0


def test_classify_chrM_carriers_are_hom_alt_not_hemizygous():
    het, hom, hemi = _classes("chrM", 100, het=[0, 2], hom=[5])
    assert (list(het), list(hom), list(hemi)) == ([], [0, 2, 5], [])

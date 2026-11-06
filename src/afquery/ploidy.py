from pyroaring import BitMap
from .constants import PAR


def is_par(chrom: str, pos: int, genome_build: str) -> bool:
    """Return True if pos is inside a PAR region on chrX or chrY. Caller must pass normalized chrom."""
    if chrom not in ("chrX", "chrY"):
        return False
    for start, end in PAR[genome_build].get(chrom, []):
        if start <= pos <= end:
            return True
    return False


def sample_ploidy(chrom: str, pos: int, sex: str, genome_build: str) -> int:
    """Number of alleles one sample carries at a position.

    Caller must pass a normalized chrom. Same rules as split_ploidy. It is 0 where
    the sample cannot be genotyped: chrY in a female, and the chrY PARs for
    everyone, because PAR loci are counted at their chrX coordinates.
    """
    if chrom == "chrM":
        return 1
    if chrom == "chrY":
        if is_par(chrom, pos, genome_build):
            return 0
        return 1 if sex == "male" else 0
    if chrom == "chrX" and not is_par(chrom, pos, genome_build):
        return 1 if sex == "male" else 2
    return 2


def compute_AN(
    eligible: BitMap,
    male_bitmap: BitMap,
    female_bitmap: BitMap,
    chrom: str,
    pos: int,
    genome_build: str,
) -> int:
    """Compute ploidy-aware AN for a set of eligible samples.

    Caller must pass a normalized chrom (e.g. 'chr1', 'chrX').

    Rules:
    - Autosomes: AN = 2 × len(eligible)
    - chrM: AN = len(eligible)
    - chrY non-PAR: only males contribute, ploidy=1 → AN = len(eligible & male_bitmap)
    - chrY PAR: AN = 0; PAR loci are counted at their chrX coordinates
    - chrX PAR: AN = 2 × len(eligible)
    - chrX non-PAR: females=2, males=1 → AN = 2×females + 1×males
    """

    if chrom == "chrM":
        return len(eligible)

    if chrom == "chrY":
        if is_par(chrom, pos, genome_build):
            return 0
        return len(eligible & male_bitmap)

    if chrom == "chrX":
        if is_par(chrom, pos, genome_build):
            return 2 * len(eligible)
        else:
            females = eligible & female_bitmap
            males = eligible & male_bitmap
            return 2 * len(females) + len(males)

    # Autosomes
    return 2 * len(eligible)


def split_ploidy(
    eligible: BitMap,
    male_bitmap: BitMap,
    female_bitmap: BitMap,
    chrom: str,
    pos: int,
    genome_build: str,
) -> tuple[BitMap, BitMap]:
    """Divide eligible into (haploid_eligible, diploid_eligible).

    Caller must pass a normalized chrom (e.g. 'chr1', 'chrX').

    Same rules as compute_AN:
    - chrM: all haploid
    - chrY non-PAR: only eligible males, all haploid
    - chrY PAR: nobody; PAR loci are counted at their chrX coordinates
    - chrX non-PAR: males haploid, females diploid
    - rest (autosomes, chrX PAR): all diploid
    """
    if chrom == "chrM":
        return eligible, BitMap()

    if chrom == "chrY":
        if is_par(chrom, pos, genome_build):
            return BitMap(), BitMap()
        return eligible & male_bitmap, BitMap()

    if chrom == "chrX" and not is_par(chrom, pos, genome_build):
        return eligible & male_bitmap, eligible & female_bitmap

    return BitMap(), eligible  # autosomes and chrX PAR: all diploid


def classify_carriers(
    chrom: str,
    haploid: BitMap,
    diploid: BitMap,
    het: BitMap,
    hom: BitMap,
) -> tuple[BitMap, BitMap, BitMap]:
    """Split the PASS carriers of one allele into (het, hom_alt, hemi).

    ``haploid`` and ``diploid`` come from split_ploidy; ``het`` and ``hom`` are
    the samples whose GT holds one or two copies of the allele.

    - Diploid samples: heterozygous or homozygous alt, as the GT says.
    - Haploid samples carry a single copy whatever the GT says (``1`` and ``1/1``
      alike). On chrX outside the PARs and on chrY they are hemizygous. chrM is
      haploid too but not hemizygous: it is a multi-copy genome, AFQuery records
      no heteroplasmy, and its carriers are counted as homozygous alt.

    Caller must pass a normalized chrom. Counts and variant-info labels both
    come from here, so they always agree.
    """
    het_dip = het & diploid
    hom_dip = hom & diploid
    hap_carriers = (het | hom) & haploid
    if chrom == "chrM":
        return het_dip, hom_dip | hap_carriers, BitMap()
    return het_dip, hom_dip, hap_carriers

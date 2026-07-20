import pytest
from afquery.constants import (
    ALL_CHROMS,
    normalize_chrom,
    is_autosome,
    is_sex_chrom,
    is_mito,
)


def test_normalize_bare_number():
    assert normalize_chrom("1") == "chr1"


def test_normalize_already_prefixed():
    assert normalize_chrom("chr1") == "chr1"


def test_normalize_X():
    assert normalize_chrom("X") == "chrX"


def test_normalize_MT():
    assert normalize_chrom("MT") == "chrM"


def test_normalize_whitespace():
    assert normalize_chrom(" chrX ") == "chrX"


def test_normalize_chrMT():
    # Regression: 'chrMT' used to survive as 'chrMT'. A capture BED using it kept
    # valid chr1..chr22 keys alongside, so the "matches no known chromosome" warning
    # never fired and mitochondrial positions were silently uncovered.
    assert normalize_chrom("chrMT") == "chrM"


def test_normalize_mito_aliases():
    for alias in ("M", "m", "mt", "chrm", "chrMt"):
        assert normalize_chrom(alias) == "chrM"


def test_normalize_lowercase_sex_chroms():
    assert normalize_chrom("x") == "chrX"
    assert normalize_chrom("chrx") == "chrX"
    assert normalize_chrom("y") == "chrY"


def test_normalize_uppercase_prefix():
    assert normalize_chrom("CHR1") == "chr1"


def test_normalize_preserves_unplaced_contigs():
    # Unplaced and alt contigs must pass through unchanged and stay outside
    # ALL_CHROMS, so the build stage keeps dropping them instead of guessing.
    assert normalize_chrom("GL000209.1") == "chrGL000209.1"
    assert normalize_chrom("chr1_KI270706v1_random") == "chr1_KI270706v1_random"
    assert normalize_chrom("chrUn_GL000220v1") == "chrUn_GL000220v1"
    assert "chrGL000209.1" not in ALL_CHROMS
    assert "chr1_KI270706v1_random" not in ALL_CHROMS


def test_normalize_is_idempotent_on_all_chroms():
    # normalize_chrom output is the on-disk Parquet partition name. A name that
    # changed under normalization would orphan already-built databases.
    for chrom in ALL_CHROMS:
        assert normalize_chrom(chrom) == chrom


def test_is_autosome_true():
    assert is_autosome("chr1") is True
    assert is_autosome("22") is True


def test_is_autosome_false():
    assert is_autosome("chrX") is False
    assert is_autosome("chrM") is False


def test_is_sex_chrom_true():
    assert is_sex_chrom("chrX") is True
    assert is_sex_chrom("Y") is True


def test_is_sex_chrom_false():
    assert is_sex_chrom("chr1") is False


def test_is_mito_true():
    assert is_mito("chrM") is True
    assert is_mito("MT") is True


def test_is_mito_accepts_chrMT():
    assert is_mito("chrMT") is True


def test_is_mito_false():
    assert is_mito("chr1") is False

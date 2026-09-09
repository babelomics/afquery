"""Tests for variant Parquet layout resolution."""
import pytest

from afquery import storage


def _touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")


@pytest.fixture
def variants(tmp_path):
    d = tmp_path / "variants"
    d.mkdir()
    return d


# ---------------------------------------------------------------------------
# bucket_id
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("pos,expected", [
    (0, 0),
    (1, 0),
    (999_999, 0),
    (1_000_000, 1),
    (1_000_001, 1),
    (2_500_000, 2),
    (887_801, 0),
])
def test_bucket_id_boundaries(pos, expected):
    assert storage.bucket_id(pos) == expected


# ---------------------------------------------------------------------------
# chrom_layout
# ---------------------------------------------------------------------------

def test_chrom_layout_partitioned(variants):
    _touch(variants / "chr1" / "bucket_0.parquet")
    assert storage.chrom_layout(variants, "chr1") == storage.PARTITIONED


def test_chrom_layout_flat(variants):
    _touch(variants / "chr1.parquet")
    assert storage.chrom_layout(variants, "chr1") == storage.FLAT


def test_chrom_layout_absent_uses_default(variants):
    assert storage.chrom_layout(variants, "chr9", storage.FLAT) == storage.FLAT
    assert storage.chrom_layout(variants, "chr9", storage.PARTITIONED) == storage.PARTITIONED


def test_chrom_layout_defaults_to_partitioned(variants):
    assert storage.chrom_layout(variants, "chr9") == storage.PARTITIONED


def test_chrom_layout_bucket_directory_wins_over_flat(variants):
    """The readers prefer the bucket directory, so the writer must too."""
    _touch(variants / "chr1.parquet")
    _touch(variants / "chr1" / "bucket_0.parquet")
    assert storage.chrom_layout(variants, "chr1") == storage.PARTITIONED


# ---------------------------------------------------------------------------
# detect_layout
# ---------------------------------------------------------------------------

def test_detect_layout_empty_dir_is_partitioned(variants):
    """A database with no variants yet must grow buckets, not a second layout."""
    assert storage.detect_layout(variants) == storage.PARTITIONED


def test_detect_layout_missing_dir_is_partitioned(tmp_path):
    assert storage.detect_layout(tmp_path / "nope") == storage.PARTITIONED


def test_detect_layout_flat(variants):
    _touch(variants / "chr1.parquet")
    _touch(variants / "chr2.parquet")
    assert storage.detect_layout(variants) == storage.FLAT


def test_detect_layout_partitioned(variants):
    _touch(variants / "chr1" / "bucket_0.parquet")
    assert storage.detect_layout(variants) == storage.PARTITIONED


def test_detect_layout_ignores_empty_chrom_dir(variants):
    """A stray empty directory must not flip the database-wide answer."""
    (variants / "chr1").mkdir()
    _touch(variants / "chr2.parquet")
    assert storage.detect_layout(variants) == storage.FLAT


# ---------------------------------------------------------------------------
# non-canonical contig names
# ---------------------------------------------------------------------------

def test_partitioned_chroms_keeps_non_canonical_contigs(variants):
    """Unplaced and alt contigs are real data; filtering them here would hide it."""
    _touch(variants / "chr1" / "bucket_0.parquet")
    _touch(variants / "chrGL000209.1" / "bucket_0.parquet")
    _touch(variants / "chr1_KI270706v1_random" / "bucket_0.parquet")
    assert storage.partitioned_chroms(variants) == {
        "chr1", "chrGL000209.1", "chr1_KI270706v1_random",
    }


def test_non_canonical_contig_resolves(variants):
    _touch(variants / "chrGL000209.1" / "bucket_0.parquet")
    assert storage.chrom_layout(variants, "chrGL000209.1") == storage.PARTITIONED
    assert storage.variant_parquet_glob(variants, "chrGL000209.1") is not None


def test_flat_chroms(variants):
    _touch(variants / "chr1.parquet")
    _touch(variants / "chrX.parquet")
    (variants / "chr2").mkdir()
    assert storage.flat_chroms(variants) == {"chr1", "chrX"}


def test_stored_chroms_spans_both_layouts(variants):
    """The Phase 2 recompute walks this set, so it must miss no chromosome."""
    _touch(variants / "chr1" / "bucket_0.parquet")
    _touch(variants / "chr2" / "bucket_3.parquet")
    _touch(variants / "chrX.parquet")
    assert storage.stored_chroms(variants) == {"chr1", "chr2", "chrX"}


def test_stored_chroms_empty_database(variants):
    assert storage.stored_chroms(variants) == set()


# ---------------------------------------------------------------------------
# mixed_layout_chroms
# ---------------------------------------------------------------------------

def test_mixed_layout_chroms_detects_both(variants):
    _touch(variants / "chr1.parquet")
    _touch(variants / "chr1" / "bucket_0.parquet")
    _touch(variants / "chr2" / "bucket_0.parquet")
    _touch(variants / "chrX.parquet")
    assert storage.mixed_layout_chroms(variants) == ["chr1"]


def test_mixed_layout_chroms_clean_db(variants):
    _touch(variants / "chr1" / "bucket_0.parquet")
    _touch(variants / "chr2" / "bucket_0.parquet")
    assert storage.mixed_layout_chroms(variants) == []


def test_mixed_layout_chroms_missing_dir(tmp_path):
    assert storage.mixed_layout_chroms(tmp_path / "nope") == []


# ---------------------------------------------------------------------------
# existing_bucket_ids
# ---------------------------------------------------------------------------

def test_existing_bucket_ids_sorted_numerically(variants):
    for b in (0, 2, 10, 1):
        _touch(variants / "chr1" / f"bucket_{b}.parquet")
    assert storage.existing_bucket_ids(variants, "chr1") == [0, 1, 2, 10]


def test_existing_bucket_ids_empty_for_flat(variants):
    _touch(variants / "chr1.parquet")
    assert storage.existing_bucket_ids(variants, "chr1") == []


# ---------------------------------------------------------------------------
# path resolution
# ---------------------------------------------------------------------------

def test_variant_parquet_for_pos_partitioned(variants):
    _touch(variants / "chr1" / "bucket_2.parquet")
    got = storage.variant_parquet_for_pos(variants, "chr1", 2_500_000)
    assert got == variants / "chr1" / "bucket_2.parquet"


def test_variant_parquet_for_pos_missing_bucket(variants):
    _touch(variants / "chr1" / "bucket_0.parquet")
    assert storage.variant_parquet_for_pos(variants, "chr1", 5_000_000) is None


def test_variant_parquet_for_pos_flat(variants):
    _touch(variants / "chr1.parquet")
    assert storage.variant_parquet_for_pos(variants, "chr1", 500) == variants / "chr1.parquet"


def test_variant_parquet_for_pos_absent(variants):
    assert storage.variant_parquet_for_pos(variants, "chr1", 500) is None


def test_variant_parquet_glob_partitioned(variants):
    _touch(variants / "chr1" / "bucket_0.parquet")
    assert storage.variant_parquet_glob(variants, "chr1").endswith("chr1/bucket_*.parquet")


def test_variant_parquet_glob_empty_chrom_dir(variants):
    (variants / "chr1").mkdir()
    assert storage.variant_parquet_glob(variants, "chr1") is None


def test_variant_parquet_glob_flat(variants):
    _touch(variants / "chr1.parquet")
    assert storage.variant_parquet_glob(variants, "chr1").endswith("chr1.parquet")


# ---------------------------------------------------------------------------
# iter_variant_parquets
# ---------------------------------------------------------------------------

def test_iter_variant_parquets_covers_both_layouts(variants):
    _touch(variants / "chr1.parquet")
    _touch(variants / "chr2" / "bucket_0.parquet")
    _touch(variants / "chr2" / "bucket_1.parquet")
    got = {p.name for p in storage.iter_variant_parquets(variants)}
    assert got == {"chr1.parquet", "bucket_0.parquet", "bucket_1.parquet"}


def test_iter_variant_parquets_skips_non_bucket_files(variants):
    _touch(variants / "chr1" / "bucket_0.parquet")
    _touch(variants / "chr1" / "notes.txt")
    _touch(variants / "chr1" / "bucket_0.parquet.tmp")
    got = [p.name for p in storage.iter_variant_parquets(variants)]
    assert got == ["bucket_0.parquet"]


def test_iter_variant_parquets_missing_dir(tmp_path):
    assert list(storage.iter_variant_parquets(tmp_path / "nope")) == []

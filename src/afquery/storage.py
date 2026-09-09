"""Resolution of on-disk paths in the variant Parquet store.

Two layouts exist. The **partitioned** layout stores one file per 1 Mbp bucket
under a per-chromosome directory::

    variants/chr1/bucket_0.parquet      # positions 0 .. 999,999
    variants/chr1/bucket_1.parquet      # positions 1,000,000 .. 1,999,999

The **flat** layout stores one file per chromosome::

    variants/chr1.parquet

``create-db`` has produced the partitioned layout since it became the default,
so every database built by the normal path is bucketed. The flat layout survives
in small hand-built databases and in test fixtures, and is still readable.

Readers have always preferred the bucket directory when a chromosome somehow has
both. That rule used to be reimplemented in every module that opened a Parquet
file — the query engine, the dump and annotate workers, the benchmark helper,
the compactor and the updater — which is how the updater came to disagree with
the readers and write rows into a file no query would ever open. Keeping the
rule in one module is what stops writers and readers from drifting apart again.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterator

BUCKET_SIZE = 1_000_000

PARTITIONED = "partitioned"
FLAT = "flat"


def bucket_id(pos: int) -> int:
    """Bucket owning `pos`."""
    return pos // BUCKET_SIZE


def bucket_path(variants_dir: Path | str, chrom: str, bucket: int) -> Path:
    """Path of one bucket file, whether or not it exists."""
    return Path(variants_dir) / chrom / f"bucket_{bucket}.parquet"


def flat_path(variants_dir: Path | str, chrom: str) -> Path:
    """Path of a chromosome's flat file, whether or not it exists."""
    return Path(variants_dir) / f"{chrom}.parquet"


def partitioned_chroms(variants_dir: Path | str) -> set[str]:
    """Chromosomes stored as a bucket directory.

    Every subdirectory counts, with no check against the canonical chromosome
    list: a database may hold unplaced or alt contigs whose names normalize to
    bodies like 'chrGL000209.1', and filtering them here would make their data
    unreadable.
    """
    variants_dir = Path(variants_dir)
    if not variants_dir.exists():
        return set()
    return {p.name for p in variants_dir.iterdir() if p.is_dir()}


def flat_chroms(variants_dir: Path | str) -> set[str]:
    """Chromosomes stored as a single flat Parquet file."""
    variants_dir = Path(variants_dir)
    if not variants_dir.exists():
        return set()
    return {p.stem for p in variants_dir.iterdir() if p.is_file() and p.suffix == ".parquet"}


def stored_chroms(variants_dir: Path | str) -> set[str]:
    """Every chromosome the database holds, in either layout."""
    return partitioned_chroms(variants_dir) | flat_chroms(variants_dir)


def detect_layout(variants_dir: Path | str) -> str:
    """Layout of the database as a whole.

    A subdirectory only votes 'partitioned' if it actually holds a bucket file,
    so a stray empty directory cannot flip the answer. An empty or missing
    variants directory answers 'partitioned', which is what create-db produces —
    that way a database with no variants yet grows buckets instead of sprouting
    a second layout the readers would ignore.
    """
    variants_dir = Path(variants_dir)
    if variants_dir.is_dir():
        for entry in variants_dir.iterdir():
            if entry.is_dir() and any(entry.glob("bucket_*.parquet")):
                return PARTITIONED
        if any(p.suffix == ".parquet" for p in variants_dir.iterdir() if p.is_file()):
            return FLAT
    return PARTITIONED


def chrom_layout(variants_dir: Path | str, chrom: str, default: str = PARTITIONED) -> str:
    """Layout of one chromosome; `default` decides for a chromosome new to the database.

    The bucket directory wins when both are present, matching what the readers do.
    """
    variants_dir = Path(variants_dir)
    if (variants_dir / chrom).is_dir():
        return PARTITIONED
    if flat_path(variants_dir, chrom).exists():
        return FLAT
    return default


def mixed_layout_chroms(variants_dir: Path | str) -> list[str]:
    """Chromosomes that have both a bucket directory and a flat file.

    A database should never be in this state. When it is, the flat file is dead
    weight: queries read only the bucket directory, so any sample whose calls
    live in the flat file is counted as homozygous reference everywhere.
    """
    variants_dir = Path(variants_dir)
    if not variants_dir.is_dir():
        return []
    return sorted(
        p.name for p in variants_dir.iterdir()
        if p.is_dir() and flat_path(variants_dir, p.name).exists()
    )


def existing_bucket_ids(variants_dir: Path | str, chrom: str) -> list[int]:
    """Bucket ids already written for `chrom`, ascending. Empty when not bucketed."""
    # The chromosome is a directory name, not part of the pattern: contig names
    # can legitimately contain glob metacharacters (GRCh38 spells HLA contigs
    # HLA-A*01:01:01:01), and this is reached with names read straight off disk.
    chrom_dir = Path(variants_dir) / chrom
    if not chrom_dir.is_dir():
        return []
    ids: list[int] = []
    for p in chrom_dir.glob("bucket_*.parquet"):
        stem = p.stem[len("bucket_"):]
        if stem.isdigit():
            ids.append(int(stem))
    return sorted(ids)


def variant_parquet_for_pos(variants_dir: Path | str, chrom: str, pos: int) -> Path | None:
    """File holding `pos` for a point query, or None when there is no data."""
    variants_dir = Path(variants_dir)
    if (variants_dir / chrom).is_dir():
        p = bucket_path(variants_dir, chrom, bucket_id(pos))
        return p if p.exists() else None
    flat = flat_path(variants_dir, chrom)
    return flat if flat.exists() else None


def variant_parquet_glob(variants_dir: Path | str, chrom: str) -> str | None:
    """Path or glob pattern covering a whole chromosome, or None when there is no data."""
    variants_dir = Path(variants_dir)
    chrom_dir = variants_dir / chrom
    if chrom_dir.is_dir():
        return str(chrom_dir / "bucket_*.parquet") if any(
            chrom_dir.glob("bucket_*.parquet")
        ) else None
    flat = flat_path(variants_dir, chrom)
    return str(flat) if flat.exists() else None


def iter_variant_parquets(variants_dir: Path | str) -> Iterator[Path]:
    """Every variant Parquet file in the store, flat files first, then buckets."""
    variants_dir = Path(variants_dir)
    if not variants_dir.is_dir():
        return
    yield from sorted(variants_dir.glob("*.parquet"))
    for chrom_dir in sorted(p for p in variants_dir.iterdir() if p.is_dir()):
        yield from sorted(chrom_dir.glob("bucket_*.parquet"))


def ensure_parent(path: Path | str) -> None:
    """Create the directory holding `path`, including a per-chromosome bucket directory."""
    os.makedirs(os.path.dirname(str(path)), exist_ok=True)

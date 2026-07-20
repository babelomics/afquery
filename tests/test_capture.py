import pickle
import tempfile
from pathlib import Path

import pytest

from afquery.capture import CaptureIndex

DATA_DIR = Path(__file__).parent / "data"
BED_A = str(DATA_DIR / "beds" / "wes_kit_a.bed")  # chr1:999-2000, chrX:999-2000
BED_B = str(DATA_DIR / "beds" / "wes_kit_b.bed")  # chr1:2999-4000
# Same regions as BED_A but with GRCh37/hs37d5-style names ('1', 'X').
BED_NOCHR = str(DATA_DIR / "beds" / "wes_kit_nochr.bed")


# --- WGS sentinel ---

def test_wgs_always_covered():
    idx = CaptureIndex.wgs()
    assert idx.covers("chr1", 1) is True
    assert idx.covers("chrX", 5_000_000) is True
    assert idx.covers("chrY", 1) is True
    assert idx.covers("chrM", 100) is True


# --- from_bed: chr1:999-2000 (0-based half-open → 1-based: 1000-2000) ---

def test_inside_region():
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr1", 1500) is True


def test_outside_region_before():
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr1", 999) is False


def test_outside_region_after():
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr1", 2001) is False


def test_exact_start_boundary():
    # BED Start=999 (0-based) → first covered 1-based pos = 1000
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr1", 1000) is True
    assert idx.covers("chr1", 999) is False


def test_exact_end_boundary():
    # BED End=2000 → last covered 1-based pos = 2000
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr1", 2000) is True
    assert idx.covers("chr1", 2001) is False


def test_wrong_chrom():
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chr2", 1500) is False


def test_chrX_covered():
    idx = CaptureIndex.from_bed(BED_A)
    assert idx.covers("chrX", 1500) is True


def test_chrX_not_covered_by_bed_b():
    idx = CaptureIndex.from_bed(BED_B)
    assert idx.covers("chrX", 1500) is False


def test_bed_b_covers_3500():
    idx = CaptureIndex.from_bed(BED_B)
    assert idx.covers("chr1", 3500) is True


def test_bed_b_does_not_cover_1500():
    idx = CaptureIndex.from_bed(BED_B)
    assert idx.covers("chr1", 1500) is False


# --- BED chrom naming: queries always arrive normalized, BEDs may not be ---

def test_nochr_bed_matches_normalized_query():
    # Regression: BEDs without the 'chr' prefix used to index under '1', so every
    # covers("chr1", ...) missed and the technology's samples silently vanished from AN.
    idx = CaptureIndex.from_bed(BED_NOCHR)
    assert idx.covers("chr1", 1500) is True
    assert idx.covers("chrX", 1500) is True


def test_nochr_bed_agrees_with_chr_bed():
    nochr = CaptureIndex.from_bed(BED_NOCHR)
    chr_ = CaptureIndex.from_bed(BED_A)
    for pos in (999, 1000, 1500, 2000, 2001):
        assert nochr.covers("chr1", pos) is chr_.covers("chr1", pos)


def test_nochr_bed_still_rejects_other_chroms():
    idx = CaptureIndex.from_bed(BED_NOCHR)
    assert idx.covers("chr2", 1500) is False


def test_index_keys_are_normalized():
    idx = CaptureIndex.from_bed(BED_NOCHR)
    assert idx.known_chroms() == {"chr1", "chrX"}


# --- public API: chrom introspection ---

def test_is_always_covered_wgs():
    assert CaptureIndex.wgs().is_always_covered is True


def test_is_always_covered_bed():
    assert CaptureIndex.from_bed(BED_A).is_always_covered is False


def test_known_chroms_excludes_unknown_contigs(tmp_path):
    bed = tmp_path / "mixed.bed"
    bed.write_text("chr1\t100\t200\ncontigZ\t100\t200\n")
    idx = CaptureIndex.from_bed(str(bed))
    assert idx.known_chroms() == {"chr1"}
    assert "chrcontigZ" in idx.indexed_chroms()


def test_is_empty_false_for_bed():
    assert CaptureIndex.from_bed(BED_A).is_empty() is False


def test_is_empty_false_for_wgs():
    # The WGS sentinel has no index but covers everything — not "empty".
    assert CaptureIndex.wgs().is_empty() is False


def test_legacy_pickle_with_unnormalized_keys_self_heals(tmp_path):
    # Simulates a database built before the fix: keys stored as '1'/'X', 2-tuple entries.
    idx = CaptureIndex.from_bed(BED_NOCHR)
    idx._index = {"1": ([999], [2000]), "X": ([999], [2000])}
    path = str(tmp_path / "legacy.pickle")
    idx.save(path)

    loaded = CaptureIndex.load(path)
    assert set(loaded._index) == {"chr1", "chrX"}
    assert loaded.covers("chr1", 1500) is True
    assert loaded.covers("chr1", 999) is False


# --- overlapping intervals: the running-max shortcut must not miss a cover ---

def test_overlapping_intervals(tmp_path):
    # An early wide interval covers pos while later, narrower ones do not.
    bed = tmp_path / "overlap.bed"
    bed.write_text("chr1\t100\t9000\nchr1\t200\t300\nchr1\t400\t500\n")
    idx = CaptureIndex.from_bed(str(bed))
    assert idx.covers("chr1", 8000) is True   # only the first interval reaches here
    assert idx.covers("chr1", 250) is True
    assert idx.covers("chr1", 9001) is False


# --- pickle save/load round-trip ---

def test_save_load_wgs(tmp_path):
    idx = CaptureIndex.wgs()
    path = str(tmp_path / "cap.pickle")
    idx.save(path)
    loaded = CaptureIndex.load(path)
    assert loaded.covers("chr1", 9999) is True


def test_save_load_bed(tmp_path):
    idx = CaptureIndex.from_bed(BED_A)
    path = str(tmp_path / "cap.pickle")
    idx.save(path)
    loaded = CaptureIndex.load(path)
    assert loaded.covers("chr1", 1500) is True
    assert loaded.covers("chr1", 999) is False

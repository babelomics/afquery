from pathlib import Path

import pytest

from afquery.capture import (
    CaptureIndex,
    describe_capture_problem,
    load_capture_indices,
)
from afquery.models import Technology

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


# --- BEDs with no intervals: report, do not crash ---

def test_from_bed_empty_file_does_not_raise(tmp_path):
    # Regression: pyranges raises IndexError on a zero-byte BED, so create-db died
    # with an opaque traceback instead of naming the offending technology.
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    idx = CaptureIndex.from_bed(str(bed))
    assert idx.is_empty() is True
    assert idx.covers("chr1", 1500) is False


def test_from_bed_comment_only_does_not_raise(tmp_path):
    # pyranges raises AssertionError here rather than IndexError — both must be caught.
    bed = tmp_path / "comment.bed"
    bed.write_text("# no regions\n")
    assert CaptureIndex.from_bed(str(bed)).is_empty() is True


# --- problem reporting ---

def test_no_problem_for_valid_bed():
    assert describe_capture_problem(CaptureIndex.from_bed(BED_A), "kit") is None


def test_no_problem_for_wgs():
    assert describe_capture_problem(CaptureIndex.wgs(), "WGS") is None


def test_problem_for_empty_bed_says_empty(tmp_path):
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    msg = describe_capture_problem(CaptureIndex.from_bed(str(bed)), "kit")
    assert "are empty" in msg
    assert "match no known chromosome" not in msg


def test_problem_for_unknown_contigs_lists_them(tmp_path):
    bed = tmp_path / "unknown.bed"
    bed.write_text("contigZ\t100\t200\n")
    msg = describe_capture_problem(CaptureIndex.from_bed(str(bed)), "kit")
    assert "match no known chromosome" in msg
    assert "chrcontigZ" in msg


def test_no_problem_when_only_some_contigs_are_unknown(tmp_path):
    # One real chromosome is enough for the index to be usable.
    bed = tmp_path / "mixed.bed"
    bed.write_text("chr1\t100\t200\ncontigZ\t100\t200\n")
    assert describe_capture_problem(CaptureIndex.from_bed(str(bed)), "kit") is None


# --- chrMT capture BEDs ---

def test_chrMT_bed_indexes_as_chrM(tmp_path):
    # Regression: a BED naming the mitochondrion 'chrMT' kept valid autosome keys, so
    # the "no known chromosome" warning never fired and chrM was silently uncovered.
    bed = tmp_path / "mito.bed"
    bed.write_text("chr1\t100\t200\nchrMT\t99\t200\n")
    idx = CaptureIndex.from_bed(str(bed))
    assert idx.known_chroms() == {"chr1", "chrM"}
    assert idx.covers("chrM", 150) is True


# --- load_capture_indices ---

def test_load_capture_indices_missing_pickle(tmp_path):
    tech = Technology(tech_id=1, tech_name="WES_kit_A", bed_path=BED_A)
    with pytest.raises(
        FileNotFoundError, match="Missing capture index for technology 'WES_kit_A'"
    ):
        load_capture_indices([tech], str(tmp_path))


def test_load_capture_indices_keys_by_tech_id(tmp_path):
    CaptureIndex.wgs().save(str(tmp_path / "tech_0.pickle"))
    CaptureIndex.from_bed(BED_A).save(str(tmp_path / "tech_1.pickle"))
    techs = [
        Technology(tech_id=0, tech_name="WGS", bed_path=None),
        Technology(tech_id=1, tech_name="WES_kit_A", bed_path=BED_A),
    ]
    loaded = load_capture_indices(techs, str(tmp_path))
    assert set(loaded) == {0, 1}
    assert loaded[0].is_always_covered is True
    assert loaded[1].covers("chr1", 1500) is True


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

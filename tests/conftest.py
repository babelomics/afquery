import json
import os
import pickle
import sqlite3
from pathlib import Path

import pytest
import pyarrow as pa
import pyarrow.parquet as pq
from pyroaring import BitMap

from afquery.bitmaps import serialize, build_sex_bitmaps, build_phenotype_bitmaps, build_tech_bitmaps
from afquery.capture import CaptureIndex
from afquery.models import Sample, Technology

GENOME_BUILD = "GRCh37"

SAMPLES = [
    (0, "S00", "male",   0),
    (1, "S01", "male",   0),
    (2, "S02", "female", 0),
    (3, "S03", "female", 0),
    (4, "S04", "male",   1),
    (5, "S05", "female", 1),
    (6, "S06", "female", 1),
    (7, "S07", "male",   2),
    (8, "S08", "male",   2),
    (9, "S09", "female", 2),
]

TECHNOLOGIES = [
    (0, "WGS",       None),
    (1, "WES_kit_A", "wes_kit_a.bed"),
    (2, "WES_kit_B", "wes_kit_b.bed"),
]

SAMPLE_PHENOTYPE = [
    (0, "E11.9"), (1, "E11.9"), (2, "E11.9"), (5, "E11.9"), (6, "E11.9"), (7, "E11.9"),
    (1, "I10"),   (3, "I10"),   (5, "I10"),   (8, "I10"),   (9, "I10"),
    (2, "C50"),   (4, "C50"),   (6, "C50"),
    (0, "J45"),   (3, "J45"),   (7, "J45"),   (8, "J45"),   (9, "J45"),
    (1, "G20"),   (2, "G20"),   (3, "G20"),   (5, "G20"),
]

# (chrom, pos, ref, alt, het_ids, hom_ids, fail_ids); as in a real build, a sample
# whose call failed FILTER is in fail_ids only
VARIANTS = [
    ("chr1",  1500,    "A", "T", [5],       [2],    [0]),      # S00 fails FILTER
    ("chr1",  3500,    "G", "C", [7],        [],     []),       # no failures
    ("chr1",  5000,    "T", "G", [0],       [],     [1, 3]),   # S01, S03 fail FILTER
    ("chrX",  5000000, "A", "G", [0, 2],    [],     []),       # no failures
    ("chrY",  500000,  "T", "C", [0, 1],    [4],    []),       # no failures
    ("chrM",  100,     "C", "A", [0, 2],    [],     [5]),      # S05 fails FILTER
]


@pytest.fixture(scope="session")
def data_dir():
    return Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def test_db(tmp_path_factory, data_dir):
    db_path = tmp_path_factory.mktemp("db")
    _build_db(db_path, data_dir)
    return str(db_path)


def _build_db(db_path: Path, data_dir: Path) -> None:
    (db_path / "variants").mkdir(exist_ok=True)
    (db_path / "capture").mkdir(exist_ok=True)

    # manifest.json
    (db_path / "manifest.json").write_text(json.dumps({
        "genome_build": GENOME_BUILD,
        "version": "0.1.0",
        "sample_count": len(SAMPLES),
        "next_sample_id": len(SAMPLES),  # Prevent ID reuse after removals
        "schema_version": "2.0",
        "pass_only_filter": True,
        "created_at": "2026-01-01T00:00:00",
    }))

    # SQLite metadata
    con = sqlite3.connect(db_path / "metadata.sqlite")
    _init_sqlite(con)
    con.close()

    # Capture indices
    _build_capture_indices(db_path, data_dir)

    # Parquet variant files
    _build_parquet(db_path)


def _init_sqlite(con: sqlite3.Connection) -> None:
    con.executescript("""
        CREATE TABLE samples (
            sample_id   INTEGER PRIMARY KEY,
            sample_name TEXT NOT NULL,
            sex         TEXT NOT NULL,
            tech_id     INTEGER NOT NULL,
            vcf_path    TEXT,
            ingested_at TEXT
        );
        CREATE TABLE technologies (
            tech_id   INTEGER PRIMARY KEY,
            tech_name TEXT NOT NULL,
            bed_path  TEXT
        );
        CREATE TABLE sample_phenotype (
            sample_id    INTEGER NOT NULL,
            phenotype_code   TEXT NOT NULL,
            PRIMARY KEY (sample_id, phenotype_code)
        );
        CREATE TABLE precomputed_bitmaps (
            bitmap_type TEXT NOT NULL,
            bitmap_key  TEXT NOT NULL,
            bitmap_data BLOB NOT NULL,
            PRIMARY KEY (bitmap_type, bitmap_key)
        );
        CREATE TABLE IF NOT EXISTS changelog (
            event_id     INTEGER PRIMARY KEY AUTOINCREMENT,
            event_type   TEXT NOT NULL,
            event_time   TEXT NOT NULL,
            sample_names TEXT,
            notes        TEXT
        );
    """)

    con.executemany(
        "INSERT INTO samples (sample_id, sample_name, sex, tech_id) VALUES (?, ?, ?, ?)",
        SAMPLES
    )
    con.executemany(
        "INSERT INTO technologies VALUES (?, ?, ?)", TECHNOLOGIES
    )
    con.executemany(
        "INSERT INTO sample_phenotype VALUES (?, ?)", SAMPLE_PHENOTYPE
    )

    samples = [Sample(s[0], s[1], s[2], s[3]) for s in SAMPLES]

    # Sex bitmaps
    sex_bms = build_sex_bitmaps(samples)
    for sex, bm in sex_bms.items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("sex", sex, serialize(bm)),
        )

    # Phenotype bitmaps
    phenotype_bms = build_phenotype_bitmaps(SAMPLE_PHENOTYPE)
    for code, bm in phenotype_bms.items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("phenotype", code, serialize(bm)),
        )

    # Tech bitmaps (key stored as string)
    tech_bms = build_tech_bitmaps(samples)
    for tech_id, bm in tech_bms.items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("tech", str(tech_id), serialize(bm)),
        )

    con.commit()


def _build_capture_indices(db_path: Path, data_dir: Path) -> None:
    beds_dir = data_dir / "beds"
    capture_dir = db_path / "capture"

    for tech_id, tech_name, bed_file in TECHNOLOGIES:
        if bed_file is None:
            idx = CaptureIndex.wgs()
        else:
            idx = CaptureIndex.from_bed(str(beds_dir / bed_file))
        idx.save(str(capture_dir / f"tech_{tech_id}.pickle"))


def _build_parquet(db_path: Path) -> None:
    schema = pa.schema([
        ("pos",                 pa.uint32()),
        ("ref",                 pa.large_utf8()),
        ("alt",                 pa.large_utf8()),
        ("het_bitmap",          pa.large_binary()),
        ("hom_bitmap",          pa.large_binary()),
        ("fail_bitmap",         pa.large_binary()),
        ("filtered_bitmap",     pa.large_binary()),
        ("quality_pass_bitmap", pa.large_binary()),
    ])

    # Group variants by chromosome
    by_chrom: dict[str, list] = {}
    for chrom, pos, ref, alt, het_ids, hom_ids, fail_ids in VARIANTS:
        by_chrom.setdefault(chrom, []).append((pos, ref, alt, het_ids, hom_ids, fail_ids))

    empty = serialize(BitMap())
    for chrom, rows in by_chrom.items():
        rows.sort(key=lambda r: r[0])  # sort by pos
        table = pa.table(
            {
                "pos":                 pa.array([r[0] for r in rows], type=pa.uint32()),
                "ref":                 pa.array([r[1] for r in rows], type=pa.large_utf8()),
                "alt":                 pa.array([r[2] for r in rows], type=pa.large_utf8()),
                "het_bitmap":          pa.array(
                    [serialize(BitMap(r[3])) for r in rows], type=pa.large_binary()
                ),
                "hom_bitmap":          pa.array(
                    [serialize(BitMap(r[4])) for r in rows], type=pa.large_binary()
                ),
                "fail_bitmap":         pa.array(
                    [serialize(BitMap(r[5])) for r in rows], type=pa.large_binary()
                ),
                "filtered_bitmap":     pa.array([empty] * len(rows), type=pa.large_binary()),
                "quality_pass_bitmap": pa.array([empty] * len(rows), type=pa.large_binary()),
            },
            schema=schema,
        )
        pq.write_table(table, db_path / "variants" / f"{chrom}.parquet")

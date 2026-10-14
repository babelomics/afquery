import glob as glob_module
import json
import logging
import os
import shutil
import sqlite3
import time
from datetime import datetime, timezone

from ..bitmaps import serialize, build_sex_bitmaps, build_phenotype_bitmaps, build_tech_bitmaps
from ..constants import VALID_GENOME_BUILDS
from ..models import Sample, Technology
from .build import build_all_parquets, consolidate_temp_files
from .ingest import ingest_all
from .manifest import parse_manifest
from .regions import build_capture_indices
from .update import add_samples, remove_samples, check_database, UpdateError, CheckResult

logger = logging.getLogger(__name__)


def run_preprocess(
    manifest_path: str,
    output_dir: str,
    genome_build: str,
    bed_dir: str | None = None,
    threads: int | None = None,
    build_threads: int | None = None,
    build_memory: str = "2GB",
    tmp_dir: str | None = None,
    force: bool = False,
    db_version: str = "1.0",
    min_dp: int = 0,
    min_gq: int = 0,
    min_qual: float = 0.0,
    min_covered: int = 0,
) -> None:
    if genome_build not in VALID_GENOME_BUILDS:
        raise ValueError(
            f"Invalid genome_build '{genome_build}', must be one of {sorted(VALID_GENOME_BUILDS)}"
        )

    # Resolve threads: None → all available CPU cores
    effective_threads = max(1, os.cpu_count() or 1) if threads is None else max(1, threads)

    samples_raw, techs_raw = parse_manifest(manifest_path, bed_dir)

    logger.info("[preprocess] Manifest: %d sample(s), %d technology/ies, build=%s",
                len(samples_raw), len(techs_raw), genome_build)

    variants_dir = os.path.join(output_dir, "variants")
    os.makedirs(variants_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "capture"), exist_ok=True)

    # Resolve tmp_dir: deterministic path so resume can find prior work
    auto_tmp = tmp_dir is None
    actual_tmp = tmp_dir if tmp_dir is not None else os.path.join(output_dir, ".tmp_preprocess")

    if force:
        logger.info("[preprocess] --force: cleaning up existing partial data and starting from scratch.")
        shutil.rmtree(actual_tmp, ignore_errors=True)
        shutil.rmtree(variants_dir, ignore_errors=True)
        os.makedirs(variants_dir, exist_ok=True)
    elif auto_tmp and os.path.exists(actual_tmp):
        logger.warning(
            "[preprocess] Existing partial run detected in %s. Resuming. "
            "Use --force to restart from scratch.", actual_tmp
        )

    os.makedirs(actual_tmp, exist_ok=True)

    # Assign tech_ids (0-indexed, first-seen order) and sample_ids (0-indexed, manifest order)
    tech_name_to_id: dict[str, int] = {}
    for idx, pt in enumerate(techs_raw):
        tech_name_to_id[pt.tech_name] = idx

    technologies = [
        Technology(
            tech_id=tech_name_to_id[pt.tech_name],
            tech_name=pt.tech_name,
            bed_path=pt.bed_path,
        )
        for pt in techs_raw
    ]

    samples = [
        Sample(
            sample_id=idx,
            sample_name=ps.sample_name,
            sex=ps.sex,
            tech_id=tech_name_to_id[ps.tech_name],
        )
        for idx, ps in enumerate(samples_raw)
    ]

    vcf_paths = [ps.vcf_path for ps in samples_raw]

    sample_phenotype_pairs: list[tuple[int, str]] = [
        (idx, code)
        for idx, ps in enumerate(samples_raw)
        for code in ps.phenotype_codes
    ]

    ingested_at = datetime.now(timezone.utc).isoformat()

    logger.debug("[preprocess] Writing SQLite metadata...")
    _write_sqlite(output_dir, samples, technologies, sample_phenotype_pairs, vcf_paths, ingested_at)

    logger.debug("[preprocess] Building capture indices...")
    capture_dir = os.path.join(output_dir, "capture")
    build_capture_indices(technologies, capture_dir)

    success = False
    try:
        consolidated_dir = os.path.join(actual_tmp, "consolidated")

        # If consolidation is already done (sample files gone, consolidated/ dir exists),
        # skip both ingest and consolidation entirely.
        sample_files = glob_module.glob(os.path.join(actual_tmp, "sample_*.parquet"))
        consolidation_done = (
            not force
            and os.path.isdir(consolidated_dir)
            and not sample_files
        )

        if consolidation_done:
            logger.info("[build] Skipping ingest + consolidation (already completed in previous run).")
            consolidated = consolidated_dir
        else:
            ingest_all(samples, vcf_paths, actual_tmp, genome_build,
                       n_workers=effective_threads, resume=(not force))

            # Re-check after ingest (edge case: consolidation was partially done)
            sample_files = glob_module.glob(os.path.join(actual_tmp, "sample_*.parquet"))
            if not force and not sample_files and os.path.isdir(consolidated_dir):
                logger.info("[build] Skipping consolidation (already done in previous run).")
                consolidated = consolidated_dir
            else:
                if os.path.isdir(consolidated_dir):
                    shutil.rmtree(consolidated_dir)
                logger.info("[build] Consolidating %d sample file(s) into per-chromosome Parquet...",
                            len(samples))
                consolidated = consolidate_temp_files(actual_tmp, threads=effective_threads)
                # Free disk space: individual sample files are now redundant
                for f in glob_module.glob(os.path.join(actual_tmp, "sample_*.parquet")):
                    os.remove(f)

        # build_threads=N → explicit cap; None → use all effective_threads
        # (build_all_parquets will further cap to n_chroms for Hive-partitioned mode)
        build_workers = (
            min(effective_threads, build_threads) if build_threads is not None
            else effective_threads
        )
        logger.info("[build] Build memory limit: %s per worker", build_memory)

        # Phase 2: pre-serialize WES tech bitmaps for picklable transfer to workers
        wes_tech_bitmaps_bytes: dict[int, bytes] = {}
        if min_covered > 0:
            for tech_id, bm in build_tech_bitmaps(samples).items():
                tech_obj = next((t for t in technologies if t.tech_id == tech_id), None)
                if tech_obj is not None and tech_obj.bed_path is not None:
                    wes_tech_bitmaps_bytes[tech_id] = serialize(bm)

        build_all_parquets(actual_tmp, variants_dir, n_workers=build_workers,
                           consolidated_path=consolidated, resume=(not force),
                           memory_limit=build_memory,
                           min_dp=min_dp, min_gq=min_gq, min_qual=min_qual,
                           min_covered=min_covered,
                           wes_tech_bitmaps_bytes=wes_tech_bitmaps_bytes or None)
        success = True
    finally:
        if auto_tmp and success:
            shutil.rmtree(actual_tmp, ignore_errors=True)
        elif not success:
            logger.warning(
                "[preprocess] Run failed. Partial results preserved in %s. "
                "Re-run the same command to resume.", actual_tmp
            )

    logger.debug("[preprocess] Writing manifest.json...")
    _write_manifest(
        output_dir, genome_build, len(samples), db_version=db_version,
        min_dp=min_dp, min_gq=min_gq, min_qual=min_qual, min_covered=min_covered,
    )

    logger.info("[preprocess] Database complete: %s", output_dir)


def _write_sqlite(
    output_dir: str,
    samples: list[Sample],
    technologies: list[Technology],
    sample_phenotype_pairs: list[tuple[int, str]],
    vcf_paths: list[str] | None = None,
    ingested_at: str | None = None,
) -> None:
    db_path = os.path.join(output_dir, "metadata.sqlite")

    # Always start fresh: preprocess is a full rebuild and SQLite creation is instantaneous.
    # This avoids schema mismatch errors when re-running after a partial failure.
    if os.path.exists(db_path):
        os.remove(db_path)

    con = sqlite3.connect(db_path)

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
        CREATE TABLE changelog (
            event_id     INTEGER PRIMARY KEY AUTOINCREMENT,
            event_type   TEXT NOT NULL,
            event_time   TEXT NOT NULL,
            sample_names TEXT,
            notes        TEXT
        );
    """)


    paths = vcf_paths if vcf_paths is not None else [None] * len(samples)
    con.executemany(
        "INSERT INTO samples (sample_id, sample_name, sex, tech_id, vcf_path, ingested_at)"
        " VALUES (?, ?, ?, ?, ?, ?)",
        [(s.sample_id, s.sample_name, s.sex, s.tech_id, paths[i], ingested_at)
         for i, s in enumerate(samples)],
    )
    con.executemany(
        "INSERT INTO technologies VALUES (?, ?, ?)",
        [(t.tech_id, t.tech_name, t.bed_path) for t in technologies],
    )
    con.executemany(
        "INSERT INTO sample_phenotype VALUES (?, ?)",
        sample_phenotype_pairs,
    )

    for sex, bm in build_sex_bitmaps(samples).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("sex", sex, serialize(bm)),
        )

    for code, bm in build_phenotype_bitmaps(sample_phenotype_pairs).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("phenotype", code, serialize(bm)),
        )

    for tech_id, bm in build_tech_bitmaps(samples).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("tech", str(tech_id), serialize(bm)),
        )

    sample_names_json = json.dumps([s.sample_name for s in samples])
    event_time = ingested_at or datetime.now(timezone.utc).isoformat()
    con.execute(
        "INSERT INTO changelog (event_type, event_time, sample_names, notes) VALUES (?, ?, ?, ?)",
        ("preprocess", event_time, sample_names_json, f"{len(samples)} samples ingested"),
    )

    con.commit()
    con.close()


def _write_manifest(
    output_dir: str,
    genome_build: str,
    sample_count: int,
    db_version: str = "1.0",
    min_dp: int = 0,
    min_gq: int = 0,
    min_qual: float = 0.0,
    min_covered: int = 0,
) -> None:
    coverage_filter = {
        "min_dp": min_dp,
        "min_gq": min_gq,
        "min_qual": min_qual,
        "min_covered": min_covered,
    }
    has_quality_filter = (min_dp > 0 or min_gq > 0 or min_qual > 0 or min_covered > 0)
    manifest = {
        "genome_build": genome_build,
        "version": "0.1.0",
        "db_version": db_version,
        "sample_count": sample_count,
        "schema_version": "3.0" if has_quality_filter else "2.0",
        "pass_only_filter": True,
        "coverage_filter": coverage_filter,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    with open(os.path.join(output_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)

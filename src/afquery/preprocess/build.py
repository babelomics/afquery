import glob as glob_module
import logging
import multiprocessing
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool

import tqdm

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from pyroaring import BitMap

from .. import storage
from ..bitmaps import serialize
from ..constants import ALL_CHROMS

logger = logging.getLogger(__name__)

BUCKET_SIZE = storage.BUCKET_SIZE

PARQUET_SCHEMA = pa.schema([
    ("pos",                 pa.uint32()),
    ("ref",                 pa.large_utf8()),
    ("alt",                 pa.large_utf8()),
    ("het_bitmap",          pa.large_binary()),
    ("hom_bitmap",          pa.large_binary()),
    ("fail_bitmap",         pa.large_binary()),
    ("filtered_bitmap",     pa.large_binary()),  # Phase 2: WES non-carriers w/ uncertain coverage
    ("quality_pass_bitmap", pa.large_binary()),  # Phase 2: carriers meeting DP/GQ/QUAL thresholds
])


def _filter_known_chroms(chroms: list[str], source: str) -> list[str]:
    """Keep only real chromosomes, reporting the dropped contigs once.

    Chrom discovery happens in three places and every one of them must drop
    unplaced and alt contigs. Centralizing the filter is what keeps the report
    from either vanishing — it used to live downstream of the filter, where it
    could never fire — or firing once per row.
    """
    known = [c for c in chroms if c in ALL_CHROMS]
    unknown = sorted({c for c in chroms if c not in ALL_CHROMS})
    if unknown:
        logger.warning(
            "[build] Skipping %d unknown contig(s) from %s: %s",
            len(unknown), source, ", ".join(unknown),
        )
    return known


def _get_chroms_from_file(parquet_path: str) -> list[str]:
    """SELECT DISTINCT chrom from a single consolidated Parquet file."""
    path = parquet_path.replace("'", "''")
    con = duckdb.connect()
    con.execute("SET enable_progress_bar=false")
    try:
        rows = con.execute(
            f"SELECT DISTINCT chrom FROM read_parquet('{path}')"
        ).fetchall()
    except Exception:
        return []
    finally:
        con.close()
    return _filter_known_chroms([r[0] for r in rows], parquet_path)


def consolidate_temp_files(
    tmp_dir: str,
    memory_limit: str = "8GB",
    threads: int | None = None,
) -> str:
    """Merge sample_*.parquet into per-chromosome Hive-partitioned directory.

    Uses DuckDB streaming COPY with PARTITION_BY to create chrom=X/ subdirectories.
    Returns path to consolidated directory (not a single file).
    """
    effective_threads = threads if threads is not None else (os.cpu_count() or 1)
    glob_pattern = os.path.join(tmp_dir, "sample_*.parquet").replace("'", "''")
    out_dir = os.path.join(tmp_dir, "consolidated")
    os.makedirs(out_dir, exist_ok=True)
    out_dir_esc = out_dir.replace("'", "''")

    con = duckdb.connect()
    con.execute(f"SET memory_limit='{memory_limit}'")
    con.execute(f"SET threads={effective_threads}")
    con.execute("SET enable_progress_bar=false")
    con.execute(
        f"COPY (SELECT * FROM read_parquet('{glob_pattern}')) "
        f"TO '{out_dir_esc}' (FORMAT PARQUET, PARTITION_BY (chrom))"
    )
    con.close()
    return out_dir


def _get_chroms_from_consolidated(path: str) -> list[str]:
    """Discover chromosomes from consolidated path (Hive-partitioned directory or single file)."""
    if os.path.isdir(path):
        # Hive-partitioned: scan chrom=X/ subdirs — fast, no SQL needed
        chroms = [
            entry.name[6:]  # strip "chrom="
            for entry in os.scandir(path)
            if entry.is_dir() and entry.name.startswith("chrom=")
        ]
        return _filter_known_chroms(chroms, path)
    else:
        # Single file (legacy) — use SQL
        return _get_chroms_from_file(path)


def get_chroms_in_temp_files(tmp_dir: str) -> list[str]:
    """SELECT DISTINCT chrom from all sample_*.parquet in tmp_dir."""
    parquet_files = glob_module.glob(os.path.join(tmp_dir, "sample_*.parquet"))
    if not parquet_files:
        return []

    glob_pattern = os.path.join(tmp_dir, "sample_*.parquet").replace("'", "''")
    con = duckdb.connect()
    con.execute("SET enable_progress_bar=false")
    try:
        rows = con.execute(
            f"SELECT DISTINCT chrom FROM read_parquet('{glob_pattern}')"
        ).fetchall()
    except Exception:
        return []
    finally:
        con.close()

    return _filter_known_chroms([r[0] for r in rows], tmp_dir)


def _make_table(
    positions, refs, alts,
    het_bitmaps, hom_bitmaps, fail_bitmaps,
    filtered_bitmaps, quality_pass_bitmaps,
) -> pa.Table:
    return pa.table(
        {
            "pos":                 pa.array(positions,            type=pa.uint32()),
            "ref":                 pa.array(refs,                 type=pa.large_utf8()),
            "alt":                 pa.array(alts,                 type=pa.large_utf8()),
            "het_bitmap":          pa.array(het_bitmaps,          type=pa.large_binary()),
            "hom_bitmap":          pa.array(hom_bitmaps,          type=pa.large_binary()),
            "fail_bitmap":         pa.array(fail_bitmaps,         type=pa.large_binary()),
            "filtered_bitmap":     pa.array(filtered_bitmaps,     type=pa.large_binary()),
            "quality_pass_bitmap": pa.array(quality_pass_bitmaps, type=pa.large_binary()),
        },
        schema=PARQUET_SCHEMA,
    )


def _write_flat(
    chrom: str,
    variants_dir: str,
    positions, refs, alts,
    het_bitmaps, hom_bitmaps, fail_bitmaps,
    filtered_bitmaps, quality_pass_bitmaps,
    row_group_size: int,
) -> None:
    table = _make_table(
        positions, refs, alts,
        het_bitmaps, hom_bitmaps, fail_bitmaps,
        filtered_bitmaps, quality_pass_bitmaps,
    )
    out_path = os.path.join(variants_dir, f"{chrom}.parquet")
    tmp_path = out_path + ".tmp"
    pq.write_table(table, tmp_path, row_group_size=row_group_size)
    os.rename(tmp_path, out_path)


def _passes_quality(dp, gq, qual, min_dp: int, min_gq: int, min_qual: float) -> bool:
    """A sample passes quality if all configured thresholds are met (None always passes)."""
    if min_dp > 0 and (dp is None or dp < min_dp):
        return False
    if min_gq > 0 and (gq is None or gq < min_gq):
        return False
    if min_qual > 0 and (qual is None or qual < min_qual):
        return False
    return True


def _compute_phase2_bitmaps(
    sample_ids: list[int],
    gt_acs: list[int],
    filter_passes: list[bool],
    dps: list[int | None],
    gqs: list[int | None],
    quals: list[float | None],
    het_ids: list[int],
    hom_ids: list[int],
    fail_ids: list[int],
    min_dp: int,
    min_gq: int,
    min_qual: float,
    min_covered: int,
    wes_tech_bitmaps: dict[int, BitMap] | None,
) -> tuple[list[int], list[int]]:
    """Compute (quality_pass_ids, filtered_ids) for one (pos, ref, alt) row.

    quality_pass_ids: PASS carriers (het|hom) meeting DP/GQ/QUAL thresholds.
    filtered_ids: non-carrier samples of WES techs whose quality_pass count < min_covered.
                  Query layer intersects this with `eligible` (BED-aware) to get the
                  effective N_NO_COVERAGE.
    """
    # Quality-pass carriers
    quality_pass_ids: list[int] = []
    for sid, ac, fp, dp, gq, qual in zip(sample_ids, gt_acs, filter_passes, dps, gqs, quals):
        if ac > 0 and fp and _passes_quality(dp, gq, qual, min_dp, min_gq, min_qual):
            quality_pass_ids.append(sid)

    # Filtered samples (non-carriers of failing WES techs)
    if not wes_tech_bitmaps or min_covered <= 0:
        return quality_pass_ids, []

    quality_pass_set = set(quality_pass_ids)
    carrier_set = set(het_ids) | set(hom_ids) | set(fail_ids)
    filtered_ids: list[int] = []
    for tech_bm in wes_tech_bitmaps.values():
        qp_count = sum(1 for sid in quality_pass_set if sid in tech_bm)
        if qp_count < min_covered:
            for sid in tech_bm:
                if sid not in carrier_set:
                    filtered_ids.append(sid)
    return quality_pass_ids, filtered_ids


def _discover_bucket_ids(
    source: str,
    base_where: str,
    params: list,
    memory_limit: str = "2GB",
) -> list[int]:
    """Discover distinct bucket IDs (pos // BUCKET_SIZE) present in source."""
    con = duckdb.connect()
    if memory_limit:
        con.execute(f"SET memory_limit='{memory_limit}'")
    con.execute("SET threads=1")
    con.execute("SET enable_progress_bar=false")
    try:
        sql = f"""
            SELECT DISTINCT CAST(pos AS BIGINT) // {BUCKET_SIZE} AS bucket_id
            FROM read_parquet('{source}')
            {base_where}
            ORDER BY 1
        """
        rows = con.execute(sql, params).fetchall()
    except Exception:
        rows = []
    finally:
        con.close()
    return [int(r[0]) for r in rows]


def _build_one_bucket_worker(
    source: str,
    base_where: str,
    params: list,
    bucket_id: int,
    out_path: str,
    row_group_size: int,
    memory_limit: str,
    min_dp: int = 0,
    min_gq: int = 0,
    min_qual: float = 0.0,
    min_covered: int = 0,
    wes_tech_bitmaps_bytes: dict[int, bytes] | None = None,
) -> tuple[int, int, float]:
    """Build one bucket Parquet file. Returns (bucket_id, variant_count, elapsed_seconds).

    Top-level function (picklable) for use with ProcessPoolExecutor.
    """
    from ..bitmaps import deserialize as _deser
    t0 = time.monotonic()
    bucket_start = bucket_id * BUCKET_SIZE
    bucket_end = bucket_start + BUCKET_SIZE
    if base_where:
        full_where = f"{base_where} AND pos >= {bucket_start} AND pos < {bucket_end}"
    else:
        full_where = f"WHERE pos >= {bucket_start} AND pos < {bucket_end}"

    wes_tech_bitmaps = (
        {tid: _deser(b) for tid, b in wes_tech_bitmaps_bytes.items()}
        if wes_tech_bitmaps_bytes else None
    )

    con = duckdb.connect()
    if memory_limit:
        con.execute(f"SET memory_limit='{memory_limit}'")
    con.execute("SET threads=1")
    con.execute("SET preserve_insertion_order=false")
    con.execute("SET enable_progress_bar=false")

    sql = f"""
        SELECT
            pos, ref, alt,
            list(sample_id   ORDER BY sample_id) AS sample_ids,
            list(gt_ac       ORDER BY sample_id) AS gt_acs,
            list(filter_pass ORDER BY sample_id) AS filter_passes,
            list(dp          ORDER BY sample_id) AS dps,
            list(gq          ORDER BY sample_id) AS gqs,
            list(qual        ORDER BY sample_id) AS quals
        FROM read_parquet('{source}')
        {full_where}
        GROUP BY pos, ref, alt
        ORDER BY pos, alt
    """
    rows = con.execute(sql, params).fetchall()
    con.close()

    if not rows:
        return bucket_id, 0, time.monotonic() - t0

    positions = []
    refs = []
    alts = []
    het_bitmaps = []
    hom_bitmaps = []
    fail_bitmaps = []
    filtered_bitmaps = []
    quality_pass_bitmaps = []

    for pos, ref, alt, sample_ids, gt_acs, filter_passes, dps, gqs, quals in rows:
        het_ids  = [sid for sid, ac, fp in zip(sample_ids, gt_acs, filter_passes) if ac == 1 and fp]
        hom_ids  = [sid for sid, ac, fp in zip(sample_ids, gt_acs, filter_passes) if ac == 2 and fp]
        fail_ids = [sid for sid, fp in zip(sample_ids, filter_passes) if not fp]
        if not het_ids and not hom_ids and not fail_ids:
            continue
        quality_pass_ids, filtered_ids = _compute_phase2_bitmaps(
            sample_ids, gt_acs, filter_passes, dps, gqs, quals,
            het_ids, hom_ids, fail_ids,
            min_dp, min_gq, min_qual, min_covered, wes_tech_bitmaps,
        )
        positions.append(pos)
        refs.append(ref)
        alts.append(alt)
        het_bitmaps.append(serialize(BitMap(het_ids)))
        hom_bitmaps.append(serialize(BitMap(hom_ids)))
        fail_bitmaps.append(serialize(BitMap(fail_ids)))
        filtered_bitmaps.append(serialize(BitMap(filtered_ids)))
        quality_pass_bitmaps.append(serialize(BitMap(quality_pass_ids)))

    del rows

    if not positions:
        return bucket_id, 0, time.monotonic() - t0

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    tmp_path = out_path + ".tmp"
    table = _make_table(
        positions, refs, alts,
        het_bitmaps, hom_bitmaps, fail_bitmaps,
        filtered_bitmaps, quality_pass_bitmaps,
    )
    pq.write_table(table, tmp_path, row_group_size=row_group_size)
    os.rename(tmp_path, out_path)

    return bucket_id, len(positions), time.monotonic() - t0


def build_chromosome_parquet(
    chrom: str,
    tmp_dir: str,
    variants_dir: str,
    row_group_size: int = 100_000,
    partitioned: bool = False,
    memory_limit: str = "2GB",
    consolidated_path: str | None = None,
    min_dp: int = 0,
    min_gq: int = 0,
    min_qual: float = 0.0,
    min_covered: int = 0,
    wes_tech_bitmaps_bytes: dict[int, bytes] | None = None,
) -> int:
    """Build Parquet for one chromosome. Returns variant count."""
    from ..bitmaps import deserialize as _deser
    params: list = []
    wes_tech_bitmaps = (
        {tid: _deser(b) for tid, b in wes_tech_bitmaps_bytes.items()}
        if wes_tech_bitmaps_bytes else None
    )
    if consolidated_path is not None and os.path.isdir(consolidated_path):
        # Hive-partitioned: read chromosome-specific subdirectory (no WHERE needed)
        chrom_dir = os.path.join(consolidated_path, f"chrom={chrom}")
        if not os.path.exists(chrom_dir):
            return 0
        source = os.path.join(chrom_dir, "*.parquet").replace("'", "''")
        where_clause = ""
    elif consolidated_path is not None:
        # Legacy: single all_samples.parquet
        source = consolidated_path.replace("'", "''")
        where_clause = "WHERE chrom = ?"
        params = [chrom]
    else:
        # Fallback: glob over sample_*.parquet
        source = os.path.join(tmp_dir, "sample_*.parquet").replace("'", "''")
        where_clause = "WHERE chrom = ?"
        params = [chrom]

    if partitioned:
        # Process bucket by bucket: avoids loading entire chromosome into memory
        bucket_ids = _discover_bucket_ids(source, where_clause, params, memory_limit)
        if not bucket_ids:
            return 0
        total_count = 0
        for bucket_id in bucket_ids:
            out_path = os.path.join(variants_dir, chrom, f"bucket_{bucket_id}.parquet")
            _, count, _ = _build_one_bucket_worker(
                source, where_clause, params, bucket_id, out_path,
                row_group_size, memory_limit,
                min_dp, min_gq, min_qual, min_covered, wes_tech_bitmaps_bytes,
            )
            total_count += count
        return total_count

    # Flat: load entire chromosome, write single file
    con = duckdb.connect()
    if memory_limit:
        con.execute(f"SET memory_limit='{memory_limit}'")
    con.execute("SET threads=1")
    con.execute("SET preserve_insertion_order=false")
    con.execute("SET enable_progress_bar=false")

    sql = f"""
        SELECT
            pos, ref, alt,
            list(sample_id   ORDER BY sample_id) AS sample_ids,
            list(gt_ac       ORDER BY sample_id) AS gt_acs,
            list(filter_pass ORDER BY sample_id) AS filter_passes,
            list(dp          ORDER BY sample_id) AS dps,
            list(gq          ORDER BY sample_id) AS gqs,
            list(qual        ORDER BY sample_id) AS quals
        FROM read_parquet('{source}')
        {where_clause}
        GROUP BY pos, ref, alt
        ORDER BY pos, alt
    """
    rows = con.execute(sql, params).fetchall()
    con.close()

    if not rows:
        return 0

    positions = []
    refs = []
    alts = []
    het_bitmaps = []
    hom_bitmaps = []
    fail_bitmaps = []
    filtered_bitmaps = []
    quality_pass_bitmaps = []

    for pos, ref, alt, sample_ids, gt_acs, filter_passes, dps, gqs, quals in rows:
        het_ids  = [sid for sid, ac, fp in zip(sample_ids, gt_acs, filter_passes) if ac == 1 and fp]
        hom_ids  = [sid for sid, ac, fp in zip(sample_ids, gt_acs, filter_passes) if ac == 2 and fp]
        fail_ids = [sid for sid, fp in zip(sample_ids, filter_passes) if not fp]
        if not het_ids and not hom_ids and not fail_ids:
            continue
        quality_pass_ids, filtered_ids = _compute_phase2_bitmaps(
            sample_ids, gt_acs, filter_passes, dps, gqs, quals,
            het_ids, hom_ids, fail_ids,
            min_dp, min_gq, min_qual, min_covered, wes_tech_bitmaps,
        )
        positions.append(pos)
        refs.append(ref)
        alts.append(alt)
        het_bitmaps.append(serialize(BitMap(het_ids)))
        hom_bitmaps.append(serialize(BitMap(hom_ids)))
        fail_bitmaps.append(serialize(BitMap(fail_ids)))
        filtered_bitmaps.append(serialize(BitMap(filtered_ids)))
        quality_pass_bitmaps.append(serialize(BitMap(quality_pass_ids)))

    del rows  # free DuckDB result before PyArrow allocation

    if not positions:
        return 0

    _write_flat(
        chrom, variants_dir,
        positions, refs, alts,
        het_bitmaps, hom_bitmaps, fail_bitmaps,
        filtered_bitmaps, quality_pass_bitmaps,
        row_group_size,
    )

    return len(positions)


def _build_chrom_worker(
    chrom: str,
    tmp_dir: str,
    variants_dir: str,
    row_group_size: int,
    partitioned: bool,
    memory_limit: str = "2GB",
    consolidated_path: str | None = None,
    min_dp: int = 0,
    min_gq: int = 0,
    min_qual: float = 0.0,
    min_covered: int = 0,
    wes_tech_bitmaps_bytes: dict[int, bytes] | None = None,
) -> tuple[str, int, float]:
    """Top-level worker for ProcessPoolExecutor (must be picklable). Returns (chrom, count, elapsed)."""
    t0 = time.monotonic()
    count = build_chromosome_parquet(
        chrom, tmp_dir, variants_dir, row_group_size, partitioned, memory_limit,
        consolidated_path,
        min_dp=min_dp, min_gq=min_gq, min_qual=min_qual, min_covered=min_covered,
        wes_tech_bitmaps_bytes=wes_tech_bitmaps_bytes,
    )
    return chrom, count, time.monotonic() - t0


def build_all_parquets(
    tmp_dir: str,
    variants_dir: str,
    row_group_size: int = 100_000,
    n_workers: int | None = None,
    partitioned: bool = True,
    memory_limit: str = "2GB",
    consolidated_path: str | None = None,
    resume: bool = True,
    min_dp: int = 0,
    min_gq: int = 0,
    min_qual: float = 0.0,
    min_covered: int = 0,
    wes_tech_bitmaps_bytes: dict[int, bytes] | None = None,
) -> dict[str, int]:
    """Build Parquet for all discovered chromosomes. Returns {chrom: count}.

    When consolidated_path is a Hive-partitioned directory and partitioned=True,
    tasks are submitted at bucket granularity (one 1 Mbp bucket per worker task),
    allowing full CPU utilization with low per-worker memory (~2 GB default).

    For legacy single-file sources, falls back to chromosome-level workers capped at 4.
    memory_limit sets the DuckDB per-worker limit.
    resume=True skips already-built Parquet files (bucket-level for Hive, chrom-level for legacy).
    """
    if consolidated_path is not None:
        chroms = _get_chroms_from_consolidated(consolidated_path)
    else:
        chroms = get_chroms_in_temp_files(tmp_dir)
    # Both discovery helpers already filter to ALL_CHROMS and report what they drop.
    valid_chroms = chroms

    if not valid_chroms:
        return {}

    effective_workers = n_workers if n_workers is not None else (os.cpu_count() or 1)

    # ── Hive-partitioned source + partitioned output ──────────────────────────
    # Use bucket-level parallelism: one task per 1 Mbp bucket across all chroms.
    # Memory per worker: ~2 GB (1 Mbp × ~20K variants × ~100 carriers → ~500 MB input,
    # 3–4× GROUP BY overhead → ~2 GB peak). All N cores can run in parallel.
    if consolidated_path is not None and os.path.isdir(consolidated_path) and partitioned:
        t0 = time.monotonic()

        # Discover bucket tasks for all chromosomes, skipping already-built buckets
        all_tasks: list[tuple[str, int, str, str, list, str]] = []
        # (chrom, bucket_id, source, base_where, params, out_path)
        skipped_full_chroms: list[str] = []
        n_skipped_buckets: int = 0

        for chrom in valid_chroms:
            chrom_src_dir = os.path.join(consolidated_path, f"chrom={chrom}")
            if not os.path.exists(chrom_src_dir):
                continue
            source = os.path.join(chrom_src_dir, "*.parquet").replace("'", "''")
            bucket_ids = _discover_bucket_ids(source, "", [], memory_limit)
            if not bucket_ids:
                continue

            chrom_out_dir = os.path.join(variants_dir, chrom)
            os.makedirs(chrom_out_dir, exist_ok=True)

            remaining = []
            for bucket_id in bucket_ids:
                out_path = os.path.join(chrom_out_dir, f"bucket_{bucket_id}.parquet")
                if resume and os.path.exists(out_path):
                    n_skipped_buckets += 1
                    continue
                remaining.append((chrom, bucket_id, source, "", [], out_path))

            if not remaining:
                skipped_full_chroms.append(chrom)
            else:
                all_tasks.extend(remaining)

        if skipped_full_chroms:
            logger.info(
                "[build] Skipping %d fully-built chrom(s): %s.",
                len(skipped_full_chroms), ", ".join(sorted(skipped_full_chroms)),
            )

        if not all_tasks:
            logger.info("[build] All chromosomes already built.")
            return {}

        n_buckets = len(all_tasks)
        n_chroms = len({t[0] for t in all_tasks})
        capped_workers = min(effective_workers, n_buckets)
        logger.info(
            "[build] Discovered %d bucket(s) across %d chromosome(s) (%d worker(s)).\n"
            "        Memory per worker: %s  |  Estimated peak: %d × %s  |"
            "  Reduce --build-threads if OOM.",
            n_buckets, n_chroms, capped_workers,
            memory_limit, capped_workers, memory_limit,
        )

        if n_buckets > 1_000_000:
            raise RuntimeError(
                f"[build] Discovered {n_buckets:,} buckets — this is abnormally high "
                f"(expected < 100,000 for a human genome). "
                f"Possible DuckDB float-division bug. Check _discover_bucket_ids."
            )

        result: dict[str, int] = {}
        total_variants = 0
        total_discovered_buckets = n_skipped_buckets + n_buckets
        try:
            _spawn_ctx = multiprocessing.get_context("spawn")
            with ProcessPoolExecutor(max_workers=capped_workers, mp_context=_spawn_ctx) as pool:
                futures = {
                    pool.submit(
                        _build_one_bucket_worker,
                        source, base_where, params, bucket_id, out_path,
                        row_group_size, memory_limit,
                        min_dp, min_gq, min_qual, min_covered,
                        wes_tech_bitmaps_bytes,
                    ): (chrom, bucket_id)
                    for chrom, bucket_id, source, base_where, params, out_path in all_tasks
                }
                with tqdm.tqdm(
                    total=total_discovered_buckets,
                    initial=n_skipped_buckets,
                    unit="bucket",
                    desc="[build]",
                    dynamic_ncols=True,
                    disable=not sys.stderr.isatty(),
                ) as bar:
                    for fut in as_completed(futures):
                        chrom, bucket_id = futures[fut]
                        try:
                            _, count, elapsed = fut.result()
                            if count > 0:
                                result[chrom] = result.get(chrom, 0) + count
                                total_variants += count
                            bar.set_postfix(variants=f"{total_variants:,}", refresh=False)
                            bar.update(1)
                            logger.debug(
                                "  [build] %s bucket_%d: %s variant(s) (%.1fs)",
                                chrom, bucket_id, f"{count:,}", elapsed,
                            )
                        except BrokenProcessPool:
                            raise RuntimeError(
                                f"A build worker was killed (likely OOM).\n"
                                f"  Current settings: --build-memory {memory_limit},"
                                f" --build-threads {capped_workers}\n"
                                f"  Try: --build-memory 4GB"
                                f"  OR  --build-threads {max(1, capped_workers // 4)}\n"
                                f"  Peak memory = --build-memory × --build-threads"
                            ) from None
        except BrokenProcessPool:
            raise RuntimeError(
                f"A build worker was killed (likely OOM).\n"
                f"  Current settings: --build-memory {memory_limit},"
                f" --build-threads {capped_workers}\n"
                f"  Try: --build-memory 4GB"
                f"  OR  --build-threads {max(1, capped_workers // 4)}\n"
                f"  Peak memory = --build-memory × --build-threads"
            ) from None

        total = sum(result.values())
        logger.info(
            "[build] Complete: %d chromosome(s), %s total variant(s) (%.1fs)",
            len(result), f"{total:,}", time.monotonic() - t0,
        )
        return result

    # ── Legacy / non-Hive or flat output: chromosome-level workers ───────────
    # Split into already-built vs to-build
    to_build = []
    skipped_chroms = []
    for chrom in valid_chroms:
        if resume:
            if partitioned:
                done = bool(storage.existing_bucket_ids(variants_dir, chrom))
            else:
                done = storage.flat_path(variants_dir, chrom).exists()
            if done:
                skipped_chroms.append(chrom)
                continue
        to_build.append(chrom)

    if skipped_chroms:
        logger.info("[build] Skipping %d already-built chrom(s): %s. %d chrom(s) remaining.",
                    len(skipped_chroms), ", ".join(sorted(skipped_chroms)), len(to_build))

    if not to_build:
        logger.info("[build] All chromosomes already built.")
        return {}

    # Hive-partitioned source (but partitioned=False output): each worker reads only its chrom subdir
    # Legacy single-file: each worker reads entire file — cap at 4 to avoid OOM
    if consolidated_path is not None and os.path.isdir(consolidated_path):
        capped_workers = min(effective_workers, len(to_build))
    else:
        capped_workers = min(effective_workers, 4)

    logger.info("[build] Building Parquet for %d chromosome(s) (%d worker(s))...",
                len(to_build), capped_workers)
    t0 = time.monotonic()

    result: dict[str, int] = {}
    try:
        _spawn_ctx = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=capped_workers, mp_context=_spawn_ctx) as pool:
            futures = {
                pool.submit(
                    _build_chrom_worker,
                    chrom, tmp_dir, variants_dir, row_group_size, partitioned, memory_limit,
                    consolidated_path,
                    min_dp, min_gq, min_qual, min_covered,
                    wes_tech_bitmaps_bytes,
                ): chrom
                for chrom in to_build
            }
            for fut in as_completed(futures):
                try:
                    chrom, count, elapsed = fut.result()  # propagates exceptions
                    if count > 0:
                        result[chrom] = count
                    logger.debug("  [build] %s: %s variant(s) (%.1fs)", chrom, f"{count:,}", elapsed)
                except BrokenProcessPool:
                    raise RuntimeError(
                        "A build worker was killed abruptly (likely out-of-memory). "
                        "Try reducing --build-threads or increasing available memory."
                    ) from None
    except BrokenProcessPool:
        raise RuntimeError(
            "A build worker was killed abruptly (likely out-of-memory). "
            "Try reducing --build-threads or increasing available memory."
        ) from None

    total = sum(result.values())
    logger.info("[build] Complete: %d chromosome(s), %s total variant(s) (%.1fs)",
                len(result), f"{total:,}", time.monotonic() - t0)

    return result

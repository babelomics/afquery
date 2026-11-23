import csv
import itertools
import logging
import os
import sys
from pathlib import Path

import duckdb

from . import storage
from .bitmaps import deserialize
from .constants import normalize_chrom, ALL_CHROMS
from .models import SampleFilter

logger = logging.getLogger(__name__)

BUCKET_SIZE = storage.BUCKET_SIZE

# Statistic columns, in output order. The whole-population block uses these
# names as they are; each group repeats them with a ``_<label>`` suffix.
STAT_FIELDS = [
    "AC", "AN", "AF", "AC_REF", "n_eligible",
    "N_HET", "N_HOM_ALT", "N_HEMI", "N_HOM_REF", "N_FAIL", "N_NO_COVERAGE",
]


def _af(stats):
    """AF of one variant's stats, or None (an empty cell) when AN is 0."""
    return stats.AC / stats.AN if stats.AN > 0 else None


def _stat_values(stats, n_eligible: int) -> dict:
    """Values of STAT_FIELDS for one variant among one set of samples."""
    return {
        "AC": stats.AC,
        "AN": stats.AN,
        "AF": _af(stats),
        "AC_REF": stats.AC_REF,
        "n_eligible": n_eligible,
        "N_HET": stats.N_HET,
        "N_HOM_ALT": stats.N_HOM_ALT,
        "N_HEMI": stats.N_HEMI,
        "N_HOM_REF": stats.N_HOM_REF,
        "N_FAIL": stats.N_FAIL,
        "N_NO_COVERAGE": stats.N_NO_COVERAGE,
    }


class DumpColumnError(ValueError):
    """Two dump columns would share a name."""


def dump_columns(group_labels: list[str]) -> list[str]:
    """CSV header for a dump with the given group labels.

    Raises DumpColumnError when two columns would share a name, which happens when
    a group label turns one column into another (a technology named ``REF``
    makes ``AC_REF`` both the base reference count and that group's AC).
    """
    columns = ["chrom", "pos", "ref", "alt"] + list(STAT_FIELDS)
    for label in group_labels:
        columns += [f"{field}_{label}" for field in STAT_FIELDS]
    seen: set[str] = set()
    duplicates = []
    for col in columns:
        if col in seen:
            duplicates.append(col)
        seen.add(col)
    if duplicates:
        raise DumpColumnError(
            "Group labels produce duplicate dump columns: "
            + ", ".join(sorted(set(duplicates)))
            + ". Rename the technology or phenotype code, or disaggregate on fewer axes."
        )
    return columns


def _build_groups(engine, base_sf, by_sex, by_tech, by_phenotype, all_groups):
    """Build list of (label, SampleFilter) for desagregation groups.

    Returns [] if no disaggregation flags are active.
    """
    if not by_sex and not by_tech and not by_phenotype and not all_groups:
        return []

    # Each axis is a list of (label_part, dimension, value)
    axes = []

    if by_sex or all_groups:
        axes.append([("male", "sex", "male"), ("female", "sex", "female")])

    if by_tech or all_groups:
        tech_names = sorted(engine._tech_name_to_id.keys())
        axes.append([(name, "tech", name) for name in tech_names])

    pheno_codes = list(dict.fromkeys(by_phenotype))
    if all_groups and not pheno_codes:
        import sqlite3
        con = sqlite3.connect(str(engine._db / "metadata.sqlite"))
        pheno_codes = [
            r[0]
            for r in con.execute(
                "SELECT DISTINCT phenotype_code FROM sample_phenotype ORDER BY phenotype_code"
            ).fetchall()
        ]
        con.close()
    if pheno_codes:
        axes.append([(code, "phenotype", code) for code in sorted(pheno_codes)])

    if not axes:
        return []

    groups = []
    for combo in itertools.product(*axes):
        label_parts = []
        sex_override = None
        tech_inc_override = None
        pheno_inc_override = None

        for label_part, dim, value in combo:
            label_parts.append(label_part)
            if dim == "sex":
                sex_override = value
            elif dim == "tech":
                tech_inc_override = [value]
            elif dim == "phenotype":
                pheno_inc_override = [value]

        label = "_".join(label_parts)
        sf = SampleFilter(
            phenotype_include=(
                pheno_inc_override
                if pheno_inc_override is not None
                else list(base_sf.phenotype_include)
            ),
            phenotype_exclude=list(base_sf.phenotype_exclude),
            tech_include=(
                tech_inc_override
                if tech_inc_override is not None
                else list(base_sf.tech_include)
            ),
            tech_exclude=list(base_sf.tech_exclude),
            sex=sex_override if sex_override is not None else base_sf.sex,
            min_pass=base_sf.min_pass,
            min_observed=base_sf.min_observed,
            min_quality_evidence=base_sf.min_quality_evidence,
            fail_as_ref=base_sf.fail_as_ref,
        )
        groups.append((label, sf))

    return groups


def _discover_flat_buckets(flat_path, pos_start, pos_end):
    """Discover distinct bucket IDs in a flat parquet file.

    Uses integer division to avoid DuckDB float-cast rounding bug.
    """
    where_parts = []
    params = [str(flat_path)]
    if pos_start is not None:
        where_parts.append("pos >= ?")
        params.append(pos_start)
    if pos_end is not None:
        where_parts.append("pos <= ?")
        params.append(pos_end)
    where_clause = ("WHERE " + " AND ".join(where_parts)) if where_parts else ""
    sql = (
        f"SELECT DISTINCT CAST(pos AS BIGINT) // {BUCKET_SIZE} AS bid"
        f" FROM read_parquet(?)"
        f" {where_clause}"
        f" ORDER BY 1"
    )
    con = duckdb.connect()
    rows = con.execute(sql, params).fetchall()
    con.close()
    return [int(r[0]) for r in rows]


def _dump_bucket_worker(
    db_path, chrom, bucket_id, base_sf, groups, pos_start, pos_end,
    include_ac_zero=False,
):
    """Compute dump rows for a single bucket.

    Top-level picklable function — no closures. Safe for ProcessPoolExecutor.
    Pass include_ac_zero=True to retain positions with AC=0 (covered, no carriers).
    """
    from .query import QueryEngine

    engine = QueryEngine(db_path)
    _db = Path(db_path)

    bucket_start = bucket_id * BUCKET_SIZE
    bucket_end = (bucket_id + 1) * BUCKET_SIZE - 1

    # Resolve parquet path and WHERE clause
    if storage.chrom_layout(_db / "variants", chrom) == storage.PARTITIONED:
        parquet_file = storage.bucket_path(_db / "variants", chrom, bucket_id)
        if not parquet_file.exists():
            return []
        where_parts = []
        params = [str(parquet_file)]
        if pos_start is not None:
            where_parts.append("pos >= ?")
            params.append(pos_start)
        if pos_end is not None:
            where_parts.append("pos <= ?")
            params.append(pos_end)
        where_clause = ("WHERE " + " AND ".join(where_parts)) if where_parts else ""
    else:
        parquet_file = storage.flat_path(_db / "variants", chrom)
        if not parquet_file.exists():
            return []
        range_start = max(bucket_start, pos_start) if pos_start is not None else bucket_start
        range_end = min(bucket_end, pos_end) if pos_end is not None else bucket_end
        where_clause = "WHERE pos BETWEEN ? AND ?"
        params = [str(parquet_file), range_start, range_end]

    cols = ", ".join(engine._bitmap_cols(with_pos=True))
    sql = (
        f"SELECT {cols}"
        f" FROM read_parquet(?)"
        f" {where_clause}"
    )

    con = duckdb.connect()
    rows = con.execute(sql, params).fetchall()
    con.close()

    if not rows:
        return []

    base_bm = engine._build_sample_bitmap(base_sf)
    group_bms = [(label, engine._build_sample_bitmap(sf)) for label, sf in groups]

    # Position-level caches (eligible/AN computed once per position per bitmap)
    pos_cache: dict[int, tuple] = {}
    group_pos_cache: dict[tuple, tuple] = {}

    unpacked = [engine._unpack_bitmaps(row[3:]) for row in rows]
    sites = engine._site_evidence_by_pos([row[0] for row in rows], unpacked)

    result_rows = []
    for row, bitmaps in zip(rows, unpacked):
        pos, ref, alt = row[0], row[1], row[2]

        # Base eligible set; a position with no eligible sample is skipped
        if pos not in pos_cache:
            pos_cache[pos] = engine._compute_eligible(chrom, pos, base_bm)
        eligible = pos_cache[pos]

        if not eligible:
            continue

        s = engine._variant_stats(
            chrom, pos, eligible, bitmaps, sites[pos],
            base_sf.min_pass, base_sf.min_observed, base_sf.min_quality_evidence,
            base_sf.fail_as_ref,
        )

        if s.AC == 0 and not include_ac_zero:
            continue  # main row filter

        out_row: dict = {"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}
        out_row.update(_stat_values(s, len(eligible)))

        # Per-group columns
        for g_idx, (label, g_bm) in enumerate(group_bms):
            g_sf = groups[g_idx][1]
            cache_key = (g_idx, pos)
            if cache_key not in group_pos_cache:
                group_pos_cache[cache_key] = engine._compute_eligible(chrom, pos, g_bm)
            g_eligible = group_pos_cache[cache_key]

            g = engine._variant_stats(
                chrom, pos, g_eligible, bitmaps, sites[pos],
                g_sf.min_pass, g_sf.min_observed, g_sf.min_quality_evidence,
                g_sf.fail_as_ref,
            )

            for field, value in _stat_values(g, len(g_eligible)).items():
                out_row[f"{field}_{label}"] = value

        result_rows.append(out_row)

    return sorted(result_rows, key=lambda r: (r["pos"], r["alt"]))


def dump_database(
    engine,
    output,
    base_sf: SampleFilter,
    groups: list,
    chrom_filter: str | None,
    pos_start: int | None,
    pos_end: int | None,
    n_workers: int | None,
    include_ac_zero: bool = False,
) -> dict:
    """Export variants to CSV.

    Args:
        engine: QueryEngine instance
        output: file-like object or path string (None = stdout)
        base_sf: base SampleFilter defining the global population
        groups: list of (label, SampleFilter) from _build_groups
        chrom_filter: restrict to this chromosome (None = all)
        pos_start: 1-based start (inclusive), requires chrom_filter
        pos_end: 1-based end (inclusive), requires chrom_filter
        n_workers: parallel workers (None = cpu_count)
        include_ac_zero: if True, include variants with AC=0

    Returns:
        {"n_rows": int, "n_buckets": int, "n_chroms": int}

    Raises:
        DumpColumnError: if the group labels produce duplicate column names
            (see dump_columns).
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    # Validate the header before any work, since a column collision is a usage error
    fieldnames = dump_columns([label for label, _ in groups])

    db_path = str(engine._db)
    variants_dir = engine._db / "variants"

    # Determine chromosomes to process in canonical order
    if chrom_filter is not None:
        chroms = [normalize_chrom(chrom_filter)]
    else:
        # All chroms that have data
        available = set()
        for chrom in ALL_CHROMS:
            if storage.variant_parquet_glob(variants_dir, chrom) is not None:
                available.add(chrom)
        chroms = [c for c in ALL_CHROMS if c in available]

    # Build work units: (chrom, bucket_id) in genomic order
    work_units: list[tuple[str, int]] = []
    for chrom in chroms:
        if storage.chrom_layout(variants_dir, chrom) == storage.PARTITIONED:
            for bid in storage.existing_bucket_ids(variants_dir, chrom):
                # Filter by region if specified
                if pos_start is not None and (bid + 1) * BUCKET_SIZE - 1 < pos_start:
                    continue
                if pos_end is not None and bid * BUCKET_SIZE > pos_end:
                    continue
                work_units.append((chrom, bid))
        else:
            flat_path = storage.flat_path(variants_dir, chrom)
            if not flat_path.exists():
                continue
            bucket_ids = _discover_flat_buckets(flat_path, pos_start, pos_end)
            for bid in bucket_ids:
                work_units.append((chrom, bid))

    if not work_units:
        logger.info("[dump] No data found for the given filters.")
        return {"n_rows": 0, "n_buckets": 0, "n_chroms": 0}

    n_units = len(work_units)
    if n_workers is None:
        effective = max(1, min(os.cpu_count() or 1, n_units))
    else:
        effective = max(1, min(n_workers, n_units))

    logger.info(
        "[dump] Processing %d bucket(s) across %d chrom(s) with %d worker(s).",
        n_units,
        len(chroms),
        effective,
    )

    # Collect results in work-unit order (preserves genomic order)
    unit_results: dict[int, list[dict]] = {}

    if effective == 1:
        for idx, (chrom, bid) in enumerate(work_units):
            unit_results[idx] = _dump_bucket_worker(
                db_path, chrom, bid, base_sf, groups, pos_start, pos_end,
                include_ac_zero,
            )
            logger.debug("  [dump] %s/bucket_%d: done", chrom, bid)
    else:
        with ProcessPoolExecutor(max_workers=effective) as executor:
            futures = {
                executor.submit(
                    _dump_bucket_worker,
                    db_path, chrom, bid, base_sf, groups,
                    pos_start, pos_end, include_ac_zero,
                ): idx
                for idx, (chrom, bid) in enumerate(work_units)
            }
            for future in as_completed(futures):
                idx = futures[future]
                chrom, bid = work_units[idx]
                unit_results[idx] = future.result()
                logger.debug("  [dump] %s/bucket_%d: done", chrom, bid)

    # Write CSV in genomic order
    close_output = False
    if output is None:
        fh = sys.stdout
    elif hasattr(output, "write"):
        fh = output
    else:
        fh = open(output, "w", newline="")
        close_output = True

    try:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()

        n_rows = 0
        n_chroms_written = set()
        for idx in range(len(work_units)):
            rows = unit_results.get(idx, [])
            chrom, _ = work_units[idx]
            for row in rows:
                # Convert None → "" for CSV
                csv_row = {
                    k: ("" if v is None else v)
                    for k, v in row.items()
                    if k in fieldnames
                }
                writer.writerow(csv_row)
                n_rows += 1
                n_chroms_written.add(chrom)
    finally:
        if close_output:
            fh.close()

    logger.info("[dump] Complete: %d row(s) written.", n_rows)
    return {"n_rows": n_rows, "n_buckets": n_units, "n_chroms": len(n_chroms_written)}

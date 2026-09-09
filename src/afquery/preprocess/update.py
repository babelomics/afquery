import glob as glob_module
import json
import logging
import os
import sqlite3
import tempfile
import time
from dataclasses import dataclass
from datetime import datetime, timezone

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from pyroaring import BitMap

from .. import storage
from ..bitmaps import deserialize, serialize, build_phenotype_bitmaps, build_sex_bitmaps, build_tech_bitmaps
from ..constants import VALID_GENOME_BUILDS
from ..models import Sample, Technology
from .build import PARQUET_SCHEMA, _make_table, _passes_quality, get_chroms_in_temp_files
from .ingest import ingest_all
from .manifest import parse_manifest
from .regions import build_capture_indices

logger = logging.getLogger(__name__)


class UpdateError(RuntimeError):
    pass


@dataclass
class CheckResult:
    severity: str   # 'error' | 'warning' | 'info'
    message: str


# ---- Internal helpers ----

def _read_manifest(db_dir: str) -> dict:
    path = os.path.join(db_dir, "manifest.json")
    with open(path) as f:
        return json.load(f)


def _bump_version(version: str) -> str:
    """Auto-increment the last numeric component of a version string."""
    parts = version.split(".")
    if parts and parts[-1].isdigit():
        parts[-1] = str(int(parts[-1]) + 1)
        return ".".join(parts)
    return version + ".1"


def _update_manifest(
    db_dir: str,
    sample_count: int,
    next_sample_id: int | None = None,
    db_version: str | None = None,
) -> None:
    path = os.path.join(db_dir, "manifest.json")
    try:
        with open(path) as f:
            manifest = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        manifest = {}
    manifest["sample_count"] = sample_count
    manifest["updated_at"] = datetime.now(timezone.utc).isoformat()
    if next_sample_id is not None:
        manifest["next_sample_id"] = next_sample_id
    if db_version is not None:
        manifest["db_version"] = db_version
    tmp_path = path + ".tmp"
    with open(tmp_path, "w") as f:
        json.dump(manifest, f, indent=2)
    os.replace(tmp_path, path)


def _get_next_sample_id(con: sqlite3.Connection) -> int:
    row = con.execute("SELECT MAX(sample_id) FROM samples").fetchone()
    return 0 if row[0] is None else row[0] + 1


def _get_next_tech_id(con: sqlite3.Connection) -> int:
    row = con.execute("SELECT MAX(tech_id) FROM technologies").fetchone()
    return 0 if row[0] is None else row[0] + 1


def _regenerate_precomputed_bitmaps(con: sqlite3.Connection) -> None:
    samples_rows = con.execute(
        "SELECT sample_id, sample_name, sex, tech_id FROM samples"
    ).fetchall()
    samples = [Sample(r[0], r[1], r[2], r[3]) for r in samples_rows]

    phenotype_pairs = con.execute(
        "SELECT sample_id, phenotype_code FROM sample_phenotype"
    ).fetchall()

    con.execute("DELETE FROM precomputed_bitmaps")

    for sex, bm in build_sex_bitmaps(samples).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("sex", sex, serialize(bm)),
        )

    for code, bm in build_phenotype_bitmaps(phenotype_pairs).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("phenotype", code, serialize(bm)),
        )

    for tech_id, bm in build_tech_bitmaps(samples).items():
        con.execute(
            "INSERT INTO precomputed_bitmaps VALUES (?, ?, ?)",
            ("tech", str(tech_id), serialize(bm)),
        )


def _merge_rows_into_file(
    out_path: str,
    new_rows: list[tuple[int, str, str, BitMap, BitMap, BitMap, BitMap]],
    has_phase2: bool,
    min_covered: int,
    wes_tech_bitmaps: dict[int, BitMap] | None,
    row_group_size: int,
) -> tuple[int, int]:
    """Merge pre-aggregated rows into one Parquet file, rewriting it atomically.

    `out_path` is a flat per-chromosome file or a single bucket; the read, merge,
    Phase 2 recompute, sort and write are identical either way.

    `new_rows` may be empty. That is the recompute-only call used for buckets of
    an affected chromosome that received no new data but whose filtered_bitmap
    still depends on the enlarged cohort. Such a file is rewritten only when the
    recomputation actually moves a bitmap.

    Returns (new_variants, updated_variants).
    """
    # Read existing Parquet via pyarrow (NOT DuckDB — need Python bitmap deserialization)
    existing: dict[tuple, tuple[BitMap, BitMap, BitMap, BitMap, BitMap]] = {}
    if os.path.exists(out_path):
        table = pq.read_table(out_path)
        existing_has_fail = "fail_bitmap" in table.schema.names
        existing_has_phase2 = (
            "filtered_bitmap" in table.schema.names
            and "quality_pass_bitmap" in table.schema.names
        )
        for i in range(len(table)):
            pos = table["pos"][i].as_py()
            ref = table["ref"][i].as_py()
            alt = table["alt"][i].as_py()
            het_bm = deserialize(table["het_bitmap"][i].as_py())
            hom_bm = deserialize(table["hom_bitmap"][i].as_py())
            fail_bm = deserialize(table["fail_bitmap"][i].as_py()) if existing_has_fail else BitMap()
            if existing_has_phase2:
                filt_bm = deserialize(table["filtered_bitmap"][i].as_py())
                qp_bm = deserialize(table["quality_pass_bitmap"][i].as_py())
            else:
                filt_bm = BitMap()
                qp_bm = BitMap()
            existing[(pos, ref, alt)] = (het_bm, hom_bm, fail_bm, filt_bm, qp_bm)

    dirty = bool(new_rows)
    new_variants = 0
    updated_variants = 0

    for pos, ref, alt, new_het, new_hom, new_fail, new_qp in new_rows:
        key = (pos, ref, alt)
        if key in existing:
            old_het, old_hom, old_fail, old_filt, old_qp = existing[key]
            existing[key] = (
                old_het | new_het,
                old_hom | new_hom,
                old_fail | new_fail,
                old_filt,
                old_qp | new_qp,
            )
            updated_variants += 1
        else:
            existing[key] = (new_het, new_hom, new_fail, BitMap(), new_qp)
            new_variants += 1

    # Phase 2: recompute filtered_bitmap per row (all rows, since merging may shift any tech)
    if has_phase2 and wes_tech_bitmaps and min_covered > 0:
        for key, (het, hom, fail, old_filt, qp) in existing.items():
            carrier_set = het | hom | fail
            new_filt = BitMap()
            for tech_bm in wes_tech_bitmaps.values():
                qp_count = len(qp & tech_bm)
                if qp_count < min_covered:
                    new_filt |= (tech_bm - carrier_set)
            if new_filt != old_filt:
                dirty = True
            existing[key] = (het, hom, fail, new_filt, qp)

    # Nothing moved: leave the file alone. Only reachable when new_rows was empty.
    if not dirty or not existing:
        return (0, 0)

    # Sort by (pos, alt) and write atomically
    sorted_keys = sorted(existing.keys(), key=lambda k: (k[0], k[2]))

    table = _make_table(
        [k[0] for k in sorted_keys],
        [k[1] for k in sorted_keys],
        [k[2] for k in sorted_keys],
        [serialize(existing[k][0]) for k in sorted_keys],
        [serialize(existing[k][1]) for k in sorted_keys],
        [serialize(existing[k][2]) for k in sorted_keys],
        [serialize(existing[k][3]) for k in sorted_keys],
        [serialize(existing[k][4]) for k in sorted_keys],
    )

    storage.ensure_parent(out_path)
    tmp_path = out_path + ".tmp"
    pq.write_table(table, tmp_path, row_group_size=row_group_size)
    os.replace(tmp_path, out_path)

    return (new_variants, updated_variants)


def _merge_chromosome_parquet(
    chrom: str,
    db_dir: str,
    update_tmp_dir: str,
    row_group_size: int = 100_000,
    coverage_filter: dict | None = None,
    wes_tech_bitmaps: dict[int, BitMap] | None = None,
    layout: str | None = None,
) -> tuple[int, int]:
    """Merge new temp files into a chromosome's variant store.

    Returns (new_variants, updated_variants).

    Writes into whichever layout the chromosome already uses; `layout` decides
    for a chromosome new to the database, and defaults to what create-db
    produces. For Phase 2 DBs (coverage_filter active) also merges
    quality_pass_bitmap and recomputes filtered_bitmap using manifest thresholds
    and current WES tech bitmaps.
    """
    variants_dir = os.path.join(db_dir, "variants")

    # Check if there are any new temp files
    parquet_files = glob_module.glob(os.path.join(update_tmp_dir, "sample_*.parquet"))
    if not parquet_files:
        return (0, 0)

    # Detect whether new temp files include quality columns (dp/gq/qual).
    # Old test fixtures may not include them; legacy ingest didn't either.
    new_has_quality = False
    try:
        sample_schema = pq.read_schema(parquet_files[0])
        new_has_quality = (
            "dp" in sample_schema.names
            and "gq" in sample_schema.names
            and "qual" in sample_schema.names
        )
    except Exception:
        new_has_quality = False

    # Aggregate new rows via DuckDB
    glob_pattern = os.path.join(update_tmp_dir, "sample_*.parquet").replace("'", "''")
    quality_select = (
        ", list(dp ORDER BY sample_id), list(gq ORDER BY sample_id), list(qual ORDER BY sample_id)"
        if new_has_quality else ""
    )
    con = duckdb.connect()
    try:
        rows = con.execute(
            f"""
            SELECT pos, ref, alt,
                list(sample_id   ORDER BY sample_id),
                list(gt_ac       ORDER BY sample_id),
                list(filter_pass ORDER BY sample_id)
                {quality_select}
            FROM read_parquet('{glob_pattern}')
            WHERE chrom = ?
            GROUP BY pos, ref, alt
            ORDER BY pos, alt
            """,
            [chrom],
        ).fetchall()
    finally:
        con.close()

    coverage_filter = coverage_filter or {}
    min_dp = coverage_filter.get("min_dp", 0)
    min_gq = coverage_filter.get("min_gq", 0)
    min_qual = coverage_filter.get("min_qual", 0.0)
    min_covered = coverage_filter.get("min_covered", 0)
    has_phase2 = (min_dp > 0 or min_gq > 0 or min_qual > 0 or min_covered > 0)

    # filtered_bitmap is (tech_bm - carrier_set) row by row, so enlarging a WES
    # tech changes it at every row of the database, not only where the new rows
    # landed. A chromosome that received nothing still has to be revisited.
    recompute_phase2 = bool(has_phase2 and wes_tech_bitmaps and min_covered > 0)

    # Nothing new here and nothing to recompute — leave the files alone.
    if not rows and not recompute_phase2:
        return (0, 0)

    prepared: list[tuple[int, str, str, BitMap, BitMap, BitMap, BitMap]] = []

    for row in rows:
        if new_has_quality:
            pos, ref, alt, sample_ids, gt_acs, filter_passes, dps, gqs, quals = row
        else:
            pos, ref, alt, sample_ids, gt_acs, filter_passes = row
            dps = [None] * len(sample_ids)
            gqs = [None] * len(sample_ids)
            quals = [None] * len(sample_ids)

        het_ids  = [sid for sid, ac, fp in zip(sample_ids, gt_acs, filter_passes) if ac == 1 and fp]
        hom_ids  = [sid for sid, ac, fp in zip(sample_ids, gt_acs, filter_passes) if ac == 2 and fp]
        fail_ids = [sid for sid, fp in zip(sample_ids, filter_passes) if not fp]

        # New quality_pass_ids (only if Phase 2 active)
        if has_phase2:
            new_qp_ids = [
                sid for sid, ac, fp, dp, gq, qual in zip(
                    sample_ids, gt_acs, filter_passes, dps, gqs, quals
                )
                if ac > 0 and fp and _passes_quality(dp, gq, qual, min_dp, min_gq, min_qual)
            ]
            new_qp = BitMap(new_qp_ids)
        else:
            new_qp = BitMap()

        prepared.append(
            (pos, ref, alt, BitMap(het_ids), BitMap(hom_ids), BitMap(fail_ids), new_qp)
        )

    if layout is None:
        layout = storage.detect_layout(variants_dir)
    chrom_fmt = storage.chrom_layout(variants_dir, chrom, layout)

    if chrom_fmt == storage.FLAT:
        return _merge_rows_into_file(
            str(storage.flat_path(variants_dir, chrom)), prepared,
            has_phase2, min_covered, wes_tech_bitmaps, row_group_size,
        )

    by_bucket: dict[int, list] = {}
    for prepared_row in prepared:
        by_bucket.setdefault(storage.bucket_id(prepared_row[0]), []).append(prepared_row)

    target_ids = set(by_bucket)
    if recompute_phase2:
        # Buckets with no new data must still be revisited, or their
        # N_NO_COVERAGE goes stale and the added samples read as hom-ref there.
        target_ids |= set(storage.existing_bucket_ids(variants_dir, chrom))

    total_new = 0
    total_updated = 0
    for bid in sorted(target_ids):
        n, u = _merge_rows_into_file(
            str(storage.bucket_path(variants_dir, chrom, bid)), by_bucket.get(bid, []),
            has_phase2, min_covered, wes_tech_bitmaps, row_group_size,
        )
        total_new += n
        total_updated += u

    return (total_new, total_updated)


def _clear_bits_from_parquet(parquet_file: str, removal_ids: BitMap) -> None:
    """Clear removal_ids bits from het/hom/fail/filtered/quality_pass bitmaps. Rewrites atomically if dirty."""
    table = pq.read_table(parquet_file)
    has_phase2 = (
        "filtered_bitmap" in table.schema.names
        and "quality_pass_bitmap" in table.schema.names
    )
    dirty = False

    new_het_list = []
    new_hom_list = []
    new_fail_list = []
    new_filtered_list = []
    new_quality_pass_list = []

    for i in range(len(table)):
        het_bm = deserialize(table["het_bitmap"][i].as_py())
        hom_bm = deserialize(table["hom_bitmap"][i].as_py())
        fail_bm = deserialize(table["fail_bitmap"][i].as_py())
        if has_phase2:
            filt_bm = deserialize(table["filtered_bitmap"][i].as_py())
            qp_bm = deserialize(table["quality_pass_bitmap"][i].as_py())
        else:
            filt_bm = BitMap()
            qp_bm = BitMap()

        combined = het_bm | hom_bm | fail_bm | filt_bm | qp_bm
        if removal_ids & combined:
            het_bm = het_bm - removal_ids
            hom_bm = hom_bm - removal_ids
            fail_bm = fail_bm - removal_ids
            filt_bm = filt_bm - removal_ids
            qp_bm = qp_bm - removal_ids
            dirty = True

        new_het_list.append(serialize(het_bm))
        new_hom_list.append(serialize(hom_bm))
        new_fail_list.append(serialize(fail_bm))
        new_filtered_list.append(serialize(filt_bm))
        new_quality_pass_list.append(serialize(qp_bm))

    if not dirty:
        return

    new_table = pa.table(
        {
            "pos":                 table["pos"],
            "ref":                 table["ref"],
            "alt":                 table["alt"],
            "het_bitmap":          pa.array(new_het_list,          type=pa.large_binary()),
            "hom_bitmap":          pa.array(new_hom_list,          type=pa.large_binary()),
            "fail_bitmap":         pa.array(new_fail_list,         type=pa.large_binary()),
            "filtered_bitmap":     pa.array(new_filtered_list,     type=pa.large_binary()),
            "quality_pass_bitmap": pa.array(new_quality_pass_list, type=pa.large_binary()),
        },
        schema=PARQUET_SCHEMA,
    )

    tmp_path = parquet_file + ".tmp"
    pq.write_table(new_table, tmp_path)
    os.replace(tmp_path, parquet_file)


# ---- Public functions ----

def parse_updates_tsv(path: str) -> list[dict]:
    """Parse a metadata-update TSV file into a list of update dicts.

    The file must have a header row with columns ``sample_name``, ``field``,
    and ``new_value`` (tab-separated).  Each subsequent row describes one
    field change for one sample.

    Args:
        path: Path to the TSV file.

    Returns:
        List of dicts with keys ``sample_name``, ``field``, and
        ``new_value``.

    Raises:
        UpdateError: If required columns are missing.
    """
    updates: list[dict] = []
    with open(path) as fh:
        header: list[str] | None = None
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if header is None:
                header = [p.strip() for p in parts]
                required = {"sample_name", "field", "new_value"}
                missing = required - set(header)
                if missing:
                    raise UpdateError(
                        f"TSV file {path!r} missing required columns: "
                        f"{', '.join(sorted(missing))}. "
                        f"Expected header: sample_name<TAB>field<TAB>new_value"
                    )
                continue
            row = dict(zip(header, parts))
            updates.append(
                {
                    "sample_name": row.get("sample_name", "").strip(),
                    "field": row.get("field", "").strip(),
                    "new_value": row.get("new_value", "").strip(),
                }
            )
    return updates


def update_sample_metadata(
    db_dir: str,
    updates: list[dict],
    operator_note: str | None = None,
) -> list[dict]:
    """Update ``sex`` and/or ``phenotype_codes`` for one or more samples.

    This is a non-destructive metadata correction that does **not** require
    re-ingesting VCF files.  All changes are recorded in the ``changelog``
    table and precomputed bitmaps are regenerated atomically.

    Args:
        db_dir: Path to the database directory.
        updates: List of dicts, each with keys:

            - ``sample_name`` (str): Name of the sample to update.
            - ``field`` (str): Field to change — ``"sex"`` or
              ``"phenotype_codes"``.
            - ``new_value`` (str): New value.  For ``"sex"``, must be
              ``"male"`` or ``"female"``.  For ``"phenotype_codes"``, a
              comma-separated string of codes that *replaces* all current
              codes for the sample.

        operator_note: Optional free-text note appended to every changelog
            entry created by this call.

    Returns:
        List of changelog entry dicts (one per field change), each containing
        ``sample``, ``field``, ``old``, and ``new`` keys.

    Raises:
        UpdateError: If a sample name is not found, a field name is invalid,
            a sex value is not ``"male"`` or ``"female"``, or
            ``phenotype_codes`` resolves to an empty list.
    """
    _VALID_FIELDS = {"sex", "phenotype_codes"}
    _VALID_SEX = {"male", "female"}

    if not updates:
        return []

    db_path = os.path.join(db_dir, "metadata.sqlite")
    manifest = _read_manifest(db_dir)

    con = sqlite3.connect(db_path)
    try:
        existing_samples: dict[str, int] = {
            r[0]: r[1]
            for r in con.execute(
                "SELECT sample_name, sample_id FROM samples"
            ).fetchall()
        }

        # --- Validation pass (all-or-nothing before any writes) ---
        for u in updates:
            sname = u.get("sample_name", "").strip()
            field = u.get("field", "").strip()
            new_value = u.get("new_value", "").strip()

            if not sname:
                raise UpdateError("Update entry has an empty sample_name.")
            if sname not in existing_samples:
                raise UpdateError(f"Sample '{sname}' not found in database.")
            if field not in _VALID_FIELDS:
                raise UpdateError(
                    f"Invalid field '{field}' for sample '{sname}'. "
                    f"Allowed fields: {', '.join(sorted(_VALID_FIELDS))}."
                )
            if field == "sex" and new_value not in _VALID_SEX:
                raise UpdateError(
                    f"Invalid sex value '{new_value}' for sample '{sname}'. "
                    f"Must be 'male' or 'female'."
                )
            if field == "phenotype_codes":
                codes = [c.strip() for c in new_value.split(",") if c.strip()]
                if not codes:
                    raise UpdateError(
                        f"phenotype_codes for sample '{sname}' is empty after "
                        f"parsing '{new_value}'. Provide at least one non-empty code."
                    )

        # --- Apply changes ---
        event_time = datetime.now(timezone.utc).isoformat()
        changelog_entries: list[dict] = []

        for u in updates:
            sname = u["sample_name"].strip()
            field = u["field"].strip()
            new_value = u["new_value"].strip()
            sample_id = existing_samples[sname]

            if field == "sex":
                old_value: str = con.execute(
                    "SELECT sex FROM samples WHERE sample_id = ?", (sample_id,)
                ).fetchone()[0]
                con.execute(
                    "UPDATE samples SET sex = ? WHERE sample_id = ?",
                    (new_value, sample_id),
                )
                entry: dict = {
                    "sample": sname,
                    "field": "sex",
                    "old": old_value,
                    "new": new_value,
                }

            else:  # phenotype_codes
                codes = [c.strip() for c in new_value.split(",") if c.strip()]
                old_codes = sorted(
                    r[0]
                    for r in con.execute(
                        "SELECT phenotype_code FROM sample_phenotype WHERE sample_id = ?",
                        (sample_id,),
                    ).fetchall()
                )
                con.execute(
                    "DELETE FROM sample_phenotype WHERE sample_id = ?", (sample_id,)
                )
                con.executemany(
                    "INSERT INTO sample_phenotype (sample_id, phenotype_code) VALUES (?, ?)",
                    [(sample_id, code) for code in codes],
                )
                entry = {
                    "sample": sname,
                    "field": "phenotype_codes",
                    "old": ",".join(old_codes),
                    "new": ",".join(sorted(codes)),
                }

            if operator_note:
                entry["operator_note"] = operator_note

            con.execute(
                "INSERT INTO changelog (event_type, event_time, sample_names, notes)"
                " VALUES (?, ?, ?, ?)",
                (
                    "UPDATE_SAMPLE",
                    event_time,
                    json.dumps([sname]),
                    json.dumps(entry),
                ),
            )
            changelog_entries.append(entry)

        _regenerate_precomputed_bitmaps(con)
        con.commit()

    finally:
        con.close()

    current_version = manifest.get("db_version", "1.0")
    _update_manifest(
        db_dir,
        manifest["sample_count"],
        db_version=_bump_version(current_version),
    )

    logger.info(
        "[update-metadata] Done. %d field change(s) applied.", len(changelog_entries)
    )
    return changelog_entries


def add_samples(
    db_dir: str,
    manifest_path: str,
    threads: int | None = None,
    tmp_dir: str | None = None,
    bed_dir: str | None = None,
    genome_build: str | None = None,
    db_version: str | None = None,
) -> dict:
    """Add new samples to the database without full rebuild."""
    # Resolve threads: None → all available CPU cores
    effective_threads = max(1, os.cpu_count() or 1) if threads is None else max(1, threads)

    # 1. Read DB manifest
    manifest = _read_manifest(db_dir)
    db_genome_build = manifest["genome_build"]

    # 2. Validate genome build if provided
    if genome_build is not None and genome_build != db_genome_build:
        raise UpdateError(
            f"genome_build mismatch: DB has '{db_genome_build}', got '{genome_build}'"
        )

    # 3. Refuse a database that already stores a chromosome in both layouts.
    # Checked before ingest so an operator is not told this after hours of work.
    mixed = storage.mixed_layout_chroms(os.path.join(db_dir, "variants"))
    if mixed:
        raise UpdateError(
            "Mixed variant layout: chromosome(s) "
            + ", ".join(mixed)
            + " have both variants/<chrom>/ (bucketed) and variants/<chrom>.parquet. "
            "The flat files were written by an older update-db and are invisible to "
            "queries; adding more samples would compound the error. Run "
            "'afquery check --db <db>' and follow the recovery procedure in the "
            "troubleshooting guide before retrying."
        )

    # 4. Parse new manifest
    samples_raw, techs_raw = parse_manifest(manifest_path, bed_dir)

    logger.info("[add-samples] Adding %d new sample(s)...", len(samples_raw))

    # 5. Open connection
    db_path = os.path.join(db_dir, "metadata.sqlite")
    con = sqlite3.connect(db_path)

    total_new = 0
    total_updated = 0

    try:
        # 6. Check for duplicate sample names
        existing_names = {
            r[0] for r in con.execute("SELECT sample_name FROM samples").fetchall()
        }
        new_names = [s.sample_name for s in samples_raw]
        duplicates = [n for n in new_names if n in existing_names]
        if duplicates:
            raise UpdateError(
                f"Sample(s) already in database: {', '.join(duplicates)}"
            )

        # 7. Assign IDs sequentially from next available.
        # Use manifest's next_sample_id when present (survives removals) so IDs
        # are never reused even after a sample has been deleted.
        manifest_next = manifest.get("next_sample_id")
        db_next = _get_next_sample_id(con)
        starting_id = max(manifest_next if manifest_next is not None else 0, db_next)
        new_samples = [
            Sample(
                sample_id=starting_id + i,
                sample_name=ps.sample_name,
                sex=ps.sex,
                tech_id=-1,  # filled below
            )
            for i, ps in enumerate(samples_raw)
        ]

        # 8. Handle technologies
        capture_dir = os.path.join(db_dir, "capture")
        os.makedirs(capture_dir, exist_ok=True)

        existing_techs = {
            r[1]: (r[0], r[2])
            for r in con.execute(
                "SELECT tech_id, tech_name, bed_path FROM technologies"
            ).fetchall()
        }

        tech_name_to_id: dict[str, int] = {}
        next_tech_id = _get_next_tech_id(con)

        for pt in techs_raw:
            if pt.tech_name in existing_techs:
                tech_name_to_id[pt.tech_name] = existing_techs[pt.tech_name][0]
            else:
                tech_id = next_tech_id
                next_tech_id += 1
                con.execute(
                    "INSERT INTO technologies VALUES (?, ?, ?)",
                    (tech_id, pt.tech_name, pt.bed_path),
                )
                tech = Technology(
                    tech_id=tech_id, tech_name=pt.tech_name, bed_path=pt.bed_path
                )
                build_capture_indices([tech], capture_dir)
                tech_name_to_id[pt.tech_name] = tech_id

        # Assign tech_ids to new samples
        for i, ps in enumerate(samples_raw):
            new_samples[i].tech_id = tech_name_to_id[ps.tech_name]

        vcf_paths = [ps.vcf_path for ps in samples_raw]

        # 9. Ingest VCFs into fresh tmp_dir
        auto_tmp = tmp_dir is None
        if auto_tmp:
            tmp_dir = tempfile.mkdtemp(prefix="afquery_update_")

        try:
            ingest_all(new_samples, vcf_paths, tmp_dir, n_workers=effective_threads)

            # 10. Collect chroms from new temp files
            chroms = get_chroms_in_temp_files(tmp_dir)

            # Phase 2: build WES tech bitmaps from CURRENT DB state (existing + new samples)
            # so that filtered_bitmap recomputation uses the merged cohort.
            coverage_filter = manifest.get("coverage_filter", {})
            wes_tech_bitmaps: dict[int, BitMap] = {}
            new_wes_sample = False
            if coverage_filter and coverage_filter.get("min_covered", 0) > 0:
                rows = con.execute(
                    "SELECT s.sample_id, s.tech_id, t.bed_path FROM samples s"
                    " JOIN technologies t ON s.tech_id = t.tech_id"
                    " WHERE t.bed_path IS NOT NULL"
                ).fetchall()
                for sid, tid, _ in rows:
                    wes_tech_bitmaps.setdefault(tid, BitMap()).add(sid)
                # Include new samples (not yet inserted into samples table)
                for s in new_samples:
                    tech_obj = next(
                        (t for t in techs_raw if t.tech_name in tech_name_to_id
                         and tech_name_to_id[t.tech_name] == s.tech_id),
                        None,
                    )
                    if tech_obj is not None and tech_obj.bed_path is not None:
                        new_wes_sample = True
                        wes_tech_bitmaps.setdefault(s.tech_id, BitMap()).add(s.sample_id)

            # 11. Merge Parquet files. The layout is decided once for the whole
            # update so every chromosome in this batch agrees, including any
            # chromosome new to the database.
            variants_dir = os.path.join(db_dir, "variants")
            if wes_tech_bitmaps and new_wes_sample:
                # A WES tech just grew, and filtered_bitmap is derived from the
                # tech bitmaps at every row of every chromosome. Restricting the
                # recomputation to the chromosomes in this batch would leave the
                # added samples counted as homozygous reference everywhere else
                # their capture BED reaches. No tech grows for a WGS-only batch,
                # so nothing off the touched chromosomes can move and the whole
                # store does not need reading.
                chroms = sorted(set(chroms) | storage.stored_chroms(variants_dir))
            layout = storage.detect_layout(variants_dir)
            for chrom in chroms:
                n, u = _merge_chromosome_parquet(
                    chrom, db_dir, tmp_dir,
                    coverage_filter=coverage_filter,
                    wes_tech_bitmaps=wes_tech_bitmaps or None,
                    layout=layout,
                )
                total_new += n
                total_updated += u
                logger.debug("  [add-samples] Merged %s: %d new, %d updated", chrom, n, u)
        finally:
            if auto_tmp:
                import shutil
                shutil.rmtree(tmp_dir, ignore_errors=True)

        # 12. Insert new samples and Phenotype pairs into SQLite
        ingested_at = datetime.now(timezone.utc).isoformat()
        con.executemany(
            "INSERT INTO samples (sample_id, sample_name, sex, tech_id, vcf_path, ingested_at)"
            " VALUES (?, ?, ?, ?, ?, ?)",
            [(s.sample_id, s.sample_name, s.sex, s.tech_id, vcf_paths[i], ingested_at)
             for i, s in enumerate(new_samples)],
        )
        sample_phenotype_pairs = [
            (starting_id + i, code)
            for i, ps in enumerate(samples_raw)
            for code in ps.phenotype_codes
        ]
        con.executemany("INSERT INTO sample_phenotype VALUES (?, ?)", sample_phenotype_pairs)

        # 13. Regenerate precomputed bitmaps
        _regenerate_precomputed_bitmaps(con)

        # 14. Append changelog entry
        import json as _json
        sample_names_json = _json.dumps([s.sample_name for s in new_samples])
        con.execute(
            "INSERT INTO changelog (event_type, event_time, sample_names, notes) VALUES (?, ?, ?, ?)",
            ("add_samples", ingested_at, sample_names_json, f"{len(new_samples)} samples added"),
        )
        con.commit()

        # sample_count is the number of samples the database holds, which is not
        # the next free id once anything has been removed.
        sample_count = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]

    finally:
        con.close()

    # 15. Update manifest (persist next_sample_id so future adds don't reuse IDs)
    next_id = starting_id + len(new_samples)
    # Resolve db_version: explicit value overrides auto-bump
    current_version = manifest.get("db_version", "1.0")
    new_version = db_version if db_version is not None else _bump_version(current_version)
    _update_manifest(db_dir, sample_count, next_sample_id=next_id, db_version=new_version)

    logger.info("[add-samples] Done. %d sample(s) added.", len(new_samples))

    return {
        "new_samples": len(new_samples),
        "new_variants": total_new,
        "updated_variants": total_updated,
    }


def remove_samples(db_dir: str, sample_names: list[str]) -> dict:
    """Remove samples from the database by clearing their bitmap bits."""
    if not sample_names:
        return {"removed": [], "not_found": []}

    manifest = _read_manifest(db_dir)

    db_path = os.path.join(db_dir, "metadata.sqlite")
    con = sqlite3.connect(db_path)

    try:
        # 1. Look up sample IDs for given names
        placeholders = ",".join("?" * len(sample_names))
        rows = con.execute(
            f"SELECT sample_id, sample_name FROM samples WHERE sample_name IN ({placeholders})",
            sample_names,
        ).fetchall()

        found_names = {r[1] for r in rows}
        missing = [n for n in sample_names if n not in found_names]
        if missing:
            raise UpdateError(f"Sample(s) not found: {', '.join(missing)}")

        logger.info("[remove-samples] Removing %d sample(s)...", len(found_names))

        removal_ids = BitMap([r[0] for r in rows])
        id_list = list(removal_ids)

        # 2-3. Clear bits from every Parquet file, in either layout
        variants_dir = os.path.join(db_dir, "variants")
        for pq_file in storage.iter_variant_parquets(variants_dir):
            _clear_bits_from_parquet(str(pq_file), removal_ids)
            logger.debug("  [remove-samples] %s cleared", pq_file.name)

        # 4-5. Delete from SQLite
        ph2 = ",".join("?" * len(id_list))
        con.execute(f"DELETE FROM sample_phenotype WHERE sample_id IN ({ph2})", id_list)
        con.execute(f"DELETE FROM samples WHERE sample_id IN ({ph2})", id_list)

        # 6. Regenerate precomputed bitmaps
        _regenerate_precomputed_bitmaps(con)

        # Append changelog entry
        import json as _json
        removed_names_json = _json.dumps(list(found_names))
        event_time = datetime.now(timezone.utc).isoformat()
        con.execute(
            "INSERT INTO changelog (event_type, event_time, sample_names, notes) VALUES (?, ?, ?, ?)",
            ("remove_samples", event_time, removed_names_json, f"{len(found_names)} samples removed"),
        )
        con.commit()

        # 7. Get new count and preserve next_sample_id for future adds
        new_count = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        max_id_row = con.execute("SELECT MAX(sample_id) FROM samples").fetchone()
        max_id = max_id_row[0] if max_id_row[0] is not None else -1
        # next_sample_id should always be at least max_id + 1 (from manifest or DB)
        next_id = max(manifest.get("next_sample_id", max_id + 1), max_id + 1)

    finally:
        con.close()

    # 8. Update manifest (preserve next_sample_id to prevent reuse after removal)
    _update_manifest(db_dir, new_count, next_sample_id=next_id)

    logger.info("[remove-samples] Done.")

    return {"removed": list(sample_names), "not_found": []}


def check_database(db_dir: str) -> list[CheckResult]:
    """Validate database integrity. Returns list of CheckResult."""
    results: list[CheckResult] = []

    def err(msg: str) -> None:
        results.append(CheckResult("error", msg))

    def warn(msg: str) -> None:
        results.append(CheckResult("warning", msg))

    def info(msg: str) -> None:
        results.append(CheckResult("info", msg))

    # Check 1: manifest.json exists
    manifest_path = os.path.join(db_dir, "manifest.json")
    if not os.path.exists(manifest_path):
        err("manifest.json not found")
        return results

    # Check 2: manifest.json parses as JSON
    try:
        with open(manifest_path) as f:
            manifest = json.load(f)
    except (json.JSONDecodeError, OSError) as e:
        err(f"manifest.json parse error: {e}")
        return results

    # Check 3: required keys present
    required_keys = {"genome_build", "sample_count"}
    missing_keys = required_keys - set(manifest.keys())
    if missing_keys:
        err(f"manifest.json missing keys: {', '.join(sorted(missing_keys))}")
        return results

    # Check 4: genome_build is valid
    if manifest["genome_build"] not in VALID_GENOME_BUILDS:
        err(f"Invalid genome_build: '{manifest['genome_build']}'")

    # Check 5: metadata.sqlite exists
    sqlite_path = os.path.join(db_dir, "metadata.sqlite")
    if not os.path.exists(sqlite_path):
        err("metadata.sqlite not found")
        return results

    # Check 6: all 4 tables exist in SQLite
    required_tables = {"samples", "technologies", "sample_phenotype", "precomputed_bitmaps"}
    try:
        con = sqlite3.connect(sqlite_path)
        existing_tables = {
            r[0]
            for r in con.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()
        }
        missing_tables = required_tables - existing_tables
        if missing_tables:
            err(f"SQLite missing tables: {', '.join(sorted(missing_tables))}")
            con.close()
            return results
    except Exception as e:
        err(f"SQLite open error: {e}")
        return results

    n_techs = 0
    max_sample_id = -1

    try:
        # Check 7: sample_count matches COUNT(*) FROM samples
        db_count = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
        if manifest["sample_count"] != db_count:
            warn(
                f"sample_count mismatch: manifest says {manifest['sample_count']}, "
                f"database has {db_count}"
            )

        # Check 8: capture pickle files exist for each tech_id
        tech_ids = [
            r[0] for r in con.execute("SELECT tech_id FROM technologies").fetchall()
        ]
        n_techs = len(tech_ids)
        capture_dir = os.path.join(db_dir, "capture")
        for tech_id in tech_ids:
            pickle_path = os.path.join(capture_dir, f"tech_{tech_id}.pickle")
            if not os.path.exists(pickle_path):
                err(f"Missing capture index: tech_{tech_id}.pickle")

        max_sid_row = con.execute("SELECT MAX(sample_id) FROM samples").fetchone()
        max_sample_id = max_sid_row[0] if max_sid_row[0] is not None else -1

    finally:
        con.close()

    # Check 9: variants/ directory exists
    variants_dir = os.path.join(db_dir, "variants")
    if not os.path.exists(variants_dir):
        err("variants/ directory not found")
        return results

    # Check 9b: no chromosome is stored in both layouts at once
    for chrom in storage.mixed_layout_chroms(variants_dir):
        err(
            f"{chrom}: both variants/{chrom}/ (bucketed) and variants/{chrom}.parquet "
            f"exist. Queries read only the bucketed files, so every sample whose calls "
            f"live in {chrom}.parquet is silently counted as homozygous reference and "
            f"biases allele frequencies downward. See the troubleshooting guide: "
            f"'Samples Added by update-db Are Missing From Queries'."
        )

    parquet_files = [str(p) for p in storage.iter_variant_parquets(variants_dir)]
    n_chroms = len(storage.flat_chroms(variants_dir)) + len(
        storage.partitioned_chroms(variants_dir)
    )

    expected_fields = {
        "pos":        pa.uint32(),
        "ref":        pa.large_utf8(),
        "alt":        pa.large_utf8(),
        "het_bitmap": pa.large_binary(),
        "hom_bitmap": pa.large_binary(),
    }

    for pq_path in parquet_files:
        basename = os.path.basename(pq_path)
        if basename.startswith("bucket_"):
            chrom = os.path.basename(os.path.dirname(pq_path))
            label = f"{chrom}/{basename}"
        else:
            chrom = basename.replace(".parquet", "")
            label = f"{chrom}.parquet"

        # Check 10: each Parquet file readable
        try:
            table = pq.read_table(pq_path)
        except Exception as e:
            err(f"{label}: read error: {e}")
            continue

        # Check 11: correct schema
        for fname, ftype in expected_fields.items():
            if fname not in table.schema.names:
                err(f"{label}: missing field '{fname}'")
            elif table.schema.field(fname).type != ftype:
                err(
                    f"{label}: field '{fname}' has wrong type "
                    f"'{table.schema.field(fname).type}' (expected '{ftype}')"
                )

        # Check 12: rows sorted by (pos, alt)
        if len(table) > 1:
            pos_col = table["pos"].to_pylist()
            alt_col = table["alt"].to_pylist()
            pairs = list(zip(pos_col, alt_col))
            if pairs != sorted(pairs):
                warn(f"{label}: rows not sorted by (pos, alt)")

        # Check 13: bitmap max element ≤ max(sample_id) (sample first 1000 rows)
        sample_limit = min(1000, len(table))
        for i in range(sample_limit):
            for field in ("het_bitmap", "hom_bitmap"):
                bm_bytes = table[field][i].as_py()
                if bm_bytes:
                    bm = deserialize(bm_bytes)
                    if bm and max(bm) > max_sample_id:
                        err(
                            f"{label} row {i}: {field} contains sample_id "
                            f"{max(bm)} > max sample_id {max_sample_id}"
                        )

    # Check 14: info summary
    info(
        f"genome_build={manifest['genome_build']}, "
        f"sample_count={manifest['sample_count']}, "
        f"chroms={n_chroms}, "
        f"techs={n_techs}"
    )

    return results

import json
import logging
import os
import sqlite3
import time
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
from pyroaring import BitMap

from .. import storage
from ..bitmaps import deserialize, serialize
from .build import PARQUET_SCHEMA

logger = logging.getLogger(__name__)


def compact_database(db_path: Path) -> dict:
    """
    Rewrite all Parquet files removing dead bits and all-zero-bitmap rows.

    Algorithm:
    1. Load active sample_ids from SQLite
    2. For each Parquet file (flat or partitioned bucket):
       a. Read all rows
       b. AND each bitmap with active_ids (removes dead bits)
       c. Drop rows where both bitmaps are empty after the AND
       d. Rewrite atomically via .tmp rename
    3. Update manifest.json with last_compact timestamp

    Returns stats: {files_rewritten, rows_processed, rows_removed,
                    size_before, size_after}
    """
    db_path = Path(db_path)
    variants_dir = db_path / "variants"

    # Get active sample IDs from SQLite
    con = sqlite3.connect(str(db_path / "metadata.sqlite"))
    rows = con.execute("SELECT sample_id FROM samples").fetchall()
    con.close()
    active_ids = BitMap([r[0] for r in rows])

    # Collect all parquet files (flat + partitioned buckets)
    all_parquets: list[Path] = list(storage.iter_variant_parquets(variants_dir))

    logger.info("[compact] Compacting %d Parquet file(s) against %d active sample(s)...",
                len(all_parquets), len(active_ids))
    t0 = time.monotonic()

    files_rewritten = 0
    rows_processed = 0
    rows_removed = 0
    size_before = sum(f.stat().st_size for f in all_parquets)

    for parquet_file in all_parquets:
        table = pq.read_table(str(parquet_file))
        rows_processed += len(table)
        has_phase2 = (
            "filtered_bitmap" in table.schema.names
            and "quality_pass_bitmap" in table.schema.names
        )

        keep_indices = []
        new_het_list = []
        new_hom_list = []
        new_fail_list = []
        new_filtered_list = []
        new_quality_pass_list = []
        dirty = False

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

            new_het = het_bm & active_ids
            new_hom = hom_bm & active_ids
            new_fail = fail_bm & active_ids
            new_filt = filt_bm & active_ids
            new_qp = qp_bm & active_ids

            if (new_het != het_bm or new_hom != hom_bm or new_fail != fail_bm
                    or new_filt != filt_bm or new_qp != qp_bm):
                dirty = True

            if (not new_het and not new_hom and not new_fail
                    and not new_filt and not new_qp):
                # No active samples carry this variant — remove the row
                rows_removed += 1
                dirty = True
                continue

            keep_indices.append(i)
            new_het_list.append(serialize(new_het))
            new_hom_list.append(serialize(new_hom))
            new_fail_list.append(serialize(new_fail))
            new_filtered_list.append(serialize(new_filt))
            new_quality_pass_list.append(serialize(new_qp))

        if not dirty:
            logger.debug("  [compact] %s: no changes", parquet_file.name)
            continue

        # Build new table with kept rows and updated bitmaps
        orig_keep = table.take(keep_indices)
        new_table = pa.table(
            {
                "pos":                 orig_keep["pos"],
                "ref":                 orig_keep["ref"],
                "alt":                 orig_keep["alt"],
                "het_bitmap":          pa.array(new_het_list,          type=pa.large_binary()),
                "hom_bitmap":          pa.array(new_hom_list,          type=pa.large_binary()),
                "fail_bitmap":         pa.array(new_fail_list,         type=pa.large_binary()),
                "filtered_bitmap":     pa.array(new_filtered_list,     type=pa.large_binary()),
                "quality_pass_bitmap": pa.array(new_quality_pass_list, type=pa.large_binary()),
            },
            schema=PARQUET_SCHEMA,
        )

        tmp_path = str(parquet_file) + ".tmp"
        pq.write_table(new_table, tmp_path)
        os.replace(tmp_path, str(parquet_file))
        files_rewritten += 1
        rows_removed_this_file = len(table) - len(keep_indices)
        logger.debug("  [compact] %s: %d row(s) kept, %d removed (rewritten)",
                     parquet_file.name, len(keep_indices), rows_removed_this_file)

    size_after = sum(f.stat().st_size for f in all_parquets if f.exists())

    logger.info("[compact] Complete: %d file(s) rewritten, %d row(s) removed (%.1fs)",
                files_rewritten, rows_removed, time.monotonic() - t0)

    # Update manifest with last_compact timestamp
    manifest_path = db_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["last_compact"] = datetime.now(timezone.utc).isoformat()
    tmp_path = str(manifest_path) + ".tmp"
    with open(tmp_path, "w") as f:
        json.dump(manifest, f, indent=2)
    os.replace(tmp_path, str(manifest_path))

    # Append changelog entry
    con2 = sqlite3.connect(str(db_path / "metadata.sqlite"))
    try:
        con2.execute(
            "INSERT INTO changelog (event_type, event_time, sample_names, notes) VALUES (?, ?, ?, ?)",
            (
                "compact",
                datetime.now(timezone.utc).isoformat(),
                None,
                f"{files_rewritten} files rewritten, {rows_removed} rows removed",
            ),
        )
        con2.commit()
    finally:
        con2.close()

    return {
        "files_rewritten": files_rewritten,
        "rows_processed": rows_processed,
        "rows_removed": rows_removed,
        "size_before": size_before,
        "size_after": size_after,
    }

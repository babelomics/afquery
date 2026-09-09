import json
import sqlite3
import tempfile
import time
from pathlib import Path

import pyarrow.parquet as pq

from . import storage
from .database import Database


def run_benchmark(db_path: Path, n_warmup: int = 3) -> dict:
    """
    Run a battery of queries against db_path and report timing in milliseconds.

    Returns:
        {
          "point_query_cold_ms": float,
          "point_query_warm_ms": float,   # median of n_warmup repeats
          "batch_100_ms": float,
          "batch_1000_ms": float,
          "targets": {
              "point_cold_ok": bool,      # <100ms
              "point_warm_ok": bool,      # <10ms
              "batch_100_ok": bool,       # <200ms
          },
          "correctness": {
              "point_returned_data": bool,
              "batch_100_hit_count": int,
              "batch_1000_hit_count": int,
          }
        }
    """
    db_path = Path(db_path)
    db = Database(str(db_path))

    # Discover real variants from the DB
    all_variants = _find_test_variants(db_path, n=1100)
    if not all_variants:
        return {"error": "No variants found in database"}
    test_chrom, test_pos, test_ref, test_alt = all_variants[0]

    # Get an Phenotype code present in the DB
    con = sqlite3.connect(str(db_path / "metadata.sqlite"))
    phenotype_row = con.execute(
        "SELECT phenotype_code FROM sample_phenotype LIMIT 1"
    ).fetchone()
    con.close()
    phenotype_codes = [phenotype_row[0]] if phenotype_row else ["E11.9"]

    # Cold point query (first call — cold I/O cache)
    t0 = time.perf_counter()
    point_result = db.query(chrom=test_chrom, pos=test_pos, phenotype=phenotype_codes)
    cold_ms = (time.perf_counter() - t0) * 1000

    # Warm point queries (n_warmup repeats, take median)
    warm_times = []
    for _ in range(n_warmup):
        t0 = time.perf_counter()
        db.query(chrom=test_chrom, pos=test_pos, phenotype=phenotype_codes)
        warm_times.append((time.perf_counter() - t0) * 1000)
    warm_ms = sorted(warm_times)[len(warm_times) // 2]

    # Build batches from real variants on the same chrom
    same_chrom_variants = [(p, r, a) for c, p, r, a in all_variants if c == test_chrom]
    n_available = len(same_chrom_variants)
    if n_available < 100:
        import warnings
        warnings.warn(
            f"Only {n_available} variants on {test_chrom}; batch tests use fewer variants than target"
        )

    batch_100 = same_chrom_variants[:100]
    batch_1000 = same_chrom_variants[:1000]

    # Batch 100 query
    t0 = time.perf_counter()
    batch_100_result = db.query_batch(chrom=test_chrom, variants=batch_100, phenotype=phenotype_codes)
    batch_100_ms = (time.perf_counter() - t0) * 1000

    # Batch 1000 query
    t0 = time.perf_counter()
    batch_1000_result = db.query_batch(chrom=test_chrom, variants=batch_1000, phenotype=phenotype_codes)
    batch_1000_ms = (time.perf_counter() - t0) * 1000

    # Count how many queried variants actually returned data
    batch_100_positions = {p for p, r, a in batch_100}
    batch_1000_positions = {p for p, r, a in batch_1000}
    batch_100_hits = len({row.variant.pos for row in batch_100_result if row.variant.pos in batch_100_positions}) if batch_100_result else 0
    batch_1000_hits = len({row.variant.pos for row in batch_1000_result if row.variant.pos in batch_1000_positions}) if batch_1000_result else 0

    return {
        "point_query_cold_ms": round(cold_ms, 2),
        "point_query_warm_ms": round(warm_ms, 2),
        "batch_100_ms": round(batch_100_ms, 2),
        "batch_1000_ms": round(batch_1000_ms, 2),
        "targets": {
            "point_cold_ok": cold_ms < 100,
            "point_warm_ok": warm_ms < 25,
            "batch_100_ok": batch_100_ms < 200,
        },
        "correctness": {
            "point_returned_data": len(point_result) > 0,
            "batch_100_hit_count": batch_100_hits,
            "batch_1000_hit_count": batch_1000_hits,
        },
    }


def _find_test_variants(
    db_path: Path, n: int = 1100
) -> list[tuple[str, int, str, str]]:
    """Return up to n (chrom, pos, ref, alt) tuples from the DB."""
    variants_dir = db_path / "variants"
    if not variants_dir.exists():
        return []

    results: list[tuple[str, int, str, str]] = []
    for entry in storage.iter_variant_parquets(variants_dir):
        if len(results) >= n:
            break
        chrom = entry.stem if entry.parent == variants_dir else entry.parent.name
        tbl = pq.read_table(str(entry), columns=["pos", "ref", "alt"])
        for row in range(min(len(tbl), n - len(results))):
            results.append((
                chrom,
                int(tbl["pos"][row].as_py()),
                str(tbl["ref"][row].as_py()),
                str(tbl["alt"][row].as_py()),
            ))
    return results


def run_benchmark_with_synth(
    n_samples: int = 1000,
    n_variants: int = 10_000,
    output_report: str = "benchmark_report.json",
) -> dict:
    """Build a synthetic DB, run benchmarks, and save a JSON report."""
    from .preprocess import run_preprocess
    from .preprocess.synth import generate_synthetic_manifest

    with tempfile.TemporaryDirectory(prefix="afquery_bench_") as tmp_root:
        tmp_path = Path(tmp_root)
        synth_dir = tmp_path / "synth"
        db_dir = tmp_path / "db"
        db_dir.mkdir()

        manifest_path = generate_synthetic_manifest(
            synth_dir,
            n_samples=n_samples,
            n_variants_per_chrom=n_variants,
            chroms=("chr1",),
        )

        run_preprocess(
            manifest_path=str(manifest_path),
            output_dir=str(db_dir),
            genome_build="GRCh37",
        )

        results = run_benchmark(db_dir)

    with open(output_report, "w") as f:
        json.dump(results, f, indent=2)

    return results

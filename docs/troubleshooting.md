# Troubleshooting

---

## create-db Runs Out of Memory

**Symptom:** `create-db` fails with `MemoryError`, `Out of memory`, or a DuckDB `Not enough memory to store` error.

**Cause:** The per-worker DuckDB memory limit is too low for your cohort size or variant density.

**Fix:**

- Lower the number of parallel build workers: `--build-threads 4`
- Increase per-worker memory: `--build-memory 4GB`
- Both together keep total RAM the same with fewer concurrent workers

```bash
# Instead of default (all CPUs × 2GB)
afquery create-db ... --build-threads 8 --build-memory 4GB
```

See [Performance Tuning](advanced/performance.md) for sizing guidance.

---

## Query Returns AN=0

**Symptom:** A query returns `AC=0, AN=0, AF=None` for a variant you know exists.

Common causes include restrictive filters, WES positions outside capture regions, or missing chromosomes. For a step-by-step diagnostic checklist, see [Debugging Results → Unexpected AN=0](advanced/debugging-results.md#1-unexpected-an0).

---

## VCF Annotation Is Slow

**Symptom:** `afquery annotate` takes much longer than expected.

**Fixes:**

- Increase thread count: `--threads 16`
- Check disk I/O: annotation reads many small Parquet files; SSDs are significantly faster than spinning disks
- For very large VCFs (1M+ variants), annotation time scales linearly with variant count

See [Performance Tuning](advanced/performance.md) for general thread and memory sizing guidance.

---

## DuckDB Errors During Build

**Symptom:** Errors like `"Not supported: Writing to Arrow IPC"` or `"Unsupported file format"` during `create-db`.

**Cause:** An older DuckDB version is attempting to write Arrow IPC format for temporary files instead of Parquet.

**Fix:** Upgrade to DuckDB ≥ 0.10:

```bash
pip install --upgrade duckdb
```

AFQuery requires DuckDB to use Parquet for all temporary files. Arrow IPC is not supported.

---

## cyvcf2 ImportError in Workers

**Symptom:** `ImportError: cannot import name 'VCF' from 'cyvcf2'` in worker processes during ingest.

**Cause:** cyvcf2 cannot be pickled for `ProcessPoolExecutor` — it must be imported inside the worker function body.

**Fix:** This is an internal invariant of AFQuery. If you see this error with the latest version, please file a bug report. Do not import cyvcf2 at module level in any code that runs in worker processes.

---

## Wrong Bucket IDs (Silent Bug)

**Symptom:** Queries return no results for known variants, but the database was built without errors.

**Cause:** A DuckDB float-division rounding bug was present in older AFQuery versions. When computing bucket IDs, `CAST(pos / 1000000 AS BIGINT)` rounds floats incorrectly (e.g., position 1,500,000 → bucket 2 instead of 1).

**Fix:** Upgrade to the latest AFQuery version. The fix uses `CAST(pos AS BIGINT) // 1000000` — the integer-division operator — for correct bucket IDs. Rebuild the database after upgrading.

---

## Compact Takes a Long Time

**Symptom:** `afquery update-db --compact` runs for many minutes or hours.

**Cause:** Compaction rewrites every Parquet file in the database. For large databases (many chromosomes × many buckets), this is expected.

**Recommendation:** Run compact during off-hours or overnight. It is safe to interrupt (resume is not supported; re-run to complete).

---

## Sample Not Found in Remove Operation

**Symptom:** `afquery update-db --remove-samples SAMP_001` fails with `"Sample not found"`.

**Cause:** The sample name is case-sensitive and must match exactly what was in the original manifest.

**Fix:** Check the exact sample name:

```bash
afquery info --db ./db/ --samples | grep SAMP
```

---

## afquery check Reports Errors

**Symptom:** `afquery check --db ./db/` exits non-zero and prints error messages.

**Common errors:**

| Error | Fix |
|-------|-----|
| `Missing Parquet for chromosome chr3` | Re-run `create-db` or investigate incomplete build |
| `Manifest mismatch: expected N samples, found M` | Database may be partially updated; re-run `update-db` |
| `Capture file missing for wes_v1` | BED file was not provided at build time; rebuild with `--bed-dir` |
| `chr1: both variants/chr1/ ... and variants/chr1.parquet exist` | Samples added by a pre-0.4.1 `update-db` are invisible to queries; see [Samples Added by update-db Are Missing From Queries](#samples-added-by-update-db-are-missing-from-queries) |

---

## WES Technology Treated as WGS

**Symptom:** AN for WES samples is much higher than expected; positions outside the capture panel return results.

**Cause:** No BED file was found for the technology. When `<tech>.bed` is missing from `--bed-dir`, AFQuery treats the technology as WGS (all positions covered) with a warning.

**Fix:** Verify BED files are present:

```bash
ls ./beds/
# Should include: wes_v1.bed, wes_v2.bed, etc.
```

Then verify the database was built correctly:
```bash
afquery check --db ./db/
```

Look for warnings like `"Capture file missing for wes_v1"`. Rebuild with the BED files if necessary.

---

## Next Steps

- [FAQ](faq.md) — common questions and answers
- [Debugging Results](advanced/debugging-results.md) — diagnostic checklist for unexpected results
- [Performance Tuning](advanced/performance.md) — memory and thread configuration for build and query phases

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

## Samples Added by update-db Are Missing From Queries

**Symptom:** After `afquery update-db --add-samples`, the new samples show up in
`afquery info --db ./db/ --samples` and `AN` grows by the expected amount, but they
never appear as carriers: `AC` does not increase at variants you know they carry, and
`afquery variant-info` does not list them. Allele frequencies across the whole database
drift downward after every update.

`afquery check --db ./db/` reports one error per affected chromosome:

```
ERROR  chr1: both variants/chr1/ (bucketed) and variants/chr1.parquet exist.
```

**Cause:** `update-db --add-samples` before 0.4.1 only understood the
single-file-per-chromosome variant layout. `create-db` produces the bucketed layout
(`variants/<chrom>/bucket_N.parquet`), so the merge found nothing to merge and wrote the
new samples to a fresh `variants/<chrom>.parquet` instead. Queries read the bucketed
files and ignore that one, so the added samples counted toward `AN` through the capture
index but never as carriers — they were treated as homozygous reference everywhere.

**Fix:** Upgrade to 0.4.1 or later, then repair the database. `afquery info --db ./db/
--changelog` lists the samples added by each `add_samples` event.

```bash
# 0. Stop all writers and back up the small files.
cp ./db/manifest.json ./db/manifest.json.bak
cp ./db/metadata.sqlite ./db/metadata.sqlite.bak

# 1. Record the affected samples BEFORE removing them: removal deletes their
#    phenotype rows, and you need them to rebuild the manifest.
sqlite3 ./db/metadata.sqlite \
  "SELECT s.sample_name, s.sex, t.tech_name, s.vcf_path,
          group_concat(p.phenotype_code)
   FROM samples s
   JOIN technologies t ON s.tech_id = t.tech_id
   LEFT JOIN sample_phenotype p ON p.sample_id = s.sample_id
   WHERE s.sample_name IN ('SAMPLE_1','SAMPLE_2') GROUP BY s.sample_id;"

# 2. List the orphan files before deleting anything.
find ./db/variants -maxdepth 1 -name '*.parquet'

# 3. Remove the affected samples. This clears their bits from both layouts and is
#    safe on a split database. Every bucket is read, so budget minutes, not
#    seconds, and do not interrupt it.
afquery update-db --db ./db/ --remove-samples SAMPLE_1 --remove-samples SAMPLE_2

# 4. Delete the orphan flat files. Each one has a sibling directory of the same name.
find ./db/variants -maxdepth 1 -name '*.parquet' -delete

# 5. Confirm the errors are gone.
afquery check --db ./db/

# 6. Re-add the samples with the fixed version.
afquery update-db --db ./db/ --add-samples repair.tsv --bed-dir ./beds/

# 7. Verify a variant you know they carry.
afquery variant-info --db ./db/ --locus chr1:887801
```

Re-added samples receive new sample IDs — IDs are never reused after a removal. That is
expected and does not affect results.

If the original VCFs are no longer available, stop after step 5. The samples are then
absent from both the metadata and `AN`, which is a correct smaller cohort rather than a
biased larger one. Rebuilding with `create-db` is always a valid fallback.

Databases that were only ever built with `create-db`, and never updated, are unaffected.

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

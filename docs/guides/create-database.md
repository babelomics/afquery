# Create a Database

`afquery create-db` builds the AFQuery database from a manifest of single-sample VCFs. This is a one-time setup step; incremental updates use `afquery update-db`.

---

## Basic Usage

```bash
afquery create-db \
  --manifest manifest.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38
```

For cohorts with WES/panel samples, provide the BED file directory:

```bash
afquery create-db \
  --manifest manifest.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38 \
  --bed-dir ./beds/
```

---

## What Happens

1. **Ingest phase** — Each VCF is parsed with cyvcf2. Genotypes and quality fields are written to a temporary per-sample Parquet file, one row per variant per sample.
2. **Build phase** — DuckDB reads the temporary Parquet files, aggregates per 1-Mbp bucket, and writes Roaring Bitmap Parquet files partitioned by chromosome and bucket.
3. **Finalize** — `manifest.json` and `metadata.sqlite` are written to the output directory.

---

## Directory Layout After Creation

```
./db/
├── manifest.json          # Build configuration (genome build, schema version, etc.)
├── metadata.sqlite        # Sample/phenotype/technology/changelog metadata
├── variants/              # Parquet files partitioned by chromosome and bucket
│   ├── chr1/
│   │   ├── bucket_0.parquet      # Positions 0–999,999
│   │   ├── bucket_1.parquet      # Positions 1,000,000–1,999,999
│   │   └── ...
│   ├── chr2/
│   └── ...
└── capture/               # Interval trees for WES technologies (pickle files)
    ├── wes_v1.pkl
    └── wes_v2.pkl
```

---

## Memory and Thread Tuning

For large cohorts, tune these options:

```bash
afquery create-db \
  --manifest manifest.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38 \
  --build-threads 32 \
  --build-memory 4GB
```

| Option | Default | Recommendation |
|--------|---------|----------------|
| `--build-threads` | all CPUs | Set to `min(cpu_count, available_RAM_GB / 2)` |
| `--build-memory` | `2GB` | Increase for dense WGS regions or large cohorts |
| `--threads` | all CPUs | Controls ingest parallelism (VCF parsing) |

The build phase uses one DuckDB process per 1-Mbp bucket. With `--build-threads 32` and `--build-memory 4GB`, peak RAM usage is approximately `32 × 4 = 128 GB`.

---

## Resume Behavior

If `create-db` is interrupted, it resumes automatically from where it left off. Individual bucket Parquet files that were already written are skipped.

To force a complete restart:

```bash
afquery create-db --manifest manifest.tsv --output-dir ./db/ --genome-build GRCh38 --force
```

!!! warning
    `--force` deletes all existing output in `--output-dir`. Use with caution.

---

## FILTER=PASS Behavior

Only `FILTER=PASS` calls (or calls with no FILTER field) contribute to AC. Calls that fail a filter are stored in `fail_bitmap`; at query time those samples are left out of AC and AN and reported as `N_FAIL`.

See [FILTER=PASS Tracking](../advanced/filter-pass-tracking.md) for details.

---

## Coverage-Evidence Filters

Four optional flags enable per-sample, quality-aware tracking of which positions
each partially-covered technology (WES, panels) actually covered. They are
fully opt-in.

| Flag | Default | Effect |
|------|---------|--------|
| `--min-dp D`     | 0   | Minimum `FORMAT/DP` for a carrier to count as quality evidence. |
| `--min-gq G`     | 0   | Minimum `FORMAT/GQ` for a carrier to count as quality evidence. |
| `--min-qual Q`   | 0.0 | Minimum VCF `QUAL` field for a carrier to count as quality evidence. |
| `--min-covered K`| 0   | Per partially-covered tech, the position is "trusted" only if at least K of its carriers pass the quality thresholds. Non-carriers of failing positions are recorded as `N_NO_COVERAGE`. |

When any of these flags is non-zero AFQuery reads `FORMAT/DP`, `FORMAT/GQ`,
and `QUAL` from each variant call during ingest. Use the bundled
`resources/normalize_vcf.sh` (which preserves these FORMAT fields) or ensure
your own preprocessing keeps them.

Example:

```bash
afquery create-db \
  --manifest samples.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38 \
  --bed-dir ./beds/ \
  --min-dp 30 --min-gq 20 --min-covered 1
```

Thresholds are fixed at creation time. `update-db --add-samples` reuses them
and re-applies them to every position whose partially-covered tech receives
new samples (see [Update Database](update-database.md)).

See [Coverage Evidence](../advanced/coverage-evidence.md) for when to reach
for each flag, how `N_NO_COVERAGE` is computed, and the query-time companion
flag `--min-quality-evidence`.

---

## Validating the Result

After creation, run:

```bash
afquery check --db ./db/
```

And inspect database metadata:

```bash
afquery info --db ./db/
```

Example `info` output:
```
Database:       ./db/
Schema version: 2.0
Genome build:   GRCh38
DB version:     1.0
Samples:        1371
Technologies:   wgs, wes_v1, wes_v2
Chromosomes:    chr1 ... chrX chrY chrM
```

These commands are available at any time after database creation — not only immediately after `create-db`. Use `afquery check` to verify database integrity (manifest consistency, Parquet file health, capture index presence) and `afquery info` to inspect sample counts, registered technologies, and phenotype codes.

---

## Full Option Reference

See [CLI Reference → create-db](../reference/cli.md#create-db).

---

## Next Steps

- [Manifest Format](manifest-format.md) — TSV manifest column reference and common mistakes
- [Query Allele Frequencies](query.md) — run your first queries against the new database
- [Performance Tuning](../advanced/performance.md) — build thread and memory configuration for large cohorts

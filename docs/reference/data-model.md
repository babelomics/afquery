# Data Model

This page documents the on-disk layout of an AFQuery database, including file formats, schemas, and schema version differences.

---

## Directory Layout

```
<db_dir>/
├── manifest.json          # Build configuration
├── metadata.sqlite        # Sample/phenotype/technology/changelog metadata
├── variants/              # Parquet variant data, partitioned by chromosome and bucket
│   ├── chr1/
│   │   ├── bucket_0.parquet      # Positions 0–999,999
│   │   ├── bucket_1.parquet      # Positions 1,000,000–1,999,999
│   │   └── ...
│   ├── chr2/
│   └── ...
└── capture/               # WES coverage interval trees
    ├── wes_v1.pkl
    └── wes_v2.pkl
```

### Storage layouts

`create-db` always writes the bucketed layout shown above. A single-file-per-chromosome
layout (`variants/chr1.parquet`, with no bucket directory) is also readable, and appears
in small hand-built databases and in test fixtures. `update-db` merges into whichever
layout a chromosome already uses, and creates buckets for a chromosome new to the
database.

A chromosome must never have both. `afquery check` reports that as an error, because
queries resolve the bucket directory first and would silently ignore the flat file — and
every sample stored only in it.

---

## manifest.json

Stores the build configuration used during `create-db`.

| Field | Type | Description |
|-------|------|-------------|
| `genome_build` | string | `"GRCh37"` or `"GRCh38"` |
| `version` | string | manifest.json format version |
| `db_version` | string | User-specified version label |
| `sample_count` | int | Number of samples at build time |
| `schema_version` | string | `"2.0"`, or `"3.0"` when the database was built with coverage-quality filters |
| `pass_only_filter` | bool | Always `true` — PASS-only counting is always enforced |
| `coverage_filter` | object | The `min_dp` / `min_gq` / `min_qual` / `min_covered` thresholds used at build time (all `0` when unused) |
| `created_at` | string | ISO 8601 timestamp |

---

## metadata.sqlite

SQLite database containing all mutable metadata. Tables:

### `samples`

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | INTEGER PRIMARY KEY | 0-indexed, monotonically increasing |
| `sample_name` | TEXT UNIQUE | Unique sample identifier |
| `sex` | TEXT | `male` or `female` |
| `tech_id` | INTEGER | Foreign key → `technologies.tech_id` |
| `active` | INTEGER | 1 = active, 0 = removed |

### `technologies`

| Column | Type | Description |
|--------|------|-------------|
| `tech_id` | INTEGER PRIMARY KEY | Auto-increment |
| `tech_name` | TEXT UNIQUE | Technology name |
| `bed_path` | TEXT \| NULL | Path to BED file; NULL = WGS |

### `sample_phenotype`

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | INTEGER | Foreign key → `samples.sample_id` |
| `phenotype_code` | TEXT | arbitrary string (e.g., ICD-10 code, HPO term, project tag) |

Many-to-many: one sample can have multiple phenotype codes.

### `precomputed_bitmaps`

| Column | Type | Description |
|--------|------|-------------|
| `key` | TEXT PRIMARY KEY | Bitmap identifier (e.g., `"all"`, `"male"`, `"female"`, `"tech:wgs"`) |
| `bitmap` | BLOB | Serialized Roaring Bitmap |

Cached bitmaps for common filter combinations, rebuilt on `update-db`.

### `changelog`

| Column | Type | Description |
|--------|------|-------------|
| `event_id` | INTEGER PRIMARY KEY | Auto-increment |
| `event_type` | TEXT | `preprocess`, `add_samples`, `remove_samples`, `compact` |
| `event_time` | TEXT | ISO 8601 timestamp |
| `sample_names` | TEXT \| NULL | JSON array of affected sample names, or NULL |
| `notes` | TEXT \| NULL | Human-readable summary |

---

## Parquet Schema

Each bucket Parquet file has this schema:

| Column | Arrow type | Description |
|--------|-----------|-------------|
| `pos` | `uint32` | 1-based genomic position |
| `ref` | `large_utf8` | Reference allele |
| `alt` | `large_utf8` | Alternate allele |
| `het_bitmap` | `large_binary` | Roaring Bitmap of heterozygous, `FILTER=PASS` sample IDs |
| `hom_bitmap` | `large_binary` | Roaring Bitmap of homozygous-alt, `FILTER=PASS` sample IDs |
| `fail_bitmap` | `large_binary` | Roaring Bitmap of sample IDs whose call has `FILTER≠PASS` at this site |
| `filtered_bitmap` | `large_binary` | Roaring Bitmap of WES non-carriers with uncertain coverage |
| `quality_pass_bitmap` | `large_binary` | Roaring Bitmap of carriers meeting the DP/GQ/QUAL thresholds |

Rows are sorted by `(pos, alt)` within each bucket.

`filtered_bitmap` and `quality_pass_bitmap` back the [coverage-evidence](../advanced/coverage-evidence.md) feature. The columns are always present, but they are empty unless the database was built with coverage-quality filters (schema version `3.0`).

!!! important "large_utf8 / large_binary"
    AFQuery uses `large_utf8` and `large_binary` (64-bit offsets) rather than `utf8` / `binary` (32-bit). This is required for compatibility with DuckDB's Parquet reader on large chromosomes.

---

## Bitmap Format

Bitmaps use the [Roaring Bitmap](https://roaringbitmap.org/) format, serialized by pyroaring's portable serialization.

- Bit position = sample ID (0-indexed integer)
- `het_bitmap`: bit set iff the sample is heterozygous with `FILTER=PASS`
- `hom_bitmap`: bit set iff the sample is homozygous alt with `FILTER=PASS`
- `fail_bitmap`: bit set iff the sample's call has `FILTER≠PASS` — either a carrier whose call failed a filter, or a missing genotype (`./.`) at a failed site
- `filtered_bitmap`, `quality_pass_bitmap`: see the Parquet schema above; both are empty unless the database was built with coverage-quality filters

To deserialize in Python:
```python
from pyroaring import BitMap

with open("bucket.parquet", "rb") as f:
    ...  # use pyarrow or duckdb to read the column

bm = BitMap.deserialize(bitmap_bytes)
sample_ids = list(bm)
```

---

## Capture Index (WES)

One pickle file per WES technology in `capture/<tech>.pkl`. Contains a `pyranges` interval tree (or dict of interval trees keyed by chromosome) for fast coverage lookup.

A position is covered by a WES technology if it falls within any interval in that technology's BED file (0-based, half-open).

---

## Partitioning

Variants are partitioned into 1-Mbp buckets:

```
bucket_id = pos // 1_000_000
```

`update-db --add-samples` writes new positions into the bucket that owns
them, creating `bucket_N.parquet` when a batch extends a chromosome past its
previous last bucket.

!!! warning "DuckDB integer arithmetic"
    When computing bucket IDs in DuckDB SQL, always use the integer-division operator:
    ```sql
    CAST(pos AS BIGINT) // 1000000
    ```
    Not `CAST(pos / 1000000 AS BIGINT)` — the `/` operator does float division, and casting the rounded result afterwards produces wrong bucket IDs.

---

See [FILTER=PASS Tracking](../advanced/filter-pass-tracking.md).

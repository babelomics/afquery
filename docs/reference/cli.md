# CLI Reference

All AFQuery commands follow the pattern `afquery <command> [OPTIONS]`.

---

## Global options

| Option | Description |
|--------|-------------|
| `--version` | Print the installed AFQuery version and exit |
| `--help` | Show help for the command and exit |

`afquery --version` prints e.g. `afquery 0.4.2`. The version is derived from the
release tag at build time.

---

## create-db

Build a new AFQuery database from a manifest of single-sample VCFs.

```
afquery create-db [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--manifest` | TEXT | **required** | Path to TSV manifest file |
| `--output-dir` | TEXT | **required** | Path to output database directory |
| `--genome-build` | `GRCh37`\|`GRCh38` | **required** | Reference genome build |
| `--threads` | INTEGER | all CPUs | Worker threads for the ingest phase (VCF parsing) |
| `--build-threads` | INTEGER | min(cpu_count, n_buckets) | Max parallel workers for the build phase (DuckDB) |
| `--build-memory` | TEXT | `2GB` | DuckDB memory limit per build worker |
| `--tmp-dir` | TEXT | `{output_dir}/.tmp_preprocess` | Temporary directory for intermediate files |
| `--bed-dir` | TEXT | None | Directory containing BED files for WES technologies |
| `--force` | flag | False | Delete any partial results and restart from scratch |
| `--db-version` | TEXT | `1.0` | Version label for this database |
| `--min-dp` | INTEGER | `0` | Minimum `FORMAT/DP` for a carrier to count as quality evidence (0 = disabled) |
| `--min-gq` | INTEGER | `0` | Minimum `FORMAT/GQ` for a carrier to count as quality evidence (0 = disabled) |
| `--min-qual` | FLOAT | `0.0` | Minimum `QUAL` for a carrier to count as quality evidence (0 = disabled) |
| `--min-covered` | INTEGER | `0` | Minimum quality-passing carriers per partially-covered tech for hom-ref to be assumed (0 = disabled) |
| `-v, --verbose` | flag | False | Verbose output with per-item progress |

---

## query

Query allele frequencies at one or more positions.

```
afquery query [OPTIONS]
```

Exactly one of `--locus`, `--region`, or `--from-file` must be provided.

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |
| `--locus` | TEXT | None | Single position as `CHROM:POS` (e.g., `chr1:925952`) |
| `--region` | TEXT | None | Genomic range as `CHROM:START-END` (e.g., `chr1:900000-1000000`) |
| `--from-file` | PATH | None | Headerless TSV with columns `chrom pos [ref [alt]]` (batch query; multi-chromosome supported) |
| `--phenotype` | TEXT | None | Phenotype filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--sex` | `male`\|`female`\|`both` | `both` | Restrict to specified sex |
| `--tech` | TEXT | None | Technology filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--ref` | TEXT | None | Filter to specific reference allele (only for `--locus`) |
| `--alt` | TEXT | None | Filter to specific alternate allele (only for `--locus`) |
| `--format` | `text`\|`json`\|`tsv` | `text` | Output format |
| `--min-pass` | INTEGER | `0` | Min PASS carriers (`het\|hom`) per partially-covered tech for hom-ref to be assumed. Non-carriers move to `N_NO_COVERAGE` if a tech falls below the threshold. (0 = disabled) |
| `--min-observed` | INTEGER | `0` | Min any-VCF entries (`het\|hom\|fail`) per partially-covered tech for hom-ref to be assumed. (0 = disabled) |
| `--min-quality-evidence` | INTEGER | `0` | Min quality-passing carriers per partially-covered tech. Requires a database built with `--min-dp`, `--min-gq`, `--min-qual`, or `--min-covered`. (0 = disabled) |
| `--no-warn` | flag | False | Suppress warnings for unknown phenotypes, technologies, and chromosomes |

---

## annotate

Annotate a VCF file with allele frequency information.

```
afquery annotate [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |
| `--input` | TEXT | **required** | Input VCF file (plain or `.gz`) |
| `--output` | TEXT | **required** | Output annotated VCF file |
| `--phenotype` | TEXT | None | Phenotype filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--sex` | `male`\|`female`\|`both` | `both` | Restrict to specified sex |
| `--tech` | TEXT | None | Technology filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--threads` | INTEGER | all CPUs | Number of worker threads for parallel annotation |
| `-v, --verbose` | flag | False | Verbose output with per-item progress |
| `--min-pass` | INTEGER | `0` | Min PASS carriers per partially-covered tech for hom-ref to be assumed (0 = disabled) |
| `--min-observed` | INTEGER | `0` | Min any-VCF entries per partially-covered tech (0 = disabled) |
| `--min-quality-evidence` | INTEGER | `0` | Min quality-passing carriers per partially-covered tech. Requires a database built with `--min-dp`, `--min-gq`, `--min-qual`, or `--min-covered`. (0 = disabled) |
| `--no-warn` | flag | False | Suppress warnings for unknown phenotypes, technologies, and chromosomes |

The annotated VCF gains `AFQUERY_AC`, `AFQUERY_AN`, `AFQUERY_AF`,
`AFQUERY_N_HET`, `AFQUERY_N_HOM_ALT`, `AFQUERY_N_HOM_REF`, `AFQUERY_N_FAIL`,
and `AFQUERY_N_NO_COVERAGE` INFO fields.

---

## dump

Export allele frequency data to CSV.

```
afquery dump [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |
| `-o, --output` | TEXT | stdout | Output CSV file path |
| `--chrom` | TEXT | None | Restrict export to this chromosome |
| `--start` | INTEGER | None | 1-based start position (inclusive). Requires `--chrom`. |
| `--end` | INTEGER | None | 1-based end position (inclusive). Requires `--chrom`. |
| `--phenotype` | TEXT | None | Phenotype filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--sex` | `male`\|`female`\|`both` | `both` | Restrict to specified sex |
| `--tech` | TEXT | None | Technology filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--by-sex` | flag | False | Disaggregate output by sex (adds `AC_male`/`AC_female` columns) |
| `--by-tech` | flag | False | Disaggregate output by technology (adds `AC_<tech>` columns) |
| `--by-phenotype` | TEXT | None | Disaggregate by specific phenotype codes. Repeatable. |
| `--all-groups` | flag | False | Disaggregate by all sexes × technologies × phenotypes (Cartesian product) |
| `--threads` | INTEGER | all CPUs | Number of worker threads for parallel export |
| `--all-variants` | flag | False | Include variants with AC=0 (covered but not observed). WARNING: may produce very large output. |
| `-v, --verbose` | flag | False | Verbose output with per-item progress |
| `--min-pass` | INTEGER | `0` | Min PASS carriers per partially-covered tech for hom-ref to be assumed (0 = disabled) |
| `--min-observed` | INTEGER | `0` | Min any-VCF entries per partially-covered tech (0 = disabled) |
| `--min-quality-evidence` | INTEGER | `0` | Min quality-passing carriers per partially-covered tech. Requires a database built with `--min-dp`, `--min-gq`, `--min-qual`, or `--min-covered`. (0 = disabled) |

CSV output adds an `N_NO_COVERAGE` column (and per-group variants
`N_NO_COVERAGE_<label>` when disaggregating).

---

## variant-info

List samples carrying a specific variant, with their metadata.

```
afquery variant-info [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |
| `--locus` | TEXT | **required** | Single position as `CHROM:POS` (e.g., `chr1:925952`) |
| `--ref` | TEXT | None | Filter to specific reference allele |
| `--alt` | TEXT | None | Filter to specific alternate allele |
| `--phenotype` | TEXT | None | Phenotype filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--sex` | `male`\|`female`\|`both` | `both` | Restrict to specified sex |
| `--tech` | TEXT | None | Technology filter. Repeatable; comma-separated or multiple flags. Use `^` prefix to exclude. |
| `--format` | `text`\|`json`\|`tsv` | `text` | Output format |
| `--min-pass` | INTEGER | `0` | Min PASS carriers per partially-covered tech (samples below threshold appear with `genotype=no_coverage`). (0 = disabled) |
| `--min-observed` | INTEGER | `0` | Min any-VCF entries per partially-covered tech. (0 = disabled) |
| `--min-quality-evidence` | INTEGER | `0` | Min quality-passing carriers per partially-covered tech. Requires a database built with `--min-dp`, `--min-gq`, `--min-qual`, or `--min-covered`. (0 = disabled) |
| `--no-warn` | flag | False | Suppress `AfqueryWarning` messages |

Returns one row per carrier sample with genotype (`het`/`hom`/`alt`/`no_coverage`) and FILTER status (`PASS`/`FAIL`). For `no_coverage` rows the FILTER column is empty in TSV, `null` in JSON, and `-` in text — these samples have no call at this position, so `PASS`/`FAIL` does not apply. Use `--ref`/`--alt` to restrict to a specific allele when multiple alleles share the same position. Samples reported as `no_coverage` are non-carriers on a partially-covered tech that has been excluded from the hom-ref assumption by one of the coverage filters.

---

## update-db

Add samples, remove samples, update sample metadata, or compact the database.

```
afquery update-db [OPTIONS]
```

At least one of `--remove-samples`, `--add-samples`, `--compact`, `--update-sample`, or `--update-samples-file` must be provided. Operations execute in order: remove → update-metadata → add → compact.

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |
| `--remove-samples` | TEXT | None | Sample name(s) to remove. Repeatable; comma-separated or multiple flags. |
| `--add-samples` | PATH | None | Manifest TSV of new samples to add. Repeatable for multiple manifests. |
| `--compact` | flag | False | Remove dead bits from removed samples to reclaim disk space |
| `--update-sample` | TEXT | None | Sample name to update (single-sample metadata mode). Requires `--set-sex` and/or `--set-phenotype`. |
| `--set-sex` | TEXT | None | New sex for `--update-sample`. Options: `male`, `female`. |
| `--set-phenotype` | TEXT | None | New phenotype codes (comma-separated) for `--update-sample`. Replaces all current codes. |
| `--update-samples-file` | PATH | None | TSV file for batch metadata update. Header: `sample_name`, `field`, `new_value`. Mutually exclusive with `--update-sample`. |
| `--operator-note` | TEXT | None | Free-text note appended to each changelog entry for this metadata update. |
| `--threads` | INTEGER | all CPUs | Number of worker threads for parallel processing |
| `--tmp-dir` | TEXT | system temp | Temporary directory for intermediate files |
| `--bed-dir` | TEXT | None | Directory containing BED files for WES technologies |
| `--db-version` | TEXT | auto-increment | New version label after update |
| `-v, --verbose` | flag | False | Verbose output with per-item progress |

---

## info

Display database metadata and statistics.

```
afquery info [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |
| `--samples` | flag | False | List all samples with metadata |
| `--changelog` | flag | False | Show full changelog history |
| `--format` | `table`\|`tsv`\|`json` | `table` | Output format |

---

## check

Validate database integrity.

```
afquery check [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |

Exits with code 0 if the database is healthy, non-zero otherwise.

---

## version show

Display the database version label.

```
afquery version show [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--db` | TEXT | **required** | Path to database directory |

---

## version set

Set a new version label on the database.

```
afquery version set [OPTIONS] NEW_VERSION
```

| Argument/Option | Type | Description |
|-----------------|------|-------------|
| `NEW_VERSION` | TEXT (positional) | New version label |
| `--db` | TEXT | **required** — Path to database directory |

---

## benchmark

Run performance benchmarks on synthetic or real data.

```
afquery benchmark [OPTIONS]
```

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--n-samples` | INTEGER | `1000` | Number of synthetic samples to generate |
| `--n-variants` | INTEGER | `10000` | Number of variants per chromosome |
| `--output` | TEXT | `benchmark_report.json` | Output path for JSON benchmark report |
| `--db` | TEXT | None | Use an existing database instead of generating synthetic data |

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (invalid arguments, file not found, database integrity failure, etc.) |
| 2 | Usage error (Click framework: missing required argument or unknown option) |

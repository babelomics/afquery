# AFQuery

[![CI](https://github.com/babelomics/afquery/actions/workflows/ci.yml/badge.svg)](https://github.com/babelomics/afquery/actions/workflows/ci.yml)
[![Coverage](https://codecov.io/gh/babelomics/afquery/graph/badge.svg)](https://codecov.io/gh/babelomics/afquery)
[![Docs](https://img.shields.io/website?url=https%3A%2F%2Fbabelomics.github.io%2Fafquery%2F&label=docs)](https://babelomics.github.io/afquery/)
<br>
[![PyPI](https://img.shields.io/pypi/v/afquery.svg?color=blue)](https://pypi.org/project/afquery/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/afquery.svg)](https://anaconda.org/bioconda/afquery)
[![Docker](https://img.shields.io/badge/ghcr.io-afquery-blue?logo=docker)](https://github.com/babelomics/afquery/pkgs/container/afquery)
[![Python](https://img.shields.io/pypi/pyversions/afquery.svg)](https://pypi.org/project/afquery/)
[![License: MIT](https://img.shields.io/github/license/babelomics/afquery)](https://github.com/babelomics/afquery/blob/master/LICENSE)

**Fast, capture-aware allele frequency queries on local genomic cohorts**

AFQuery is a bitmap-indexed engine that recomputes AC/AN/AF for arbitrary subcohorts (by phenotype, sex, sequencing technology, or any combination) in tens of milliseconds, independently of cohort size. It accounts for capture-kit heterogeneity, ploidy on sex chromosomes, and FILTER/coverage evidence, and runs locally as a file-based system (Parquet + SQLite) with no server or cluster required.

[Full documentation →](https://babelomics.github.io/afquery/)

---

## Headline results

- **~14 ms point queries**, constant from 1,000 to 50,000 samples (O(1) scaling).
- **33× faster than bcftools** for full-chromosome bulk export at 2,504 samples; **R² > 0.99999** AF concordance over 1.1M common variants.
- **~25,000 variants/s** VCF annotation on 4 cores — a typical exome annotates in **~1 second**.
- **Up to 45× reduction** in toward-pathogenic ACMG classification errors vs. naive AN on mixed-capture-kit cohorts.
- **9.5×–13.4× storage compression** vs. input single-sample VCFs.

## Why AFQuery

Local allele frequencies are a core input to ACMG/AMP variant classification (criteria BA1, BS1, PM2), and global resources like gnomAD systematically underrepresent local ancestries and disease-enriched institutional cohorts. Computing AF on your own cohort sounds simple — until the cohort mixes WGS, WES, and panels, or several versions of the same capture kit. Successive Agilent SureSelect versions (v5/v6/v7), for example, only share ~57% of their targets on chromosome 22. Naive AN counting then *inflates* AN at positions outside some kits' targets, *deflates* AF, and systematically shifts variants toward "pathogenic" by ACMG criteria.

AFQuery solves this with per-position, per-technology capture-aware AN, ploidy-aware sex chromosome handling, and an explicit `N_NO_COVERAGE` channel that separates trusted hom-ref from "we cannot tell". Queries on subcohorts are answered by intersecting precomputed Roaring Bitmaps, so latency is constant in cohort size.

## When to use AFQuery

- You need allele frequencies for phenotype-defined or arbitrarily filtered subcohorts.
- Your cohort mixes sequencing technologies (WGS, WES, panels, multiple capture-kit versions).
- You want fast, repeated, interactive queries instead of one-off VCF re-scans.
- You need a local, reproducible workflow — no cloud, no Spark cluster.

## Features

- **Constant-time subcohort queries** — bitmap intersections at query time; no per-query VCF re-scan.
- **Capture-aware AN** — per-position eligibility from each technology's BED, eliminating systematic AF bias when mixing WGS / WES / panels and kit versions.
- **Ploidy-aware sex chromosomes** — correct AN on chrX PAR / non-PAR, chrY, chrM, by sample sex.
- **Coverage evidence model** — `N_NO_COVERAGE` separates trusted hom-ref from samples lacking sufficient evidence; query-time gates (`--min-pass`, `--min-observed`, `--min-quality-evidence`) keep AF conservative.
- **ACMG-compatible AC/AN/AF** — per-standard definitions, exposed as both query output and VCF INFO fields.
- **Flexible metadata filtering** — arbitrary phenotype labels (ICD-10, HPO, OMIM, custom tags), inclusion or exclusion (`^` prefix), combined with sex and technology.
- **Parallel VCF annotation** — multi-threaded; adds `AFQUERY_AC/AN/AF/N_HET/N_HOM_ALT/N_HOM_REF/N_FAIL/N_NO_COVERAGE` INFO fields.
- **Bulk CSV export** — per-variant frequencies with optional disaggregation by sex, technology, or phenotype.
- **Incremental updates** — add or remove samples, edit phenotype/sex metadata, compact storage, without full rebuilds.
- **Audit changelog** — every database operation is logged with timestamps and operator notes.
- **Database validation** — `afquery check` with scripted exit codes.
- **Serverless** — Parquet + SQLite on disk; no daemon, no Java, no Spark.

## Installation

```bash
# PyPI
pip install afquery

# Bioconda
conda install -c bioconda afquery

# Docker (linux/amd64, linux/arm64)
docker pull ghcr.io/babelomics/afquery:latest

# From source
git clone https://github.com/babelomics/afquery.git
cd afquery
pip install -e .
```

Requires Python ≥ 3.10. Core dependencies: `pyroaring`, `pyarrow`, `duckdb`, `pyranges`, `cyvcf2`, `click`, `tqdm`.

## Quickstart

### 1. Prepare a manifest

One row per sample (TSV, header required):

```tsv
sample_name	vcf_path	sex	tech_name	phenotype_codes
SAMP_001	/data/vcfs/SAMP_001.vcf.gz	female	wgs	E11.9,I10
SAMP_002	/data/vcfs/SAMP_002.vcf.gz	male	wes_v6	E11.9
SAMP_003	/data/vcfs/SAMP_003.vcf.gz	female	panel_card	I42.0
```

For every non-WGS technology, place a 3-column BED in `--bed-dir` named `<tech_name>.bed`.

### 2. Build the database

```bash
afquery create-db \
  --manifest manifest.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38 \
  --bed-dir ./beds/
```

### 3. Query

```bash
# Single locus
afquery query --db ./db/ --locus chr1:925952

# Locus filtered to a phenotype-and-sex subcohort
afquery query --db ./db/ --locus chr1:925952 --phenotype E11.9 --sex female

# Genomic region
afquery query --db ./db/ --region chr1:900000-1000000

# Batch from file (chrom pos [ref [alt]] per line)
afquery query --db ./db/ --from-file variants.tsv
```

### 4. Inspect carriers of a variant

```bash
afquery variant-info --db ./db/ --locus chr17:43093454
```

### 5. Annotate a VCF

```bash
afquery annotate \
  --db ./db/ \
  --input patient.vcf \
  --output patient.annotated.vcf \
  --threads 4
```

### 6. Export

```bash
# Region export, disaggregated by sex
afquery dump --db ./db/ --chrom chr17 --start 43044292 --end 43170327 \
  --output brca1.csv --by-sex
```

### 7. Update

```bash
afquery update-db --db ./db/ --add-samples new_batch.tsv
afquery update-db --db ./db/ --update-sample SAMP_007 --set-phenotype I42.0
```

## Output fields

Every query and annotated VCF reports:

| Field | Meaning |
|---|---|
| `AC` | Alt-allele count over eligible samples (FILTER=PASS) |
| `AN` | Total alleles considered, ploidy- and capture-aware |
| `AF` | `AC / AN` |
| `N_HET` | Heterozygous PASS carriers |
| `N_HOM_ALT` | Homozygous-alt PASS carriers |
| `N_HOM_REF` | Samples trusted as homozygous reference |
| `N_FAIL` | Samples with a failed call (FILTER ≠ PASS) at the position (excluded from AC/AN) |
| `N_NO_COVERAGE` | Non-carriers on a partially-covered tech without sufficient evidence to call hom-ref |

VCF INFO field names use the `AFQUERY_` prefix (e.g. `AFQUERY_AF`, `AFQUERY_N_NO_COVERAGE`).

## CLI commands

| Command | Purpose |
|---|---|
| `create-db` | Build a database from a manifest of single-sample VCFs |
| `query` | Point / region / batch AF queries |
| `variant-info` | List samples carrying a variant, with metadata |
| `annotate` | Annotate a VCF with cohort `AFQUERY_*` INFO fields |
| `dump` | Bulk CSV export, optionally disaggregated |
| `update-db` | Add / remove samples, edit metadata, compact |
| `info` | Show database metadata, sample list, changelog |
| `check` | Validate database integrity (scripted exit code) |
| `version show` / `version set` | Inspect or set the database version label |
| `benchmark` | Run synthetic or on-database performance benchmarks |

See the [CLI reference](https://babelomics.github.io/afquery/reference/cli/) for all options.

## How it works

AFQuery indexes each variant as three Roaring Bitmaps — heterozygous PASS carriers, homozygous-alt PASS carriers, and FILTER≠PASS carriers — stored in Apache Parquet, partitioned by chromosome and 1-Mbp positional buckets. Sample metadata (sex, technology, phenotype) is precomputed as bitmaps in SQLite. A query resolves its sample filter into a single candidate bitmap by intersection/difference in microseconds, then intersects it against each variant's genotype bitmaps to compute AC. AN is computed per position from the same candidate bitmap, restricted to samples whose technology actually covers the position (via BED-derived capture indices) and adjusted for ploidy on sex chromosomes. See the [data model reference](https://babelomics.github.io/afquery/reference/data-model/) for details.

## How AFQuery compares

| | AFQuery | bcftools | GATK GenomicsDB | Hail |
|---|---|---|---|---|
| Capture-aware AN | **Yes** | No | No | No |
| Metadata filtering | **Arbitrary labels** | No | No | Custom code |
| Ploidy-aware sex chromosomes | **Yes** | Manual | No | Manual |
| Dynamic subcohort queries | **Yes** | No | Limited | Requires code |
| FILTER / coverage tracking | **Per variant** | Manual | No | Manual |
| Incremental updates | **Yes** | No | **Yes** | No |
| Infrastructure required | **None** | **None** | Java/server | Spark cluster |

### Benchmarks vs. bcftools (1000 Genomes Phase 3, n = 2,504, chr22)

| Workload | AFQuery | bcftools | Speedup |
|---|---|---|---|
| Full-chromosome AC/AN/AF export | ~7.0 s | ~3.8 min | **~33×** |

Point-query latency on AFQuery is ~14 ms and constant from 1K to 50K samples (median over 50 replicates, warm cache).

## Documentation

- [Full documentation](https://babelomics.github.io/afquery/)
- [Getting started](https://babelomics.github.io/afquery/getting-started/)
- [Why local allele frequencies?](https://babelomics.github.io/afquery/getting-started/motivation/)
- [Manifest format](https://babelomics.github.io/afquery/guides/manifest-format/)
- [CLI reference](https://babelomics.github.io/afquery/reference/cli/)
- [FAQ](https://babelomics.github.io/afquery/faq/)

## Citation

If you use AFQuery in your work, please cite:

> AFQuery: fast, capture-aware allele frequency queries on local genomic cohorts.
> *(manuscript in preparation)*

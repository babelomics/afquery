# AFQuery

**Fast, file-based genomic allele frequency queries for large cohorts. No server, no cloud — just files.**

AFQuery stores genotype data as Roaring Bitmaps in Parquet files and answers allele frequency queries in under 100 ms across large cohorts, with flexible filtering by sex, phenotype codes (arbitrary sample labels), and sequencing technology.

---

## When to Use AFQuery

- You need allele frequencies for **phenotype-defined subgroups** (not just whole-population AF)
- You mix **WGS, WES, and panels** in one cohort and need technology-aware AN
- You require **reproducible local AF** computed on your own samples — not just public reference databases
- You run **repeated queries** on the same dataset (annotation, clinical interpretation, research)
- You need **sub-100 ms query latency** without database servers or cloud infrastructure

---

## How It Works

AFQuery pre-indexes genotypes as [Roaring Bitmaps](https://roaringbitmap.org/) in Parquet files. At query time, it intersects carrier bitmaps with eligible-sample bitmaps (determined by sex, phenotype, and capture filters) and counts bits — reducing each query to microsecond-scale bitmap operations with sub-100 ms end-to-end latency.

---

## Features

- **Sub-100 ms queries** — bitmap operations, not VCF scanning. Latency is independent of cohort size.
- **Incremental updates** — add or remove samples without rebuilding the database.
- **Multi-dimensional filtering** — filter by sex, phenotype codes, and sequencing technology with include/exclude semantics.
- **Server-less** — a directory of Parquet files + SQLite. Copy to share, no daemon required.
- **Ploidy-aware** — correct AN on chrX PAR/non-PAR, chrY, and chrM, with hemizygous carriers counted on their own (`N_HEMI`).
- **Technology-aware AN** — per-position capture BED intersection across WGS, WES kits, and panels.
- **Carrier lookup** — list samples carrying any variant with full metadata (sex, tech, phenotypes, genotype, FILTER status).
- **VCF annotation** — add `AFQUERY_AC/AN/AF/N_HET/N_HOM_ALT/N_HEMI/N_HOM_REF/N_FAIL/N_NO_COVERAGE` INFO fields from any sample subset.
- **Audit changelog** — every database operation is recorded for reproducibility.

---

## Architecture

```mermaid
graph TD
    A["🔍 Input VCFs<br/>single-sample"]
    B["📥 Ingest<br/>cyvcf2 reads VCFs →<br/>per-sample Parquet files emitted"]
    C["🏗️ Build<br/>DuckDB aggregates per bucket →<br/>Roaring Bitmaps → Parquet"]
    D["💾 Database on Disk<br/>variants/chr*/bucket_*.parquet<br/>capture/*.pkl<br/>metadata.sqlite<br/>manifest.json"]
    E["⚡ Query Engine<br/>Load bitmap → filter samples →<br/>compute AC/AN/AF<br/>~10-100ms"]
    F["📊 Annotate"]
    G["🔄 Update"]

    A --> B
    B --> C
    C --> D
    D --> E
    E -.->|VCF| F
    E -.->|Batch| G
    G -.->|Updated| D

    style A fill:#e1f5ff
    style B fill:#fff3e0
    style C fill:#f3e5f5
    style D fill:#e8f5e9
    style E fill:#fce4ec
    style F fill:#fce4ec
    style G fill:#fce4ec
```

---


## Where to Start

```mermaid
graph TD
    A["What do you want to do?"]
    B["Build a database<br/>from VCFs"]
    C["Query allele<br/>frequencies"]
    D["Annotate a<br/>patient VCF"]
    E["Classify variants<br/>using ACMG criteria"]
    F["Compare AF across<br/>groups"]

    M["Find carriers of<br/>a variant"]

    A -->|First time| G["5-Min Quickstart"]
    A -->|Build| B
    A -->|Query| C
    A -->|Annotate| D
    A -->|Classify| E
    A -->|Compare| F
    A -->|Carriers| M

    B --> H["Create a Database"]
    C --> I["Query Guide"]
    D --> J["Annotate a VCF"]
    E --> K["ACMG Criteria"]
    F --> L["Cohort Stratification"]
    M --> N["Variant Info"]

    click G "getting-started/quickstart/"
    click H "guides/create-database/"
    click I "guides/query/"
    click J "guides/annotate-vcf/"
    click K "use-cases/acmg-use-cases/"
    click L "use-cases/cohort-stratification/"
    click N "guides/variant-info/"

    style A fill:#e3f2fd
    style G fill:#e8f5e9
```

---

## Next Steps

- [Why Local AF Matters](getting-started/motivation.md) — population gaps, technology bias, and the research behind AFQuery
- [Installation](getting-started/installation.md) — pip, conda, from source
- [Quickstart](getting-started/quickstart.md) — 5-minute end-to-end tutorial
- [Key Concepts](getting-started/concepts.md) — bitmaps, Parquet, manifest, metadata model

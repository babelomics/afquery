# Quickstart

This tutorial walks through building a small AFQuery database and running your first queries. It takes about 5 minutes.

!!! tip "New to AFQuery?"
    If you want to understand how bitmaps, Parquet storage, and metadata filtering work together before diving in, read [Key Concepts](concepts.md) first. Otherwise, follow along — you can always come back to the theory later.

!!! tip "Normalize your VCFs first"
    AFQuery works best with normalized, left-aligned VCFs with ploidy-corrected sex chromosome calls. See [VCF Preprocessing](preprocessing.md) for a reference normalization pipeline using bcftools.

---

## 1. Create a Manifest

The manifest is a TSV file describing your initial cohort. Create `manifest.tsv`:

```tsv
sample_name	vcf_path	sex	tech_name	phenotype_codes
SAMPLE_001	/data/vcfs/sample001.vcf.gz	female	wgs	E11.9,I10
SAMPLE_002	/data/vcfs/sample002.vcf.gz	male	wgs	E11.9
SAMPLE_003	/data/vcfs/sample003.vcf.gz	female	wes_v1	I10
```

Fields:

- `sample_name`: unique identifier
- `vcf_path`: path to single-sample VCF (plain or `.gz`)
- `sex`: `male` or `female`
- `tech_name`: sequencing technology name. Use `WGS` for whole genome sequencing
- `phenotype_codes`: comma-separated phenotype codes (arbitrary strings)

See [Manifest Format](../guides/manifest-format.md) for full details.

---

## 2. Create the Database

```bash
afquery create-db \
  --manifest manifest.tsv \
  --output-dir ./my_db/ \
  --genome-build GRCh38 \
  --bed-dir ./beds/ \
  --threads 12
```

BED files are required for WES and panel technologies. Place one file per technology in the `--bed-dir` directory, named `<tech_name>.bed` (e.g., `wes_v1.bed`). WGS technologies do not need a BED file.


The command will:
1. Ingest all VCFs
2. Build Roaring Bitmap Parquet files per chromosome/bucket
3. Write `manifest.json` and `metadata.sqlite`

---

## 3. Inspect the Database (optional)

```bash
afquery info --db ./my_db/
```

---

## 4. Query a Single Position

```bash
afquery query --db ./my_db/ --locus chr1:925952
```

Example output:
```
chr1:925952 G>A  AC=2  AN=6  AF=0.3333  n_eligible=3  N_HET=2  N_HOM_ALT=0  N_HOM_REF=1  N_FAIL=0
```

Filter to female samples with a specific phenotype:

```bash
afquery query \
  --db ./my_db/ \
  --locus chr1:925952 \
  --sex female \
  --phenotype E11.9
```

See [Sample Filtering](../guides/sample-filtering.md) for the full include/exclude syntax.

!!! note "Warnings for missing data"
    If a phenotype code, technology name, or chromosome is not found in the database, afquery prints a warning to stderr and returns empty results. Use `--no-warn` to suppress these warnings.

---

## 5. Inspect Carriers (optional)

See which samples carry the variant you just queried:

```bash
afquery variant-info --db ./my_db/ --locus chr1:925952
```

This lists each carrier with their sex, technology, phenotype codes, genotype (het/hom), and FILTER status. See [Variant Info](../guides/variant-info.md) for details.

---

## 6. Query a Region

```bash
afquery query \
  --db ./my_db/ \
  --region chr1:900000-1000000
```

---

## 7. Annotate a VCF

Given a VCF with variants you want to annotate:

```bash
afquery annotate \
  --db ./my_db/ \
  --input variants.vcf \
  --output annotated.vcf \
  --threads 12
```

The output VCF gains INFO fields (see [Annotate a VCF](../guides/annotate-vcf.md) for parallelism options and downstream usage):

| Field | Number | Description |
|-------|--------|-------------|
| `AFQUERY_AC` | A (per ALT) | Allele count |
| `AFQUERY_AN` | 1 (per site) | Allele number |
| `AFQUERY_AF` | A (per ALT) | Allele frequency |
| `AFQUERY_N_HET` | A (per ALT) | Heterozygous sample count |
| `AFQUERY_N_HOM_ALT` | A (per ALT) | Homozygous alt sample count |
| `AFQUERY_N_HOM_REF` | A (per ALT) | Homozygous ref sample count |
| `AFQUERY_N_FAIL` | 1 (per site) | Samples with a failed call (FILTER≠PASS), left out of AN |
| `AFQUERY_N_NO_COVERAGE` | A (per ALT) | Eligible samples lacking coverage evidence (0 unless a coverage-evidence filter is active) |

---

## Next Steps

- [Key Concepts](concepts.md) — understand how bitmaps, Parquet, and metadata filtering work together
- [Sample Filtering](../guides/sample-filtering.md) — full syntax for phenotype, sex, and technology filters
- [Variant Info](../guides/variant-info.md) — list carriers of any variant with metadata
- [Annotate a VCF](../guides/annotate-vcf.md) — annotation options, parallelism, and downstream usage
- [ACMG Criteria](../use-cases/acmg-use-cases.md) — applying local AF to BA1, PM2, and PS4

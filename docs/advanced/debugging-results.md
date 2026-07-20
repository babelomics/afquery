# Debugging Results

When AFQuery returns unexpected results — AN=0, surprising AF values, or missing variants — use this diagnostic checklist to identify the root cause.

---

## Diagnostic Checklist

### 1. Unexpected AN=0

AN=0 means no eligible samples at the queried position. Work through these checks in order:

| Check | Command | What to look for |
|-------|---------|-----------------|
| Chromosome normalization | `afquery query --db ./db/ --chrom 1 ...` vs `--chrom chr1` | Database may use `chr1` while you're querying `1` (or vice versa). Check `manifest.json` for `genome_build`. |
| Position exists in database | `afquery query --db ./db/ --locus chr1:12345678` | If no result at all, the variant was not observed in any sample during ingestion. |
| BED coverage (WES) | `afquery info --db ./db/` | If all eligible samples are WES and the position is outside capture regions, AN=0 is correct. |
| Sample filter too restrictive | Remove `--phenotype` and `--sex` filters | Query with no filters first. If AN>0 without filters, the filter is excluding all samples. |
| WES samples missing from AN | Compare `variant-info` carriers against `query` counts | If `variant-info` shows many WES carriers but AN reflects only the WGS samples, that technology's capture index is not matching the position. `afquery` warns on open when a capture BED is empty or matches no known chromosome, naming the offending contigs. |
| AN=0 only on chrM | `afquery query --db ./db/ --locus chrM:3243` | Mitochondrial coverage is easy to lose on its own: a BED naming it `chrMT` still matches on the autosomes, so nothing else looks wrong. Names are normalized (`MT`, `M`, `chrMT` → `chrM`) since v0.4.0 — databases built before that may need the capture BED re-checked. |
| Technology filter | Remove `--tech` filter | Check if any samples match the requested technology. |

### 2. Unexpected AF Value

| Symptom | Likely Cause | Resolution |
|---------|-------------|------------|
| AF higher than expected | Cohort enriched for disease with this variant | Compare with `--phenotype ^disease_code` to get control AF |
| AF lower than expected | Many WES samples without coverage at this position → diluted AN | Filter by `--tech wgs` to check WGS-only AF |
| AF=1.0 | All eligible samples carry the variant | Check if AN is very small (e.g., AN=2 means just 1 sample) |
| AF differs from gnomAD | Expected — local cohort AF ≠ global population AF | This is by design; see [Why Cohort-Specific AF Matters](../getting-started/concepts.md#why-cohort-specific-af-matters) |

### 3. Missing Variants

A variant you expect to find is not in the database:

| Check | Details |
|-------|---------|
| Was it in the source VCFs? | AFQuery only stores variants present in ingested VCFs |
| Was the alt allele ever called? | A variant is stored only if at least one sample carries the alt allele or has a failed call there. A site seen only as `0/0` in every sample is genuinely absent. Note that non-PASS calls *are* kept (in `fail_bitmap`) — such a variant is not missing, it just shows `AC=0`. |
| Multiallelic sites | AFQuery stores each ALT separately. Query the specific ALT allele, not just position |
| Chromosome naming | Ensure consistent `chr` prefix usage |

### 4. Unexpected N_FAIL > 0

`N_FAIL > 0` means some eligible samples had a call with `FILTER≠PASS` at this position. These samples are excluded from AC, but they remain eligible and still count in AN. This is usually benign (1–2 samples), but a high N_FAIL warrants investigation:

| N_FAIL relative to n_eligible | Likely cause |
|---|---|
| 1–2 samples | Isolated low-quality calls — not concerning |
| > 5% of n_eligible | Systematic sequencing artifact at this site |
| All eligible samples | Site-wide QC failure — AF=0 but variant is present in source VCFs |

To inspect a site with high N_FAIL, query with `--format json` to see all fields:

```bash
afquery query --db ./db/ --locus chr1:12345678 --format json
```

!!! tip "Identify failing samples"
    Use `afquery variant-info --db ./db/ --locus chr1:12345678` to see exactly which samples have FAIL status and their metadata (technology, phenotype codes). This helps determine if failures cluster in a specific technology or sample subset. See [Variant Info](../guides/variant-info.md).

If N_FAIL is consistently high across many sites, check the variant calling pipeline and FILTER field settings in your VCFs.

---

## Diagnostic Commands

### Database Info

```bash
afquery info --db ./db/
```

Shows: sample count, technology list, schema version, genome build, PASS-only status.

### Check Database Integrity

```bash
afquery check --db ./db/
```

Validates: manifest consistency, Parquet file integrity, capture index presence.

### Query with Full Output

```bash
# Point query with all fields visible
afquery query --db ./db/ --locus chr1:12345678 --format json
```

JSON format shows all fields including N_HET, N_HOM_ALT, N_HOM_REF, and N_FAIL — useful for understanding the composition of the result.


### Direct Metadata Inspection

```python
import sqlite3

conn = sqlite3.connect("./db/metadata.sqlite")

# List all phenotype codes and sample counts
cursor = conn.execute("""
    SELECT code, COUNT(*) as n_samples
    FROM sample_phenotypes
    GROUP BY code
    ORDER BY n_samples DESC
""")
for row in cursor:
    print(f"{row[0]:30s}  {row[1]} samples")

# List technologies
cursor = conn.execute("""
    SELECT technology, COUNT(*) as n_samples
    FROM samples
    GROUP BY technology
""")
for row in cursor:
    print(f"{row[0]:20s}  {row[1]} samples")

conn.close()
```

---

## Common Root Causes

| Symptom | Root Cause | Fix |
|---------|-----------|-----|
| AN=0 for all queries | Wrong `--db` path or empty database | Verify path; run `afquery info` |
| AN=0 for specific region | WES-only cohort, position outside capture | Check BED file coverage |
| AN much lower than sample count | Mixed WGS/WES, position outside WES capture | Filter by `--tech wgs` to isolate |
| AF=None in output | AN=0 (division by zero) | See AN=0 diagnosis above |
| Different AF between `query` and `annotate` | Different default filters or phenotype context | Ensure same `--phenotype`, `--sex`, `--tech` flags |
| N_FAIL high at a site | Systematic QC failure in source VCFs at this position | Inspect site with `--format json`; check VCF FILTER annotations |

---

## Next Steps

- [Understanding Output](../getting-started/understanding-output.md) — what each field means
- [FAQ](../faq.md) — common questions and answers
- [Troubleshooting](../troubleshooting.md) — error messages and solutions
- [FILTER=PASS Tracking](filter-pass-tracking.md) — how PASS filtering works

# Population-Specific Allele Frequency

## Scenario

You manage a rare disease registry for a specific regional population. A variant in a candidate gene is reported at AF=0.001 in gnomAD (classified as very rare, supporting pathogenicity). But your cohort shows this variant at much higher frequency locally. Correct classification requires a locally calibrated AF.

## Why Standard Databases Fall Short

gnomAD aggregates data from wide ancestry groups worldwide. A variant at AF=0.001 at European (non finish) group may be at AF=0.01 in a specific Iberian, French, or German population — a 10× difference that directly impacts ACMG criteria:

- BA1 (benign standalone): AF > 0.05 in any population
- BS1 (benign strong): AF > 0.01 in a matched population
- PM2 (pathogenic moderate): AF < 0.001

Without a locally calibrated database, these thresholds are applied to population-averaged frequencies, potentially misclassifying variants enriched in your specific ancestry.

## Solution with AFQuery

Build your AFQuery database from samples representative of your population. Queries then reflect the actual frequency distribution in your cohort — the most relevant reference for your patients.

## Querying Population-Specific Subsets

### 1. Build a population-specific database

```bash
afquery create-db \
  --manifest cohort_manifest.tsv \
  --output-dir ./local_db/ \
  --genome-build GRCh38
```

### 2. Query a candidate variant

```bash
afquery query --db ./local_db/ --locus chr1:925952
```

```
chr1:925952 G>A  AC=142  AN=2742  AF=0.0518  n_eligible=1371  N_HET=138  N_HOM_ALT=2  N_HOM_REF=1231  N_FAIL=0  N_NO_COVERAGE=0
```

AF=0.052 locally, compared to gnomAD NFE AF=0.005 (10× higher). Under ACMG BS1, this variant is likely benign in this population.

### 3. Compare with a specific subgroup

If your cohort includes samples from multiple ancestry groups tagged in the manifest:

```bash
# Samples tagged with population label
afquery query --db ./local_db/ --locus chr1:925952 --phenotype iberian_population
```

### 4. Python: systematic comparison across many variants

```python
from afquery import Database

db = Database("./local_db/")
variants = [(925952, "G", "A"), (1014541, "C", "T"), (3456789, "A", "G")]

results = db.query_batch("chr1", variants=variants)
for r in results:
    print(f"{r.variant.pos}: local_AF={r.AF:.4f}  AN={r.AN}")
```

## Biological Interpretation

| Variant | gnomAD NFE AF | Local AF | Fold difference | ACMG impact |
|---------|--------------|----------|-----------------|-------------|
| chr1:925952 G>A | 0.005 | 0.052 | 10× | BS1 (benign strong) |
| chr1:1014541 C>T | 0.0001 | 0.0001 | 1× | PM2 supported |
| chr1:3456789 A>G | 0.0008 | 0.009 | 11× | BS1 in local pop |

For the third variant, gnomAD would support PM2 (very rare), but local frequency shows it is common in this population (BS1 level). Using gnomAD alone would misclassify this variant.

## Next Steps

- [Cohort Stratification](cohort-stratification.md) — compare AF across subgroups within one database
- [Sample Filtering](../guides/sample-filtering.md) — filter by population tags
- [Clinical Prioritization](clinical-prioritization.md) — apply AF filters in VCF annotation

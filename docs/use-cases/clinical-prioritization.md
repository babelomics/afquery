# Clinical Variant Prioritization

## Scenario

A rare disease patient has undergone whole-exome sequencing. After standard filtering, 500,000 candidate variants remain. You want to annotate each variant with cohort-specific allele frequency and filter to those that are genuinely rare in your population — removing variants that are common locally but might appear rare in gnomAD.

## Why Standard Databases Fall Short

gnomAD provides an excellent first filter, but:

1. **Population mismatch**: A variant at AF=0.001 in gnomAD may be at AF=0.02 in your local cohort — common locally but appearing rare globally. This discrepancy often reflects natural allele frequency variation between subpopulations driven by genetic drift, founder effects, and historical bottlenecks: alleles that are rare on a global scale may have reached appreciable frequencies in geographically or ethnically isolated groups.

2. **Fine-grained Control Cohort Selection**: Unlike resources such as gnomAD, where allele frequencies are derived from largely phenotype-agnostic populations, AFQuery allows the dynamic inclusion or exclusion of samples based on any annotated feature. This is particularly valuable in rare disease studies, where overlapping genetic architectures may confound analyses: for example, samples associated with a related condition can be selectively excluded to avoid bias. Because phenotypes are treated as flexible annotations, this control extends to any variable of interest, enabling more precise and context-aware frequency estimation.


3. **Local artifacts**: Systematic sequencing artifacts specific to your pipeline or capture kit manifest as recurrent variants that appear rare in gnomAD but accumulate high frequency in your cohort. These are best identified by elevated AF in your local database paired with low allele number (AN) or high N_FAIL counts, indicating poor genotype quality at the site.

AFQuery lets you apply cohort-specific AF as an additional filter layer on top of gnomAD, removing locally common variants that standard databases miss.

## Step-by-Step Example

### 1. Annotate patient VCF with cohort AF

```bash
afquery annotate \
  --db ./db/ \
  --input patient.vcf.gz \
  --output patient_annotated.vcf.gz \
  --threads 16
```

This adds to each variant:

- `AFQUERY_AC`: allele count in cohort
- `AFQUERY_AN`: allele number (eligible samples at this position)
- `AFQUERY_AF`: allele frequency
- `AFQUERY_N_HET`, `AFQUERY_N_HOM_ALT`, `AFQUERY_N_HEMI`: genotype counts (heterozygous, homozygous, and hemizygous males on chrX/chrY)

### 2. Filter for rare variants with reliable AN

```bash
bcftools filter \
  -i 'AFQUERY_AF < 0.001 && AFQUERY_AN >= 1000' \
  patient_annotated.vcf.gz \
  -o patient_rare.vcf.gz
```

The `AFQUERY_AN >= 1000` threshold ensures the AF estimate is based on sufficient data. An AF estimate from AN=10 is meaningless — with AN=10, even a single carrier gives AF=0.1.

### 3. Handle variants not in the cohort

Variants absent from the database have `AFQUERY_AN=0` (no eligible samples, or variant not observed). These require separate treatment:

```bash
# Variants NOT in cohort (AN=0): treat as novel
bcftools filter -i 'AFQUERY_AN == 0' patient_annotated.vcf.gz

# Variants in cohort with sufficient coverage
bcftools filter -i 'AFQUERY_AN >= 1000 && AFQUERY_AF < 0.001' patient_annotated.vcf.gz
```

### 4. Annotation with subgroup AF (optional)

If you want AF relative to a matched control group:

```bash
afquery annotate \
  --db ./db/ \
  --input patient.vcf.gz \
  --output patient_control_af.vcf.gz \
  --phenotype ^rare_disease     # AF in non-rare-disease samples
```

### 5. Python workflow

```python
import cyvcf2
from afquery import Database

db = Database("./db/")

# Annotate and filter in memory
vcf = cyvcf2.VCF("patient.vcf.gz")
rare_candidates = []

for variant in vcf:
    results = db.query(variant.CHROM, pos=variant.POS, alt=variant.ALT[0])
    if not results:
        rare_candidates.append(variant)  # Not in cohort → novel
        continue
    r = results[0]
    if r.AN >= 1000 and r.AF < 0.001:
        rare_candidates.append(variant)

print(f"Rare candidates: {len(rare_candidates)}")
```

## Biological Interpretation

| Filter | Retained | Removed | Reason |
|--------|----------|---------|--------|
| AFQUERY_AN >= 1000 | 45,000 | 455,000 | Insufficient cohort coverage |
| AFQUERY_AF < 0.001 | 1,200 | 43,800 | Locally common variants |
| Novel (AN=0) | 300 | — | Not observed in cohort |

A typical clinical pipeline retains ~1,500 rare/novel candidates after cohort AF filtering, compared to 500,000 before.

**AN threshold guidance:**

- `AN >= 100`: minimum for any AF interpretation
- `AN >= 500`: recommended for rare variant filtering
- `AN >= 1000`: conservative threshold for robust AF estimates

For detailed ACMG workflows with worked examples and AN threshold guidance, see [ACMG Criteria (BA1/PM2/PS4)](acmg-use-cases.md).


### 6. Inspect carriers of rare candidates

After identifying rare variants, use `variant-info` to see who carries each one — useful for confirming case/control enrichment or checking sample quality:

```bash
afquery variant-info --db ./db/ --locus chr2:166845670 --ref G --alt A
```

This lists each carrier with their phenotype codes, sequencing technology, genotype, and FILTER status. See [Variant Info](../guides/variant-info.md) for filtering options.

## Next Steps

- [Variant Info](../guides/variant-info.md) — list carriers of any variant with metadata
- [Annotate a VCF](../guides/annotate-vcf.md) — full annotation options
- [FILTER=PASS Tracking](../advanced/filter-pass-tracking.md) — using N_FAIL in filtering
- [Population-Specific AF](population-specific-af.md) — local vs. gnomAD comparison

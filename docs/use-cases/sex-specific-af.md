# Sex-Specific Allele Frequency

## Scenario

You are analyzing a candidate X-linked recessive gene. A hemizygous variant in males may be pathogenic, but the same variant in females acts as a heterozygous carrier state. Allele frequency means something different for each sex due to ploidy, and a combined AF obscures the clinically relevant frequencies.

## Why Standard Databases Fall Short

Most population databases report a single AF for chrX variants, combining diploid female calls and haploid male calls. This average is appropriate for some purposes but misleading for X-linked recessive analysis, where the male hemizygous frequency is the clinically relevant metric.

For a variant at chrX non-PAR with AC=10 in 500 males (AN=500) and AC=5 in 1000 females (AN=2000):

- Male hemizygous rate: 10/500 = 2%
- Female carrier rate: 5/2000 = 0.25%
- Combined AF: 15/2500 = 0.6%

The combined AF (0.6%) obscures the fact that 2% of males are hemizygous carriers — a crucial number for X-linked recessive risk assessment.

## Solution with AFQuery

AFQuery applies ploidy-aware AN computation automatically. For chrX non-PAR positions:

- Males contribute AN=1 (haploid)
- Females contribute AN=2 (diploid)

Using `--sex male` returns AF over a purely haploid denominator — the true hemizygous frequency. Using `--sex female` returns AF over the diploid denominator.

## Computing Sex-Specific AF at chrX

### 1. Query full cohort (combined)

```bash
afquery query --db ./db/ --locus chrX:153296777
```

```
chrX:153296777 C>T  AC=15  AN=2500  AF=0.0060  n_eligible=1500  N_HET=5  N_HOM_ALT=0  N_HEMI=10  N_HOM_REF=1485  N_FAIL=0  N_NO_COVERAGE=0
```

AN=2500: 500 males × 1 + 1000 females × 2 = 2500 (mixed ploidy). The 10 hemizygous male carriers are counted in N_HEMI, and the 5 heterozygous females in N_HET; N_HOM_ALT counts only females homozygous for the allele (none here). N_HOM_REF holds 490 males hemizygous for the reference and 995 homozygous-reference females.

### 2. Query males only (hemizygous frequency)

```bash
afquery query --db ./db/ --locus chrX:153296777 --sex male
```

```
chrX:153296777 C>T  AC=10  AN=500  AF=0.0200  n_eligible=500  N_HET=0  N_HOM_ALT=0  N_HEMI=10  N_HOM_REF=490  N_FAIL=0  N_NO_COVERAGE=0
```

Hemizygous rate: 2%. Males are haploid at chrX non-PAR, so every carrier is hemizygous: N_HEMI = AC = 10, and N_HET and N_HOM_ALT are 0. A male whose call there is heterozygous (`0/1`) is counted in N_FAIL instead, since a single copy cannot hold both alleles.

### 3. Query females only (carrier frequency)

```bash
afquery query --db ./db/ --locus chrX:153296777 --sex female
```

```
chrX:153296777 C>T  AC=5  AN=2000  AF=0.0025  n_eligible=1000  N_HET=5  N_HOM_ALT=0  N_HEMI=0  N_HOM_REF=995  N_FAIL=0  N_NO_COVERAGE=0
```

Carrier rate: 0.25%

### 4. Python API comparison

```python
from afquery import Database

db = Database("./db/")

male = db.query("chrX", pos=153296777, sex="male")[0]
female = db.query("chrX", pos=153296777, sex="female")[0]
combined = db.query("chrX", pos=153296777)[0]

print(f"Male hemizygous rate:  {male.AF:.4f}  (AN={male.AN}, hemizygotes={male.N_HEMI})")
print(f"Female carrier rate:   {female.AF:.4f}  (AN={female.AN})")
print(f"Combined AF:           {combined.AF:.4f}  (AN={combined.AN})")
```

### 5. Annotating with sex-stratified AF

```bash
# Annotate using only male samples (for X-linked recessive analysis)
afquery annotate \
  --db ./db/ \
  --input patient.vcf.gz \
  --output annotated_males.vcf.gz \
  --sex male
```

## Biological Interpretation

| Metric | Value | Clinical relevance |
|--------|-------|-------------------|
| Male hemizygous rate | 2% | ACMG BS1 for X-linked recessive (exceeds 0.01 threshold) |
| Female carrier rate | 0.25% | Heterozygous carriers; does not inform X-linked recessive pathogenicity |
| Combined AF | 0.6% | Misleading for X-linked analysis — neither the hemizygous nor carrier rate |

For X-linked recessive disorders, always use `--sex male` when evaluating the hemizygous carrier frequency.

## PAR Regions

The Pseudoautosomal Regions (PAR1/PAR2) on chrX behave like autosomes. At PAR positions, males contribute AN=2 and sex-stratified queries are less meaningful. AFQuery applies the correct PAR/non-PAR boundary automatically.

See [Ploidy & Special Chromosomes](../advanced/ploidy-and-sex-chroms.md) for exact PAR coordinates and genotype counting rules for haploid regions.

## Next Steps

- [Ploidy & Special Chromosomes](../advanced/ploidy-and-sex-chroms.md) — PAR regions, ploidy rules
- [Sample Filtering](../guides/sample-filtering.md) — sex filter syntax
- [Clinical Prioritization](clinical-prioritization.md) — VCF annotation with sex filter

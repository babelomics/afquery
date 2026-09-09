# Ploidy & Special Chromosomes

AFQuery computes ploidy-aware AN for sex chromosomes (chrX, chrY) and the mitochondrial chromosome (chrM). This ensures that allele frequencies are correct when querying these chromosomes, where the number of alleles per sample differs from the diploid autosomes.

!!! note "Chromosome name normalization"
    AFQuery accepts `MT`, `chrMT`, and `chrM` as input; output always uses `chrM`.

---

## Ploidy Rules

| Chromosome | Female AN contribution | Male AN contribution |
|------------|----------------------|---------------------|
| Autosomes (chr1–22) | 2 | 2 |
| chrX (non-PAR) | 2 | 1 |
| chrX (PAR1, PAR2) | 2 | 2 |
| chrY | 0 | 1 |
| chrM | 1 | 1 |

For each eligible sample at a given position, AFQuery adds the appropriate ploidy count to AN based on the sample's sex and the chromosome/position.

A sample contributing ploidy 0 is not eligible at that position at all: a female has no chrY to genotype, so she is neither a carrier nor homozygous reference there. On chrY, `n_eligible`, `N_HOM_REF` and `AN` therefore count males only.

---

## Pseudoautosomal Regions (PAR)

The pseudoautosomal regions on chrX and chrY behave like autosomes — both males and females contribute AN=2 on chrX PAR. PAR coordinates by genome build:

<!-- Source of truth: src/afquery/constants.py PAR dict -->

### GRCh38

**chrX:**

| Region | Start | End |
|--------|-------|-----|
| PAR1 | 10,001 | 2,781,479 |
| PAR2 | 155,701,383 | 156,030,895 |

**chrY:**

| Region | Start | End |
|--------|-------|-----|
| PAR1 | 10,001 | 2,781,479 |
| PAR2 | 56,887,903 | 57,217,415 |

### GRCh37 / hg19

**chrX:**

| Region | Start | End |
|--------|-------|-----|
| PAR1 | 60,001 | 2,699,520 |
| PAR2 | 154,931,044 | 155,260,560 |

**chrY:**

| Region | Start | End |
|--------|-------|-----|
| PAR1 | 10,001 | 2,649,520 |
| PAR2 | 59,034,050 | 59,363,566 |

Positions within PAR1 or PAR2 on chrX are treated as diploid for all samples.

---

## Effect on AF Queries

### chrY

Querying chrY with `--sex female` returns `AN=0` (females have no Y chromosome):

```bash
afquery query --db ./db/ --locus chrY:2787758 --sex female
# chrY:2787758 — no results (AN=0 for all variants)

afquery query --db ./db/ --locus chrY:2787758 --sex male
# chrY:2787758 C>T  AC=3  AN=856  AF=0.0035  n_eligible=856  N_HET=0  N_HOM_ALT=3  N_HOM_REF=853  N_FAIL=0
```

### chrX non-PAR

Male samples contribute AN=1, female samples contribute AN=2. This means a cohort of 500 females and 500 males has AN = 500×2 + 500×1 = 1500 at a non-PAR X position.

```bash
afquery query --db ./db/ --locus chrX:100000000
# N_total = 1000 samples, AN = 1500 (not 2000)
```

### chrM

All samples are haploid at mitochondrial loci:

```bash
afquery query --db ./db/ --locus chrM:3243
# AN = n_samples (one allele per sample)
```

---

## Genotype Counting

### Counting Identity

For every query result, the following identity holds:

**N_HET + N_HOM_ALT + N_HOM_REF + N_FAIL + N_NO_COVERAGE = n_eligible**

This can be used to validate results. N_HOM_REF is the number of eligible samples that are homozygous reference (i.e., do not carry the alt allele and passed quality filters). N_NO_COVERAGE is 0 unless a coverage-evidence filter is active — see [Coverage Evidence](coverage-evidence.md).

!!! note "Mutual exclusivity"
    N_HET, N_HOM_ALT, N_HOM_REF, N_FAIL, and N_NO_COVERAGE are mutually exclusive. A sample with a non-ref allele but FILTER≠PASS is counted in N_FAIL only — it does not appear in N_HET or N_HOM_ALT. Likewise, N_HOM_REF counts only PASS-filtered samples.

### chrX non-PAR

- A male with GT=`1` contributes AC=1, AN=1
- A female with GT=`0/1` contributes AC=1, AN=2
- A female with GT=`1/1` contributes AC=2, AN=2

N_HET and N_HOM_ALT are counted per sample (not per allele):

- Males at chrX non-PAR (haploid positions) are counted in **N_HOM_ALT** when GT=1, because all alleles at that position are alternate. N_HET is reserved for diploid positions where both reference and alternate alleles are present.
- Females at chrX with GT=`0/1` are counted in N_HET; with GT=`1/1` in N_HOM_ALT.

### chrY

chrY is fully haploid (males only, females contribute AN=0):

- All carriers are counted in **N_HOM_ALT** (never N_HET)
- N_HET is always 0 on chrY

### chrM

chrM is haploid for all samples (both sexes contribute AN=1):

- All carriers are counted in **N_HOM_ALT** (never N_HET)
- N_HET is always 0 on chrM

!!! warning "N_HET is always 0 on haploid regions"
    On chrY, chrM, and chrX non-PAR for males, all carriers are counted in N_HOM_ALT. N_HET is 0 because haploid samples have only one allele copy — there is no heterozygous state. This is correct behavior, not a bug.

---

## Sex Filter Interaction

When `--sex female` is used on chrX (non-PAR), AN is purely diploid:

- Each eligible female contributes AN=2
- AF is computed over a fully diploid denominator

When `--sex male` is used on chrX (non-PAR), AN is purely haploid:

- Each eligible male contributes AN=1
- AF reflects the observed allele frequency in haploid male calls

This makes it straightforward to compare X-linked variant frequencies between sexes without manual ploidy adjustment.

---

## Next Steps

- [Sex-Specific AF](../use-cases/sex-specific-af.md) — X-linked variant analysis using sex-stratified queries
- [Key Concepts](../getting-started/concepts.md) — AC/AN/AF overview and the ploidy table
- [Sample Filtering](../guides/sample-filtering.md) — `--sex` filter syntax

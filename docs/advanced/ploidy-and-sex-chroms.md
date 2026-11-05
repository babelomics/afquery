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
| chrY (non-PAR) | 0 | 1 |
| chrY (PAR1, PAR2) | 0 | 0 — counted on chrX |
| chrM | 1 | 1 |

For each eligible sample at a given position, AFQuery adds the appropriate ploidy count to AN based on the sample's sex and the chromosome/position.

A sample contributing ploidy 0 is not eligible at that position at all: a female has no chrY to genotype, so she is neither a carrier nor homozygous reference there. On chrY, `n_eligible`, `N_HOM_REF` and `AN` therefore count males only. In the chrY PARs nobody is eligible (see below).

---

## Pseudoautosomal Regions (PAR)

PAR1 and PAR2 are the same sequence at the two ends of chrX and chrY, and they recombine like an autosome: a male carries one copy on his X and one on his Y, so at a PAR locus every sample is diploid. The reference genome holds that sequence twice, once on each chromosome, so AFQuery counts every PAR locus once, at its chrX coordinates:

- **chrX PAR**: diploid for both sexes — males and females contribute AN=2.
- **chrY PAR**: not genotyped. Nobody is eligible there, so a query returns no results and AN is 0. This keeps a PAR variant from being counted twice, and keeps a male carrier from being taken for a hemizygote.

Most analysis references hard-mask the chrY PARs (the GRCh38 analysis sets do), so variant callers report PAR variants on chrX anyway. Calls at chrY PAR coordinates only appear when the reads were aligned to an unmasked reference; they are stored, but not counted. A point query (`query --locus`, `variant-info`) at such a position warns that it lies in a PAR; region and batch queries, `dump` and `annotate` skip it silently (`annotate` writes `AFQUERY_AN=0`).

PAR coordinates by genome build:

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

Positions within PAR1 or PAR2 on chrX are treated as diploid for all samples; positions within PAR1 or PAR2 on chrY are not genotyped.

---

## Effect on AF Queries

### chrY

Querying chrY with `--sex female` returns `AN=0` (females have no Y chromosome):

```bash
afquery query --db ./db/ --locus chrY:2787758 --sex female
# chrY:2787758 — no results (AN=0 for all variants)

afquery query --db ./db/ --locus chrY:2787758 --sex male
# chrY:2787758 C>T  AC=3  AN=856  AF=0.0035  n_eligible=856  N_HET=0  N_HOM_ALT=3  N_HOM_REF=853  N_FAIL=0  N_NO_COVERAGE=0
```

Inside a chrY PAR (here GRCh38 PAR1), nothing is reported, whatever the filters:

```bash
afquery query --db ./db/ --locus chrY:1500000
# AfqueryWarning: chrY:1500000 lies in a pseudoautosomal region (PAR); PAR variants
# are counted at their chrX coordinates, so nobody is genotyped here.
# No variants found for the given filters.
```

Query the same locus at its chrX coordinates instead (`chrX:1500000`).

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

For every query result at a biallelic site, the following identity holds:

**N_HET + N_HOM_ALT + N_HOM_REF + N_FAIL + N_NO_COVERAGE = n_eligible**

This can be used to validate results. At a multi-allelic site, eligible samples that carry only another ALT allele at the position fall in none of these categories, so the sum is lower by exactly that number of samples (see [Multi-allelic sites](../getting-started/understanding-output.md#multi-allelic-sites)). N_HOM_REF is the number of eligible samples that are homozygous reference (i.e., carry no ALT allele at the position and have no failed call there). N_NO_COVERAGE is 0 unless a coverage-evidence filter is active — see [Coverage Evidence](coverage-evidence.md).

!!! note "Mutual exclusivity"
    N_HET, N_HOM_ALT, N_HOM_REF, N_FAIL, and N_NO_COVERAGE are mutually exclusive. A sample with a failed call (FILTER≠PASS) for any allele at the position is counted in N_FAIL only — it does not appear in N_HET, N_HOM_ALT or N_HOM_REF, and is left out of AC and AN.

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

### Genotypes inconsistent with ploidy

Ploidy comes from each sample's sex and the position, not from how the VCF writes the genotype. A male carrier on chrX non-PAR or chrY counts as one allele whether its GT is written haploid (`1`) or diploid (`1/1`): both name a single allele. A call that names more alleles than the sample has cannot be counted that way:

- on chrX non-PAR or chrY in a male, a call naming two different alleles: a heterozygous call (`0/1`) or two ALT alleles (`1/2`). A male has a single copy there, so REF and ALT cannot both be his genotype; such a call usually comes from mis-mapped reads or contamination, and a 47,XXY sample registered as male shows it genuinely;
- on chrM, two different ALT alleles (`1/2`). A `0/1` call on chrM is not inconsistent: it is heteroplasmy, a mix of mtDNA copies, and the sample is a carrier (see [chrM](#chrm_1));
- on a diploid position, more than two copies across the records at that position (e.g. `1/1` for one ALT and `0/1` for another).

`create-db` and `update-db --add-samples` record every call of the sample at such a position as a failed call, so it counts in `N_FAIL` and is left out of AC and AN (counted as hom-ref with `--fail-as-ref`), and log a warning naming the sample and the positions.

A cohort called in diploid mode on the sex chromosomes therefore shows male `N_FAIL` at chrX positions where the caller saw both alleles in a male. That is the intended reading, the one gnomAD applies too: those genotypes are unknown, not hom-ref and not carriers. Do not rewrite them with `bcftools +fixploidy` beforehand; it keeps the first allele of each call (`0/1` becomes `0`), which turns them into hom-ref without a trace (see [Preprocessing](../getting-started/preprocessing.md#ploidy-of-sex-chromosomes-and-chrm)).

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

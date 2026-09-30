# Understanding Output

This page explains what each field in AFQuery output means and how to interpret special cases.

---

## Output Fields

| Field | Type | Description |
|-------|------|-------------|
| **AC** | int | Allele count — number of alt allele copies in eligible samples |
| **AN** | int | Allele number — total alleles examined (adjusted for ploidy and eligible samples) |
| **AF** | float | Allele frequency — `AC / AN`. `None` when AN=0 |
| **N_HET** | int | Number of eligible samples heterozygous for the alt allele (GT=0/1) |
| **N_HOM_ALT** | int | Number of eligible samples homozygous for the alt allele (GT=1/1 or GT=1). Includes haploid carriers on sex chromosomes and chrM. See [Ploidy](../advanced/ploidy-and-sex-chroms.md#genotype-counting). |
| **N_HOM_REF** | int | Number of eligible samples homozygous reference (GT=0/0 or GT=0). Samples carrying another ALT allele at the same position are not counted (see [Multi-allelic sites](#multi-allelic-sites)). |
| **n_eligible** | int | Number of eligible samples — those passing the sex/phenotype/tech filters *and* covered at this position |
| **N_FAIL** | int | Number of eligible samples whose call at this position had FILTER≠PASS. These samples are counted *only* in N_FAIL — not in N_HET, N_HOM_ALT, or N_HOM_REF — but they stay eligible and still count toward AN. |
| **N_NO_COVERAGE** | int | Number of eligible samples whose tech lacks coverage evidence at this position. Excluded from `N_HOM_REF` to keep AC/AN conservative. Always `0` unless a coverage-evidence filter is active. See [Coverage Evidence](../advanced/coverage-evidence.md). |


---

## Output Formats

### Text (default)

```bash
afquery query --db ./db/ --locus chr1:925952 --ref G --alt A
```

```
chr1:925952 G>A  AC=3  AN=120  AF=0.0250  n_eligible=60  N_HET=1  N_HOM_ALT=1  N_HOM_REF=57  N_FAIL=0  N_NO_COVERAGE=0
```

### TSV

```bash
afquery query --db ./db/ --locus chr1:925952 --ref G --alt A --format tsv
```

```
chrom	pos	ref	alt	AC	AN	AF	n_eligible	N_HET	N_HOM_ALT	N_HOM_REF	N_FAIL	N_NO_COVERAGE
chr1	925952	G	A	3	120	0.025000	60	1	1	57	0	0
```

### JSON

```bash
afquery query --db ./db/ --locus chr1:925952 --ref G --alt A --format json
```

```json
{
  "chrom": "chr1",
  "pos": 925952,
  "ref": "G",
  "alt": "A",
  "AC": 3,
  "AN": 120,
  "AF": 0.025,
  "n_eligible": 60,
  "N_HET": 1,
  "N_HOM_ALT": 1,
  "N_HOM_REF": 57,
  "N_FAIL": 0,
  "N_NO_COVERAGE": 0
}
```

---

## Special Cases

### AN=0 and AF=None

AN=0 means no eligible samples have coverage at this position. This happens when:

- All eligible samples are WES or panels and the position is outside capture regions
- The phenotype/sex/technology filter excludes all samples
- The chromosome is not in the database


!!! warning "AN=0 does not mean the variant is absent"
    AN=0 means AFQuery has no data to compute frequency. It is not evidence of rarity.

### Multi-allelic sites

Each ALT allele at a position is reported on its own line, and its counts refer to
that allele only. A sample carrying a different ALT allele at the same position
(for example `G/T` on the `G>A` line) is neither a carrier nor homozygous
reference for this allele, so it appears in none of `N_HET`, `N_HOM_ALT`,
`N_HOM_REF`, `N_FAIL` or `N_NO_COVERAGE`. The five counts add up to
`n_eligible` at biallelic sites; at a multi-allelic site they fall short by the
number of such samples. AC, AN and AF are not affected.

A sample with two different ALT alleles (`GT=1/2`) is counted as heterozygous on
both lines.

### Warnings

afquery emits a `AfqueryWarning` to stderr when a query may silently return fewer or no results. Common causes:

| Situation | Warning message |
|-----------|----------------|
| Chromosome not in database | `Chromosome 'chrXX' has no data in this database. Available: [...]` |
| Unknown phenotype code (include) | `Phenotype 'CODE' not in database — include will match 0 samples.` |
| Unknown phenotype code (exclude) | `Phenotype 'CODE' not in database — exclude has no effect.` |
| Unknown technology name (include) | `Technology 'NAME' not in database — include will match 0 samples.` |
| Unknown technology name (exclude) | `Technology 'NAME' not in database — exclude has no effect.` |
| Contradictory filters (e.g. include + exclude same code) | `Sample filter produces an empty eligible set — all queries will return AN=0.` |

Use `--no-warn` to suppress these warnings:

```bash
afquery query --db ./my_db/ --locus chr22:1000 --no-warn
afquery annotate --db ./my_db/ --input in.vcf --output out.vcf --no-warn
```

### AC=0 with High AN

AC=0 with a high AN (e.g., AN=4000) means the variant was **genuinely not observed** in a well-covered cohort. This is meaningful evidence that the variant is rare or absent in your population.


### N_FAIL > 0

When N_FAIL > 0, some eligible samples had a non-PASS filter at this position in their source VCF. A high N_FAIL relative to AN may indicate a problematic site (e.g., systematic sequencing artifacts). Consider filtering positions with more than 10% of failing (`N_FAIL / (AN / 2) > 0.1`).

---

## VCF Annotation Fields

When using `afquery annotate`, the following INFO fields are added to each variant:

| INFO field | Number | Description |
|-----------|--------|-------------|
| `AFQUERY_AC` | A (per ALT) | Allele count — one value per ALT allele |
| `AFQUERY_AN` | 1 (per site) | Allele number — shared across all ALT alleles |
| `AFQUERY_AF` | A (per ALT) | Allele frequency — one value per ALT allele |
| `AFQUERY_N_HET` | A (per ALT) | Heterozygous sample count per ALT allele |
| `AFQUERY_N_HOM_ALT` | A (per ALT) | Homozygous alt sample count per ALT allele |
| `AFQUERY_N_HOM_REF` | A (per ALT) | Homozygous ref sample count per ALT allele, excluding carriers of other ALT alleles at the position |
| `AFQUERY_N_FAIL` | 1 (per site) | Fail sample count — shared across all ALT alleles |
| `AFQUERY_N_NO_COVERAGE` | A (per ALT) | Eligible samples whose tech lacks coverage evidence at this position. Always `0` unless a coverage-evidence filter is active. See [Coverage Evidence](../advanced/coverage-evidence.md). |

!!! note "Multi-allelic sites"
    Number=A fields have one value per ALT allele (comma-separated for multi-allelic sites). Number=1 fields are shared across all ALT alleles at the same position.

These fields can be used directly in downstream filtering with `bcftools filter`:

```bash
# Keep variants rare in cohort with sufficient coverage
bcftools filter -i 'AFQUERY_AF < 0.001 && AFQUERY_AN >= 1000' annotated.vcf.gz
```

---

## Next Steps

- [Key Concepts](concepts.md) — how AC, AN, and AF are computed
- [Sample Filtering Guide](../guides/sample-filtering.md) — phenotype, sex, and technology filters
- [Annotate a VCF](../guides/annotate-vcf.md) — add AFQUERY_AF/AC/AN fields to a patient VCF
- [Clinical ACMG Use Cases](../use-cases/acmg-use-cases.md) — applying local AF to ACMG BS1/PM2 criteria

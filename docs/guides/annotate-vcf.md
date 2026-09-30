# Annotate a VCF

`afquery annotate` adds allele frequency information to an existing VCF file as INFO fields. It queries each variant in the VCF against the database and writes results inline.

---

## Basic Usage


```bash
afquery annotate \
  --db ./db/ \
  --input variants.vcf.gz \
  --output annotated.vcf.gz
```

---

## Added INFO Fields

| Field | Type | Number | Description |
|-------|------|--------|-------------|
| `AFQUERY_AC` | Integer | A (per ALT) | Allele count in eligible samples |
| `AFQUERY_AN` | Integer | 1 (per site) | Allele number (total alleles examined) |
| `AFQUERY_AF` | Float | A (per ALT) | Allele frequency (`AC / AN`) |
| `AFQUERY_N_HET` | Integer | A (per ALT) | Heterozygous sample count |
| `AFQUERY_N_HOM_ALT` | Integer | A (per ALT) | Homozygous alt sample count |
| `AFQUERY_N_HOM_REF` | Integer | A (per ALT) | Homozygous ref sample count, excluding carriers of other ALT alleles at the position (including alleles stored in the database but absent from the input record) |
| `AFQUERY_N_FAIL` | Integer | 1 (per site) | Eligible samples whose call had FILTER≠PASS. Excluded from AC, but still counted in AN. Mutually exclusive with N_HET/N_HOM_ALT/N_HOM_REF. |
| `AFQUERY_N_NO_COVERAGE` | Integer | A (per ALT) | Eligible samples whose tech lacks coverage evidence at this position. Excluded from `N_HOM_REF` to keep AC/AN conservative. Always `0` unless a coverage-evidence filter is active. See [Coverage Evidence](../advanced/coverage-evidence.md). |

!!! note "Multi-allelic sites"
    Number=A fields have one value per ALT allele (comma-separated for multi-allelic sites). Number=1 fields are shared across all ALT alleles at the same position.

---

## Sample Filtering

Annotate with a specific subgroup:

```bash
afquery annotate \
  --db ./db/ \
  --input variants.vcf \
  --output annotated.vcf \
  --phenotype E11.9 \
  --sex female \
  --tech wgs
```

The INFO fields will reflect AF computed only over the filtered sample set. This allows generating population-specific frequency tracks.

---

## Parallelism

Annotation runs in parallel across variants. By default, all available CPU cores are used:

```bash
afquery annotate \
  --db ./db/ \
  --input variants.vcf \
  --output annotated.vcf \
  --threads 8
```

For large VCFs (100K+ variants), set `--threads` to the number of available cores.

---


## Using Annotated Output

### BCFtools

Filter variants with high AF:

```bash
bcftools filter -i 'AFQUERY_AF > 0.01' annotated.vcf
```

Extract specific fields:

```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AFQUERY_AC\t%AFQUERY_AN\t%AFQUERY_AF\n' annotated.vcf
```

### Python (pysam / cyvcf2)

```python
import cyvcf2

vcf = cyvcf2.VCF("annotated.vcf")
for variant in vcf:
    ac = variant.INFO.get("AFQUERY_AC")
    an = variant.INFO.get("AFQUERY_AN")
    af = variant.INFO.get("AFQUERY_AF")
    print(f"{variant.CHROM}:{variant.POS} AC={ac} AN={an} AF={af}")
```

### R (VariantAnnotation)

```r
library(VariantAnnotation)
vcf <- readVcf("annotated.vcf")
info(vcf)$AFQUERY_AF
```

---

## Full Option Reference

See [CLI Reference → annotate](../reference/cli.md#annotate).

---

## Next Steps

- [Clinical Prioritization](../use-cases/clinical-prioritization.md) — filter annotated VCFs to retain only rare local variants
- [Sample Filtering](sample-filtering.md) — annotate using subgroup-specific AF (phenotype, sex, tech)
- [Understanding Output](../getting-started/understanding-output.md) — what AFQUERY_AF/AC/AN/N_FAIL fields mean

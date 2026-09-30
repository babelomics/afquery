# Query Allele Frequencies

`afquery query` retrieves allele frequencies from the database. Three query modes are available: point, region, and batch.

---

## Point Query

Query a single genomic position:

```bash
afquery query --db ./db/ --locus chr1:925952
```

Filter to a specific alt allele (useful at multi-allelic sites):

```bash
afquery query --db ./db/ --locus chr1:925952 --alt A
```

---

## Region Query

Query all variants in a genomic range:

```bash
afquery query --db ./db/ --region chr1:900000-1000000
```

The range is 1-based, inclusive on both ends.

---

### Python API — multi-chromosome regions

To query variants across multiple regions (including different chromosomes)
in a single call, use `query_region_multi`:

```python
from afquery import Database

db = Database("./db/")
regions = [
    ("chr1",  900000,   1000000),
    ("chr17", 41196311, 41277500),
]
results = db.query_region_multi(regions, phenotype=["E11.9"])
```

Results are returned in **genomic order** (chr1, chr2, …, chr22, chrX, chrY,
chrM). Overlapping regions are automatically deduplicated — each variant
appears at most once. Chromosome names are normalized, so `"1"` and `"chr1"`
are equivalent.

For querying specific variants across chromosomes, use `query_batch_multi`:

```python
variants = [
    ("chr1",  925952,  "G", "A"),
    ("chrX",  5000000, "A", "G"),
]
results = db.query_batch_multi(variants)
```

Results are returned in **input order** (by original index). Duplicate entries
are deduplicated per chromosome — if the same `(chrom, pos, ref, alt)` appears
more than once, only the first occurrence is included. Chromosome names are
normalized, so `"1"` and `"chr1"` are equivalent.

---

## Batch Query

Query multiple positions at once from a file:

```bash
afquery query --db ./db/ --from-file variants.tsv
```

The input file is a **headerless TSV** with columns `chrom pos [ref [alt]]` (ref and alt are optional):

```tsv
chr1	925952	G	A
chr1	1014541	C	T
chrX	5000000	A	G
```

Batch queries support variants across multiple chromosomes in a single file.

---

## Output Formats

### text (default)

Human-readable, one block per variant:

```
chr1:925952 G>A  AC=142  AN=2742  AF=0.0518  n_eligible=1371  N_HET=138  N_HOM_ALT=2  N_HOM_REF=1231  N_FAIL=0  N_NO_COVERAGE=0
```

### tsv

Tab-separated, one row per variant, suitable for downstream processing:

```bash
afquery query --db ./db/ --region chr1:900000-1000000 --format tsv
```

```
chrom	pos	ref	alt	AC	AN	AF	n_eligible	N_HET	N_HOM_ALT	N_HOM_REF	N_FAIL	N_NO_COVERAGE
chr1	925952	G	A	142	2742	0.051782	1371	138	2	1231	0	0
```


### json

JSON array, one object per variant:

```bash
afquery query --db ./db/ --locus chr1:925952 --format json
```

```json
[
  {
    "chrom": "chr1",
    "pos": 925952,
    "ref": "G",
    "alt": "A",
    "AC": 142,
    "AN": 2742,
    "AF": 0.05178,
    "n_eligible": 1371,
    "N_HET": 138,
    "N_HOM_ALT": 2,
    "N_HOM_REF": 1231,
    "N_FAIL": 0,
    "N_NO_COVERAGE": 0
  }
]
```

---

## Coverage-Evidence Filters (no_coverage)

By default AFQuery counts every BED-covered sample without a variant call as
hom-ref. With standard variant-only VCFs that assumption can be wrong: a missing
position may simply mean the sample was not sequenced deeply enough at that locus.
Three optional flags let you trade hom-ref aggressiveness for confidence. Samples
that fall below a threshold are reported in **N_NO_COVERAGE** instead of N_HOM_REF
(they remain in `eligible` and `AN`, like `N_FAIL`).

| Flag | Meaning |
|------|---------|
| `--min-pass K` | A partially-covered tech is valid for hom-ref at a position only if it has ≥K PASS carriers (het\|hom) of any ALT allele there. Otherwise its non-carrier samples move to `N_NO_COVERAGE`. |
| `--min-observed K` | Same as `--min-pass`, but counts any VCF entry (`het\|hom\|fail`). Useful when you want to include calls that failed FILTER as evidence the position was sequenced. |
| `--min-quality-evidence K` | Requires ≥K quality-passing carriers per partially-covered tech. Requires a database built with `--min-dp`, `--min-gq`, `--min-qual`, or `--min-covered`. |

`--min-pass` and `--min-observed` combine with AND (both must hold). Both
default to `0`, which disables the gate.

```bash
afquery query --db ./db/ --locus chr1:925952 --min-pass 1
afquery query --db ./db/ --region chr1:900000-1000000 --min-observed 2 --min-pass 1
```

At a biallelic site the genotype invariant becomes:
`N_HET + N_HOM_ALT + N_HOM_REF + N_FAIL + N_NO_COVERAGE = n_eligible`.
At a multi-allelic site, eligible samples that carry only another ALT allele at
the position fall in none of these categories for this allele, so the sum is
lower than `n_eligible` by exactly that number of samples.

Fully-covered samples (those whose tech was registered without a BED) are
never affected. Carrier samples (het/hom/fail) are never moved to
`N_NO_COVERAGE`. See [Coverage Evidence](../advanced/coverage-evidence.md)
for when to reach for each flag.

---

## Sample Filtering

All query modes support the same filter options:

```bash
afquery query \
  --db ./db/ \
  --locus chr1:925952 \
  --phenotype E11.9 \
  --sex female \
  --tech wgs
```

Filters compose with AND:

- `--phenotype E11.9 --sex female` = female samples with E11.9

Multiple values for the same filter compose with OR:

- `--phenotype E11.9 --phenotype I10` = samples with E11.9 OR I10

Exclude with `^` prefix:

- `--tech ^wes_v1` = all technologies except `wes_v1`

See [Sample Filtering](sample-filtering.md) for full syntax.

---

## Results When AN=0

If all samples are excluded by your filters, the result will have `AC=0`, `AN=0`, and `AF=None`. This is expected behavior — it means no samples in the selected subgroup were eligible at that position (e.g., all WES samples and the position is not in their capture regions).

---

## Comparing AF Across Subgroups

Run two queries and compare:

```bash
# Diabetic patients
afquery query --db ./db/ --locus chr1:925952 --phenotype E11.9 --format json

# Healthy controls (exclude diabetic)
afquery query --db ./db/ --locus chr1:925952 --phenotype ^E11.9 --format json
```

For systematic comparison across many variants, consider [Bulk Export](dump-export.md) with `--by-phenotype`, or see [Cohort Stratification](../use-cases/cohort-stratification.md) for a worked multi-group comparison.

---

## Full Option Reference

See [CLI Reference → query](../reference/cli.md#query).

---

## Next Steps

- [Sample Filtering](sample-filtering.md) — full filter syntax for phenotype, sex, and technology
- [Understanding Output](../getting-started/understanding-output.md) — field definitions and special cases (AN=0, N_FAIL)
- [Cohort Stratification](../use-cases/cohort-stratification.md) — comparing AF across multiple groups systematically

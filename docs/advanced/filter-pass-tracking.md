# FILTER=PASS Tracking

AFQuery keeps the calls that fail quality filters (`FILTER≠PASS`) instead of discarding them. A failed call tells you the position was looked at but the sample's genotype there is unknown, which is different both from "variant not present" and from "variant present".

---

## What counts as a failed call

A sample has a failed call at a position when its VCF holds, for any ALT allele there:

- a call carrying the ALT allele with `FILTER` other than `PASS` or `.`;
- a missing genotype (`./.`) at a site that failed a filter; or
- calls that carry more alleles than the sample has at the position, such as two ALT alleles on chrX in a male. These are detected when the VCF is ingested, which logs a warning (see [Genotypes inconsistent with ploidy](ploidy-and-sex-chroms.md#genotypes-inconsistent-with-ploidy)).

These samples are stored in the variant's `fail_bitmap` and reported as `N_FAIL`. A homozygous-reference call (`0/0`) is never stored, so at a failed site it would count as hom-ref; [`normalize_vcf.sh`](../getting-started/preprocessing.md) avoids this by masking every genotype at a failed site as missing.

---

## How failed calls are counted

A failed call is judged per **position**, not per allele: a sample with a failed call for any ALT allele at a position has no reliable genotype there. For every allele at that position, the sample:

- adds nothing to **AC**, even if another record at the position called it with `FILTER=PASS`;
- is left out of **AN**, which counts only the alleles of samples with a usable call;
- is counted in **N_FAIL**, which therefore has the same value on every line of a multi-allelic site;
- stays in **n_eligible**, so `N_HET + N_HOM_ALT + N_HOM_REF + N_FAIL + N_NO_COVERAGE` still accounts for every eligible sample (see [Multi-allelic sites](../getting-started/understanding-output.md#multi-allelic-sites) for the one exception).

Leaving failed calls out of AN keeps AF an estimate over the genotypes that were actually observed. At every position, AN equals the reference alleles plus the AC of all its ALT alleles.

```bash
afquery query --db ./db/ --locus chr1:925952
```

```
chr1:925952 G>A  AC=142  AN=2728  AF=0.0521  n_eligible=1371  N_HET=138  N_HOM_ALT=2  N_HOM_REF=1224  N_FAIL=7  N_NO_COVERAGE=0
```

The 7 samples with a failed call are part of `n_eligible` (138 + 2 + 1224 + 7 = 1371) but not of AN: 2 × (1371 − 7) = 2728.

### When every eligible sample failed

If all eligible samples have a failed call at a position, AN is 0 and AF is undefined: `NA` in text and TSV output, `null` in JSON and the Python API, an empty cell in `dump`, and `.` in `AFQUERY_AF`. `N_FAIL` still reports how many samples failed.

---

## Counting failed calls as reference

`--fail-as-ref` (`fail_as_ref=True` in the Python API) keeps samples with a failed call in AN, counting them as homozygous reference. It is available on `query`, `annotate` and `dump`, and changes only AN and AF: the failed samples still add nothing to AC and are still reported in `N_FAIL`.

```bash
afquery query --db ./db/ --locus chr1:925952 --fail-as-ref
```

```
chr1:925952 G>A  AC=142  AN=2742  AF=0.0518  n_eligible=1371  N_HET=138  N_HOM_ALT=2  N_HOM_REF=1224  N_FAIL=7  N_NO_COVERAGE=0
```

Because every failed sample is then treated as a non-carrier, the AF is never higher than the default one, which makes it a conservative lower estimate when a site has many failures. AN then no longer equals the reference alleles plus the AC of all ALT alleles.

---

## Finding the failed samples

`afquery variant-info` lists each sample with a failed call at the position as genotype `alt` with filter `FAIL`, along with its technology and phenotypes. This shows whether failures cluster in one technology or sample group. See [Variant Info](../guides/variant-info.md).

```bash
afquery variant-info --db ./db/ --locus chr1:925952
```

---

## VCF annotation

`afquery annotate` writes the same per-position values: `AFQUERY_AN` excludes failed calls and `AFQUERY_N_FAIL` counts each failed sample once per site. See [Understanding Output](../getting-started/understanding-output.md#vcf-annotation-fields) for the full set of `AFQUERY_*` fields.

---

## Next Steps

- [Understanding Output](../getting-started/understanding-output.md) — interpreting N_FAIL in query and annotate results
- [ACMG Criteria](../use-cases/acmg-use-cases.md) — using N_FAIL to assess site quality before applying BA1/PM2
- [Debugging Results](debugging-results.md) — diagnosing unexpectedly high N_FAIL values

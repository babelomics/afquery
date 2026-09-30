# Coverage Evidence

Standard variant-only VCFs do not record hom-ref calls. When a sample has no
entry at a position, AFQuery has to decide whether the sample is genuinely
homozygous reference or simply was not sequenced there. For **fully-covered
techs** — those registered without a BED capture file in the manifest, so
every position is assumed to be sequenced — the answer is unambiguous. For
**partially-covered techs** (whole-exome kits, gene panels), the BED proves a
position was *targeted* by the assay, not that *this* sample was sequenced
deeply enough to call a confident hom-ref.

`N_NO_COVERAGE` lets you label that uncertain subset instead of forcing it
into `N_HOM_REF`. The flags below decide *which* samples land there.

---

## What `N_NO_COVERAGE` represents

`N_NO_COVERAGE` counts eligible samples whose hom-ref status is not trusted
under the active criteria. At a biallelic site the genotype invariant becomes:

```
N_HET + N_HOM_ALT + N_HOM_REF + N_FAIL + N_NO_COVERAGE = n_eligible
```

At a multi-allelic site, eligible samples that carry only another ALT allele at
the position fall in none of these categories for this allele, so the sum is
lower than `n_eligible` by exactly that number of samples.

Samples in `N_NO_COVERAGE` remain in `eligible` and contribute to `AN` (just
like `N_FAIL`), so AC/AN/AF stay conservative — the field never inflates
allele frequencies. Two rules always hold:

- **Carriers are never reclassified.** A sample with a `het`, `hom`, or
  `fail` call at the position stays in its category. `N_NO_COVERAGE` only
  draws from non-carriers.
- **Fully-covered samples are never gated.** Every coverage flag is a
  *per-tech* decision evaluated only on partially-covered techs. Samples on
  fully-covered techs are always treated as hom-ref when they have no
  carrier call.

---

## Cohort-evidence gates at query time

These flags use only the carriers already present in your cohort to decide
whether each partially-covered tech has enough evidence to trust hom-ref at a
position. They run at query time, so no database rebuild is needed.

| Flag | Effect |
|------|--------|
| `--min-pass K`     | A partially-covered tech must have ≥K PASS carriers (`het ∪ hom`) at the position. Carriers of any ALT allele there count, since each shows the position was sequenced. If it falls short, all of its non-carrier samples move from `N_HOM_REF` to `N_NO_COVERAGE`. |
| `--min-observed K` | Same shape, but counts every recorded carrier (`het ∪ hom ∪ fail`). Useful when a non-PASS call still proves the position was sequenced. |

When both flags are >0, both must hold (AND). The default `0` disables the
gate.

!!! tip
    If your VCFs do not carry `FORMAT/DP` or `FORMAT/GQ`, these are the
    flags you want. They are the cheapest option and apply to any database.

### Worked example

The numbers below are illustrative; concrete values depend on your cohort.

Default query — every BED-covered non-carrier counts as hom-ref:

```bash
afquery query --db ./db/ --locus chr1:925952
```

```
chr1:925952 G>A  AC=142  AN=2742  AF=0.0518  n_eligible=1371  N_HET=138  N_HOM_ALT=2  N_HOM_REF=1231  N_FAIL=0  N_NO_COVERAGE=0
```

Now require at least one PASS carrier per partially-covered tech:

```bash
afquery query --db ./db/ --locus chr1:925952 --min-pass 1
```

```
chr1:925952 G>A  AC=142  AN=2742  AF=0.0518  n_eligible=1371  N_HET=138  N_HOM_ALT=2  N_HOM_REF=1108  N_FAIL=0  N_NO_COVERAGE=123
```

Samples on partially-covered techs that did not contribute a single PASS
carrier at this position have moved out of `N_HOM_REF` and into
`N_NO_COVERAGE`. `AC`, `AN`, and `AF` are unchanged: the samples are still
eligible, they just no longer count as confident hom-refs.

---

## Quality-aware filtering at database creation

If your VCFs carry `FORMAT/DP`, `FORMAT/GQ`, or you trust the `QUAL` column,
you can demand that carriers meet quality thresholds before they count as
evidence for hom-ref. These flags apply when you create the database, so the
coverage decision is baked in.

| Flag (`create-db`) | Effect |
|--------------------|--------|
| `--min-dp D`     | Minimum `FORMAT/DP` per carrier. |
| `--min-gq G`     | Minimum `FORMAT/GQ` per carrier. |
| `--min-qual Q`   | Minimum VCF `QUAL` per carrier. |
| `--min-covered K`| Per partially-covered tech, the position is "trusted" only if at least K of its carriers pass the quality thresholds. Non-carriers of failing positions are recorded as `N_NO_COVERAGE`. The gate is evaluated at build time for each ALT allele separately; samples carrying another allele at the position are never counted as `N_NO_COVERAGE`. |

A carrier counts as quality-passing only if **all** active thresholds hold
(unset thresholds are simply ignored). At least one of these flags must be
non-zero to enable quality-aware coverage filtering — without that, queries
fall back to the cohort-evidence gates above.

```bash
afquery create-db \
  --manifest samples.tsv \
  --output-dir ./db/ \
  --genome-build GRCh38 \
  --bed-dir ./beds/ \
  --min-dp 30 --min-gq 20 --min-covered 1
```

!!! note
    The chosen thresholds are recorded with the database and re-applied
    automatically when you grow it via `update-db --add-samples`. You do
    not re-pass them on each update.

Enabling quality-aware filtering requires creating (or re-creating) the
database; existing databases without quality data must be rebuilt.

### Tightening at query time — `--min-quality-evidence`

Once a database has been built with at least one of `--min-dp`,
`--min-gq`, `--min-qual`, or `--min-covered`, you can tighten the gate at
query time without rebuilding:

```bash
afquery query --db ./db/ --locus chr1:925952 --min-quality-evidence 5
```

`--min-quality-evidence K` requires each partially-covered tech to have ≥K
quality-passing carriers of any ALT allele at the position. Non-carriers of failing techs
(other than those already filtered at build time) move to `N_NO_COVERAGE`.

Running the flag against a database that was not built with quality data
exits with a clear error:

```
This database was not built with coverage quality data.
Re-create with --min-dp / --min-gq to use --min-quality-evidence.
```

---

## Choosing thresholds

Three concrete profiles, ordered from cheapest to strictest:

- **Pure-genotype cohorts** (no `FORMAT/DP` / `FORMAT/GQ` / reliable `QUAL`)
  Use `--min-pass 1` at query time. Or `--min-observed 1` if you want
  failed calls to also count as evidence the position was sequenced.
  No rebuild needed; conservative — positions where your cohort happens to
  have zero PASS calls flip to `N_NO_COVERAGE`.

- **Cohorts with `FORMAT/DP` and `FORMAT/GQ`**
  Build with `--min-dp 20 --min-gq 20 --min-covered 1`. Carriers with low
  confidence stop validating positions, and the decision is stored in the
  database — every query benefits without further flags.

- **High-stakes clinical interpretation**
  Layer `--min-quality-evidence 3` (or higher) on top of a quality-aware
  database to demand multiple independent quality-passing carriers per tech
  before trusting hom-ref.

---

## How the filters combine

`N_NO_COVERAGE` is the union of:

1. samples whose tech failed the build-time `--min-covered` gate;
2. samples whose tech failed `--min-pass` / `--min-observed` at query time;
3. samples whose tech failed `--min-quality-evidence`.

Carriers are never included; the same sample is never counted twice.

---

## Next Steps

- [Understanding Output](../getting-started/understanding-output.md) —
  field definitions for `N_HOM_REF`, `N_FAIL`, and `N_NO_COVERAGE`
- [FILTER=PASS Tracking](filter-pass-tracking.md) —
  the related `N_FAIL` field for failed-quality carrier calls
- [Technology Integration](../use-cases/technology-integration.md) —
  mixing whole-genome, whole-exome, and panel data in one cohort
- [Debugging Results](debugging-results.md) —
  diagnosing unexpected `N_NO_COVERAGE` or AN values

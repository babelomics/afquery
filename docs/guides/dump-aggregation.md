# Summing Dump Groups

A dump disaggregated by group (`--by-sex`, `--by-tech`, `--by-phenotype` or
`--all-groups`) reports the counts of every group in its own columns. Adding
up the counts of several groups gives the frequency of a variant in their
union, for example all the samples sequenced with any exome kit, or the
females of two technologies, without running a new query. This page explains
which columns can be added, how to derive frequencies from the sums and when
the result is exact.

A small worked example, with a database, its dump and a script that sums
groups, is available in
[`examples/dump_aggregation/`](https://github.com/babelomics/afquery/tree/master/examples/dump_aggregation).

---

## Columns that can be summed

| Column | Summable | Notes |
|---|---|---|
| `AC` | Yes | Copies of the ALT allele of the row |
| `AN` | Yes | Alleles with a usable call at the position |
| `AC_REF` | Yes | Reference alleles at the position |
| `n_eligible` | Yes | Samples that can be genotyped at the position |
| `N_HET`, `N_HOM_ALT`, `N_HEMI`, `N_HOM_REF`, `N_FAIL`, `N_NO_COVERAGE` | Yes | Samples by genotype |
| `AF` | No | Recompute it from the summed `AC` and `AN` |

`AN`, `AC_REF`, `n_eligible` and `N_FAIL` describe the position, so they take
the same value on every row of a multi-allelic site. Add them up across groups
within a row, never across the rows of a position, which would count the same
samples more than once.

Ploidy, capture coverage and failed calls are already taken into account in
the counts of each group, so the sums need no further correction. See
[Understanding Output](../getting-started/understanding-output.md) for the
definition of each column.

---

## Frequencies

For a selection `S` of disjoint groups:

```
AC(S)     = Σ AC_g
AN(S)     = Σ AN_g
AC_REF(S) = Σ AC_REF_g

AF(S)     = AC(S) / AN(S)          alternate allele frequency
REF_AF(S) = AC_REF(S) / AN(S)      reference allele frequency
```

The genotype counts are added in the same way. Both frequencies are undefined
when `AN(S)` is 0, which means that no sample in the selection has a usable
call at the position, and should be reported as unavailable rather than as 0.

The identity `AN = AC_REF + ΣAC`, with `AC` summed over every allele stored at
the position (including alleles that have no row in the dump), holds in every
group and in any selection of groups.

!!! warning "Common mistakes"
    - **Summing `AF`.** Frequencies are not additive. Sum `AC` and `AN`, then divide.
    - **Using `1 − AF` as the reference frequency.** At a multi-allelic
      position the other ALT alleles also take part of `AN`. With `AN` = 16, a
      deletion with `AC` = 4 and an insertion with `AC` = 2 at the same
      position, `AC_REF` is 10 and the reference frequency is 0.625, whereas
      `1 − AF` would give 0.75 and 0.875.

---

## When groups are disjoint

Sums are exact only when no sample belongs to two of the selected groups.

| Axis | Disjoint |
|---|---|
| Sex (`--by-sex`) | Always, since each sample has one sex |
| Technology (`--by-tech`) | Always, since each sample has one technology |
| Phenotype (`--by-phenotype`) | Only if no sample has more than one of the selected codes |
| Combinations (`--all-groups`, or several `--by-*` options) | When every axis is disjoint |

A sample with the codes `E11.9,I10` is counted in both the `E11.9` and the
`I10` groups, so adding those two groups counts it twice. This command lists
the samples that have more than one phenotype code:

```bash
afquery info --db ./db/ --samples --format tsv | awk -F'\t' 'NR > 1 && $5 ~ /,/'
```

---

## Producing a dump for summing

```bash
afquery dump --db ./db/ --all-groups --output groups.csv
```

Some options change which samples each group contains.

- **Base filters.** The axis of a group replaces the matching inclusion
  filter of the command. With `--sex male --by-sex` the `_female` columns
  still count the female samples, while `--tech wgs --by-sex` restricts both
  sex groups to WGS. Exclusions (`--phenotype ^CODE`, `--tech ^NAME`) apply to
  every group. With no base filters and disjoint groups that together include
  every sample, the base columns are the sum of all the groups. `--by-sex`,
  `--by-tech` and `--all-groups` always include every sample, whereas
  `--by-phenotype` leaves out the samples that have none of the requested
  codes. Groups that involve phenotypes are disjoint only under the condition
  described in the previous section.
- **Coverage gates.** `--min-pass`, `--min-observed` and
  `--min-quality-evidence` are evaluated on the samples of each group, so a
  group can fail a gate that the whole cohort passes. Under a gate
  `N_HOM_REF` and `N_NO_COVERAGE` stop adding up to the base row, while the
  other counts still do. A database built with `--min-covered` checks the
  evidence over the whole cohort, which keeps every sum exact.
- **`--fail-as-ref`** applies to every group in the same way and keeps every
  count additive.
- **Rows with `AC` = 0** in the base columns are written only with
  `--all-variants`. Without base filters, an allele missing from the dump has
  no carrier in any group.

Group columns are named `<FIELD>_<label>`. How labels are built is described
in [Bulk Export](dump-export.md#disaggregate-output).

---

## Limitations

- Sites are defined by position. An insertion, a deletion and a SNV starting
  at the same position share `AN` and `AC_REF`, but a deletion is not carried
  over to the positions it spans, and its carriers count as reference at a
  SNV inside the deleted bases.
- A missing genotype (`./.`) with `FILTER=PASS` is not stored, so the sample
  counts as reference. If those samples should be left out of `AN`, give the
  record a failing `FILTER` value before ingestion so that the call is stored
  as failed.

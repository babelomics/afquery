# Bulk Export (dump)

`afquery dump` exports allele frequency data to CSV. It supports filtering by region and disaggregating output by sex, technology, or phenotype group.

---

## Basic Usage

Export all variants to stdout:

```bash
afquery dump --db ./db/
```

Write to a file:

```bash
afquery dump --db ./db/ --output all_variants.csv
```

!!! note "AC > 0 filter"
    By default `dump` exports only variants with AC > 0. Variants at covered positions with no carriers are omitted. Use `--all-variants` to include AC=0 rows, or `afquery query --locus` to verify coverage at a specific position.

    ```bash
    afquery dump --db ./db/ --all-variants --output all_covered.csv
    ```

---

## Filter by Region

Export a single chromosome:

```bash
afquery dump --db ./db/ --chrom chr1 --output chr1.csv
```

Export a specific region:

```bash
afquery dump --db ./db/ --chrom chr1 --start 900000 --end 1000000 --output region.csv
```

Positions are 1-based, inclusive on both ends.

---

## Sample Filtering

Apply the same sample filters as `query`:

```bash
afquery dump --db ./db/ \
  --phenotype E11.9 \
  --sex female \
  --tech wgs \
  --output diabetic_female_wgs.csv
```

---

## Disaggregate Output

All three disaggregation modes work on the same principle: add stratified columns alongside the totals. They can be combined in a single command. Each group adds the same columns as the base block, suffixed with the group label; `n_eligible_<label>` counts the samples of that group that can be genotyped at the position, and is `0` where none can.

Base columns (always present):
```
chrom,pos,ref,alt,AC,AN,AF,n_eligible,N_HET,N_HOM_ALT,N_HOM_REF,N_FAIL,N_NO_COVERAGE
```

`n_eligible` is the number of samples that can be genotyped at the position:
they pass the sample filters, their technology covers the position, and they
carry at least one copy of the locus. It counts samples, not alleles, so on
chrX non-PAR a cohort of 2 males and 2 females has `n_eligible=4` and `AN=6`.
It is the same value [`query`](query.md) reports.

`N_NO_COVERAGE` is always emitted but is `0` unless a coverage-evidence
filter is active (see [Coverage Evidence](../advanced/coverage-evidence.md)).
`AF` and every `AF_<group>` column are left empty when the matching `AN` is 0,
i.e. when no eligible sample of that group has a usable call at the position.

=== "--by-sex"

    Add separate columns for male and female:

    ```bash
    afquery dump --db ./db/ --by-sex --output by_sex.csv
    ```

    Output columns:
    ```
    chrom,pos,ref,alt,AC,AN,AF,n_eligible,N_HET,N_HOM_ALT,N_HOM_REF,N_FAIL,N_NO_COVERAGE,AC_male,AN_male,AF_male,n_eligible_male,N_HET_male,N_HOM_ALT_male,N_HOM_REF_male,N_FAIL_male,N_NO_COVERAGE_male,AC_female,AN_female,AF_female,n_eligible_female,N_HET_female,N_HOM_ALT_female,N_HOM_REF_female,N_FAIL_female,N_NO_COVERAGE_female
    ```

=== "--by-tech"

    Add separate columns per sequencing technology:

    ```bash
    afquery dump --db ./db/ --by-tech --output by_tech.csv
    ```

    Output columns include `AC_wgs`, `AN_wgs`, `AF_wgs`, `n_eligible_wgs`, `N_HET_wgs`, `N_HOM_ALT_wgs`, `N_HOM_REF_wgs`, `N_FAIL_wgs`, `N_NO_COVERAGE_wgs`, `AC_wes_v1`, `AN_wes_v1`, etc. (one group of nine columns per registered technology).

=== "--by-phenotype"

    Add separate columns for specific phenotype groups:

    ```bash
    afquery dump --db ./db/ \
      --by-phenotype E11.9 \
      --by-phenotype I10 \
      --output by_phenotype.csv
    ```

    Output includes `AC_E11.9`, `AN_E11.9`, `AF_E11.9`, `n_eligible_E11.9`, `N_HET_E11.9`, `N_HOM_ALT_E11.9`, `N_HOM_REF_E11.9`, `N_FAIL_E11.9`, `N_NO_COVERAGE_E11.9`, `AC_I10`, etc.

=== "--all-groups"

    Disaggregate by all combinations of sex × technology × phenotype simultaneously:

    ```bash
    afquery dump --db ./db/ --all-groups --output all_groups.csv
    ```

    !!! warning
        The number of columns grows as a Cartesian product. Disaggregating by all combinations may produce a very wide CSV file.

---

## Parallelism

Use multiple threads for faster export:

```bash
afquery dump --db ./db/ --threads 16 --output all_variants.csv
```

---


## Full Option Reference

See [CLI Reference → dump](../reference/cli.md#dump).

---

## Next Steps

- [Cohort Stratification](../use-cases/cohort-stratification.md) — systematic cross-group AF comparison on specific loci
- [Sample Filtering](sample-filtering.md) — filter syntax shared across query, annotate, and dump
- [Performance Tuning](../advanced/performance.md) — thread tuning for large export jobs

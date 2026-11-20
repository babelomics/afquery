# Summing Dump Groups

This example builds a small database, dumps it with `--all-groups` and adds up
the counts of a selection of groups to obtain their allele and reference
frequencies. The rules behind it are explained in
[Summing Dump Groups](../../docs/guides/dump-aggregation.md).

## Usage

```bash
python create_example_data.py            # builds example_output/
python aggregate_groups.py example_output/dump_groups.csv
python aggregate_groups.py example_output/dump_groups.csv \
    female_WGS_E11 male_WGS_E11 female_WES_E11 male_WES_E11
```

`create_example_data.py` writes the VCFs, BED, manifest and database to
`example_output/` (or to the directory given as its argument), plus two dumps:
`dump_total.csv` without groups and `dump_groups.csv` with `--all-groups`.
`aggregate_groups.py` needs only the Python standard library.

## Cohort

There is one sample per sex × technology × phenotype cell, so the eight groups
are disjoint and together make up the whole cohort. Each VCF is already
normalized, with one ALT per record, so a sample carrying two ALT alleles at a
position has two records there.

| Sample | Sex | Tech | Phenotype | Genotypes |
|---|---|---|---|---|
| S01 | female | WGS | E11 | chr1:1000 A>AT 0/1; chr1:5000 G>A 0/1; chrX:5000000 G>A 0/1 |
| S02 | male | WGS | E11 | chr1:1000 AC>A 0/1; chr1:5000 G>A 1/1; chrX:5000000 G>A 1 |
| S03 | female | WGS | I10 | chr1:1000 A>AT 0/1 and AC>A 0/1; chr1:5000 G>T 0/1 |
| S04 | male | WGS | I10 | chr1:5000 G>A 0/1 and G>T 0/1; chr1:7000 TAGG>T 0/1 |
| S05 | female | WES | E11 | chr1:1000 AC>A 1/1; chr1:7002 G>C 0/1; chrX:5000000 G>A 1/1 |
| S06 | male | WES | E11 | chr1:7000 TAGG>T 0/1; chr1:7002 G>C 0/1 |
| S07 | female | WES | I10 | chr1:5000 G>T 0/1 |
| S08 | male | WES | I10 | chr1:9000 C>G 0/1 |

The WES capture covers every variant, so all eight samples are eligible
everywhere.

## Whole cohort

`dump_total.csv` holds the whole cohort. Summing all eight groups with
`aggregate_groups.py` gives the same counts.

```
chrom,pos,ref,alt,AC,AN,AF,AC_REF,n_eligible,N_HET,N_HOM_ALT,N_HEMI,N_HOM_REF,N_FAIL,N_NO_COVERAGE
chr1,1000,AC,A,4,16,0.25,10,8,2,1,0,4,0,0
chr1,1000,A,AT,2,16,0.125,10,8,2,0,0,4,0,0
chr1,5000,G,A,4,16,0.25,9,8,2,1,0,3,0,0
chr1,5000,G,T,3,16,0.1875,9,8,3,0,0,3,0,0
chr1,7000,TAGG,T,2,16,0.125,14,8,2,0,0,6,0,0
chr1,7002,G,C,2,16,0.125,14,8,2,0,0,6,0,0
chr1,9000,C,G,1,16,0.0625,15,8,1,0,0,7,0,0
chrX,5000000,G,A,4,12,0.3333333333333333,8,8,1,1,1,5,0,0
```

Each position illustrates a different case.

- **chr1:1000.** The deletion `AC>A` and the insertion `A>AT` have different
  REF but start at the same position, so they form one site and share `AN`
  and `AC_REF`. The reference frequency is 10/16 = 0.625, whereas `1 − AF`
  would give 0.75 and 0.875. S03 carries both alleles and is counted in
  `N_HET` on both rows. S01 carries only the insertion, so it appears in none
  of the genotype counts of the deletion row, where `N_HOM_REF` is 4 rather
  than 8 − 3.
- **chr1:5000.** Two SNVs, with S04 carrying both. `AN = AC_REF + ΣAC` gives
  16 = 9 + 4 + 3.
- **chr1:7000 and chr1:7002.** The deletion `TAGG>T` removes positions
  7001–7003, but the SNV at 7002 is a separate site with its own `AC_REF`.
  S04, which carries only the deletion, counts as reference at 7002, because
  afquery does not carry a deletion over to the positions it spans. S06 also
  carries the deletion but has its own `G>C` call at 7002.
- **chr1:9000.** A biallelic site, where `AC_REF = AN − AC`.
- **chrX:5000000.** Outside the PARs males are haploid, so `AN` is 12
  (4 females × 2 + 4 males × 1) rather than `2 × n_eligible`. S02 is
  hemizygous and counted in `N_HEMI`.

## A selection of groups

Summing the four E11 groups (`female_WGS_E11 male_WGS_E11 female_WES_E11 male_WES_E11`)
gives this result.

```
chrom,pos,ref,alt,AC,AN,AC_REF,n_eligible,N_HET,N_HOM_ALT,N_HEMI,N_HOM_REF,N_FAIL,N_NO_COVERAGE,AF,REF_AF
chr1,1000,AC,A,3,8,4,4,1,1,0,1,0,0,0.375,0.5
chr1,1000,A,AT,1,8,4,4,1,0,0,1,0,0,0.125,0.5
chr1,5000,G,A,3,8,5,4,1,1,0,2,0,0,0.375,0.625
chr1,5000,G,T,0,8,5,4,0,0,0,2,0,0,0.0,0.625
chr1,7000,TAGG,T,1,8,7,4,1,0,0,3,0,0,0.125,0.875
chr1,7002,G,C,2,8,6,4,2,0,0,2,0,0,0.25,0.75
chr1,9000,C,G,0,8,8,4,0,0,0,4,0,0,0.0,1.0
chrX,5000000,G,A,4,6,2,4,1,1,1,1,0,0,0.6666666666666666,0.3333333333333333
```

Every value comes from a single row of `dump_groups.csv`. `chr1:5000 G>T` has
no carrier among the E11 samples, yet its row still gives the `AN` and `AC_REF`
of the selection. `chr1:7000` and `chr1:7002` now differ in `AC_REF` (7 and 6)
because they are separate sites. In a dump without base filters, a row is
missing only when no sample in the cohort carries the allele, so it would add
nothing to any selection.

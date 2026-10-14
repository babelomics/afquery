# Glossary

## AC (Allele Count)

Number of copies of the alternate allele observed across eligible samples. A heterozygous carrier contributes 1; a homozygous alt carrier contributes 2. See [Key Concepts](../getting-started/concepts.md#allele-frequency-ac-an-af).

## AF (Allele Frequency)

Ratio of AC to AN (`AC / AN`). Represents the proportion of alleles that are the alternate allele among eligible samples. Returns `None` when AN=0. See [Key Concepts](../getting-started/concepts.md#allele-frequency-ac-an-af).

## AN (Allele Number)

Total number of alleles examined across eligible samples, leaving out samples with a failed call at the position ([N_FAIL](#n_fail)). For autosomes, AN = 2 × (eligible samples − N_FAIL). For sex chromosomes, AN depends on ploidy. See [Ploidy & Sex Chromosomes](../advanced/ploidy-and-sex-chroms.md).

## Bitmap

See [Roaring Bitmap](#roaring-bitmap).

## Bucket

A 1 Mbp (1,000,000 base pair) genomic interval used to partition variant data within each chromosome. Bucket 0 covers positions 0–999,999; bucket 1 covers 1,000,000–1,999,999; and so on. See [Data Model](data-model.md).

## Carrier

A sample that has at least one copy of the alternate allele at a given position (heterozygous or homozygous alt). Use `variant-info` to list carriers with their metadata. See [Variant Info](../guides/variant-info.md).

## Capture Index

An interval tree (stored as a `.pkl` file) built from a BED file that defines the genomic regions covered by a whole-exome sequencing (WES) technology. Used to determine which WES samples are eligible at a given position. WGS samples do not require a capture index. See [Key Concepts](../getting-started/concepts.md#capture-index-wes).

## Eligible Sample

A sample that passes all query filters (sex, phenotype, technology) and has coverage at the queried position. Only eligible samples contribute to AC and AN, and only if their call there did not fail. See [Sample Filtering](../guides/sample-filtering.md).

## N_FAIL

The count of eligible samples with a failed call (FILTER≠PASS) for any allele at a position. Their genotype there is unknown, so they are left out of AC and AN for every allele at the position. See [FILTER=PASS Tracking](../advanced/filter-pass-tracking.md).

## N_NO_COVERAGE

The count of eligible samples whose tech lacks coverage evidence at this position. Excluded from `N_HOM_REF`; they still count in AN. Always `0` unless a coverage-evidence filter (`--min-pass`, `--min-observed`, `--min-quality-evidence`, or build-time `--min-covered`) is active. See [Coverage Evidence](../advanced/coverage-evidence.md).

## Manifest

A TSV file that maps each sample to its VCF path, sex, sequencing technology, and phenotype codes. The manifest drives database creation and is parsed into `metadata.sqlite`. See [Manifest Format](../guides/manifest-format.md).

## Metadata Code

Legacy term. See [Phenotype Code](#phenotype-code).

## Parquet

A columnar storage format used by AFQuery to store variant bitmaps on disk. Each Parquet file covers one bucket within one chromosome and contains rows sorted by `(pos, alt)`. See [Data Model](data-model.md).

## Phenotype Code

An arbitrary string label assigned to a sample in the manifest (e.g., ICD-10 codes, HPO terms, project tags). Multiple codes per sample are supported. Phenotype codes enable cohort stratification at query time. See [Key Concepts](../getting-started/concepts.md#the-metadata-model).

## Ploidy

The number of allele copies at a given locus, which depends on chromosome type and sex. Autosomes are diploid (ploidy=2) for all samples. chrX non-PAR is haploid for males, diploid for females. See [Ploidy & Sex Chromosomes](../advanced/ploidy-and-sex-chroms.md).

## Roaring Bitmap

A compressed bitset data structure that efficiently stores sets of integers. AFQuery uses Roaring Bitmaps (via the `pyroaring` library) to represent which samples carry a given genotype. Typical compression: ~2 bytes per sample per variant. See [roaringbitmap.org](https://roaringbitmap.org/).

## Schema Version

The AFQuery database format version stored in `manifest.json` — `2.0`, or `3.0` when the database was built with coverage-quality filters. See [Data Model](data-model.md).

## Technology

The sequencing platform or assay type assigned to a sample (e.g., `wgs`, `wes_v1`, `panel_cardiac`). Determines capture index eligibility and enables technology-stratified queries. See [Technology Integration](../use-cases/technology-integration.md).

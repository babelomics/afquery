# Technology-Aware AN: Managing Mixed Sequencing Cohorts

## Why This Is Hard Without AFQuery

### Panel proliferation in hospitals

Many clinical genomics laboratories do not perform WGS or WES for all patients. Instead, they sequence targeted gene panels tailored to specific disease groups: a cardiac panel for arrhythmias, a retinal panel for inherited eye diseases, a metabolic panel for inborn errors of metabolism. Over time, a hospital's cohort may contain samples sequenced with dozens of distinct panels, each covering a different genomic footprint.

Even when a laboratory standardizes on WES, **multiple capture kit versions are used over time**. Kit versions from the same vendor may differ by hundreds to thousands of base pairs at exon boundaries. If kit version is ignored and all WES samples are treated identically, AN is miscounted at positions that differ between versions — inflating AF at some sites and deflating it at others.

### The manual calculation problem

To compute correct AN at any given position in a mixed-technology cohort, you must:

1. For each sample, determine whether the queried position falls within its capture BED file
2. Count only the samples with coverage as the AN denominator
3. Repeat this per-position, per-sample calculation across every variant in every query

Doing this correctly for a cohort with 10+ technologies and thousands of samples is prohibitively complex using standard VCF tools. **No general-purpose bioinformatics tool automates technology-aware AN computation.** Tools like bcftools, GATK, and VCFtools operate on static VCF files and have no concept of per-sample capture regions.

### Consequences of ignoring technology

When AN is computed naively (all samples included regardless of coverage):

| Scenario | Effect |
|----------|--------|
| Panel sample at position outside panel | AN inflated → AF underestimated |
| WES_v1 sample at position unique to WES_v2 | AN inflated → AF underestimated |
| WGS-only position with WES samples included | AN correct only if WES samples happen to have the position |

In a cohort with 60% WES and 40% WGS, a WGS-only position would have its AN inflated by up to 60% if WES samples are incorrectly counted — making a variant appear 1.6× rarer than it actually is.

## How AFQuery Solves It

When building the database, you provide one BED file per non-WGS technology:

```bash
afquery create-db \
  --manifest manifest.tsv \
  --bed-dir ./beds/ \
  --output-dir ./db/ \
  --genome-build GRCh38
```

AFQuery builds an **interval tree** per technology from each BED file. At query time:

- A WGS sample is always eligible (no BED constraint)
- A WES or panel sample is eligible at a position **only if** that position falls within its capture interval tree

This is computed automatically for every query — no manual BED intersection required.

## Setup Guide

### Manifest with multiple technologies

```tsv
sample_name	vcf_path	sex	tech_name	phenotype_codes
SAMP_001	vcfs/SAMP_001.vcf.gz	female	wgs	I42
SAMP_002	vcfs/SAMP_002.vcf.gz	male	wgs	I42
SAMP_003	vcfs/SAMP_003.vcf.gz	female	wes_v1	rare_disease
SAMP_004	vcfs/SAMP_004.vcf.gz	male	wes_v1	rare_disease
SAMP_005	vcfs/SAMP_005.vcf.gz	female	wes_v2	rare_disease
SAMP_006	vcfs/SAMP_006.vcf.gz	male	panel_cardio	I42
SAMP_007	vcfs/SAMP_007.vcf.gz	female	panel_retina	H35
```

### BED directory structure

```
beds/
├── wes_v1.bed          # Agilent SureSelect v5 (or whichever kit)
├── wes_v2.bed          # Agilent SureSelect v6
├── panel_cardio.bed    # Cardiac gene panel BED
└── panel_retina.bed    # Retinal gene panel BED
```

!!! note
    No BED file is needed for `wgs` — WGS samples are always eligible at every position. The `tech_name` field for WGS samples must be exactly `wgs` (case-insensitive).

### Naming conventions

- BED files must be named `<tech_name>.bed` and placed in the `--bed-dir` directory
- `tech_name` in the manifest must match the BED filename exactly (minus the `.bed` extension)
- Technology names are case-sensitive: `WES_v1` ≠ `wes_v1`

## Worked Examples

=== "Query at WES-covered position"

    A position within the WES_v1 capture region — WGS and WES_v1 samples are eligible; panel and WES_v2 samples are not:

    ```bash
    afquery query --db ./db/ --locus chr1:925952
    # AN reflects WGS + WES_v1 samples only
    ```

    ```
    chr1:925952 G>A  AC=8  AN=200  AF=0.0400  n_eligible=100  N_HET=8  N_HOM_ALT=0  N_HOM_REF=92  N_FAIL=0  N_NO_COVERAGE=0
    ```

    Compare WGS-only vs WES_v1-only to check for technology-specific bias:

    ```bash
    afquery query --db ./db/ --locus chr1:925952 --tech wgs
    afquery query --db ./db/ --locus chr1:925952 --tech wes_v1
    ```

=== "Query at WGS-only position"

    A position outside all WES capture regions — only WGS samples are eligible:

    ```bash
    afquery query --db ./db/ --locus chr1:12345678
    # AN reflects WGS samples only
    ```

    ```
    chr1:12345678 C>T  AC=2  AN=80  AF=0.0250  n_eligible=40  N_HET=2  N_HOM_ALT=0  N_HOM_REF=38  N_FAIL=0  N_NO_COVERAGE=0
    ```

    This is correct. WES samples are automatically excluded because their BED does not cover this position.

=== "Compare technologies"

    Scan multiple technologies to detect capture artifacts:

    ```bash
    for tech in wgs wes_v1 wes_v2 panel_cardio; do
      echo -n "$tech: "
      afquery query --db ./db/ --locus chr1:925952 --tech $tech --format tsv | \
        awk 'NR>1 {printf "AF=%.4f AN=%s\n", $7, $6}'
    done
    ```

    If AF differs substantially between WGS and WES at the same position, this may indicate:

    - A capture artifact (systematic over- or under-calling near exon boundaries)
    - A WES coverage bias at this specific site
    - True biological difference if WES/WGS cohorts have different phenotypic composition

## Detecting Technology Artifacts

Technology-stratified AF comparison is a useful quality check:

```bash
# Flag positions with large WGS-WES AF discrepancy
afquery query --db ./db/ --from-file variants.tsv --tech wgs   --format tsv > af_wgs.tsv
afquery query --db ./db/ --from-file variants.tsv --tech wes_v1 --format tsv > af_wes.tsv
```

Sites with AF ratio > 5× between WGS and WES may reflect systematic sequencing artifacts rather than true biological signal.

## Best Practices

- Always provide BED files for **all** non-WGS technologies, including each kit version separately
- Verify registration with `afquery check --db ./db/` — it validates BED file presence and interval tree integrity
- Use `--tech wgs` for cohort-wide analyses where capture bias is a concern
- Report `AN` alongside `AF`: AN < 100 indicates insufficient coverage to interpret AF reliably
- For variant interpretation, prefer positions where both WGS and WES samples are covered (high AN from both tech groups)

## Next Steps

- [Manifest Format](../guides/manifest-format.md) — tech_name and BED file setup
- [Create a Database](../guides/create-database.md) — `--bed-dir` option
- [Sample Filtering](../guides/sample-filtering.md) — `--tech` filter syntax
- [Cohort Stratification](cohort-stratification.md) — compare WGS vs. WES AF

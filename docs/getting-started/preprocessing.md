# VCF Preprocessing

AFQuery ingests single-sample, normalized VCF files. While AFQuery itself does not perform VCF normalization, the accuracy of your allele frequency estimates depends on the quality and consistency of your input VCFs. This page explains why normalization matters and provides a reference pipeline.

A ready-to-use normalization script is provided at [`resources/normalize_vcf.sh`](https://github.com/babelomics/afquery/blob/master/resources/normalize_vcf.sh). It requires **bcftools 1.20 or later**: older releases refuse to split multi-allelic records and remove duplicates in the same `bcftools norm` call, and the script exits with an error if it finds one.

!!! warning "VCFs normalized with the script from afquery 0.4.2 or earlier"
    Those versions of the script did not split multi-allelic records, and removed
    one of two different indels at the same position as a duplicate. Indels from
    multi-allelic sites could therefore be stored untrimmed (e.g. `CAAAAAAA>CAAAAAAAAAAAA`
    instead of `C>CAAAAA`), so the same variant appeared under two representations, or
    were dropped. Re-normalize those VCFs with the current script and rebuild the
    database to correct them.


---

## Why Normalize?

### Left-alignment and decomposition

Multi-allelic variants and complex indels can be represented in multiple equivalent ways in VCF format. Without normalization:

- The same deletion may appear as different ALT alleles depending on the caller
- Multi-allelic sites may not be decomposed into biallelic records
- Duplicate records can inflate AC

`bcftools norm -m -both` decomposes multi-allelic sites into biallelic records, and left-aligns and trims each allele against the reference genome, ensuring consistent representation across samples. `-d exact` then drops only records whose alleles are identical, so two different alleles at the same position are both kept.

### Ploidy correction for sex chromosomes

Male samples should have haploid genotype calls at chrX non-PAR regions (GT=`1`, not GT=`0/1` or `1/1`). Many variant callers output diploid genotypes for males at chrX by default. Without ploidy correction:

- Males at chrX are counted as diploid → AN is inflated by up to 50%
- AF estimates for X-linked variants are systematically underestimated

`bcftools +fixploidy` corrects genotype ploidy based on a sex file.

### Filtering and stripping

- **Non-PASS genotypes**: Masking them as missing ensures that low-quality calls do not inflate AC. AFQuery tracks these as N_FAIL.
- **Homozygous reference calls**: Removing ref/ref genotypes reduces file size and speeds ingestion; they contribute AC=0 and are not needed
- **INFO fields**: Stripping INFO reduces file size and speeds ingestion. Additionally, malformed or non-standard INFO fields produced by some variant callers can break downstream parsing; stripping them pre-emptively prevents these errors.
- **FORMAT fields**: Only `GT`, `DP`, and `GQ` are preserved. `GT` is the genotype required by AFQuery for all queries. `DP` and `GQ` are read by the [coverage-evidence](../advanced/coverage-evidence.md) quality flags (`afquery create-db --min-dp / --min-gq / --min-qual / --min-covered`). VCFs without `DP`/`GQ` are still valid; their carriers simply contribute no quality evidence to the cohort. All other FORMAT fields (`PL`, `AD`, etc.) are dropped to reduce file size.

---

## Next Steps

- [Manifest Format](../guides/manifest-format.md) — describe your cohort samples for ingestion
- [Create a Database](../guides/create-database.md) — build the AFQuery database from normalized VCFs
- [5-Min Quickstart](quickstart.md) — end-to-end tutorial from manifest to first query

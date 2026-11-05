# VCF Preprocessing

AFQuery ingests single-sample, normalized VCF files. While AFQuery itself does not perform VCF normalization, the accuracy of your allele frequency estimates depends on the quality and consistency of your input VCFs. This page explains why normalization matters and provides a reference pipeline.

A ready-to-use normalization script is provided at [`resources/normalize_vcf.sh`](https://github.com/babelomics/afquery/blob/master/resources/normalize_vcf.sh). It requires **bcftools 1.20 or later**: older releases refuse to split multi-allelic records and remove duplicates in the same `bcftools norm` call, and the script exits with an error if it finds one.

```bash
resources/normalize_vcf.sh <vcf_id> <vcf> <ref_genome> [threads]
# e.g.
resources/normalize_vcf.sh sample1 input.vcf.gz GRCh38.fa 8   # writes sample1_norm.vcf.gz
```

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

### Ploidy of sex chromosomes and chrM

No ploidy correction is needed, and the script makes none: AFQuery takes each sample's ploidy from its sex in the manifest and the position (see [Ploidy & Sex Chromosomes](../advanced/ploidy-and-sex-chroms.md)). Genotypes are read as follows, whether the caller wrote them haploid or diploid:

- **Male on chrX non-PAR or chrY**: `1` and `1/1` name a single ALT allele, so he is a hemizygous carrier. A heterozygous call (`0/1`) names two alleles where he has one copy; it is recorded as a failed call (`N_FAIL`), with a warning during ingestion (see [Genotypes inconsistent with ploidy](../advanced/ploidy-and-sex-chroms.md#genotypes-inconsistent-with-ploidy)). In a cohort called in diploid mode on the sex chromosomes this makes male `N_FAIL` visible at chrX positions where the caller saw both alleles, rather than hiding those genotypes.
- **chrM**: any call carrying the ALT allele (`1`, `0/1`, `1/1`) is a carrier. A `0/1` there is heteroplasmy, not an error.
- **chrY PARs**: not genotyped; PAR variants are counted at their chrX coordinates.

Do not run `bcftools +fixploidy` before ingestion. It keeps the first allele of each call when it lowers the ploidy (`0/1` becomes `0`), so a male's heterozygous call on chrX would be stored as hom-ref, and a heteroplasmic chrM call would be dropped (chrM is haploid in its rules for both sexes).

### Filtering and stripping

- **Non-PASS genotypes**: Masking every genotype at a failed site as missing ensures that low-quality calls do not inflate AC, and that no sample there is taken for hom-ref. AFQuery tracks these samples as N_FAIL and leaves them out of AN.
- **Homozygous reference calls**: Removing ref/ref genotypes reduces file size and speeds ingestion; they contribute AC=0 and are not needed
- **INFO fields**: Stripping INFO reduces file size and speeds ingestion. Additionally, malformed or non-standard INFO fields produced by some variant callers can break downstream parsing; stripping them pre-emptively prevents these errors.
- **FORMAT fields**: Only `GT`, `DP`, and `GQ` are preserved. `GT` is the genotype required by AFQuery for all queries. `DP` and `GQ` are read by the [coverage-evidence](../advanced/coverage-evidence.md) quality flags (`afquery create-db --min-dp / --min-gq / --min-qual / --min-covered`). VCFs without `DP`/`GQ` are still valid; their carriers simply contribute no quality evidence to the cohort. All other FORMAT fields (`PL`, `AD`, etc.) are dropped to reduce file size.

---

## Next Steps

- [Manifest Format](../guides/manifest-format.md) — describe your cohort samples for ingestion
- [Create a Database](../guides/create-database.md) — build the AFQuery database from normalized VCFs
- [5-Min Quickstart](quickstart.md) — end-to-end tutorial from manifest to first query

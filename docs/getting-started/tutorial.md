# Tutorial: End-to-End Walkthrough

This tutorial walks through every major AFQuery feature using a synthetic demo dataset. By the end, you will have built a database, queried variants, filtered by metadata, annotated a VCF, and exported results.

!!! tip "Prerequisites"
    Make sure AFQuery is installed. See [Installation](installation.md) if needed.

---

## 1. Generate Demo Data

AFQuery ships with a script that creates 10 synthetic VCFs, a manifest, and BED files for two WES technologies:

```bash
python examples/demo/create_demo_data.py
```

This creates `examples/demo/demo_output/` with:

- `vcfs/` — 10 single-sample VCFs (DEMO_001 through DEMO_010)
- `beds/` — capture BED files for `wes_v1` and `wes_v2`
- `manifest.tsv` — sample metadata

The demo cohort has 4 WGS samples, 3 wes_v1 samples, and 3 wes_v2 samples, with phenotype codes `E11.9`, `I10`, and `control`. It includes a set of variants at fixed positions so all query examples in this tutorial produce reproducible output.

---

## 2. Create the Database

```bash
afquery create-db \
  --manifest examples/demo/demo_output/manifest.tsv \
  --output-dir ./demo_db/ \
  --genome-build GRCh38 \
  --bed-dir examples/demo/demo_output/beds/
```

This ingests all VCFs, builds Roaring Bitmap Parquet files, and writes `manifest.json` and `metadata.sqlite`.

---

## 3. Inspect the Database

```bash
afquery info --db ./demo_db/
```

Example output:

```
Database:     ./demo_db/
Version:      1.0
Genome build: GRCh38   Schema: 2.0
Created: ...    Updated: N/A

Samples: 10 total
  By sex:       female=5  male=5
  By tech:      wes_v1=3  wes_v2=3  wgs=4
  By phenotype: E11.9=4  I10=5  control=4
```

The 10 samples span 3 technologies. Phenotype counts reflect samples tagged with each code; a sample may have multiple codes.

---

## 4. Query a Single Variant

Query the variant at `chr1:925952`:

```bash
afquery query --db ./demo_db/ --locus chr1:925952
```

Example output:

```
chr1:925952 G>A  AC=6  AN=20  AF=0.3000  n_eligible=10  N_HET=4  N_HOM_ALT=1  N_HOM_REF=5  N_FAIL=0  N_NO_COVERAGE=0
```

### Reading the output

Every query result includes the following fields:

| Field | Description |
|-------|-------------|
| **AC** | Allele count — alt allele copies in eligible samples |
| **AN** | Allele number — total alleles examined (adjusted for coverage and ploidy) |
| **AF** | Allele frequency — `AC / AN` |
| **n_eligible** | Eligible samples at this position — those passing the metadata filters *and* covered here. Ploidy is applied afterwards, to AN. Here, 10 = no filter applied and every sample is covered. |
| **N_HET** | Eligible samples heterozygous for the alt allele |
| **N_HOM_ALT** | Eligible samples homozygous for the alt allele |
| **N_HOM_REF** | Eligible samples homozygous reference |
| **N_FAIL** | Eligible samples with a non-ref allele called but `FILTER≠PASS` in the source VCF |

For `chr1:925952`: AN=20 = 10 diploid samples × 2 alleles; AC=6 = 4 het samples (1 copy each) + 1 hom-alt (2 copies).

**Accounting identity:** `n_eligible = N_HET + N_HOM_ALT + N_HOM_REF + N_FAIL + N_NO_COVERAGE` always holds. `N_NO_COVERAGE` is 0 here (it only becomes non-zero with a coverage-evidence filter), so 10 = 4+1+5+0+0 ✓.

### N_FAIL in practice

Some variants have samples whose call passed genotyping but failed a quality filter in the source VCF. Query the variant at `chr1:946000`:

```bash
afquery query --db ./demo_db/ --locus chr1:946000
```

Example output:

```
chr1:946000 T>C  AC=2  AN=16  AF=0.1250  n_eligible=10  N_HET=2  N_HOM_ALT=0  N_HOM_REF=6  N_FAIL=2  N_NO_COVERAGE=0
```

Here, 2 samples were flagged as `LowQual` in their source VCF at this position. Their genotype there is unknown, so they are left out of both AC and AN: `AN=2×(10−2)=16` and `AF=2/16=0.1250`. They remain eligible, so the accounting identity still holds: 10 = 2+0+6+2+0 ✓. See [FILTER=PASS Tracking](../advanced/filter-pass-tracking.md).

!!! tip
    N_FAIL > 0 is a signal to inspect source VCF quality at this site. A high N_FAIL relative to n_eligible may indicate a systematic sequencing artifact. See [Understanding Output](understanding-output.md) for thresholds.

### Region query

To find all variants in the WES capture region:

```bash
afquery query --db ./demo_db/ --region chr1:900000-950000
```

Example output:

```
chr1:925100 C>T  AC=4  AN=20  AF=0.2000  n_eligible=10  N_HET=4  N_HOM_ALT=0  N_HOM_REF=6  N_FAIL=0  N_NO_COVERAGE=0
chr1:925952 G>A  AC=6  AN=20  AF=0.3000  n_eligible=10  N_HET=4  N_HOM_ALT=1  N_HOM_REF=5  N_FAIL=0  N_NO_COVERAGE=0
chr1:946000 T>C  AC=2  AN=16  AF=0.1250  n_eligible=10  N_HET=2  N_HOM_ALT=0  N_HOM_REF=6  N_FAIL=2  N_NO_COVERAGE=0
```

---

## 5. Inspect Variant Carriers

After finding a variant of interest, use `variant-info` to see which specific samples carry it:

```bash
afquery variant-info --db ./demo_db/ --locus chr1:925952
```

Example output:

```
sample_id  sample_name  sex     tech    phenotypes       genotype  filter
---------  -----------  ------  ------  ---------------  --------  ------
0          DEMO_001     female  wgs     E11.9,I10        het       PASS
2          DEMO_003     male    wgs     E11.9            het       PASS
4          DEMO_005     female  wes_v1  E11.9,control    het       PASS
6          DEMO_007     male    wes_v2  control          het       PASS
8          DEMO_009     female  wgs     E11.9            hom       PASS
```

Each row is one carrier. The `genotype` column shows `het` (heterozygous), `hom` (homozygous alt), or `alt` (non-ref with FILTER≠PASS). The `filter` column indicates whether the call passed quality filters in the source VCF.

For machine-readable output, use `--format tsv`:

```bash
afquery variant-info --db ./demo_db/ --locus chr1:925952 --format tsv > carriers.tsv
```

See [Variant Info](../guides/variant-info.md) for full options including allele-specific queries and sample filtering.

---

## 6. Filter by Sex

Query only female samples:

```bash
afquery query --db ./demo_db/ --locus chr1:925952 --sex female
```

Example output:

```
chr1:925952 G>A  AC=3  AN=10  AF=0.3000  n_eligible=5  N_HET=3  N_HOM_ALT=0  N_HOM_REF=2  N_FAIL=0  N_NO_COVERAGE=0
```

Now query male samples:

```bash
afquery query --db ./demo_db/ --locus chr1:925952 --sex male
```

Example output:

```
chr1:925952 G>A  AC=3  AN=10  AF=0.3000  n_eligible=5  N_HET=1  N_HOM_ALT=1  N_HOM_REF=3  N_FAIL=0  N_NO_COVERAGE=0
```

AN drops from 20 to 10 in both cases because only 5 samples are eligible. The AF happens to be identical here (0.3000), but the genotype distributions differ: the male group has one homozygous-alt carrier while the female group has three heterozygous carriers.

---

## 7. Filter by Phenotype

Query samples tagged with `E11.9`:

```bash
afquery query --db ./demo_db/ --locus chr1:925952 --phenotype E11.9
```

Example output:

```
chr1:925952 G>A  AC=5  AN=8  AF=0.6250  n_eligible=4  N_HET=3  N_HOM_ALT=1  N_HOM_REF=0  N_FAIL=0  N_NO_COVERAGE=0
```

Exclude control samples:

```bash
afquery query --db ./demo_db/ --locus chr1:925952 --phenotype ^control
```

Example output:

```
chr1:925952 G>A  AC=6  AN=12  AF=0.5000  n_eligible=6  N_HET=4  N_HOM_ALT=1  N_HOM_REF=1  N_FAIL=0  N_NO_COVERAGE=0
```

The `^` prefix means "exclude". Excluding controls removes 4 samples, leaving 6. See [Sample Filtering](../guides/sample-filtering.md) for the full syntax.

---

## 8. Filter by Technology

Restrict to WGS samples only:

```bash
afquery query --db ./demo_db/ --locus chr1:925952 --tech wgs
```

Example output:

```
chr1:925952 G>A  AC=4  AN=8  AF=0.5000  n_eligible=4  N_HET=2  N_HOM_ALT=1  N_HOM_REF=1  N_FAIL=0  N_NO_COVERAGE=0
```

Now query each WES technology separately:

```bash
afquery query --db ./demo_db/ --locus chr1:925952 --tech wes_v1
```

```
chr1:925952 G>A  AC=1  AN=6  AF=0.1667  n_eligible=3  N_HET=1  N_HOM_ALT=0  N_HOM_REF=2  N_FAIL=0  N_NO_COVERAGE=0
```

```bash
afquery query --db ./demo_db/ --locus chr1:925952 --tech wes_v2
```

```
chr1:925952 G>A  AC=1  AN=6  AF=0.1667  n_eligible=3  N_HET=1  N_HOM_ALT=0  N_HOM_REF=2  N_FAIL=0  N_NO_COVERAGE=0
```

`chr1:925952` falls inside the `chr1:900000-950000` capture region shared by both WES kits. All 3 wes_v1 samples and all 3 wes_v2 samples are covered, so AN=6 (3 samples × 2 alleles) for both.

### The critical case: a position outside all WES capture regions

Now query a position that is not in any WES BED file:

```bash
afquery query --db ./demo_db/ --locus chr1:1399914 --tech wgs
```

```
chr1:1399914 G>T  AC=2  AN=8  AF=0.2500  n_eligible=4  N_HET=2  N_HOM_ALT=0  N_HOM_REF=2  N_FAIL=0  N_NO_COVERAGE=0
```

```bash
afquery query --db ./demo_db/ --locus chr1:1399914 --tech wes_v1
```

```
No variants found for the given filters.
```

`chr1:1399914` lies between the two chr1 WES capture regions (`900000–950000` and `1200000–1250000`). When you restrict to `--tech wes_v1`, AFQuery excludes all wes_v1 samples from AN because none of them have coverage at this position — AN=0, so no results are returned.

!!! note "No variants found ≠ rare"
    "No variants found" here means AFQuery has no coverage data for the requested technology at this position. The variant exists in the database (visible with `--tech wgs`), but frequency cannot be estimated for wes_v1. Always check that AN is sufficient before interpreting AF.

---

## 9. Combine Filters

All filter dimensions compose with AND:

```bash
afquery query \
  --db ./demo_db/ \
  --locus chr1:925952 \
  --sex female \
  --phenotype E11.9 \
  --tech wgs
```

Example output:

```
chr1:925952 G>A  AC=1  AN=2  AF=0.5000  n_eligible=1  N_HET=1  N_HOM_ALT=0  N_HOM_REF=0  N_FAIL=0  N_NO_COVERAGE=0
```

Only one sample meets all three criteria (DEMO_001: female, wgs, E11.9). With n_eligible=1, AN=2 (one diploid sample on an autosome) and AC=1 (heterozygous carrier).

---

## 10. Annotate a VCF

Use one of the demo VCFs as input:

```bash
afquery annotate \
  --db ./demo_db/ \
  --input examples/demo/demo_output/vcfs/DEMO_001.vcf.gz \
  --output ./annotated_demo.vcf
```

The output VCF gains `AFQUERY_AC`, `AFQUERY_AN`, `AFQUERY_AF`, `AFQUERY_N_HET`, `AFQUERY_N_HOM_ALT`, `AFQUERY_N_HOM_REF`, `AFQUERY_N_FAIL`, and `AFQUERY_N_NO_COVERAGE` INFO fields. Inspect the chr1 variants:

```bash
grep -v "^##" ./annotated_demo.vcf | grep "^chr1" | head -3
```

Example output (columns: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, DEMO_001):

```
chr1	925100	.	C	T	50	PASS	AFQUERY_AC=4;AFQUERY_AN=20;AFQUERY_AF=0.2;AFQUERY_N_HET=4;AFQUERY_N_HOM_ALT=0;AFQUERY_N_HOM_REF=6;AFQUERY_N_FAIL=0;AFQUERY_N_NO_COVERAGE=0	GT	0/1
chr1	925952	.	G	A	50	PASS	AFQUERY_AC=6;AFQUERY_AN=20;AFQUERY_AF=0.3;AFQUERY_N_HET=4;AFQUERY_N_HOM_ALT=1;AFQUERY_N_HOM_REF=5;AFQUERY_N_FAIL=0;AFQUERY_N_NO_COVERAGE=0	GT	0/1
chr1	946000	.	T	C	50	PASS	AFQUERY_AC=2;AFQUERY_AN=16;AFQUERY_AF=0.125;AFQUERY_N_HET=2;AFQUERY_N_HOM_ALT=0;AFQUERY_N_HOM_REF=6;AFQUERY_N_FAIL=2;AFQUERY_N_NO_COVERAGE=0	GT	0/1
```

`AFQUERY_N_FAIL=2` for `chr1:946000` matches the query output from Step 4. The annotated values reflect cohort-wide frequencies across all 10 samples, not just DEMO_001.

See [Annotate a VCF](../guides/annotate-vcf.md) for filtering and downstream usage.

---

## 11. Bulk Export with Dump

Export all variant frequencies to CSV:

```bash
afquery dump --db ./demo_db/ --output demo_dump.csv
```

Preview the first rows:

```bash
head -6 demo_dump.csv
```

Example output:

```
chrom,pos,ref,alt,AC,AN,AF,N_HET,N_HOM_ALT,N_HOM_REF,N_FAIL,N_NO_COVERAGE
chr1,925100,C,T,4,20,0.2,4,0,6,0,0
chr1,925952,G,A,6,20,0.3,4,1,5,0,0
chr1,946000,T,C,2,16,0.125,2,0,6,2,0
chr1,1225000,A,G,5,20,0.25,5,0,5,0,0
chr1,1399914,G,T,2,8,0.25,2,0,2,0,0
```

Note that `chr1:1399914` has AN=8 even in the full-cohort dump — technology-aware AN is applied automatically (WES samples are excluded for positions outside their capture regions).

Disaggregate by sex and technology:

```bash
afquery dump --db ./demo_db/ --output demo_dump_stratified.csv --by-sex --by-tech
```

This adds columns following the pattern `AC_{sex}_{tech}`, `AN_{sex}_{tech}`, `AF_{sex}_{tech}` (and genotype counts) for each sex-technology combination. The stratified dump is useful for case-control comparisons and population-specific frequency analysis.

---

## 12. Interpret Results with ACMG Criteria

With the annotated VCF or query results, you can apply ACMG criteria:

- **BA1**: Is AF > 5% in your cohort? → Variant is benign.
- **PM2_Supporting**: Is the variant absent (AC=0) with high AN? → Supporting pathogenic evidence.
- **PS4**: Is the variant enriched in cases vs. controls? Use `--phenotype` and `--phenotype ^disease_code` to compare.

For detailed guidance, see [ACMG Criteria (BA1/PM2/PS4)](../use-cases/acmg-use-cases.md).

---

## Next Steps

- [Key Concepts](concepts.md) — understand how bitmaps, Parquet, and metadata filtering work
- [Understanding Output](understanding-output.md) — field definitions and interpretation of N_FAIL, AN=0, and other special cases
- [Manifest Format](../guides/manifest-format.md) — prepare your own cohort manifest
- [Create a Database](../guides/create-database.md) — build a database from your real data
- [ACMG Criteria](../use-cases/acmg-use-cases.md) — apply local AF to variant classification

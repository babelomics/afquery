# Variant Info

`afquery variant-info` returns the list of samples carrying a specific variant, together with each sample's metadata: sex, sequencing technology, phenotype codes, genotype (het/hom), and FILTER status (PASS/FAIL).

This avoids the need to re-query raw VCF files when inspecting individual variant carriers.

---

## Basic usage

```bash
afquery variant-info --db ./db/ --locus chr1:925952
```

!!! tip
    `variant-info` is the natural next step after `query` — once you find a variant of interest, use it to see which specific samples carry it.

!!! note "Carrier counts may exceed the counts reported by `query`"
    `variant-info` lists every carrier present in the source VCFs. `query` additionally
    restricts to samples whose capture regions cover the position, so a WES sample with a
    call just outside its own BED appears here but is not counted in `AN`, `AC` or the
    genotype tallies. The two commands answer different questions — "who carries it?"
    versus "what is the frequency among samples that could have been called?" — so a
    modest difference is expected.

    A *large* gap is worth investigating: if `variant-info` shows many carriers from a
    technology that contributes nothing to `AN`, check that technology's capture index
    (see [Debugging Results](../advanced/debugging-results.md)).

By default all samples are queried and results are printed as an aligned text table:

```
sample_id  sample_name  sex     tech       phenotypes    genotype  filter
---------  -----------  ------  ---------  ------------  --------  ------
3          P003         male    WGS        E11.9,J45     het       PASS
17         P017         female  WES_kit_A  E11.9         hom       PASS
42         P042         male    WGS        I10           alt       FAIL
```

---

## Filtering to a specific allele

When multiple alleles exist at the same position, use `--ref` and `--alt` to restrict to one:

```bash
afquery variant-info --db ./db/ --locus chr17:41245466 --ref A --alt T
```

!!! note
    Without `--ref`/`--alt`, carriers for all alleles at the locus are returned and a warning is emitted if more than one allele is found. Specify both flags to disambiguate at multi-allelic sites.

---

## Sample filters

`variant-info` accepts the same sample filters as `query`:

```bash
# Only female carriers
afquery variant-info --db ./db/ --locus chr1:925952 --sex female

# Only carriers with phenotype E11.9
afquery variant-info --db ./db/ --locus chr1:925952 --phenotype E11.9

# Exclude phenotype I10
afquery variant-info --db ./db/ --locus chr1:925952 --phenotype ^I10

# Restrict to WGS samples
afquery variant-info --db ./db/ --locus chr1:925952 --tech WGS

# Combine filters
afquery variant-info --db ./db/ --locus chr1:925952 \
  --sex female --phenotype E11.9 --tech WGS,WES_kit_A
```

See [Sample Filtering](sample-filtering.md) for the full filter syntax.

---

## Output formats

### TSV

Machine-readable tab-separated output, suitable for downstream processing:

```bash
afquery variant-info --db ./db/ --locus chr1:925952 --format tsv > carriers.tsv
```

```
sample_id	sample_name	sex	tech	phenotypes	genotype	filter
3	P003	male	WGS	E11.9,J45	het	PASS
17	P017	female	WES_kit_A	E11.9	hom	PASS
42	P042	male	WGS	I10	alt	FAIL
```

### JSON

Structured output with variant metadata and a sample list:

```bash
afquery variant-info --db ./db/ --locus chr1:925952 --format json
```

```json
{
  "variant": {
    "chrom": "chr1",
    "pos": 925952,
    "ref": ".",
    "alt": "."
  },
  "samples": [
    {
      "sample_id": 3,
      "sample_name": "P003",
      "sex": "male",
      "tech": "WGS",
      "phenotypes": ["E11.9", "J45"],
      "genotype": "het",
      "filter": "PASS"
    },
    {
      "sample_id": 42,
      "sample_name": "P042",
      "sex": "male",
      "tech": "WGS",
      "phenotypes": ["I10"],
      "genotype": "alt",
      "filter": "FAIL"
    }
  ]
}
```

When `--ref` and `--alt` are specified, the `variant` block contains the actual alleles. Otherwise, `"."` is used as a placeholder.

---

## Genotype values

| Value | Meaning |
|---|---|
| `het` | Heterozygous carrier, FILTER=PASS |
| `hom` | Homozygous alt carrier, FILTER=PASS |
| `alt` | Non-ref carrier with FILTER≠PASS (ploidy unknown) |
| `no_coverage` | Sample's tech lacks coverage evidence at this position; not a carrier. Only appears when a coverage-evidence filter (`--min-pass`, `--min-observed`, `--min-quality-evidence`, or build-time `--min-covered`) is active. The FILTER column is empty (text/tsv) or `null` (JSON) — `PASS`/`FAIL` does not apply because there is no call. See [Coverage Evidence](../advanced/coverage-evidence.md). |

---

## All options

| Option | Default | Description |
|---|---|---|
| `--db` | required | Path to database directory |
| `--locus` | required | `CHROM:POS` (e.g. `chr1:925952`) |
| `--ref` | — | Filter to specific reference allele |
| `--alt` | — | Filter to specific alternate allele |
| `--phenotype` | all | Include phenotype (repeatable; `^CODE` excludes) |
| `--sex` | `both` | `male`, `female`, or `both` |
| `--tech` | all | Include technology (repeatable; `^NAME` excludes) |
| `--format` | `text` | `text`, `tsv`, or `json` |
| `--no-warn` | off | Suppress `AfqueryWarning` messages |

See also [CLI Reference → variant-info](../reference/cli.md#variant-info).

---

## Next Steps

- [Sample Filtering](sample-filtering.md) — full filter syntax for phenotype, sex, and technology
- [Understanding Output](../getting-started/understanding-output.md) — field definitions and special cases
- [FILTER=PASS Tracking](../advanced/filter-pass-tracking.md) — understanding FAIL genotypes
- [Python API → variant_info](../reference/python-api.md#variant_info) — programmatic access
- [ACMG Criteria](../use-cases/acmg-use-cases.md) — using carrier info for variant classification

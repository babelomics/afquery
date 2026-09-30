# Python API

AFQuery exposes a clean Python API through the `Database` class. All CLI functionality is available programmatically.

---

## Installation

```python
from afquery import Database
```

---

## Database

### Constructor

```python
db = Database(db_path: str)
```

Opens an AFQuery database at the given path. The manifest and sample metadata are loaded on initialization.

```python
db = Database("./my_db/")
```

---

### query

```python
db.query(
    chrom: str,
    pos: int,
    phenotype: list[str] | None = None,
    sex: str = "both",
    ref: str | None = None,
    alt: str | None = None,
    tech: list[str] | None = None,
) -> list[QueryResult]
```

Query allele frequencies at a single genomic position.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `chrom` | str | Chromosome name (e.g., `"chr1"`, `"chrX"`) |
| `pos` | int | 1-based genomic position |
| `phenotype` | list[str] \| None | Phenotype filter codes. Use `"^CODE"` prefix to exclude. |
| `sex` | str | `"both"` (default), `"male"`, or `"female"` |
| `ref` | str \| None | Filter to specific reference allele |
| `alt` | str \| None | Filter to specific alternate allele |
| `tech` | list[str] \| None | Technology filter. Use `"^TECH"` prefix to exclude. |

**Returns:** List of `QueryResult` objects, one per variant at the position (sorted by `(pos, alt)`).

**Example:**

```python
results = db.query("chr1", pos=925952)
for r in results:
    print(f"{r.variant.ref}>{r.variant.alt}  AC={r.AC}  AN={r.AN}  AF={r.AF:.4f}")

# With filters
results = db.query(
    chrom="chr1",
    pos=925952,
    phenotype=["E11.9"],
    sex="female",
    tech=["wgs"],
)
```

---

### query_region

```python
db.query_region(
    chrom: str,
    start: int,
    end: int,
    phenotype: list[str] | None = None,
    sex: str = "both",
    tech: list[str] | None = None,
) -> list[QueryResult]
```

Query allele frequencies for all variants in a genomic range.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `chrom` | str | Chromosome name |
| `start` | int | 1-based start position (inclusive) |
| `end` | int | 1-based end position (inclusive) |
| `phenotype` | list[str] \| None | Phenotype filter |
| `sex` | str | Sex filter |
| `tech` | list[str] \| None | Technology filter |

**Example:**

```python
results = db.query_region("chr1", start=900000, end=1000000)
print(f"Found {len(results)} variants")
```

---

### query_region_multi

```python
db.query_region_multi(
    regions: list[tuple[str, int, int]],
    phenotype: list[str] | None = None,
    sex: str = "both",
    tech: list[str] | None = None,
) -> list[QueryResult]
```

Query allele frequencies across multiple genomic regions, which may span
different chromosomes. Overlapping regions are deduplicated automatically.
Chromosome names are normalized (`"1"` and `"chr1"` are equivalent).

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `regions` | list[tuple[str, int, int]] | List of `(chrom, start, end)` tuples, 1-based inclusive |
| `phenotype` | list[str] \| None | Phenotype filter |
| `sex` | str | Sex filter |
| `tech` | list[str] \| None | Technology filter |

**Returns:** List of `QueryResult` objects sorted in genomic order
(chr1, chr2, …, chr22, chrX, chrY, chrM).

**Example:**

```python
# Gene panel spanning multiple chromosomes
regions = [
    ("chr1",  925000,   1000000),
    ("chr17", 41196311, 41277500),  # BRCA1
    ("chr13", 32315086, 32400266),  # BRCA2
]
results = db.query_region_multi(regions, phenotype=["C50"])
for r in results:
    print(f"{r.variant.chrom}:{r.variant.pos}  AF={r.AF:.4f}")
```

---

### query_batch

```python
db.query_batch(
    chrom: str,
    variants: list[tuple[int, str, str]],
    phenotype: list[str] | None = None,
    sex: str = "both",
    tech: list[str] | None = None,
) -> list[QueryResult]
```

Query allele frequencies for a list of specific variants.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `chrom` | str | Chromosome name |
| `variants` | list[tuple[int, str, str]] | List of `(pos, ref, alt)` tuples |
| `phenotype` | list[str] \| None | Phenotype filter |
| `sex` | str | Sex filter |
| `tech` | list[str] \| None | Technology filter |

**Example:**

```python
variants = [(925952, "G", "A"), (1014541, "C", "T"), (1020172, "A", "G")]
results = db.query_batch("chr1", variants=variants)
```

---

### query_batch_multi

```python
db.query_batch_multi(
    variants: list[tuple[str, int, str, str]],
    phenotype: list[str] | None = None,
    sex: str = "both",
    tech: list[str] | None = None,
) -> list[QueryResult]
```

Query allele frequencies for a list of specific variants across multiple
chromosomes. Chromosome names are normalized (`"1"` and `"chr1"` are
equivalent). Duplicate input entries are deduplicated per chromosome — if the
same `(chrom, pos, ref, alt)` appears more than once, only the first
occurrence is included.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `variants` | list[tuple[str, int, str, str]] | List of `(chrom, pos, ref, alt)` tuples |
| `phenotype` | list[str] \| None | Phenotype filter |
| `sex` | str | Sex filter |
| `tech` | list[str] \| None | Technology filter |

**Returns:** List of `QueryResult` objects in input order (by original index).
Variants not found in the database are omitted.

**Example:**

```python
variants = [
    ("chr1",  925952,   "G", "A"),
    ("chrX",  5000000,  "A", "G"),
    ("chr17", 41223094, "C", "T"),
]
results = db.query_batch_multi(variants, phenotype=["E11.9"])
```

---

### annotate_vcf

```python
db.annotate_vcf(
    input_vcf: str,
    output_vcf: str,
    phenotype: list[str] | None = None,
    sex: str = "both",
    n_workers: int | None = None,
    tech: list[str] | None = None,
) -> dict
```

Annotate a VCF file with allele frequency INFO fields.

**Returns:** Stats dict:
```python
{
    "n_variants": int,    # total variants in input VCF
    "n_annotated": int,   # variants with at least one allele found in DB
    "n_uncovered": int,   # variants with no allele found in DB
}
```

---

### dump

```python
db.dump(
    output: str | None = None,
    phenotype: list[str] | None = None,
    sex: str = "both",
    tech: list[str] | None = None,
    by_sex: bool = False,
    by_tech: bool = False,
    by_phenotype: list[str] | None = None,
    all_groups: bool = False,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    n_workers: int | None = None,
    include_ac_zero: bool = False,
) -> dict
```

Export allele frequency data to CSV. If `output` is None, writes to stdout.
`include_ac_zero=True` includes positions with `AC=0` (equivalent to `--all-variants` in the CLI).

---

### variant_info

```python
db.variant_info(
    chrom: str,
    pos: int,
    ref: str | None = None,
    alt: str | None = None,
    phenotype: list[str] | None = None,
    sex: str = "both",
    tech: list[str] | None = None,
) -> list[SampleCarrier]
```

Return all samples carrying the variant at the given position, with their metadata.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `chrom` | str | Chromosome name (e.g., `"chr1"`, `"chrX"`) |
| `pos` | int | 1-based genomic position |
| `ref` | str \| None | Filter to specific reference allele |
| `alt` | str \| None | Filter to specific alternate allele |
| `phenotype` | list[str] \| None | Phenotype filter codes. Use `"^CODE"` prefix to exclude. |
| `sex` | str | `"both"` (default), `"male"`, or `"female"` |
| `tech` | list[str] \| None | Technology filter. Use `"^TECH"` prefix to exclude. |

**Returns:** List of `SampleCarrier` objects sorted by `sample_id`. Empty list if no eligible carrier exists.

**Example:**

```python
carriers = db.variant_info("chr1", pos=925952)
for c in carriers:
    print(f"{c.sample_name}  {c.genotype}  {c.tech_name}  {c.phenotypes}")

# With allele and sample filters
carriers = db.variant_info(
    chrom="chr17", pos=41245466, ref="A", alt="T",
    phenotype=["E11.9"], sex="female", tech=["WGS"],
)
```

The `variant_info` function is also available at the top level for one-off use:

```python
from afquery import variant_info

carriers = variant_info("./db/", "chr1", 925952, ref="G", alt="A")
```

---

### add_samples

```python
db.add_samples(
    manifest_path: str,
    threads: int = 8,
    tmp_dir: str | None = None,
    bed_dir: str | None = None,
    genome_build: str | None = None,
) -> dict
```

Add samples from a manifest TSV. Returns a stats dict.

---

### remove_samples

```python
db.remove_samples(sample_names: list[str]) -> dict
```

Remove samples by name. Returns a stats dict with `n_removed`.

---

### compact

```python
db.compact() -> dict
```

Compact the database to reclaim space from removed samples.

---

### info

```python
db.info() -> dict
```

Return database metadata as a dict.

---

### list_samples

```python
db.list_samples() -> list[dict]
```

Return a list of all samples with their metadata (name, sex, tech, phenotypes).

---

### check

```python
db.check() -> list
```

Validate database integrity. Returns `list[CheckResult]`. Each item has `.severity` (`"error"`, `"warning"`, or `"info"`) and `.message` (str). An empty list means the database is healthy.

---

### changelog

```python
db.changelog(limit: int | None = None) -> list[dict]
```

Return changelog history. Each entry is a dict:
```python
{
    "event_id": int,
    "event_type": str,        # "preprocess", "add_samples", "remove_samples", "compact"
    "event_time": str,        # ISO datetime string
    "sample_names": list[str] | None,
    "notes": str | None,
}
```

---

### set_db_version

```python
db.set_db_version(version: str) -> None
```

Set the database version label.

---

### get_all_phenotypes

```python
db.get_all_phenotypes() -> list[str]
```

Return all distinct phenotype codes present in the database.

---

## QueryResult

```python
@dataclass
class QueryResult:
    variant: VariantKey      # chrom, pos, ref, alt
    AC: int                  # Allele count
    AN: int                  # Allele number
    AF: float | None         # Allele frequency (None if AN == 0)
    n_samples_eligible: int  # Number of eligible samples at this position
    N_HET: int               # Heterozygous count
    N_HOM_ALT: int           # Homozygous alt count
    N_HOM_REF: int           # Homozygous ref count
    N_FAIL: int              # Eligible samples whose call had FILTER≠PASS (excluded from AC, kept in AN)
    N_NO_COVERAGE: int       # Eligible samples whose tech lacks evidence (excluded from N_HOM_REF)
```

At a biallelic site the genotype categories partition the eligible samples:
`N_HET + N_HOM_ALT + N_HOM_REF + N_FAIL + N_NO_COVERAGE == n_samples_eligible`.
At a multi-allelic site, eligible samples that carry only another ALT allele at
the position fall in none of these categories for this allele, so the sum is
lower than `n_samples_eligible` by exactly that number of samples.
See [Coverage Evidence](../advanced/coverage-evidence.md) for details on
`N_NO_COVERAGE`.

### VariantKey

```python
@dataclass
class VariantKey:
    chrom: str   # Canonical form: 'chr1', 'chrX', etc.
    pos: int     # 1-based
    ref: str
    alt: str
```

---

## SampleCarrier

```python
@dataclass
class SampleCarrier:
    sample_id: int        # 0-based internal ID
    sample_name: str      # Name from manifest
    sex: str              # 'male' | 'female'
    tech_name: str        # Sequencing technology
    phenotypes: list[str] # Sorted phenotype codes
    genotype: str            # 'het' | 'hom' | 'alt' | 'no_coverage'
    filter_pass: bool | None # True=PASS, False=FILTER≠PASS, None=no call (no_coverage)
```

Returned by `Database.variant_info()`. Each instance represents one sample carrying the queried variant.

| Field | Type | Description |
|-------|------|-------------|
| `sample_id` | int | 0-based internal sample ID |
| `sample_name` | str | Sample name from manifest |
| `sex` | str | `"male"` or `"female"` |
| `tech_name` | str | Sequencing technology name |
| `phenotypes` | list[str] | Sorted list of phenotype codes |
| `genotype` | str | `"het"` (heterozygous, PASS), `"hom"` (homozygous alt, PASS), `"alt"` (non-ref, FILTER≠PASS), or `"no_coverage"` (WES sample whose tech lacks evidence — see [Coverage Evidence](../advanced/coverage-evidence.md)) |
| `filter_pass` | bool \| None | `True` if FILTER=PASS, `False` if FILTER≠PASS, `None` when `genotype == "no_coverage"` (the sample has no call at this position) |

---

## SampleFilter

```python
@dataclass
class SampleFilter:
    phenotype_include: list[str] = []     # Empty = all samples
    phenotype_exclude: list[str] = []
    tech_include: list[str] = []          # Empty = all samples
    tech_exclude: list[str] = []
    sex: str = "both"                     # 'male' | 'female' | 'both'
    min_pass: int = 0                     # partially-covered tech needs ≥K PASS carriers
    min_observed: int = 0                 # partially-covered tech needs ≥K any-VCF entries
    min_quality_evidence: int = 0         # partially-covered tech needs ≥K quality_pass carriers (DB built with --min-dp/--min-gq)

    @staticmethod
    def parse(
        phenotype_tokens: list[str],
        tech_tokens: list[str],
        sex: str = "both",
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
    ) -> SampleFilter
```

`min_pass`, `min_observed`, and `min_quality_evidence` control how
non-carrier WES samples are classified at query time. See
[Coverage Evidence](../advanced/coverage-evidence.md) for the model and the
build-time companions (`--min-dp`, `--min-gq`, `--min-qual`, `--min-covered`).

`SampleFilter.parse` handles the `^` prefix exclusion syntax:

```python
from afquery.models import SampleFilter

sf = SampleFilter.parse(
    phenotype_tokens=["E11.9", "^I10"],
    tech_tokens=["wgs"],
    sex="female",
)
# sf.phenotype_include = ["E11.9"]
# sf.phenotype_exclude = ["I10"]
# sf.tech_include = ["wgs"]
# sf.sex = "female"
```

---

## Full Example

```python
from afquery import Database

db = Database("./my_db/")

# Point query
results = db.query("chr1", pos=925952, phenotype=["E11.9"], sex="female")
for r in results:
    print(f"AF={r.AF:.4f}  AC={r.AC}/{r.AN}  HET={r.N_HET}  HOM={r.N_HOM_ALT}")

# Region query
region_results = db.query_region("chrX", start=154931044, end=155270560)
print(f"Found {len(region_results)} variants in PAR2")

# Batch query
variants = [(925952, "G", "A"), (1014541, "C", "T")]
batch_results = db.query_batch("chr1", variants=variants)

# Database info
meta = db.info()
print(f"Samples: {meta['n_samples']}, Build: {meta['genome_build']}")

# List all phenotypes
phenotypes = db.get_all_phenotypes()
print("Phenotypes:", phenotypes)

# Multi-chromosome batch query
variants = [("chr1", 925952, "G", "A"), ("chrX", 5000000, "A", "G")]
batch_results = db.query_batch_multi(variants)

# Multi-region query
regions = [("chr1", 900000, 1000000), ("chrX", 5000000, 6000000)]
region_results = db.query_region_multi(regions)
print(f"Found {len(region_results)} variants across {len(regions)} regions")

# Carrier lookup
carriers = db.variant_info("chr1", pos=925952)
for c in carriers:
    print(f"{c.sample_name}  {c.genotype}  {c.tech_name}  PASS={c.filter_pass}")
```

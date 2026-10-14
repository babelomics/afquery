from dataclasses import dataclass, field


class AfqueryWarning(UserWarning):
    """Warning emitted when a query may silently return fewer or no results."""


@dataclass
class SampleFilter:
    """Filtro unificado de muestras. Todos los campos vacíos = sin restricción."""
    phenotype_include: list[str] = field(default_factory=list)  # [] = todas
    phenotype_exclude: list[str] = field(default_factory=list)
    tech_include: list[str] = field(default_factory=list)       # [] = todas
    tech_exclude: list[str] = field(default_factory=list)
    sex: str = "both"  # 'male' | 'female' | 'both'
    min_pass: int = 0      # WES tech needs ≥K PASS carriers (het|hom) at position
    min_observed: int = 0  # WES tech needs ≥K any-VCF entries (het|hom|fail) at position
    min_quality_evidence: int = 0  # WES tech needs ≥K quality_pass carriers (Phase 2 DB only)
    fail_as_ref: bool = False  # count samples with a failed call as hom-ref in AN

    @staticmethod
    def parse(
        phenotype_tokens: list[str],
        tech_tokens: list[str],
        sex: str = "both",
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
        fail_as_ref: bool = False,
    ) -> "SampleFilter":
        """Parsea tokens con prefijo ^ (exclusión) estilo bcftools."""
        return SampleFilter(
            phenotype_include=[t for t in phenotype_tokens if not t.startswith("^")],
            phenotype_exclude=[t[1:] for t in phenotype_tokens if t.startswith("^")],
            tech_include=[t for t in tech_tokens if not t.startswith("^")],
            tech_exclude=[t[1:] for t in tech_tokens if t.startswith("^")],
            sex=sex,
            min_pass=min_pass,
            min_observed=min_observed,
            min_quality_evidence=min_quality_evidence,
            fail_as_ref=fail_as_ref,
        )


@dataclass
class Sample:
    sample_id: int      # 0-indexed; used as bit position in Roaring Bitmaps
    sample_name: str
    sex: str            # 'male' | 'female'
    tech_id: int


@dataclass
class Technology:
    tech_id: int
    tech_name: str
    bed_path: str | None  # None means WGS (always covered)


@dataclass
class VariantKey:
    chrom: str   # canonical form: 'chr1', 'chrX', etc.
    pos: int     # 1-based
    ref: str
    alt: str


@dataclass
class QueryParams:
    chrom: str
    pos: int               # 1-based
    filter: SampleFilter = field(default_factory=SampleFilter)
    ref: str | None = None   # optional filter to specific ref allele
    alt: str | None = None   # optional filter to specific alt allele


@dataclass
class QueryResult:
    variant: VariantKey
    AC: int
    AN: int
    AF: float | None       # None if AN == 0
    n_samples_eligible: int
    N_HET: int = 0
    N_HOM_ALT: int = 0
    N_HOM_REF: int = 0
    N_FAIL: int = 0
    N_NO_COVERAGE: int = 0


@dataclass
class SampleCarrier:
    """A sample carrying a given variant, with its associated metadata."""
    sample_id: int
    sample_name: str
    sex: str            # 'male' | 'female'
    tech_name: str
    phenotypes: list[str]
    genotype: str              # 'het' | 'hom' | 'alt' (FILTER≠PASS) | 'no_coverage' (uncertain hom-ref)
    filter_pass: bool | None   # True = PASS, False = FILTER≠PASS, None = no call (no_coverage)

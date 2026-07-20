# PAR1 and PAR2 for GRCh37 and GRCh38
PAR = {
    "GRCh38": {
        "chrX": [(10_001, 2_781_479), (155_701_383, 156_030_895)],
        "chrY": [(10_001, 2_781_479), (56_887_903, 57_217_415)],
    },
    "GRCh37": {
        "chrX": [(60_001, 2_699_520), (154_931_044, 155_260_560)],
        "chrY": [(10_001, 2_649_520), (59_034_050, 59_363_566)],
    },
}

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
SEX_CHROMS = ["chrX", "chrY"]
MITO_CHROM = "chrM"
# Mitochondrion aliases seen in the wild (Ensembl 'MT', UCSC 'chrM', hybrid 'chrMT').
_MITO_ALIASES = frozenset(["m", "mt", "chrm", "chrmt"])
# Chromosomes whose canonical spelling is uppercase after the 'chr' prefix.
_UPPERCASE_SUFFIXES = {"x": "X", "y": "Y"}
ALL_CHROMS = AUTOSOMES + SEX_CHROMS + [MITO_CHROM]
CHROM_ORDER: dict[str, int] = {c: i for i, c in enumerate(ALL_CHROMS)}
VALID_GENOME_BUILDS = frozenset(["GRCh37", "GRCh38"])
VALID_SEX = frozenset(["male", "female"])


def normalize_chrom(chrom: str) -> str:
    """Normalize chromosome name to 'chr'-prefixed canonical form.

    '1' -> 'chr1', 'chr1' -> 'chr1', 'CHR1' -> 'chr1'
    'X' -> 'chrX', 'x' -> 'chrX', 'MT' -> 'chrM', 'chrMT' -> 'chrM'

    Only the 'chr' prefix and the known chromosome suffixes are case-normalized;
    anything else keeps its body byte-for-byte ('GL000209.1' -> 'chrGL000209.1',
    'chr1_KI270706v1_random' unchanged), so unplaced and alt contigs stay outside
    ALL_CHROMS and are still dropped by the build stage rather than guessed at.

    Idempotent on every member of ALL_CHROMS: the output is what gets persisted as
    the Parquet partition name, so a name that changed under normalization would
    orphan the partitions of already-built databases.
    """
    chrom = chrom.strip()
    if chrom.lower() in _MITO_ALIASES:
        return MITO_CHROM
    body = chrom[3:] if chrom[:3].lower() == "chr" else chrom
    if body.lower() in _UPPERCASE_SUFFIXES:
        return "chr" + _UPPERCASE_SUFFIXES[body.lower()]
    return "chr" + body


def is_autosome(chrom: str) -> bool:
    return normalize_chrom(chrom) in AUTOSOMES


def is_sex_chrom(chrom: str) -> bool:
    return normalize_chrom(chrom) in SEX_CHROMS


def is_mito(chrom: str) -> bool:
    return normalize_chrom(chrom) == MITO_CHROM

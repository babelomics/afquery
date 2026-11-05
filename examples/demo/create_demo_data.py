#!/usr/bin/env python3
"""Generate synthetic demo data for the AFQuery tutorial.

Creates:
- 10 single-sample VCFs (~100 variants each) in demo_output/vcfs/
- A manifest.tsv in demo_output/
- BED files for 2 WES technologies in demo_output/beds/

All data is synthetic and not based on real patient genomes.
"""

import gzip
import os
import random
from pathlib import Path

SEED = 42
OUTPUT_DIR = Path(__file__).resolve().parent / "demo_output"

# --- Sample definitions -------------------------------------------------------

SAMPLES = [
    # (name, sex, tech, phenotype_codes)
    ("DEMO_001", "female", "wgs",      "E11.9,I10"),
    ("DEMO_002", "male",   "wgs",      "E11.9"),
    ("DEMO_003", "female", "wgs",      "I10"),
    ("DEMO_004", "male",   "wgs",      "control"),
    ("DEMO_005", "female", "wes_v1",   "E11.9,I10"),
    ("DEMO_006", "male",   "wes_v1",   "control"),
    ("DEMO_007", "female", "wes_v1",   "I10"),
    ("DEMO_008", "male",   "wes_v2",   "E11.9"),
    ("DEMO_009", "female", "wes_v2",   "control"),
    ("DEMO_010", "male",   "wes_v2",   "I10,control"),
]

# --- Variant pool -------------------------------------------------------------
# Shared pool of positions; each sample gets a random subset with random GTs.

CHROMS_AUTOSOME = ["chr1", "chr2", "chr7", "chr12", "chr15"]
BASES = ["A", "C", "G", "T"]

# WES capture regions (shared by both kits with slight differences)
WES_V1_REGIONS = [
    # (chrom, start, end)
    ("chr1",  900000,  950000),
    ("chr1",  1200000, 1250000),
    ("chr2",  1000000, 1080000),
    ("chr7",  5000000, 5100000),
    ("chr12", 2000000, 2050000),
    ("chr15", 4000000, 4100000),
]

WES_V2_REGIONS = WES_V1_REGIONS + [
    # V2 adds extra regions
    ("chr1",  1800000, 1850000),
    ("chr7",  7000000, 7050000),
    ("chr12", 3000000, 3050000),
]

# --- Pinned variants -----------------------------------------------------------
# These variants are injected into every sample with fixed genotypes and filters.
# They exist at tutorial-specific positions so query examples are reproducible.
#
# Format: (chrom, pos, ref, alt, gt_by_sample_name)
# gt_by_sample_name values:
#   str          → GT string, FILTER=PASS implied
#   (str, str)   → (GT, FILTER) — used for the N_FAIL teaching variant

PINNED_VARIANTS = [
    # chr1:925952 — inside wes_v1 + wes_v2 capture region (chr1:900000-950000)
    # Primary teaching variant used throughout steps 4–8 of the tutorial.
    # Expected (no filter): AC=6, AN=20, AF=0.3000, N_HET=4, N_HOM_ALT=1, N_HOM_REF=5, N_FAIL=0
    ("chr1", 925952, "G", "A", {
        "DEMO_001": "0/1",   # female, wgs,    E11.9+I10
        "DEMO_002": "1/1",   # male,   wgs,    E11.9
        "DEMO_003": "0/1",   # female, wgs,    I10
        "DEMO_004": "0/0",   # male,   wgs,    control
        "DEMO_005": "0/1",   # female, wes_v1, E11.9+I10
        "DEMO_006": "0/0",   # male,   wes_v1, control
        "DEMO_007": "0/0",   # female, wes_v1, I10
        "DEMO_008": "0/1",   # male,   wes_v2, E11.9
        "DEMO_009": "0/0",   # female, wes_v2, control
        "DEMO_010": "0/0",   # male,   wes_v2, I10+control
    }),
    # chr1:925100 — second variant in the same WES region, adds a second row to
    # region queries on chr1:900000-950000.
    # Expected (no filter): AC=4, AN=20, AF=0.2000, N_HET=4, N_HOM_ALT=0, N_HOM_REF=6, N_FAIL=0
    ("chr1", 925100, "C", "T", {
        "DEMO_001": "0/1",
        "DEMO_002": "0/0",
        "DEMO_003": "0/0",
        "DEMO_004": "0/1",
        "DEMO_005": "0/0",
        "DEMO_006": "0/1",
        "DEMO_007": "0/0",
        "DEMO_008": "0/0",
        "DEMO_009": "0/1",
        "DEMO_010": "0/0",
    }),
    # chr1:946000 — inside wes_v1 + wes_v2 capture region (chr1:900000-950000)
    # N_FAIL teaching variant: DEMO_002 and DEMO_004 carry the alt allele with
    # a LowQual filter → they contribute to N_FAIL, not AC/AN.
    # Expected (no filter): AC=2, AN=16, AF=0.1250, N_HET=2, N_HOM_ALT=0, N_HOM_REF=6, N_FAIL=2
    ("chr1", 946000, "T", "C", {
        "DEMO_001": ("0/1", "PASS"),
        "DEMO_002": ("0/1", "LowQual"),
        "DEMO_003": ("0/0", "PASS"),
        "DEMO_004": ("0/1", "LowQual"),
        "DEMO_005": ("0/0", "PASS"),
        "DEMO_006": ("0/0", "PASS"),
        "DEMO_007": ("0/1", "PASS"),
        "DEMO_008": ("0/0", "PASS"),
        "DEMO_009": ("0/0", "PASS"),
        "DEMO_010": ("0/0", "PASS"),
    }),
    # chr1:1399914 — OUTSIDE all WES capture regions (between chr1:950000–1200000)
    # Used in step 7 to show that WES queries return no results at non-captured sites.
    # Expected (--tech wgs): AC=2, AN=8, AF=0.2500, N_HET=2, N_HOM_REF=2, N_FAIL=0
    # Expected (--tech wes_v1): "No variants found" (AN=0 → engine returns [])
    ("chr1", 1399914, "G", "T", {
        "DEMO_001": "0/1",
        "DEMO_002": "0/0",
        "DEMO_003": "0/1",
        "DEMO_004": "0/0",
        "DEMO_005": "0/1",
        "DEMO_006": "0/0",
        "DEMO_007": "0/0",
        "DEMO_008": "0/1",
        "DEMO_009": "0/0",
        "DEMO_010": "0/0",
    }),
    # chr1:1225000 — inside wes_v1 + wes_v2 capture region (chr1:1200000-1250000)
    # Adds coverage in a second WES region for completeness.
    ("chr1", 1225000, "A", "G", {
        "DEMO_001": "0/1",
        "DEMO_002": "0/0",
        "DEMO_003": "0/0",
        "DEMO_004": "0/1",
        "DEMO_005": "0/1",
        "DEMO_006": "0/0",
        "DEMO_007": "0/1",
        "DEMO_008": "0/0",
        "DEMO_009": "0/0",
        "DEMO_010": "0/1",
    }),
]


def _alt_for(ref: str) -> str:
    """Pick a random ALT allele different from REF."""
    choices = [b for b in BASES if b != ref]
    return random.choice(choices)


def _generate_variant_pool(n_per_chrom: int = 20) -> list[tuple[str, int, str, str]]:
    """Generate a pool of (chrom, pos, ref, alt) variants."""
    pool = []

    for chrom in CHROMS_AUTOSOME:
        for _ in range(n_per_chrom):
            pos = random.randint(900000, 10000000)
            ref = random.choice(BASES)
            alt = _alt_for(ref)
            pool.append((chrom, pos, ref, alt))

    # Add chrX, chrY, chrM variants
    for _ in range(10):
        pos = random.randint(5000000, 100000000)
        ref = random.choice(BASES)
        pool.append(("chrX", pos, ref, _alt_for(ref)))

    for _ in range(5):
        pos = random.randint(2700000, 10000000)
        # Stay outside the chrY PAR1 of GRCh38 (ends at 2,781,479): the chrY PARs
        # are not genotyped. Shifting instead of redrawing keeps the random
        # sequence, and so every other variant and genotype, unchanged.
        if pos <= 2781479:
            pos += 100000
        ref = random.choice(BASES)
        pool.append(("chrY", pos, ref, _alt_for(ref)))

    for _ in range(5):
        pos = random.randint(100, 16500)
        ref = random.choice(BASES)
        pool.append(("chrM", pos, ref, _alt_for(ref)))

    # Sort by chrom, pos
    chrom_order = {c: i for i, c in enumerate(
        CHROMS_AUTOSOME + ["chrX", "chrY", "chrM"]
    )}
    pool.sort(key=lambda v: (chrom_order.get(v[0], 99), v[1]))

    return pool


def _write_vcf(
    path: Path,
    sample_name: str,
    sex: str,
    variants: list[tuple[str, int, str, str]],
    pinned: list | None = None,
) -> None:
    """Write a minimal single-sample VCF.

    Random variants get random GTs and a 5% chance of LowQual filter.
    Pinned variants use fixed GTs and filters from PINNED_VARIANTS.
    """
    chrom_order = {c: i for i, c in enumerate(
        CHROMS_AUTOSOME + ["chrX", "chrY", "chrM"]
    )}

    header = (
        "##fileformat=VCFv4.2\n"
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##FILTER=<ID=LowQual,Description="Low quality">\n'
    )

    # Add contig headers
    contigs = sorted({v[0] for v in variants})
    if pinned:
        contigs = sorted(set(contigs) | {v[0] for v in pinned})
    for c in contigs:
        header += f"##contig=<ID={c}>\n"

    header += f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"

    lines = []
    for chrom, pos, ref, alt in variants:
        # Random genotype
        r = random.random()
        if r < 0.60:
            gt = "0/0"  # hom ref (most common)
        elif r < 0.85:
            gt = "0/1"  # het
        elif r < 0.95:
            gt = "1/1"  # hom alt
        else:
            continue  # skip (not called)

        # Males are haploid on chrX outside the PARs (all chrX pool positions
        # are) and on chrY; females have no chrY. A haploid genotype is written
        # with one allele, as a sex-aware caller writes it.
        if chrom == "chrY" and sex == "female":
            continue
        if chrom in ("chrX", "chrY") and sex == "male":
            gt = "1" if "1" in gt else "0"

        if chrom == "chrM":
            gt = "1" if "1" in gt else "0"

        # Small chance of LowQual filter
        filt = "LowQual" if random.random() < 0.05 else "PASS"

        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{filt}\t.\tGT\t{gt}\n")

    # Inject pinned variants with controlled GTs and filters
    if pinned:
        for chrom, pos, ref, alt, gt_map in pinned:
            entry = gt_map.get(sample_name)
            if entry is None:
                continue
            if isinstance(entry, tuple):
                gt, filt = entry
            else:
                gt, filt = entry, "PASS"
            lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{filt}\t.\tGT\t{gt}\n")

    # Sort all lines by (chrom, pos) so the VCF is in genomic order
    lines.sort(key=lambda l: (chrom_order.get(l.split("\t")[0], 99), int(l.split("\t")[1])))

    with gzip.open(path, "wt") as f:
        f.write(header)
        f.writelines(lines)


def _write_bed(path: Path, regions: list[tuple[str, int, int]]) -> None:
    """Write a BED file."""
    with open(path, "w") as f:
        for chrom, start, end in sorted(regions):
            f.write(f"{chrom}\t{start}\t{end}\n")


def main() -> None:
    random.seed(SEED)

    vcf_dir = OUTPUT_DIR / "vcfs"
    bed_dir = OUTPUT_DIR / "beds"
    vcf_dir.mkdir(parents=True, exist_ok=True)
    bed_dir.mkdir(parents=True, exist_ok=True)

    # Generate variant pool, excluding any positions already pinned
    pool = _generate_variant_pool(n_per_chrom=20)
    pinned_positions = {(chrom, pos) for chrom, pos, *_ in PINNED_VARIANTS}
    pool = [v for v in pool if (v[0], v[1]) not in pinned_positions]

    # Write VCFs
    for name, sex, tech, pheno in SAMPLES:
        # Each sample gets ~80% of the pool (random subset)
        sample_variants = [v for v in pool if random.random() < 0.80]
        vcf_path = vcf_dir / f"{name}.vcf.gz"
        _write_vcf(vcf_path, name, sex, sample_variants, pinned=PINNED_VARIANTS)

    # Write BED files
    _write_bed(bed_dir / "wes_v1.bed", WES_V1_REGIONS)
    _write_bed(bed_dir / "wes_v2.bed", WES_V2_REGIONS)

    # Write manifest
    manifest_path = OUTPUT_DIR / "manifest.tsv"
    with open(manifest_path, "w") as f:
        f.write("sample_name\tvcf_path\tsex\ttech_name\tphenotype_codes\n")
        for name, sex, tech, pheno in SAMPLES:
            vcf_abs = str((vcf_dir / f"{name}.vcf.gz").resolve())
            f.write(f"{name}\t{vcf_abs}\t{sex}\t{tech}\t{pheno}\n")

    print(f"Demo data created in: {OUTPUT_DIR}")
    print(f"  VCFs:     {vcf_dir} ({len(SAMPLES)} files)")
    print(f"  BED files: {bed_dir} (wes_v1.bed, wes_v2.bed)")
    print(f"  Manifest: {manifest_path}")
    print()
    print("Next steps:")
    print(f"  afquery create-db --manifest {manifest_path} --output-dir ./demo_db/ "
          f"--genome-build GRCh38 --bed-dir {bed_dir}")


if __name__ == "__main__":
    main()

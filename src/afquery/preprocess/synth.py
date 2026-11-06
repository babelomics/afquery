import random
from pathlib import Path

_NUCLEOTIDES = ("A", "C", "G", "T")


def generate_synthetic_manifest(
    output_dir: Path,
    n_samples: int = 10_000,
    n_variants_per_chrom: int = 100_000,
    chroms: tuple[str, ...] = ("chr1", "chr2", "chr3"),
    seed: int = 42,
) -> Path:
    """
    Generate synthetic VCF files and a manifest.tsv for benchmark testing.

    All samples share the same variant positions (high bitmap overlap).
    Genotype distribution: 95% 0/0, 3% 0/1, 2% 1/1 (sparse / rare-variant model).
    Sex alternates male/female. Tech is WGS (no BED file required).

    Genotypes are diploid, so ``chroms`` should be autosomes: on chrX outside
    the PARs or on chrY, a male's 0/1 is inconsistent with his single copy and
    would be ingested as a failed call.

    Returns path to manifest.tsv.
    """
    output_dir = Path(output_dir)
    vcf_dir = output_dir / "vcfs"
    vcf_dir.mkdir(parents=True, exist_ok=True)

    rng = random.Random(seed)

    # Generate fixed variant positions for each chrom (shared across all samples)
    chrom_variants: dict[str, list[tuple[int, str, str]]] = {}
    for chrom in chroms:
        raw: list[tuple[int, str, str]] = []
        for _ in range(n_variants_per_chrom):
            pos = rng.randint(1, 250_000_000)
            ref = rng.choice(_NUCLEOTIDES)
            alt = rng.choice([n for n in _NUCLEOTIDES if n != ref])
            raw.append((pos, ref, alt))
        # Deduplicate by position (keep first occurrence after sort)
        seen: set[int] = set()
        unique: list[tuple[int, str, str]] = []
        for pos, ref, alt in sorted(raw):
            if pos not in seen:
                seen.add(pos)
                unique.append((pos, ref, alt))
        chrom_variants[chrom] = unique

    manifest_lines = ["sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes"]

    for i in range(n_samples):
        sample_name = f"synth_{i:06d}"
        sex = "male" if i % 2 == 0 else "female"
        vcf_path = vcf_dir / f"{sample_name}.vcf"

        _write_synthetic_vcf(vcf_path, sample_name, chrom_variants, rng)

        manifest_lines.append(
            f"{sample_name}\t{sex}\tWGS\t{vcf_path.resolve()}\tE11.9"
        )

    manifest_path = output_dir / "manifest.tsv"
    manifest_path.write_text("\n".join(manifest_lines) + "\n")
    return manifest_path


def _write_synthetic_vcf(
    vcf_path: Path,
    sample_name: str,
    chrom_variants: dict[str, list[tuple[int, str, str]]],
    rng: random.Random,
) -> None:
    """Write a sparse single-sample VCF with realistic genotype distribution."""
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        for chrom in sorted(chrom_variants):
            f.write(f"##contig=<ID={chrom}>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write(
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n"
        )
        for chrom in sorted(chrom_variants):
            for pos, ref, alt in chrom_variants[chrom]:
                r = rng.random()
                if r < 0.95:
                    continue       # 95% hom-ref — skip (no row)
                elif r < 0.98:
                    gt = "0/1"    # 3% het
                else:
                    gt = "1/1"   # 2% hom-alt
                f.write(
                    f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n"
                )

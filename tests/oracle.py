"""An independent reimplementation of what a query is supposed to return.

This module deliberately imports nothing from afquery. It reads the same raw
inputs the database was built from — VCF text, BED text, the manifest TSV — and
works out the expected counts from first principles. Reusing the code under test
to predict its own output would only assert that the code is self-consistent,
which is exactly what the existing suite already did while add-samples was
silently dropping samples on the floor.

The coverage half is the half that matters. AC, N_HET and N_HOM_ALT come from
the VCFs, but AN, N_HOM_REF and N_NO_COVERAGE come from capture-BED coverage
crossed with sex-dependent ploidy, and every silent counting bug this project
has had lived on that side.
"""

from __future__ import annotations

from pathlib import Path

# PAR coordinates, 1-based inclusive, written out here rather than imported.
PAR = {
    "GRCh37": {
        "chrX": [(60_001, 2_699_520), (154_931_044, 155_260_560)],
        "chrY": [(10_001, 2_649_520), (59_034_050, 59_363_566)],
    },
    "GRCh38": {
        "chrX": [(10_001, 2_781_479), (155_701_383, 156_030_895)],
        "chrY": [(10_001, 2_781_479), (56_887_903, 57_217_415)],
    },
}


def _canon(chrom: str) -> str:
    """Fold the naming conventions a BED or VCF may use onto one spelling."""
    c = chrom.strip()
    if c[:3].lower() == "chr":
        c = c[3:]
    if c.upper() in ("M", "MT"):
        return "chrM"
    if c.upper() in ("X", "Y"):
        return "chr" + c.upper()
    return "chr" + c


class Cohort:
    """Samples, their genotypes and their capture regions, read from raw files."""

    def __init__(self, manifest_path, bed_dir=None, genome_build="GRCh37"):
        self.genome_build = genome_build
        self.samples: list[str] = []
        self.sex: dict[str, str] = {}
        self.tech: dict[str, str] = {}
        self.phenotypes: dict[str, list[str]] = {}
        # (chrom, pos, ref, alt) -> {sample: "het" | "hom" | "fail"}
        self.calls: dict[tuple[str, int, str, str], dict[str, str]] = {}
        # tech -> list of (chrom, start_1based, end_1based); None means whole genome
        self.regions: dict[str, list[tuple[str, int, int]] | None] = {}

        manifest_path = Path(manifest_path)
        base = manifest_path.parent
        rows = _read_tsv(manifest_path)
        for row in rows:
            name = row["sample_name"]
            self.samples.append(name)
            self.sex[name] = row["sex"]
            self.tech[name] = row["tech_name"]
            self.phenotypes[name] = [
                c.strip() for c in row.get("phenotype_codes", "").split(",") if c.strip()
            ]
            vcf = Path(row["vcf_path"])
            if not vcf.is_absolute():
                vcf = base / vcf
            self._read_vcf(vcf, name)

        for tech in set(self.tech.values()):
            bed = None if bed_dir is None else Path(bed_dir) / f"{tech}.bed"
            if bed is not None and bed.exists():
                self.regions[tech] = _read_bed(bed)
            else:
                self.regions[tech] = None  # no BED: covered everywhere

    def _read_vcf(self, path: Path, sample: str) -> None:
        for line in Path(path).read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            chrom, pos, ref, alt, filt, fmt, call = (
                _canon(f[0]), int(f[1]), f[3], f[4], f[6], f[8], f[9]
            )
            gt = call.split(":")[fmt.split(":").index("GT")]
            alleles = [a for a in gt.replace("|", "/").split("/") if a != "."]
            passed = filt in ("PASS", ".", "")
            for i, one_alt in enumerate(alt.split(","), start=1):
                if one_alt == "*":
                    continue
                n = alleles.count(str(i))
                key = (chrom, pos, ref, one_alt)
                if not passed and not alleles:
                    self.calls.setdefault(key, {})[sample] = "fail"
                elif n == 0:
                    continue
                elif not passed:
                    self.calls.setdefault(key, {})[sample] = "fail"
                else:
                    self.calls.setdefault(key, {})[sample] = "hom" if n >= 2 else "het"

    def covers(self, sample: str, chrom: str, pos: int) -> bool:
        regions = self.regions[self.tech[sample]]
        if regions is None:
            return True
        return any(c == chrom and start <= pos <= end for c, start, end in regions)

    def is_haploid(self, sample: str, chrom: str, pos: int) -> bool:
        if chrom == "chrM":
            return True
        if chrom == "chrY":
            return True
        if chrom == "chrX" and not self._in_par(chrom, pos):
            return self.sex[sample] == "male"
        return False

    def _in_par(self, chrom: str, pos: int) -> bool:
        return any(s <= pos <= e for s, e in PAR[self.genome_build].get(chrom, []))

    def eligible(self, chrom: str, pos: int, samples=None) -> list[str]:
        pool = self.samples if samples is None else samples
        out = []
        for s in pool:
            if not self.covers(s, chrom, pos):
                continue
            # chrY carries no alleles for females, so they are not eligible there
            if chrom == "chrY" and self.sex[s] != "male":
                continue
            out.append(s)
        return out

    def expect(self, chrom, pos, ref, alt, samples=None) -> dict:
        """Expected query counts for one variant."""
        elig = self.eligible(chrom, pos, samples)
        calls = self.calls.get((chrom, pos, ref, alt), {})

        AN = 0
        AC = 0
        n_het = n_hom_alt = n_fail = 0
        for s in elig:
            haploid = self.is_haploid(s, chrom, pos)
            AN += 1 if haploid else 2
            state = calls.get(s)
            if state == "fail":
                n_fail += 1
            elif state == "het":
                if haploid:
                    AC += 1
                    n_hom_alt += 1
                else:
                    AC += 1
                    n_het += 1
            elif state == "hom":
                if haploid:
                    AC += 1
                else:
                    AC += 2
                n_hom_alt += 1

        return {
            "AC": AC,
            "AN": AN,
            "N_HET": n_het,
            "N_HOM_ALT": n_hom_alt,
            "N_FAIL": n_fail,
            "N_HOM_REF": len(elig) - n_het - n_hom_alt - n_fail,
        }

    def variants(self) -> list[tuple[str, int, str, str]]:
        return sorted(self.calls)

    def with_phenotype(self, code: str) -> list[str]:
        return [s for s in self.samples if code in self.phenotypes[s]]


def _read_tsv(path: Path) -> list[dict]:
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


def _read_bed(path: Path) -> list[tuple[str, int, int]]:
    """BED is 0-based half-open; convert to 1-based inclusive."""
    out = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split("\t")
        out.append((_canon(f[0]), int(f[1]) + 1, int(f[2])))
    return out

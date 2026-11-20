#!/usr/bin/env python3
"""Build a small database and dump it disaggregated by every group.

The cohort has one sample per sex x technology x phenotype cell, so the eight
groups of the dump are disjoint and cover the whole cohort. The variants
include these cases:

- chr1:1000  an insertion and a deletion starting at the same position
- chr1:5000  two SNVs at one position, with samples carrying both
- chr1:7000  a deletion, and chr1:7002 a SNV inside the deleted bases
- chr1:9000  a plain biallelic SNV
- chrX:5000000  outside the PARs, where males are haploid

Creates, next to this script (or in the directory given as argument):

- example_output/vcfs/        one VCF per sample, already normalized
- example_output/beds/        the capture BED of the WES technology
- example_output/manifest.tsv
- example_output/db/          the afquery database
- example_output/dump_total.csv    the dump without groups
- example_output/dump_groups.csv   the dump with --all-groups

All data is synthetic.
"""

import shutil
import sys
from pathlib import Path

from afquery import Database
from afquery.preprocess import run_preprocess

DEFAULT_OUTPUT = Path(__file__).resolve().parent / "example_output"

# (name, sex, technology, phenotype, [(chrom, pos, ref, alt, GT), ...])
# Records are written as normalize_vcf.sh would write them, with one ALT per
# record, so a sample carrying two ALT alleles at a position has two records.
SAMPLES = [
    ("S01", "female", "WGS", "E11", [("chr1", 1000, "A", "AT", "0/1"),
                                     ("chr1", 5000, "G", "A", "0/1"),
                                     ("chrX", 5000000, "G", "A", "0/1")]),
    ("S02", "male",   "WGS", "E11", [("chr1", 1000, "AC", "A", "0/1"),
                                     ("chr1", 5000, "G", "A", "1/1"),
                                     ("chrX", 5000000, "G", "A", "1")]),
    ("S03", "female", "WGS", "I10", [("chr1", 1000, "A", "AT", "1/0"),
                                     ("chr1", 1000, "AC", "A", "0/1"),
                                     ("chr1", 5000, "G", "T", "0/1")]),
    ("S04", "male",   "WGS", "I10", [("chr1", 5000, "G", "A", "1/0"),
                                     ("chr1", 5000, "G", "T", "0/1"),
                                     ("chr1", 7000, "TAGG", "T", "0/1")]),
    ("S05", "female", "WES", "E11", [("chr1", 1000, "AC", "A", "1/1"),
                                     ("chr1", 7002, "G", "C", "0/1"),
                                     ("chrX", 5000000, "G", "A", "1/1")]),
    ("S06", "male",   "WES", "E11", [("chr1", 7000, "TAGG", "T", "0/1"),
                                     ("chr1", 7002, "G", "C", "0/1")]),
    ("S07", "female", "WES", "I10", [("chr1", 5000, "G", "T", "0/1")]),
    ("S08", "male",   "WES", "I10", [("chr1", 9000, "C", "G", "0/1")]),
]

# The WES capture covers every variant of the example
WES_BED = "chr1\t900\t9100\nchrX\t4999000\t5001000\n"

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##FILTER=<ID=PASS,Description="All filters passed">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "##contig=<ID=chr1>\n"
    "##contig=<ID=chrX>\n"
)


def write_inputs(out: Path) -> Path:
    vcf_dir = out / "vcfs"
    bed_dir = out / "beds"
    vcf_dir.mkdir(parents=True)
    bed_dir.mkdir()
    (bed_dir / "WES.bed").write_text(WES_BED)

    rows = ["sample_name\tsex\ttech_name\tvcf_path\tphenotype_codes"]
    for name, sex, tech, phenotype, records in SAMPLES:
        vcf = vcf_dir / f"{name}.vcf"
        lines = [VCF_HEADER + f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}\n"]
        for chrom, pos, ref, alt, gt in sorted(records, key=lambda r: (r[0], r[1])):
            lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
        vcf.write_text("".join(lines))
        rows.append(f"{name}\t{sex}\t{tech}\t{vcf}\t{phenotype}")
    manifest = out / "manifest.tsv"
    manifest.write_text("\n".join(rows) + "\n")
    return manifest


def main() -> None:
    out = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else DEFAULT_OUTPUT
    shutil.rmtree(out, ignore_errors=True)
    manifest = write_inputs(out)

    db_dir = out / "db"
    db_dir.mkdir()
    run_preprocess(
        manifest_path=str(manifest), output_dir=str(db_dir),
        genome_build="GRCh37", threads=1, bed_dir=str(out / "beds"),
    )
    db = Database(str(db_dir))
    db.dump(output=str(out / "dump_total.csv"))
    db.dump(output=str(out / "dump_groups.csv"), all_groups=True)
    print(f"Wrote {out / 'dump_total.csv'} and {out / 'dump_groups.csv'}")


if __name__ == "__main__":
    main()

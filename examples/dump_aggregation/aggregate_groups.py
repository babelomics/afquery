#!/usr/bin/env python3
"""Sum the groups of an afquery dump and derive frequencies from the sums.

Usage:
    python aggregate_groups.py DUMP_CSV [GROUP_LABEL ...]

With no labels, every group in the dump is used. Writes CSV to stdout with,
for each row of the dump, the summed counts of the selected groups and two
frequencies:

    AF       = sum(AC) / sum(AN)          alternate allele frequency
    REF_AF   = sum(AC_REF) / sum(AN)      reference allele frequency

Each row is computed on its own, without looking at other rows of the same
position. The selected groups must be disjoint, with no sample in two of them.
Sex and technology groups always are, and phenotype groups are when every
sample has a single phenotype code. Frequencies are left empty when the summed
AN is 0.

Only the standard library is used.
"""

import csv
import sys

# Counts that add up over disjoint groups. AF is never summed.
SUMMED = ["AC", "AN", "AC_REF", "n_eligible",
          "N_HET", "N_HOM_ALT", "N_HEMI", "N_HOM_REF", "N_FAIL", "N_NO_COVERAGE"]


def group_labels(header: list[str]) -> list[str]:
    """Group labels of a dump, in column order.

    Labels are read from the n_eligible_<label> columns, because no other
    column name starts with "n_eligible_" and a label may itself contain "_"
    (female_WES_E11).
    """
    prefix = "n_eligible_"
    return [col[len(prefix):] for col in header if col.startswith(prefix)]


def aggregate(row: dict, labels: list[str]) -> dict:
    out = {key: row[key] for key in ("chrom", "pos", "ref", "alt")}
    for field in SUMMED:
        out[field] = sum(int(row[f"{field}_{label}"]) for label in labels)
    an = out["AN"]
    if an > 0:
        af = out["AC"] / an
        ref_af = out["AC_REF"] / an
        out.update(AF=af, REF_AF=ref_af)
    else:
        out.update(AF="", REF_AF="")
    return out


def main() -> None:
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    with open(sys.argv[1], newline="") as fh:
        reader = csv.DictReader(fh)
        available = group_labels(reader.fieldnames or [])
        labels = sys.argv[2:] or available
        unknown = [label for label in labels if label not in available]
        if unknown:
            sys.exit(f"Unknown group label(s): {', '.join(unknown)}. "
                     f"Available: {', '.join(available)}")
        writer = csv.DictWriter(
            sys.stdout,
            fieldnames=["chrom", "pos", "ref", "alt"] + SUMMED + ["AF", "REF_AF"],
            lineterminator="\n",
        )
        writer.writeheader()
        for row in reader:
            writer.writerow(aggregate(row, labels))


if __name__ == "__main__":
    main()

"""The dump aggregation example runs and gives the values its README shows."""
import csv
import importlib.util
import io
import sys
from contextlib import redirect_stdout
from pathlib import Path

EXAMPLE = Path(__file__).resolve().parent.parent / "examples" / "dump_aggregation"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, EXAMPLE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _aggregate(dump_csv, labels):
    agg = _load("aggregate_groups")
    argv, out = sys.argv, io.StringIO()
    sys.argv = ["aggregate_groups.py", str(dump_csv), *labels]
    try:
        with redirect_stdout(out):
            agg.main()
    finally:
        sys.argv = argv
    return list(csv.DictReader(io.StringIO(out.getvalue())))


def test_dump_aggregation_example(tmp_path):
    out = tmp_path / "example"
    argv = sys.argv
    sys.argv = ["create_example_data.py", str(out)]
    try:
        _load("create_example_data").main()
    finally:
        sys.argv = argv

    total = list(csv.DictReader(open(out / "dump_total.csv")))
    ac_ref = {(r["chrom"], r["pos"], r["alt"]): int(r["AC_REF"]) for r in total}
    assert ac_ref == {
        ("chr1", "1000", "A"): 10, ("chr1", "1000", "AT"): 10,
        ("chr1", "5000", "A"): 9, ("chr1", "5000", "T"): 9,
        ("chr1", "7000", "T"): 14, ("chr1", "7002", "C"): 14,
        ("chr1", "9000", "G"): 15, ("chrX", "5000000", "A"): 8,
    }

    # Summing every group gives back the whole-cohort dump
    summed = _aggregate(out / "dump_groups.csv", [])
    for want, got in zip(total, summed):
        for field in ("AC", "AN", "AC_REF", "n_eligible", "N_HET", "N_HOM_ALT",
                      "N_HEMI", "N_HOM_REF", "N_FAIL", "N_NO_COVERAGE"):
            assert got[field] == want[field], (want["pos"], want["alt"], field)

    e11 = _aggregate(out / "dump_groups.csv",
                     ["female_WGS_E11", "male_WGS_E11", "female_WES_E11", "male_WES_E11"])
    row = next(r for r in e11 if (r["pos"], r["alt"]) == ("1000", "A"))
    assert (row["AN"], row["AC_REF"], row["REF_AF"]) == ("8", "4", "0.5")

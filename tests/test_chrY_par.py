"""The chrY pseudoautosomal regions are not genotyped.

PAR1 and PAR2 are the same sequence on chrX and chrY, and a male carries one copy
on each, so a PAR locus is diploid in both sexes. AFQuery counts it once, at its
chrX coordinates. A call at chrY PAR coordinates (only possible against a
reference that does not mask the chrY PARs) is stored but nobody is eligible
there: it is not reported, and a male carrier is never taken for a hemizygote.
"""
import io
import warnings

import cyvcf2
import pytest

from afquery.models import AfqueryWarning

from test_synthetic_stats import _build_db

PAR_POS = 500_000        # GRCh37 chrY PAR1: 10,001-2,649,520
NON_PAR_POS = 2_700_000  # GRCh37 chrY outside the PARs


@pytest.fixture
def db(tmp_path):
    return _build_db(tmp_path, [
        ("M1", "male",   [("chrY", PAR_POS, "A", "G", "1"),
                          ("chrY", NON_PAR_POS, "C", "T", "1")]),
        ("M2", "male",   [("chrY", PAR_POS, "A", "G", "1/1")]),
        ("F1", "female", [("chr1", 1000, "A", "T", "0/1")]),
    ])


def test_point_query_in_par_warns_and_returns_nothing(db):
    with pytest.warns(AfqueryWarning, match="pseudoautosomal"):
        assert db.query(chrom="chrY", pos=PAR_POS) == []


def test_variant_info_in_par_warns_and_returns_nothing(db):
    with pytest.warns(AfqueryWarning, match="pseudoautosomal"):
        assert db.variant_info("chrY", PAR_POS) == []


def test_outside_par_is_reported_without_warning(db):
    with warnings.catch_warnings():
        warnings.simplefilter("error", AfqueryWarning)
        [r] = db.query(chrom="chrY", pos=NON_PAR_POS)
    assert (r.AC, r.AN, r.n_samples_eligible) == (1, 2, 2)


def test_region_and_batch_skip_the_par(db):
    region = db.query_region(chrom="chrY", start=1, end=10_000_000)
    assert [r.variant.pos for r in region] == [NON_PAR_POS]
    batch = db.query_batch(chrom="chrY", variants=[(PAR_POS, "A", "G"), (NON_PAR_POS, "C", "T")])
    assert [r.variant.pos for r in batch] == [NON_PAR_POS]


def test_dump_has_no_par_row(db):
    buf = io.StringIO()
    db.dump(output=buf, chrom="chrY", include_ac_zero=True)
    positions = [line.split(",")[1] for line in buf.getvalue().splitlines()[1:]]
    assert positions == [str(NON_PAR_POS)]


def test_annotate_gives_no_allele_number_in_par(db, tmp_path):
    vcf_in = tmp_path / "in.vcf"
    vcf_in.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chrY>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chrY\t{PAR_POS}\t.\tA\tG\t.\tPASS\t.\n"
        f"chrY\t{NON_PAR_POS}\t.\tC\tT\t.\tPASS\t.\n"
    )
    out = tmp_path / "out.vcf"
    stats = db.annotate_vcf(str(vcf_in), str(out))
    by_pos = {v.POS: v for v in cyvcf2.VCF(str(out))}
    assert by_pos[PAR_POS].INFO.get("AFQUERY_AN") == 0
    assert by_pos[PAR_POS].INFO.get("AFQUERY_AC") is None
    assert by_pos[NON_PAR_POS].INFO.get("AFQUERY_AN") == 2
    assert stats["n_uncovered"] == 1

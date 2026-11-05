"""Checks on resources/normalize_vcf.sh that need no bcftools."""
import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "resources" / "normalize_vcf.sh"

pytestmark = pytest.mark.skipif(shutil.which("bash") is None, reason="bash not available")


def test_script_parses():
    subprocess.run(["bash", "-n", str(SCRIPT)], check=True)


@pytest.mark.parametrize("args", [[], ["s1", "in.vcf.gz"], ["s1", "in.vcf.gz", "ref.fa", "8", "extra"]])
def test_wrong_argument_count_prints_usage(args):
    result = subprocess.run(["bash", str(SCRIPT), *args], capture_output=True, text=True)
    assert result.returncode == 1
    assert "<vcf_id> <vcf> <ref_genome> [threads" in result.stdout


def test_genotype_ploidy_is_left_to_afquery():
    """AFQuery takes ploidy from the sample's sex and the position. fixploidy
    would keep only the first allele of a call (0/1 -> 0), turning a male's
    heterozygous call on chrX into hom-ref and dropping heteroplasmic chrM calls."""
    assert "fixploidy" not in SCRIPT.read_text()

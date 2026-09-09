"""Properties that must hold however a database was assembled.

These compare a database against another database rather than against a stored
expectation, so they keep holding as the fixtures change. The build-equivalence
property is the one that would have caught add-samples writing into a layout
the query engine does not read.
"""
import shutil

import pytest

from afquery import storage
from afquery.database import Database
from afquery.preprocess import run_preprocess
from afquery.preprocess.build import build_all_parquets
from afquery.preprocess.ingest import ingest_all
from afquery.models import Sample

from test_update import write_manifest, write_vcf

# (name, sex, tech, [(chrom, pos, ref, alt, gt), ...])
GROUP_A = [
    ("A0", "male",   "wgs", [("chr1", 1000, "A", "T", "0/1"),
                             ("chr1", 2_400_000, "G", "C", "1/1"),
                             ("chrX", 4_000_000, "T", "A", "0/1")]),
    ("A1", "female", "wgs", [("chr1", 1000, "A", "T", "1/1"),
                             ("chrM", 300, "C", "T", "0/1")]),
]
GROUP_B = [
    ("B0", "female", "wgs", [("chr1", 1000, "A", "T", "0/1"),
                             ("chr1", 5_600_000, "A", "G", "0/1"),
                             ("chr2", 700, "T", "C", "1/1")]),
    ("B1", "male",   "wgs", [("chr1", 2_400_000, "G", "C", "0/1"),
                             ("chrY", 700_000, "A", "T", "1/1"),
                             ("chrX", 4_000_000, "T", "A", "1/1")]),
]

ALL_SITES = sorted({
    (chrom, pos, ref, alt)
    for _n, _s, _t, variants in GROUP_A + GROUP_B
    for chrom, pos, ref, alt, _gt in variants
})


def _write_group(tmp_path, group):
    entries = []
    for name, sex, tech, variants in group:
        vcf = tmp_path / f"{name}.vcf"
        write_vcf(str(vcf), name, variants)
        entries.append((name, sex, tech, str(vcf), "COHORT"))
    return entries


def _manifest(tmp_path, entries, filename):
    path = tmp_path / filename
    write_manifest(str(path), entries)
    return str(path)


def _snapshot(db_dir):
    """Every site's counts, keyed by variant — comparable across sample id orders."""
    db = Database(str(db_dir))
    out = {}
    for chrom, pos, ref, alt in ALL_SITES:
        for r in db.query(chrom=chrom, pos=pos, sex="both"):
            out[(chrom, pos, r.variant.ref, r.variant.alt)] = (
                r.AC, r.AN, r.N_HET, r.N_HOM_ALT, r.N_HOM_REF, r.N_FAIL,
            )
    return out


@pytest.fixture
def groups(tmp_path):
    a = _write_group(tmp_path, GROUP_A)
    b = _write_group(tmp_path, GROUP_B)
    return a, b


def test_add_samples_equals_building_the_union(tmp_path, groups):
    """build(A + B) and build(A) then add(B) must answer identically.

    Compared at query level rather than bitmap level: the two routes assign
    sample ids in a different order, so the bytes legitimately differ while
    every answer must not.
    """
    a, b = groups

    together = tmp_path / "together"
    run_preprocess(manifest_path=_manifest(tmp_path, a + b, "ab.tsv"),
                   output_dir=str(together), genome_build="GRCh37", threads=1)

    incremental = tmp_path / "incremental"
    run_preprocess(manifest_path=_manifest(tmp_path, a, "a.tsv"),
                   output_dir=str(incremental), genome_build="GRCh37", threads=1)
    from afquery.preprocess.update import add_samples
    add_samples(str(incremental), _manifest(tmp_path, b, "b.tsv"), threads=1)

    assert _snapshot(incremental) == _snapshot(together)


def test_flat_and_partitioned_answer_identically(tmp_path, groups):
    """The storage layout must not be observable through any query."""
    a, b = groups
    manifest = _manifest(tmp_path, a + b, "all.tsv")

    part = tmp_path / "part"
    run_preprocess(manifest_path=manifest, output_dir=str(part),
                   genome_build="GRCh37", threads=1)

    flat = tmp_path / "flat"
    shutil.copytree(part, flat)
    shutil.rmtree(flat / "variants")
    (flat / "variants").mkdir()
    samples = [Sample(i, name, sex, 0)
               for i, (name, sex, _t, _v, _p) in enumerate(a + b)]
    vcfs = [str(tmp_path / f"{name}.vcf") for name, _s, _t, _v, _p in a + b]
    ingest_tmp = tmp_path / "ingest"
    ingest_tmp.mkdir()
    ingest_all(samples, vcfs, str(ingest_tmp), n_workers=1)
    build_all_parquets(str(ingest_tmp), str(flat / "variants"),
                       n_workers=1, partitioned=False)

    assert storage.detect_layout(part / "variants") == storage.PARTITIONED
    assert storage.detect_layout(flat / "variants") == storage.FLAT
    assert _snapshot(flat) == _snapshot(part)


def test_remove_undoes_add(tmp_path, groups):
    """build(A) then add(B) then remove(B) must answer like build(A).

    Removal clears bits but leaves the rows behind — that is what compact is
    for — so sites only group B carried survive as AC=0 rows. Every site the
    smaller database knows must answer identically, and the leftovers must
    carry nobody.
    """
    a, b = groups

    alone = tmp_path / "alone"
    run_preprocess(manifest_path=_manifest(tmp_path, a, "a_only.tsv"),
                   output_dir=str(alone), genome_build="GRCh37", threads=1)
    before = _snapshot(alone)

    grown = tmp_path / "grown"
    run_preprocess(manifest_path=_manifest(tmp_path, a, "a_dup.tsv"),
                   output_dir=str(grown), genome_build="GRCh37", threads=1)
    from afquery.preprocess.update import add_samples, remove_samples
    add_samples(str(grown), _manifest(tmp_path, b, "b_dup.tsv"), threads=1)
    remove_samples(str(grown), [name for name, _s, _t, _v, _p in b])

    after = _snapshot(grown)
    for key, counts in before.items():
        assert after[key] == counts, f"{key} changed across add-then-remove"
    for key, counts in after.items():
        if key not in before:
            AC, _AN, n_het, n_hom_alt, _n_hom_ref, n_fail = counts
            assert (AC, n_het, n_hom_alt, n_fail) == (0, 0, 0, 0), (
                f"{key} still has carriers after their only samples were removed"
            )


def test_compact_preserves_carried_variants_and_is_idempotent(tmp_path, groups):
    """Compact may drop rows nobody carries, but must not change any that remain.

    B1 is the only carrier on chrY:700000, so removing it empties that bucket
    entirely — the case where compact used to crash on an empty row selection.
    """
    a, b = groups
    db = tmp_path / "compactme"
    run_preprocess(manifest_path=_manifest(tmp_path, a + b, "cab.tsv"),
                   output_dir=str(db), genome_build="GRCh37", threads=1)
    from afquery.preprocess.compact import compact_database
    from afquery.preprocess.update import remove_samples

    remove_samples(str(db), ["B1"])
    before = _snapshot(db)

    compact_database(db)
    once = _snapshot(db)

    for key, counts in once.items():
        assert counts == before[key], f"{key} changed across compact"
    for key in set(before) - set(once):
        assert before[key][0] == 0, f"compact dropped {key}, which still had carriers"

    compact_database(db)
    assert _snapshot(db) == once

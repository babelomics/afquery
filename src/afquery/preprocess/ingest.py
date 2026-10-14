import logging
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import tqdm

import pyarrow as pa
import pyarrow.parquet as pq

from ..constants import normalize_chrom
from ..models import Sample
from ..ploidy import sample_ploidy

logger = logging.getLogger(__name__)

INGEST_SCHEMA = pa.schema([
    ("chrom",       pa.utf8()),
    ("pos",         pa.uint32()),
    ("ref",         pa.utf8()),
    ("alt",         pa.utf8()),
    ("gt_ac",       pa.uint8()),
    ("sample_id",   pa.uint32()),
    ("filter_pass", pa.bool_()),
    ("dp",          pa.int32()),    # FORMAT/DP — None when absent
    ("gq",          pa.int32()),    # FORMAT/GQ — None when absent
    ("qual",        pa.float32()),  # VCF QUAL — None when '.'
])


class IngestError(RuntimeError):
    pass


def _read_format_int(variant, field: str) -> int | None:
    """Read an integer FORMAT field for the (single) sample.

    Returns None when the field is absent, missing (cyvcf2 sentinel), or non-numeric.
    """
    try:
        arr = variant.format(field)
    except Exception:
        return None
    if arr is None:
        return None
    try:
        v = arr[0][0]
    except (IndexError, TypeError):
        return None
    if v is None:
        return None
    try:
        v_int = int(v)
    except (TypeError, ValueError):
        return None
    # cyvcf2 returns -2147483648 for missing integers
    if v_int < 0:
        return None
    return v_int


def _exceeds_ploidy(alts: list[str], gt_acs: list[int], ploidy: int) -> bool:
    """True if one sample's calls at a position carry more alleles than it has.

    A haploid sample can carry one ALT allele at most; a diploid one, two copies.
    """
    carried = [(alt, ac) for alt, ac in zip(alts, gt_acs) if ac > 0]
    if ploidy == 1:
        return len({alt for alt, _ac in carried}) > 1
    if ploidy == 2:
        return sum(ac for _alt, ac in carried) > 2
    return False  # no allele here (chrY in a female): never eligible


def ingest_sample(
    sample_id: int, vcf_path: str, tmp_dir: str, sex: str, genome_build: str,
) -> tuple[str, float, list[str]]:
    """Parse one VCF, write Parquet file.

    A position where the sample's calls carry more alleles than its ploidy
    allows (e.g. a 1/2 call on chrX in a male) has no usable genotype: all its
    calls there are recorded as failed.

    Returns (output_path, elapsed_seconds, positions recorded as failed for
    that reason, as "chrom:pos").
    """
    import cyvcf2  # local import for worker-safety

    t0 = time.monotonic()
    out_path = os.path.join(tmp_dir, f"sample_{sample_id}.parquet")

    chroms: list[str] = []
    positions: list[int] = []
    refs: list[str] = []
    alts: list[str] = []
    gt_acs: list[int] = []
    sample_ids: list[int] = []
    filter_passes: list[bool] = []
    dps: list[int | None] = []
    gqs: list[int | None] = []
    quals: list[float | None] = []
    inconsistent: list[str] = []

    # Rows of the current position, which are contiguous in a sorted VCF
    site: tuple[str, int] | None = None
    site_start = 0

    def close_site() -> None:
        if site is None:
            return
        ploidy = sample_ploidy(site[0], site[1], sex, genome_build)
        if _exceeds_ploidy(alts[site_start:], gt_acs[site_start:], ploidy):
            for i in range(site_start, len(filter_passes)):
                filter_passes[i] = False
            inconsistent.append(f"{site[0]}:{site[1]}")

    vcf = cyvcf2.VCF(vcf_path)
    for variant in vcf:
        if not (variant.is_snp or variant.is_indel):
            continue

        gt = variant.genotypes[0]
        alleles = [a for a in gt[:-1] if a >= 0]

        chrom = normalize_chrom(variant.CHROM)
        pos = variant.POS
        if (chrom, pos) != site:
            close_site()
            site, site_start = (chrom, pos), len(chroms)
        ref = variant.REF
        # cyvcf2: FILTER is None for PASS/missing, string for others (e.g. "LowQual")
        fp = variant.FILTER is None or variant.FILTER == "PASS"

        # Quality fields (None when missing/absent)
        dp_val = _read_format_int(variant, "DP")
        gq_val = _read_format_int(variant, "GQ")
        qual_val = variant.QUAL  # float or None

        # Missing GT (./.) at a failed site: track as N_FAIL for all ALTs
        if not alleles:
            if not fp:
                for alt_str in variant.ALT:
                    if alt_str == "*":
                        continue
                    chroms.append(chrom)
                    positions.append(pos)
                    refs.append(ref)
                    alts.append(alt_str)
                    gt_acs.append(0)
                    sample_ids.append(sample_id)
                    filter_passes.append(False)
                    dps.append(dp_val)
                    gqs.append(gq_val)
                    quals.append(qual_val)
            continue

        for idx, alt_str in enumerate(variant.ALT):
            if alt_str == "*":
                continue
            alt_allele_num = idx + 1
            ac = alleles.count(alt_allele_num)
            if ac == 0:
                continue

            chroms.append(chrom)
            positions.append(pos)
            refs.append(ref)
            alts.append(alt_str)
            gt_acs.append(ac)
            sample_ids.append(sample_id)
            filter_passes.append(fp)
            dps.append(dp_val)
            gqs.append(gq_val)
            quals.append(qual_val)

    close_site()
    vcf.close()

    table = pa.table(
        {
            "chrom":       pa.array(chroms,        type=pa.utf8()),
            "pos":         pa.array(positions,     type=pa.uint32()),
            "ref":         pa.array(refs,          type=pa.utf8()),
            "alt":         pa.array(alts,          type=pa.utf8()),
            "gt_ac":       pa.array(gt_acs,        type=pa.uint8()),
            "sample_id":   pa.array(sample_ids,    type=pa.uint32()),
            "filter_pass": pa.array(filter_passes, type=pa.bool_()),
            "dp":          pa.array(dps,           type=pa.int32()),
            "gq":          pa.array(gqs,           type=pa.int32()),
            "qual":        pa.array(quals,         type=pa.float32()),
        },
        schema=INGEST_SCHEMA,
    )

    pq.write_table(table, out_path)
    elapsed = time.monotonic() - t0
    return out_path, elapsed, inconsistent


def _ingest_worker(args: tuple) -> tuple[str, float, list[str]]:
    return ingest_sample(*args)


def _warn_inconsistent(sample_name: str, positions: list[str]) -> None:
    shown = ", ".join(positions[:10])
    more = f" and {len(positions) - 10} more" if len(positions) > 10 else ""
    logger.warning(
        "[ingest] %s: %d position(s) with a genotype inconsistent with the sample's "
        "ploidy, recorded as failed calls: %s%s",
        sample_name, len(positions), shown, more,
    )


def ingest_all(
    samples: list[Sample],
    vcf_paths: list[str],
    tmp_dir: str,
    genome_build: str,
    n_workers: int = 8,
    resume: bool = True,
) -> list[str]:
    """Parallel ingestion via ProcessPoolExecutor. Collects all errors before raising.

    Warns, per sample, about positions recorded as failed because the genotype
    exceeds the sample's ploidy. Samples skipped on resume were reported by the
    run that ingested them.
    """
    args_list = [
        (s.sample_id, vcf_paths[i], tmp_dir, s.sex, genome_build)
        for i, s in enumerate(samples)
    ]
    names = {s.sample_id: s.sample_name for s in samples}

    if resume:
        pending = []
        n_skip = 0
        for a in args_list:
            parquet_path = os.path.join(tmp_dir, f"sample_{a[0]}.parquet")
            if os.path.exists(parquet_path):
                n_skip += 1
            else:
                pending.append(a)
        if n_skip:
            logger.info("[ingest] Skipping %d already-ingested sample(s); ingesting %d remaining.",
                        n_skip, len(pending))
        args_list = pending

    if not args_list:
        logger.info("[ingest] All samples already ingested, nothing to do.")
        return []

    logger.info("[ingest] Ingesting %d VCF(s) with %d worker(s)...", len(args_list), n_workers)
    t0 = time.monotonic()

    results: list[str] = []
    errors: list[str] = []

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        future_to_args = {
            executor.submit(_ingest_worker, args): args for args in args_list
        }
        with tqdm.tqdm(
            total=len(args_list),
            unit="VCF",
            desc="[ingest]",
            dynamic_ncols=True,
            disable=not sys.stderr.isatty(),
        ) as bar:
            for future in as_completed(future_to_args):
                args = future_to_args[future]
                try:
                    out_path, elapsed, inconsistent = future.result()
                    results.append(out_path)
                    if inconsistent:
                        _warn_inconsistent(names[args[0]], inconsistent)
                    bar.update(1)
                    logger.debug("  [ingest] sample_id=%d (%s) done in %.2fs",
                                 args[0], os.path.basename(args[1]), elapsed)
                except Exception as e:
                    errors.append(f"Sample {args[0]} ({args[1]}): {e}")

    logger.info("[ingest] Complete: %d ok, %d error(s) (%.1fs)",
                len(results), len(errors), time.monotonic() - t0)

    if errors:
        raise IngestError("VCF ingestion failed:\n" + "\n".join(errors))

    return results

import json
import sqlite3
import warnings
from pathlib import Path

import duckdb
from pyroaring import BitMap

from .bitmaps import deserialize
from .capture import CaptureIndex, describe_capture_problem, load_capture_indices
from .constants import normalize_chrom, ALL_CHROMS, CHROM_ORDER
from .models import AfqueryWarning, QueryParams, QueryResult, SampleCarrier, VariantKey, Sample, Technology, SampleFilter
from .ploidy import compute_AN, split_ploidy

BATCH_IN_THRESHOLD = 10_000


def _parse_schema_version(s: str) -> tuple[int, ...]:
    try:
        return tuple(int(p) for p in s.split("."))
    except (ValueError, AttributeError):
        return (0,)


class QueryEngine:
    def __init__(self, db_path: str):
        self._db = Path(db_path)
        self._manifest = json.loads((self._db / "manifest.json").read_text())
        self._genome_build = self._manifest["genome_build"]
        self._schema_version = self._manifest.get("schema_version", "1.0")
        self._has_coverage_data = _parse_schema_version(self._schema_version) >= (3, 0)
        self._coverage_filter = self._manifest.get("coverage_filter", {})

        con = sqlite3.connect(self._db / "metadata.sqlite")

        # Load precomputed bitmaps
        self._bitmaps: dict[str, dict[str, BitMap]] = {}
        for btype, bkey, bdata in con.execute(
            "SELECT bitmap_type, bitmap_key, bitmap_data FROM precomputed_bitmaps"
        ).fetchall():
            self._bitmaps.setdefault(btype, {})[bkey] = deserialize(bytes(bdata))

        # Load samples and technologies
        samples = [
            Sample(r[0], r[1], r[2], r[3])
            for r in con.execute(
                "SELECT sample_id, sample_name, sex, tech_id FROM samples"
            ).fetchall()
        ]
        techs = [
            Technology(r[0], r[1], r[2])
            for r in con.execute(
                "SELECT tech_id, tech_name, bed_path FROM technologies"
            ).fetchall()
        ]

        # Phenotype cache: sample_id → list of phenotype codes
        self._sample_phenotypes: dict[int, list[str]] = {}
        for row in con.execute(
            "SELECT sample_id, phenotype_code FROM sample_phenotype"
        ).fetchall():
            self._sample_phenotypes.setdefault(row[0], []).append(row[1])

        con.close()

        self._samples = samples
        self._samples_by_id: dict[int, Sample] = {s.sample_id: s for s in samples}
        self._tech_map = {t.tech_id: t for t in techs}
        sex_bms = self._bitmaps.get("sex", {})
        self._male_bm = sex_bms.get("male", BitMap())
        self._female_bm = sex_bms.get("female", BitMap())
        self._capture = load_capture_indices(techs, str(self._db / "capture"))
        for tech_id, capture_idx in self._capture.items():
            problem = describe_capture_problem(
                capture_idx, self._tech_map[tech_id].tech_name
            )
            if problem is not None:
                warnings.warn(problem, AfqueryWarning, stacklevel=2)
        self._tech_bitmaps = self._bitmaps.get("tech", {})
        self._all_samples_bm = BitMap(s.sample_id for s in self._samples)
        self._tech_name_to_id: dict[str, str] = {
            t.tech_name: str(t.tech_id) for t in self._tech_map.values()
        }

        # Cache which chroms use partitioned storage (variants/{chrom}/ directory)
        self._partitioned_chroms: set[str] = set()
        variants_dir = self._db / "variants"
        if variants_dir.exists():
            for p in variants_dir.iterdir():
                if p.is_dir():
                    self._partitioned_chroms.add(p.name)

        self._flat_chroms: set[str] = set()
        if variants_dir.exists():
            for p in variants_dir.iterdir():
                if p.is_file() and p.suffix == ".parquet":
                    self._flat_chroms.add(p.stem)
        self._all_known_chroms: set[str] = self._partitioned_chroms | self._flat_chroms

        # Precompute covered bitmap for positions where ALL technologies cover (common case)
        # For each tech, store its bitmap; also precompute the union of all WGS tech bitmaps
        self._covered_all_bm: BitMap = BitMap()
        for tech_id, capture_idx in self._capture.items():
            if capture_idx.is_always_covered:
                bm = self._tech_bitmaps.get(str(tech_id), BitMap())
                self._covered_all_bm |= bm

        # Targeted (non-WGS) technologies only. The per-position loops below run once
        # per queried position, so the WGS ones are skipped here rather than re-tested
        # on every iteration — they are already folded into _covered_all_bm.
        self._targeted_capture: dict[int, CaptureIndex] = {
            tech_id: capture_idx
            for tech_id, capture_idx in self._capture.items()
            if not capture_idx.is_always_covered
        }

        # Cache parquet glob patterns per chrom (avoids filesystem scan on every query)
        self._glob_cache: dict[str, str | None] = {}

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _validate_sample_filter(self, sf: SampleFilter) -> None:
        ph_bms = self._bitmaps.get("phenotype", {})
        for code in sf.phenotype_include:
            if code not in ph_bms:
                warnings.warn(
                    f"Phenotype {code!r} not in database — include will match 0 samples.",
                    AfqueryWarning, stacklevel=4,
                )
        for code in sf.phenotype_exclude:
            if code not in ph_bms:
                warnings.warn(
                    f"Phenotype {code!r} not in database — exclude has no effect.",
                    AfqueryWarning, stacklevel=4,
                )
        for name in sf.tech_include:
            if name not in self._tech_name_to_id:
                warnings.warn(
                    f"Technology {name!r} not in database — include will match 0 samples.",
                    AfqueryWarning, stacklevel=4,
                )
        for name in sf.tech_exclude:
            if name not in self._tech_name_to_id:
                warnings.warn(
                    f"Technology {name!r} not in database — exclude has no effect.",
                    AfqueryWarning, stacklevel=4,
                )
        if sf.sex not in ("male", "female", "both"):
            raise ValueError(f"Invalid sex {sf.sex!r}. Must be 'male', 'female', or 'both'.")

    def _build_sample_bitmap(self, sf: SampleFilter) -> BitMap:
        """Construye bitmap de samples elegibles según filtros de fenotipo, sexo y tecnología."""
        self._validate_sample_filter(sf)
        # --- Fenotipo ---
        ph_bms = self._bitmaps.get("phenotype", {})
        if sf.phenotype_include:
            ph_bm = BitMap()
            for code in sf.phenotype_include:
                if code in ph_bms:
                    ph_bm |= ph_bms[code]
        else:
            ph_bm = BitMap(self._all_samples_bm)  # copia, todas las muestras
        for code in sf.phenotype_exclude:
            if code in ph_bms:
                ph_bm -= ph_bms[code]

        # --- Sexo ---
        sex_bms = self._bitmaps.get("sex", {})
        if sf.sex == "both":
            sex_bm = sex_bms.get("male", BitMap()) | sex_bms.get("female", BitMap())
        else:
            sex_bm = sex_bms.get(sf.sex, BitMap())

        # --- Tecnología (nivel de sample; cobertura por posición sigue en _compute_eligible) ---
        tech_bms = self._bitmaps.get("tech", {})
        if sf.tech_include or sf.tech_exclude:
            if sf.tech_include:
                tech_bm = BitMap()
                for name in sf.tech_include:
                    tid = self._tech_name_to_id.get(name)
                    if tid and tid in tech_bms:
                        tech_bm |= tech_bms[tid]
            else:
                tech_bm = BitMap(self._all_samples_bm)
            for name in sf.tech_exclude:
                tid = self._tech_name_to_id.get(name)
                if tid and tid in tech_bms:
                    tech_bm -= tech_bms[tid]
        else:
            tech_bm = self._all_samples_bm  # sin copia — solo lectura en la intersección

        result_bm = ph_bm & sex_bm & tech_bm
        if len(result_bm) == 0:
            warnings.warn(
                "Sample filter produces an empty eligible set — all queries will return AN=0. "
                "Check that include/exclude filters are not contradictory.",
                AfqueryWarning, stacklevel=3,
            )
        return result_bm

    def _compute_no_coverage_bm(
        self,
        eligible: BitMap,
        het_bm: BitMap,
        hom_bm: BitMap,
        fail_bm: BitMap,
        min_pass: int,
        min_observed: int,
        filtered_bm: "BitMap | None" = None,
        quality_pass_bm: "BitMap | None" = None,
        min_quality_evidence: int = 0,
    ) -> BitMap:
        """Return bitmap of WES non-carrier samples that should not be assumed hom-ref.

        Phase 1 (query-time): count-based per-tech gate using existing bitmaps.
        Phase 2 (build-time): stored filtered_bitmap + optional quality_pass gate.
        Results are unioned; carriers (het/hom/fail) are never included.
        """
        if min_quality_evidence > 0 and not self._has_coverage_data:
            raise ValueError(
                "This database was not built with coverage quality data. "
                "Re-create with --min-dp / --min-gq to use --min-quality-evidence."
            )
        no_cov = BitMap()
        all_carriers = het_bm | hom_bm | fail_bm

        # Phase 1: count-based gate
        if min_pass > 0 or min_observed > 0:
            for tech_id in self._targeted_capture:
                tech_bm = self._tech_bitmaps.get(str(tech_id), BitMap())
                tech_eligible = eligible & tech_bm
                if len(tech_eligible) == 0:
                    continue
                pass_count     = len((het_bm | hom_bm) & tech_eligible)
                observed_count = len(all_carriers & tech_eligible)
                if pass_count < min_pass or observed_count < min_observed:
                    no_cov |= tech_eligible - (all_carriers & tech_eligible)

        # Phase 2: stored filtered_bitmap
        if filtered_bm is not None:
            no_cov |= filtered_bm & eligible

        # Phase 2: quality_pass gate (--min-quality-evidence K)
        if quality_pass_bm is not None and min_quality_evidence > 0:
            already_filtered = filtered_bm if filtered_bm is not None else BitMap()
            for tech_id in self._targeted_capture:
                tech_bm = self._tech_bitmaps.get(str(tech_id), BitMap())
                tech_eligible = eligible & tech_bm
                if len(tech_eligible) == 0:
                    continue
                quality_count = len(quality_pass_bm & tech_eligible)
                if quality_count < min_quality_evidence:
                    # additionally filter non-carrier, not-already-filtered samples
                    extra = tech_eligible - (all_carriers & tech_eligible) - already_filtered
                    no_cov |= extra

        return no_cov

    def _bitmap_cols(self, with_pos: bool = False) -> list[str]:
        """List of SELECT columns for variant Parquet reads.

        Returns 5 entries for legacy DBs and 7 for Phase 2 DBs.
        """
        cols = (["pos"] if with_pos else []) + [
            "ref", "alt", "het_bitmap", "hom_bitmap", "fail_bitmap",
        ]
        if self._has_coverage_data:
            cols += ["filtered_bitmap", "quality_pass_bitmap"]
        return cols

    def _unpack_bitmaps(
        self, bitmap_cols: tuple,
    ) -> tuple[BitMap, BitMap, BitMap, "BitMap | None", "BitMap | None"]:
        """Deserialize the bitmap suffix of a Parquet row.

        Accepts the trailing 3 cols (legacy) or 5 cols (Phase 2) of a SELECT row.
        Returns: (het, hom, fail, filtered, quality_pass).
        filtered and quality_pass are None for legacy DBs.
        """
        het_bm = deserialize(bytes(bitmap_cols[0]))
        hom_bm = deserialize(bytes(bitmap_cols[1]))
        fail_bm = deserialize(bytes(bitmap_cols[2]))
        if self._has_coverage_data and len(bitmap_cols) >= 5:
            filtered_bm = deserialize(bytes(bitmap_cols[3]))
            quality_pass_bm = deserialize(bytes(bitmap_cols[4]))
        else:
            filtered_bm = None
            quality_pass_bm = None
        return het_bm, hom_bm, fail_bm, filtered_bm, quality_pass_bm

    def _compute_eligible(
        self,
        chrom: str,
        pos: int,
        sample_bitmap: BitMap,
    ) -> tuple[BitMap, int]:
        """Return (eligible_bitmap, AN) for a single position."""
        covered = BitMap(self._covered_all_bm)  # start with WGS samples (always covered)
        for tech_id, capture_idx in self._targeted_capture.items():
            if capture_idx.covers(chrom, pos):
                covered |= self._tech_bitmaps.get(str(tech_id), BitMap())
        eligible = sample_bitmap & covered
        AN = compute_AN(
            eligible,
            self._male_bm,
            self._female_bm,
            chrom, pos, self._genome_build,
        )
        return eligible, AN

    def _parquet_path(self, chrom: str, pos: int) -> Path | None:
        """Resolve Parquet path for a point query. Partitioned takes priority over flat."""
        if chrom in self._partitioned_chroms:
            bucket = pos // 1_000_000
            p = self._db / "variants" / chrom / f"bucket_{bucket}.parquet"
            return p if p.exists() else None
        flat = self._db / "variants" / f"{chrom}.parquet"
        return flat if flat.exists() else None

    def _parquet_glob(self, chrom: str) -> str | None:
        """Return path/glob pattern for batch/region queries. None if no data for chrom."""
        if chrom in self._glob_cache:
            return self._glob_cache[chrom]
        if chrom in self._partitioned_chroms:
            chrom_dir = self._db / "variants" / chrom
            result = str(chrom_dir / "bucket_*.parquet") if any(chrom_dir.glob("bucket_*.parquet")) else None
        else:
            flat = self._db / "variants" / f"{chrom}.parquet"
            result = str(flat) if flat.exists() else None
        self._glob_cache[chrom] = result
        return result

    # ------------------------------------------------------------------
    # Public query methods
    # ------------------------------------------------------------------

    def query(self, params: QueryParams) -> list[QueryResult]:
        chrom = normalize_chrom(params.chrom)
        pos = params.pos

        if chrom not in self._all_known_chroms:
            warnings.warn(
                f"Chromosome {chrom!r} has no data in this database. "
                f"Available: {sorted(self._all_known_chroms)!r}",
                AfqueryWarning, stacklevel=3,
            )
            return []

        sample_bm = self._build_sample_bitmap(params.filter)
        eligible, AN = self._compute_eligible(chrom, pos, sample_bm)

        if AN == 0:
            return []

        parquet_path = self._parquet_path(chrom, pos)
        if parquet_path is None:
            return []

        cols = ", ".join(self._bitmap_cols(with_pos=False))
        con = duckdb.connect(config={"temp_directory": "/tmp"})
        rows = con.execute(
            f"SELECT {cols} FROM read_parquet(?) WHERE pos = ?",
            [str(parquet_path), pos],
        ).fetchall()
        con.close()

        if not rows:
            return []

        results = []
        sf = params.filter
        for row in rows:
            ref, alt = row[0], row[1]
            het_bm, hom_bm, fail_bm, filtered_bm, quality_pass_bm = self._unpack_bitmaps(row[2:])
            haploid_elig, diploid_elig = split_ploidy(
                eligible, self._male_bm, self._female_bm, chrom, pos, self._genome_build
            )
            het_elig = het_bm & eligible
            hom_elig = hom_bm & eligible
            AC = (len((het_elig | hom_elig) & haploid_elig)
                  + len(het_elig & diploid_elig)
                  + 2 * len(hom_elig & diploid_elig))
            N_HET = len(het_elig & diploid_elig)
            N_HOM_ALT = len(hom_elig & diploid_elig) + len((het_elig | hom_elig) & haploid_elig)
            AF = AC / AN if AN > 0 else None
            N_FAIL = len(fail_bm & eligible)
            no_cov_bm = self._compute_no_coverage_bm(
                eligible, het_bm, hom_bm, fail_bm,
                sf.min_pass, sf.min_observed,
                filtered_bm=filtered_bm,
                quality_pass_bm=quality_pass_bm,
                min_quality_evidence=sf.min_quality_evidence,
            )
            N_NO_COVERAGE = len(no_cov_bm)
            N_HOM_REF = len(eligible) - N_HET - N_HOM_ALT - N_FAIL - N_NO_COVERAGE
            results.append(QueryResult(
                variant=VariantKey(chrom=chrom, pos=pos, ref=ref, alt=alt),
                AC=AC, AN=AN, AF=AF, n_samples_eligible=len(eligible),
                N_HET=N_HET, N_HOM_ALT=N_HOM_ALT, N_HOM_REF=N_HOM_REF,
                N_FAIL=N_FAIL, N_NO_COVERAGE=N_NO_COVERAGE,
            ))
        if params.ref is not None:
            results = [r for r in results if r.variant.ref == params.ref]
        if params.alt is not None:
            results = [r for r in results if r.variant.alt == params.alt]
        return results

    def query_batch(
        self,
        chrom: str,
        variants: list[tuple[int, str, str]],
        sf: SampleFilter,
    ) -> list[QueryResult]:
        """Query multiple variants (pos, ref, alt) in a single Parquet read."""
        chrom = normalize_chrom(chrom)

        if chrom not in self._all_known_chroms:
            warnings.warn(
                f"Chromosome {chrom!r} has no data in this database. "
                f"Available: {sorted(self._all_known_chroms)!r}",
                AfqueryWarning, stacklevel=3,
            )
            return []

        unique_variants = list(dict.fromkeys(variants))
        sample_bm = self._build_sample_bitmap(sf)
        return self._query_batch_inner(
            chrom, unique_variants, sample_bm,
            min_pass=sf.min_pass,
            min_observed=sf.min_observed,
            min_quality_evidence=sf.min_quality_evidence,
        )

    def query_region_multi(
        self,
        regions: list[tuple[str, int, int]],
        sf: SampleFilter,
    ) -> list[QueryResult]:
        """Query variants across multiple genomic regions (may span chromosomes).

        The sample bitmap is resolved once for all regions.
        Overlapping regions are deduplicated by variant key.
        Results are sorted in genomic order (chr1, chr2, …, chr22, chrX, chrY, chrM).

        Args:
            regions: List of ``(chrom, start, end)`` tuples, 1-based inclusive.
            sf: Sample filter.

        Returns:
            List of :class:`~afquery.models.QueryResult` objects.
        """
        if not regions:
            return []
        sample_bm = self._build_sample_bitmap(sf)
        seen: set[tuple[str, int, str, str]] = set()
        results: list[QueryResult] = []
        for chrom, start, end in regions:
            chrom = normalize_chrom(chrom)
            if chrom not in self._all_known_chroms:
                warnings.warn(
                    f"Chromosome {chrom!r} has no data in this database.",
                    AfqueryWarning, stacklevel=3,
                )
                continue
            for r in self._query_region_inner(
                chrom, start, end, sample_bm,
                min_pass=sf.min_pass,
                min_observed=sf.min_observed,
                min_quality_evidence=sf.min_quality_evidence,
            ):
                key = (r.variant.chrom, r.variant.pos, r.variant.ref, r.variant.alt)
                if key not in seen:
                    seen.add(key)
                    results.append(r)
        results.sort(key=lambda r: (
            CHROM_ORDER.get(r.variant.chrom, len(ALL_CHROMS)),
            r.variant.pos,
            r.variant.alt,
        ))
        return results

    def query_batch_multi(
        self,
        variants: list[tuple[str, int, str, str]],
        sf: SampleFilter,
    ) -> list[QueryResult]:
        """Query variants across multiple chromosomes, preserving input order.

        The sample bitmap is resolved once for all chromosomes.
        Duplicate input variants are deduplicated per chromosome (same as query_batch).

        Args:
            variants: List of ``(chrom, pos, ref, alt)`` tuples.
            sf: Sample filter.

        Returns:
            List of :class:`~afquery.models.QueryResult` objects in input order.
            Variants not found in the database are omitted.
        """
        if not variants:
            return []
        from collections import defaultdict
        sample_bm = self._build_sample_bitmap(sf)
        by_chrom: dict[str, list[tuple[int, tuple[int, str, str]]]] = defaultdict(list)
        for idx, (chrom, pos, ref, alt) in enumerate(variants):
            by_chrom[normalize_chrom(chrom)].append((idx, (pos, ref, alt)))
        tagged: list[tuple[int, QueryResult]] = []
        for chrom, idx_variants in by_chrom.items():
            if chrom not in self._all_known_chroms:
                warnings.warn(
                    f"Chromosome {chrom!r} has no data in this database.",
                    AfqueryWarning, stacklevel=3,
                )
                continue
            idxs = [i for i, _ in idx_variants]
            per_chrom = [v for _, v in idx_variants]
            unique = list(dict.fromkeys(per_chrom))
            results = self._query_batch_inner(
                chrom, unique, sample_bm,
                min_pass=sf.min_pass,
                min_observed=sf.min_observed,
                min_quality_evidence=sf.min_quality_evidence,
            )
            result_map = {(r.variant.pos, r.variant.ref, r.variant.alt): r for r in results}
            seen_on_chrom: set[tuple[int, str, str]] = set()
            for i, (pos, ref, alt) in zip(idxs, per_chrom):
                if (pos, ref, alt) in result_map and (pos, ref, alt) not in seen_on_chrom:
                    seen_on_chrom.add((pos, ref, alt))
                    tagged.append((i, result_map[(pos, ref, alt)]))
        tagged.sort(key=lambda x: x[0])
        return [r for _, r in tagged]

    def _query_region_inner(
        self,
        chrom: str,
        start: int,
        end: int,
        sample_bm: BitMap,
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
    ) -> list[QueryResult]:
        """Run a region query with a pre-built sample bitmap.

        Caller is responsible for normalization and chrom existence check.
        """
        parquet_glob = self._parquet_glob(chrom)
        if parquet_glob is None:
            return []

        escaped_glob = parquet_glob.replace("'", "''")

        cols = ", ".join(self._bitmap_cols(with_pos=True))
        con = duckdb.connect(config={"temp_directory": "/tmp"})
        rows = con.execute(
            f"SELECT {cols}"
            f" FROM read_parquet('{escaped_glob}') WHERE pos BETWEEN ? AND ?",
            [start, end],
        ).fetchall()
        con.close()

        if not rows:
            return []

        pos_data: dict[int, tuple[BitMap, int]] = {}
        for row in rows:
            pos = row[0]
            if pos not in pos_data:
                eligible, AN = self._compute_eligible(chrom, pos, sample_bm)
                pos_data[pos] = (eligible, AN)

        results = []
        for row in rows:
            pos, ref, alt = row[0], row[1], row[2]
            eligible, AN = pos_data[pos]
            if AN == 0:
                continue
            het_bm, hom_bm, fail_bm, filtered_bm, quality_pass_bm = self._unpack_bitmaps(row[3:])
            haploid_elig, diploid_elig = split_ploidy(
                eligible, self._male_bm, self._female_bm, chrom, pos, self._genome_build
            )
            het_elig = het_bm & eligible
            hom_elig = hom_bm & eligible
            AC = (len((het_elig | hom_elig) & haploid_elig)
                  + len(het_elig & diploid_elig)
                  + 2 * len(hom_elig & diploid_elig))
            N_HET = len(het_elig & diploid_elig)
            N_HOM_ALT = len(hom_elig & diploid_elig) + len((het_elig | hom_elig) & haploid_elig)
            AF = AC / AN if AN > 0 else None
            N_FAIL = len(fail_bm & eligible)
            no_cov_bm = self._compute_no_coverage_bm(
                eligible, het_bm, hom_bm, fail_bm,
                min_pass, min_observed,
                filtered_bm=filtered_bm,
                quality_pass_bm=quality_pass_bm,
                min_quality_evidence=min_quality_evidence,
            )
            N_NO_COVERAGE = len(no_cov_bm)
            N_HOM_REF = len(eligible) - N_HET - N_HOM_ALT - N_FAIL - N_NO_COVERAGE
            results.append(QueryResult(
                variant=VariantKey(chrom=chrom, pos=pos, ref=ref, alt=alt),
                AC=AC, AN=AN, AF=AF, n_samples_eligible=len(eligible),
                N_HET=N_HET, N_HOM_ALT=N_HOM_ALT, N_HOM_REF=N_HOM_REF,
                N_FAIL=N_FAIL, N_NO_COVERAGE=N_NO_COVERAGE,
            ))
        results.sort(key=lambda r: (r.variant.pos, r.variant.alt))
        return results

    def _query_batch_inner(
        self,
        chrom: str,
        unique_variants: list[tuple[int, str, str]],
        sample_bm: BitMap,
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
    ) -> list[QueryResult]:
        """Run a batch query with a pre-built sample bitmap.

        Caller is responsible for normalization, chrom existence check, and deduplication.
        """
        unique_positions = list(dict.fromkeys(pos for pos, _ref, _alt in unique_variants))

        pos_data: dict[int, tuple[BitMap, int]] = {}
        for pos in unique_positions:
            eligible, AN = self._compute_eligible(chrom, pos, sample_bm)
            pos_data[pos] = (eligible, AN)

        valid_positions = [p for p in unique_positions if pos_data[p][1] > 0]
        if not valid_positions:
            return []

        requested_variants: set[tuple[int, str, str]] = {
            (pos, ref, alt) for pos, ref, alt in unique_variants if pos in set(valid_positions)
        }

        parquet_glob = self._parquet_glob(chrom)
        if parquet_glob is None:
            return []

        escaped_glob = parquet_glob.replace("'", "''")

        col_list = self._bitmap_cols(with_pos=True)
        cols = ", ".join(col_list)
        cols_aliased = ", ".join(f"v.{c}" for c in col_list)

        con = duckdb.connect(config={"temp_directory": "/tmp"})
        if len(valid_positions) < BATCH_IN_THRESHOLD:
            placeholders = ", ".join("?" * len(valid_positions))
            rows = con.execute(
                f"SELECT {cols}"
                f" FROM read_parquet('{escaped_glob}') WHERE pos IN ({placeholders})",
                valid_positions,
            ).fetchall()
        else:
            con.execute("CREATE TEMP TABLE pos_filter (pos UINTEGER)")
            con.executemany("INSERT INTO pos_filter VALUES (?)", [(p,) for p in valid_positions])
            rows = con.execute(
                f"SELECT {cols_aliased}"
                f" FROM read_parquet('{escaped_glob}') v JOIN pos_filter pf ON v.pos = pf.pos"
            ).fetchall()
        con.close()

        results = []
        for row in rows:
            pos, ref, alt = row[0], row[1], row[2]
            if (pos, ref, alt) not in requested_variants:
                continue
            eligible, AN = pos_data[pos]
            het_bm, hom_bm, fail_bm, filtered_bm, quality_pass_bm = self._unpack_bitmaps(row[3:])
            haploid_elig, diploid_elig = split_ploidy(
                eligible, self._male_bm, self._female_bm, chrom, pos, self._genome_build
            )
            het_elig = het_bm & eligible
            hom_elig = hom_bm & eligible
            AC = (len((het_elig | hom_elig) & haploid_elig)
                  + len(het_elig & diploid_elig)
                  + 2 * len(hom_elig & diploid_elig))
            N_HET = len(het_elig & diploid_elig)
            N_HOM_ALT = len(hom_elig & diploid_elig) + len((het_elig | hom_elig) & haploid_elig)
            AF = AC / AN if AN > 0 else None
            N_FAIL = len(fail_bm & eligible)
            no_cov_bm = self._compute_no_coverage_bm(
                eligible, het_bm, hom_bm, fail_bm,
                min_pass, min_observed,
                filtered_bm=filtered_bm,
                quality_pass_bm=quality_pass_bm,
                min_quality_evidence=min_quality_evidence,
            )
            N_NO_COVERAGE = len(no_cov_bm)
            N_HOM_REF = len(eligible) - N_HET - N_HOM_ALT - N_FAIL - N_NO_COVERAGE
            results.append(QueryResult(
                variant=VariantKey(chrom=chrom, pos=pos, ref=ref, alt=alt),
                AC=AC, AN=AN, AF=AF, n_samples_eligible=len(eligible),
                N_HET=N_HET, N_HOM_ALT=N_HOM_ALT, N_HOM_REF=N_HOM_REF,
                N_FAIL=N_FAIL, N_NO_COVERAGE=N_NO_COVERAGE,
            ))
        results.sort(key=lambda r: (r.variant.pos, r.variant.alt))
        return results

    def query_region(
        self,
        chrom: str,
        start: int,
        end: int,
        sf: SampleFilter,
    ) -> list[QueryResult]:
        """Query all variants in a 1-based inclusive range [start, end]."""
        chrom = normalize_chrom(chrom)

        if chrom not in self._all_known_chroms:
            warnings.warn(
                f"Chromosome {chrom!r} has no data in this database. "
                f"Available: {sorted(self._all_known_chroms)!r}",
                AfqueryWarning, stacklevel=3,
            )
            return []

        sample_bm = self._build_sample_bitmap(sf)
        return self._query_region_inner(
            chrom, start, end, sample_bm,
            min_pass=sf.min_pass,
            min_observed=sf.min_observed,
            min_quality_evidence=sf.min_quality_evidence,
        )

    def variant_info(self, params: QueryParams) -> list[SampleCarrier]:
        """Return all samples carrying the variant at params.chrom:params.pos.

        Each returned :class:`~afquery.models.SampleCarrier` includes the
        sample's name, sex, technology, phenotype codes, genotype (het/hom/alt),
        and FILTER status.

        Args:
            params: Query parameters (chrom, pos, optional ref/alt, sample filter).

        Returns:
            List of :class:`~afquery.models.SampleCarrier`, sorted by
            ``sample_id``.  Empty if the variant is absent or no eligible
            carrier exists.
        """
        chrom = normalize_chrom(params.chrom)
        pos = params.pos

        if chrom not in self._all_known_chroms:
            warnings.warn(
                f"Chromosome {chrom!r} has no data in this database. "
                f"Available: {sorted(self._all_known_chroms)!r}",
                AfqueryWarning, stacklevel=3,
            )
            return []

        sample_bm = self._build_sample_bitmap(params.filter)
        parquet_path = self._parquet_path(chrom, pos)
        if parquet_path is None:
            return []

        cols = ", ".join(self._bitmap_cols(with_pos=True))
        con = duckdb.connect(config={"temp_directory": "/tmp"})
        rows = con.execute(
            f"SELECT {cols}"
            " FROM read_parquet(?) WHERE pos = ?",
            [str(parquet_path), pos],
        ).fetchall()
        con.close()

        if not rows:
            return []

        if params.ref is not None:
            rows = [r for r in rows if r[1] == params.ref]
        if params.alt is not None:
            rows = [r for r in rows if r[2] == params.alt]
        if not rows:
            return []

        # Warn if multiple alleles found without explicit ref/alt filter
        if len(rows) > 1 and params.ref is None and params.alt is None:
            alts = ", ".join(r[2] for r in rows)
            warnings.warn(
                f"Multiple alleles found at {chrom}:{pos} ({alts}). "
                "Use --ref/--alt to restrict to a specific variant.",
                AfqueryWarning, stacklevel=3,
            )

        # Compute eligible (BED-aware) for no_coverage assessment
        eligible, _AN = self._compute_eligible(chrom, pos, sample_bm)
        sf = params.filter

        carriers: list[SampleCarrier] = []
        for row in rows:
            row_pos, ref, alt = row[0], row[1], row[2]
            het_bm, hom_bm, fail_bm, filtered_bm, quality_pass_bm = self._unpack_bitmaps(row[3:])
            het_elig = het_bm & sample_bm
            hom_elig = hom_bm & sample_bm
            fail_elig = fail_bm & sample_bm
            no_cov_bm = self._compute_no_coverage_bm(
                eligible, het_bm, hom_bm, fail_bm,
                sf.min_pass, sf.min_observed,
                filtered_bm=filtered_bm,
                quality_pass_bm=quality_pass_bm,
                min_quality_evidence=sf.min_quality_evidence,
            )
            no_cov_elig = no_cov_bm & sample_bm

            seen: set[int] = set()
            for sid in sorted(hom_elig):
                seen.add(sid)
                is_fail = sid in fail_elig
                carriers.append(self._make_carrier(sid, "alt" if is_fail else "hom", not is_fail))
            for sid in sorted(het_elig):
                if sid not in seen:
                    seen.add(sid)
                    is_fail = sid in fail_elig
                    carriers.append(self._make_carrier(sid, "alt" if is_fail else "het", not is_fail))
            for sid in sorted(fail_elig):
                if sid not in seen:
                    seen.add(sid)
                    carriers.append(self._make_carrier(sid, "alt", False))
            for sid in sorted(no_cov_elig):
                if sid not in seen:
                    seen.add(sid)
                    carriers.append(self._make_carrier(sid, "no_coverage", None))

        carriers.sort(key=lambda c: c.sample_id)
        return carriers

    def _make_carrier(
        self,
        sample_id: int,
        genotype: str,
        filter_pass: bool | None,
    ) -> SampleCarrier:
        s = self._samples_by_id[sample_id]
        tech_name = self._tech_map[s.tech_id].tech_name
        phenotypes = self._sample_phenotypes.get(sample_id, [])
        return SampleCarrier(
            sample_id=sample_id,
            sample_name=s.sample_name,
            sex=s.sex,
            tech_name=tech_name,
            phenotypes=sorted(phenotypes),
            genotype=genotype,
            filter_pass=filter_pass,
        )

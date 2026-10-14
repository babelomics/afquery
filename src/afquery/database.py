import json
import sqlite3
from pathlib import Path

from .models import QueryParams, QueryResult, SampleCarrier, SampleFilter
from .query import QueryEngine


class Database:
    def __init__(self, db_path: str):
        self._path = Path(db_path)
        self._engine = QueryEngine(db_path)
        self._manifest = json.loads((self._path / "manifest.json").read_text())

    def _reload(self) -> None:
        """Reinitialize engine and manifest after a mutating operation."""
        self._engine = QueryEngine(str(self._path))
        self._manifest = json.loads((self._path / "manifest.json").read_text())

    def _make_filter(
        self, phenotype, sex, tech=None,
        min_pass: int = 0, min_observed: int = 0, min_quality_evidence: int = 0,
        fail_as_ref: bool = False,
    ) -> SampleFilter:
        return SampleFilter.parse(
            phenotype_tokens=phenotype or [],
            tech_tokens=tech or [],
            sex=sex,
            min_pass=min_pass,
            min_observed=min_observed,
            min_quality_evidence=min_quality_evidence,
            fail_as_ref=fail_as_ref,
        )

    def query(
        self,
        chrom: str,
        pos: int,
        phenotype: list[str] | None = None,
        sex: str = "both",
        ref: str | None = None,
        alt: str | None = None,
        tech: list[str] | None = None,
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
        fail_as_ref: bool = False,
    ) -> list[QueryResult]:
        sf = self._make_filter(
            phenotype, sex, tech, min_pass, min_observed, min_quality_evidence, fail_as_ref,
        )
        params = QueryParams(chrom=chrom, pos=pos, filter=sf, ref=ref, alt=alt)
        return self._engine.query(params)

    def query_batch(
        self,
        chrom: str,
        variants: list[tuple[int, str, str]],
        phenotype: list[str] | None = None,
        sex: str = "both",
        tech: list[str] | None = None,
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
        fail_as_ref: bool = False,
    ) -> list[QueryResult]:
        sf = self._make_filter(
            phenotype, sex, tech, min_pass, min_observed, min_quality_evidence, fail_as_ref,
        )
        return self._engine.query_batch(chrom, variants, sf)

    def query_batch_multi(
        self,
        variants: list[tuple[str, int, str, str]],
        phenotype: list[str] | None = None,
        sex: str = "both",
        tech: list[str] | None = None,
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
        fail_as_ref: bool = False,
    ) -> list[QueryResult]:
        """Query variants across multiple chromosomes, preserving input order.

        Duplicate input variants are deduplicated per chromosome.
        Chromosome names are normalized (``"1"`` and ``"chr1"`` are equivalent).

        Args:
            variants: List of ``(chrom, pos, ref, alt)`` tuples.
            phenotype: Phenotype filter tokens. Use ``"^CODE"`` prefix to exclude.
            sex: ``"both"`` (default), ``"male"``, or ``"female"``.
            tech: Technology filter tokens. Use ``"^TECH"`` prefix to exclude.

        Returns:
            List of :class:`~afquery.models.QueryResult` objects in input order
            (by original index). Variants not found in the database are omitted.
        """
        sf = self._make_filter(
            phenotype, sex, tech, min_pass, min_observed, min_quality_evidence, fail_as_ref,
        )
        return self._engine.query_batch_multi(variants, sf)

    def query_region(
        self,
        chrom: str,
        start: int,
        end: int,
        phenotype: list[str] | None = None,
        sex: str = "both",
        tech: list[str] | None = None,
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
        fail_as_ref: bool = False,
    ) -> list[QueryResult]:
        sf = self._make_filter(
            phenotype, sex, tech, min_pass, min_observed, min_quality_evidence, fail_as_ref,
        )
        return self._engine.query_region(chrom, start, end, sf)

    def variant_info(
        self,
        chrom: str,
        pos: int,
        ref: str | None = None,
        alt: str | None = None,
        phenotype: list[str] | None = None,
        sex: str = "both",
        tech: list[str] | None = None,
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
    ) -> list[SampleCarrier]:
        """Return all samples carrying a variant, with their metadata.

        Each element includes sample name, sex, technology, phenotype codes,
        genotype (``"het"``, ``"hom"``, or ``"alt"`` for FILTER≠PASS samples),
        and FILTER status.

        Args:
            chrom: Chromosome (e.g. ``"chr1"`` or ``"1"``).
            pos: Position, 1-based.
            ref: Reference allele filter. If omitted, all alleles at *pos* are returned.
            alt: Alternate allele filter. If omitted, all alleles at *pos* are returned.
            phenotype: Phenotype filter tokens. Use ``"^CODE"`` prefix to exclude.
            sex: ``"both"`` (default), ``"male"``, or ``"female"``.
            tech: Technology filter tokens. Use ``"^TECH"`` prefix to exclude.

        Returns:
            List of :class:`~afquery.models.SampleCarrier` sorted by
            ``sample_id``.  Empty list if the variant is absent or no eligible
            carrier exists.
        """
        sf = self._make_filter(phenotype, sex, tech, min_pass, min_observed, min_quality_evidence)
        params = QueryParams(chrom=chrom, pos=pos, filter=sf, ref=ref, alt=alt)
        return self._engine.variant_info(params)

    def query_region_multi(
        self,
        regions: list[tuple[str, int, int]],
        phenotype: list[str] | None = None,
        sex: str = "both",
        tech: list[str] | None = None,
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
        fail_as_ref: bool = False,
    ) -> list[QueryResult]:
        """Query variants across multiple genomic regions (may span chromosomes).

        Args:
            regions: List of ``(chrom, start, end)`` tuples, 1-based inclusive.
            phenotype: Phenotype filter tokens. Use ``"^CODE"`` prefix to exclude.
            sex: ``"both"`` (default), ``"male"``, or ``"female"``.
            tech: Technology filter tokens. Use ``"^TECH"`` prefix to exclude.

        Returns:
            List of :class:`~afquery.models.QueryResult` objects sorted in
            genomic order (chr1, chr2, …, chr22, chrX, chrY, chrM).
            Overlapping regions are deduplicated.
        """
        sf = self._make_filter(
            phenotype, sex, tech, min_pass, min_observed, min_quality_evidence, fail_as_ref,
        )
        return self._engine.query_region_multi(regions, sf)

    def dump(
        self,
        output=None,
        phenotype: list[str] | None = None,
        sex: str = "both",
        tech: list[str] | None = None,
        by_sex: bool = False,
        by_tech: bool = False,
        by_phenotype: list[str] | None = None,
        all_groups: bool = False,
        chrom: str | None = None,
        start: int | None = None,
        end: int | None = None,
        n_workers: int | None = None,
        include_ac_zero: bool = False,
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
        fail_as_ref: bool = False,
    ) -> dict:
        from .dump import dump_database, _build_groups
        base_sf = self._make_filter(
            phenotype, sex, tech, min_pass, min_observed, min_quality_evidence, fail_as_ref,
        )
        groups = _build_groups(
            self._engine, base_sf, by_sex, by_tech,
            by_phenotype or [], all_groups,
        )
        return dump_database(
            self._engine, output, base_sf, groups,
            chrom, start, end, n_workers, include_ac_zero,
        )

    def annotate_vcf(
        self,
        input_vcf: str,
        output_vcf: str,
        phenotype: list[str] | None = None,
        sex: str = "both",
        n_workers: int | None = None,
        tech: list[str] | None = None,
        min_pass: int = 0,
        min_observed: int = 0,
        min_quality_evidence: int = 0,
        fail_as_ref: bool = False,
    ) -> dict:
        sf = self._make_filter(
            phenotype, sex, tech, min_pass, min_observed, min_quality_evidence, fail_as_ref,
        )
        from .annotate import annotate_vcf as _annotate
        return _annotate(self._engine, input_vcf, output_vcf, sf, n_workers=n_workers)

    def add_samples(
        self,
        manifest_path: str,
        threads: int = 8,
        tmp_dir: str | None = None,
        bed_dir: str | None = None,
        genome_build: str | None = None,
    ) -> dict:
        from .preprocess.update import add_samples as _add
        result = _add(
            str(self._path), manifest_path, threads=threads,
            tmp_dir=tmp_dir, bed_dir=bed_dir, genome_build=genome_build,
        )
        self._reload()
        return result

    def remove_samples(self, sample_names: list[str]) -> dict:
        from .preprocess.update import remove_samples as _remove
        result = _remove(str(self._path), sample_names)
        self._reload()
        return result

    def update_sample_metadata(
        self,
        updates: list[dict],
        operator_note: str | None = None,
    ) -> list[dict]:
        """Update ``sex`` and/or ``phenotype_codes`` for one or more samples.

        Args:
            updates: List of dicts with keys ``sample_name``, ``field``, and
                ``new_value``.  See
                :func:`~afquery.preprocess.update.update_sample_metadata` for
                full details.
            operator_note: Optional free-text note appended to each changelog
                entry.

        Returns:
            List of changelog entry dicts created (one per field change).
        """
        from .preprocess.update import update_sample_metadata as _update
        result = _update(str(self._path), updates, operator_note=operator_note)
        self._reload()
        return result

    def compact(self) -> dict:
        from .preprocess.compact import compact_database
        result = compact_database(self._path)
        self._reload()
        return result

    def check(self) -> list:
        from .preprocess.update import check_database
        return check_database(str(self._path))

    def info(self) -> dict:
        con = sqlite3.connect(str(self._path / "metadata.sqlite"))
        try:
            total = con.execute("SELECT COUNT(*) FROM samples").fetchone()[0]
            by_sex: dict[str, int] = {}
            for row in con.execute(
                "SELECT sex, COUNT(*) FROM samples GROUP BY sex"
            ).fetchall():
                by_sex[row[0]] = row[1]
            by_tech: dict[str, int] = {}
            for row in con.execute(
                "SELECT t.tech_name, COUNT(*) FROM samples s"
                " JOIN technologies t ON s.tech_id = t.tech_id GROUP BY t.tech_name"
            ).fetchall():
                by_tech[row[0]] = row[1]
            by_phenotype: dict[str, int] = {}
            for row in con.execute(
                "SELECT phenotype_code, COUNT(*) FROM sample_phenotype GROUP BY phenotype_code"
            ).fetchall():
                by_phenotype[row[0]] = row[1]
            tables = {r[0] for r in con.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()}
            changelog_recent: list[dict] = []
            if "changelog" in tables:
                for r in con.execute(
                    "SELECT event_id, event_type, event_time, sample_names, notes"
                    " FROM changelog ORDER BY event_id DESC LIMIT 5"
                ).fetchall():
                    changelog_recent.append({
                        "event_id": r[0],
                        "event_type": r[1],
                        "event_time": r[2],
                        "sample_names": json.loads(r[3]) if r[3] else None,
                        "notes": r[4],
                    })
        finally:
            con.close()

        m = dict(self._manifest)
        return {
            "db_path": str(self._path),
            "db_version": m.get("db_version", "unknown"),
            "genome_build": m["genome_build"],
            "schema_version": m["schema_version"],
            "created_at": m.get("created_at"),
            "updated_at": m.get("updated_at"),
            "sample_count": total,
            "by_sex": by_sex,
            "by_tech": by_tech,
            "by_phenotype": by_phenotype,
            "changelog_recent": changelog_recent,
        }

    def list_samples(self) -> list[dict]:
        """Return all samples with full metadata."""
        con = sqlite3.connect(str(self._path / "metadata.sqlite"))
        try:
            rows = con.execute(
                "SELECT s.sample_id, s.sample_name, s.sex, t.tech_name,"
                " s.vcf_path, s.ingested_at"
                " FROM samples s JOIN technologies t ON s.tech_id = t.tech_id"
                " ORDER BY s.sample_id"
            ).fetchall()
            pheno_map: dict[int, list[str]] = {}
            for sid, code in con.execute(
                "SELECT sample_id, phenotype_code FROM sample_phenotype ORDER BY sample_id"
            ).fetchall():
                pheno_map.setdefault(sid, []).append(code)
        finally:
            con.close()

        return [
            {
                "sample_id": r[0],
                "sample_name": r[1],
                "sex": r[2],
                "tech": r[3],
                "phenotypes": pheno_map.get(r[0], []),
                "vcf_path": r[4],
                "ingested_at": r[5],
            }
            for r in rows
        ]

    def changelog(self, limit: int | None = None) -> list[dict]:
        """Return changelog entries ordered by event_id descending."""
        con = sqlite3.connect(str(self._path / "metadata.sqlite"))
        try:
            q = "SELECT event_id, event_type, event_time, sample_names, notes FROM changelog ORDER BY event_id DESC"
            if limit is not None:
                q += f" LIMIT {int(limit)}"
            rows = con.execute(q).fetchall()
        finally:
            con.close()

        return [
            {
                "event_id": r[0],
                "event_type": r[1],
                "event_time": r[2],
                "sample_names": json.loads(r[3]) if r[3] else None,
                "notes": r[4],
            }
            for r in rows
        ]

    def set_db_version(self, version: str) -> None:
        """Set the db_version label in manifest.json atomically."""
        import os
        from datetime import datetime, timezone
        manifest_path = self._path / "manifest.json"
        try:
            manifest = json.loads(manifest_path.read_text())
        except (FileNotFoundError, json.JSONDecodeError):
            manifest = {}
        old_version = manifest.get("db_version", "unknown")
        manifest["db_version"] = version
        tmp_path = str(manifest_path) + ".tmp"
        with open(tmp_path, "w") as f:
            json.dump(manifest, f, indent=2)
        os.replace(tmp_path, str(manifest_path))
        self._reload()
        # Log to changelog
        con = sqlite3.connect(str(self._path / "metadata.sqlite"))
        try:
            con.execute(
                "INSERT INTO changelog (event_type, event_time, sample_names, notes) VALUES (?, ?, ?, ?)",
                (
                    "version_set",
                    datetime.now(timezone.utc).isoformat(),
                    None,
                    f"version changed from {old_version} to {version}",
                ),
            )
            con.commit()
        finally:
            con.close()

    def get_all_phenotypes(self) -> list[str]:
        """Return all available phenotype codes in the database."""
        import sqlite3
        con = sqlite3.connect(str(self._path / "metadata.sqlite"))
        phenotypes = [
            row[0] for row in con.execute(
                "SELECT DISTINCT phenotype_code FROM sample_phenotype ORDER BY phenotype_code"
            ).fetchall()
        ]
        con.close()
        return phenotypes

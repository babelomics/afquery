import json
import logging
import sys

import click

from afquery import __version__

from .database import Database


def _expand_tokens(values: tuple[str, ...]) -> list[str]:
    """Expand comma-separated tokens into a flat list."""
    result = []
    for v in values:
        result.extend(t.strip() for t in v.split(",") if t.strip())
    return result


def _configure_logging(verbose: bool) -> None:
    """Configure logging to stderr with optional DEBUG verbosity."""
    level = logging.DEBUG if verbose else logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter("%(message)s"))
    root = logging.getLogger("afquery")
    root.setLevel(level)
    root.handlers.clear()
    root.addHandler(handler)
    root.propagate = False


def _parse_locus(locus: str) -> tuple[str, int]:
    """Parse 'chrom:pos' → (chrom, pos). Raises UsageError on malformed input."""
    try:
        chrom, pos_str = locus.rsplit(":", 1)
        return chrom, int(pos_str)
    except ValueError:
        raise click.UsageError(f"--locus must be CHROM:POS (e.g., chr1:925952), got: {locus!r}")


def _parse_region(region: str) -> tuple[str, int, int]:
    """Parse 'chrom:start-end' → (chrom, start, end). Raises UsageError on malformed input."""
    try:
        chrom, span = region.rsplit(":", 1)
        start, end = map(int, span.split("-", 1))
        return chrom, start, end
    except ValueError:
        raise click.UsageError(f"--region must be CHROM:START-END (e.g., chr1:900000-1000000), got: {region!r}")


def _parse_variants_file(path: str) -> list[tuple[str, int, str, str]]:
    """Parse a TSV file with columns: chrom pos [ref [alt]] (no header, whitespace-separated)."""
    variants = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            chrom = parts[0]
            pos = int(parts[1])
            ref = parts[2] if len(parts) > 2 else ""
            alt = parts[3] if len(parts) > 3 else ""
            variants.append((chrom, pos, ref, alt))
    return variants


def _print_results(results, fmt: str) -> None:
    """Print query results in the requested format."""
    if fmt == "json":
        out = []
        for r in results:
            entry = {
                "chrom": r.variant.chrom, "pos": r.variant.pos,
                "ref": r.variant.ref, "alt": r.variant.alt,
                "AC": r.AC, "AN": r.AN, "AF": r.AF, "n_eligible": r.n_samples_eligible,
                "N_HET": r.N_HET, "N_HOM_ALT": r.N_HOM_ALT, "N_HOM_REF": r.N_HOM_REF,
                "N_FAIL": r.N_FAIL, "N_NO_COVERAGE": r.N_NO_COVERAGE,
            }
            out.append(entry)
        click.echo(json.dumps(out, indent=2))
    elif fmt == "tsv":
        header = (
            "chrom\tpos\tref\talt\tAC\tAN\tAF\tn_eligible\t"
            "N_HET\tN_HOM_ALT\tN_HOM_REF\tN_FAIL\tN_NO_COVERAGE"
        )
        click.echo(header)
        for r in results:
            af = f"{r.AF:.6f}" if r.AF is not None else "NA"
            line = (
                f"{r.variant.chrom}\t{r.variant.pos}\t{r.variant.ref}\t"
                f"{r.variant.alt}\t{r.AC}\t{r.AN}\t{af}\t{r.n_samples_eligible}\t"
                f"{r.N_HET}\t{r.N_HOM_ALT}\t{r.N_HOM_REF}\t{r.N_FAIL}\t{r.N_NO_COVERAGE}"
            )
            click.echo(line)
    else:  # text
        if not results:
            click.echo("No variants found for the given filters.")
            return
        for r in results:
            af = f"{r.AF:.4f}" if r.AF is not None else "NA"
            click.echo(
                f"{r.variant.chrom}:{r.variant.pos} {r.variant.ref}>{r.variant.alt}  "
                f"AC={r.AC}  AN={r.AN}  AF={af}  n_eligible={r.n_samples_eligible}  "
                f"N_HET={r.N_HET}  N_HOM_ALT={r.N_HOM_ALT}  N_HOM_REF={r.N_HOM_REF}  "
                f"N_FAIL={r.N_FAIL}  N_NO_COVERAGE={r.N_NO_COVERAGE}"
            )


def _print_carriers(carriers, variant_key, fmt: str) -> None:
    """Print variant_info results in the requested format."""
    if fmt == "json":
        out = {
            "variant": {
                "chrom": variant_key[0], "pos": variant_key[1],
                "ref": variant_key[2], "alt": variant_key[3],
            },
            "samples": [
                {
                    "sample_id": c.sample_id,
                    "sample_name": c.sample_name,
                    "sex": c.sex,
                    "tech": c.tech_name,
                    "phenotypes": c.phenotypes,
                    "genotype": c.genotype,
                    "filter": (
                        None if c.filter_pass is None
                        else ("PASS" if c.filter_pass else "FAIL")
                    ),
                }
                for c in carriers
            ],
        }
        click.echo(json.dumps(out, indent=2))
    elif fmt == "tsv":
        click.echo("sample_id\tsample_name\tsex\ttech\tphenotypes\tgenotype\tfilter")
        for c in carriers:
            filter_cell = (
                "" if c.filter_pass is None
                else ("PASS" if c.filter_pass else "FAIL")
            )
            click.echo(
                f"{c.sample_id}\t{c.sample_name}\t{c.sex}\t{c.tech_name}\t"
                f"{','.join(c.phenotypes)}\t{c.genotype}\t"
                f"{filter_cell}"
            )
    else:  # text
        if not carriers:
            click.echo("No carriers found for the given filters.")
            return
        headers = ["sample_id", "sample_name", "sex", "tech", "phenotypes", "genotype", "filter"]
        rows = [
            [
                str(c.sample_id), c.sample_name, c.sex, c.tech_name,
                ",".join(c.phenotypes) if c.phenotypes else ".",
                c.genotype,
                (
                    "-" if c.filter_pass is None
                    else ("PASS" if c.filter_pass else "FAIL")
                ),
            ]
            for c in carriers
        ]
        widths = [max(len(h), max((len(r[i]) for r in rows), default=0)) for i, h in enumerate(headers)]
        fmt_row = "  ".join(f"{{:<{w}}}" for w in widths)
        click.echo(fmt_row.format(*headers))
        click.echo("  ".join("-" * w for w in widths))
        for row in rows:
            click.echo(fmt_row.format(*row))


_CLI_HELP = f"""\
AFQuery v{__version__} — bitmap-indexed allele frequency engine for local genomic cohorts.

Enables fast AC/AN/AF queries on user-defined subcohorts (phenotype, sex,
technology) without rescanning VCFs.
"""


@click.group(help=_CLI_HELP)
@click.version_option(version=__version__, prog_name="afquery", message="%(prog)s %(version)s")
def cli():
    pass


@cli.command()
@click.option("--db",       required=True, help="Path to database directory.")
@click.option("--locus",     default=None, help="Single position as CHROM:POS (e.g., chr1:925952).")
@click.option("--region",    default=None, help="Genomic region as CHROM:START-END (e.g., chr1:900000-1000000).")
@click.option("--from-file", default=None, type=click.Path(exists=True), help="TSV file with columns: chrom pos [ref [alt]] (batch query, no header).")
@click.option("--phenotype", multiple=True, help="Phenotype to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--sex",   default="both", type=click.Choice(["male", "female", "both"]), help="Restrict to specified sex. Options: male, female, both. (default: both)")
@click.option("--ref",   default=None, help="Filter to specific reference allele (only for --locus).")
@click.option("--alt",   default=None, help="Filter to specific alternate allele (only for --locus).")
@click.option("--tech",  multiple=True, help="Technology filter to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--format", "fmt", default="text", type=click.Choice(["text", "json", "tsv"]), help="Output format. Options: text, json, tsv. (default: text)")
@click.option("--min-pass", default=0, type=int, help="Min PASS carriers (het|hom) per WES tech for hom-ref to be assumed. Non-carriers move to N_NO_COVERAGE if tech falls below the threshold. (default: 0 = disabled)")
@click.option("--min-observed", default=0, type=int, help="Min any-VCF entries (het|hom|fail) per WES tech for hom-ref to be assumed. (default: 0 = disabled)")
@click.option("--min-quality-evidence", default=0, type=int, help="Min quality-passing carriers per WES tech (Phase 2 DBs only — created with --min-dp/--min-gq). (default: 0 = disabled)")
@click.option("--no-warn", is_flag=True, default=False, help="Suppress AfqueryWarning messages.")
def query(db, locus, region, from_file, phenotype, sex, ref, alt, tech, fmt,
          min_pass, min_observed, min_quality_evidence, no_warn):
    """Query allele frequencies.

    Exactly one of --locus, --region, or --from-file must be provided:

    \b
      --locus chr1:925952            single position
      --region chr1:900000-1000000   all variants in [start, end]
      --from-file variants.tsv       batch from TSV (chrom pos ref alt, no header)
    """
    if no_warn:
        import warnings
        from .models import AfqueryWarning
        warnings.filterwarnings("ignore", category=AfqueryWarning)

    modes = sum(x is not None for x in [locus, region, from_file])
    if modes == 0:
        raise click.UsageError("One of --locus, --region, or --from-file is required.")
    if modes > 1:
        raise click.UsageError("--locus, --region, and --from-file are mutually exclusive.")

    database = Database(db)

    if locus is not None:
        chrom, pos = _parse_locus(locus)
        results = database.query(
            chrom=chrom, pos=pos,
            phenotype=_expand_tokens(phenotype), sex=sex,
            ref=ref, alt=alt, tech=_expand_tokens(tech),
            min_pass=min_pass, min_observed=min_observed,
            min_quality_evidence=min_quality_evidence,
        )
    elif region is not None:
        chrom, start, end = _parse_region(region)
        results = database.query_region(
            chrom=chrom, start=start, end=end,
            phenotype=_expand_tokens(phenotype), sex=sex, tech=_expand_tokens(tech),
            min_pass=min_pass, min_observed=min_observed,
            min_quality_evidence=min_quality_evidence,
        )
    else:
        variants = _parse_variants_file(from_file)
        results = database.query_batch_multi(
            variants=variants,
            phenotype=_expand_tokens(phenotype), sex=sex, tech=_expand_tokens(tech),
            min_pass=min_pass, min_observed=min_observed,
            min_quality_evidence=min_quality_evidence,
        )

    _print_results(results, fmt)


@cli.command("variant-info")
@click.option("--db",       required=True, help="Path to database directory.")
@click.option("--locus",    required=True, help="Position as CHROM:POS (e.g., chr1:925952).")
@click.option("--ref",      default=None, help="Filter to specific reference allele.")
@click.option("--alt",      default=None, help="Filter to specific alternate allele.")
@click.option("--phenotype", multiple=True, help="Phenotype to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--sex",   default="both", type=click.Choice(["male", "female", "both"]), help="Restrict to specified sex. Options: male, female, both. (default: both)")
@click.option("--tech",  multiple=True, help="Technology filter. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--format", "fmt", default="text", type=click.Choice(["text", "json", "tsv"]), help="Output format. Options: text, json, tsv. (default: text)")
@click.option("--min-pass", default=0, type=int, help="Min PASS carriers per WES tech for hom-ref to be assumed (samples below threshold show as 'no_coverage' genotype). (default: 0 = disabled)")
@click.option("--min-observed", default=0, type=int, help="Min any-VCF entries per WES tech (default: 0 = disabled).")
@click.option("--min-quality-evidence", default=0, type=int, help="Min quality-passing carriers per WES tech (Phase 2 DBs only). (default: 0 = disabled)")
@click.option("--no-warn", is_flag=True, default=False, help="Suppress AfqueryWarning messages.")
def variant_info_cmd(db, locus, ref, alt, phenotype, sex, tech, fmt,
                     min_pass, min_observed, min_quality_evidence, no_warn):
    """List samples carrying a specific variant, with their metadata.

    Returns one row per carrier sample with genotype (het/hom/alt) and FILTER
    status (PASS/FAIL).  Use --ref/--alt to restrict to a specific allele when
    multiple alleles share the same position.

    \b
    Examples:
      afquery variant-info --db ./db --locus chr1:925952
      afquery variant-info --db ./db --locus chr1:925952 --ref A --alt T
      afquery variant-info --db ./db --locus chrX:2700000 --sex female --format tsv
    """
    if no_warn:
        import warnings
        from .models import AfqueryWarning
        warnings.filterwarnings("ignore", category=AfqueryWarning)

    chrom, pos = _parse_locus(locus)
    database = Database(db)
    carriers = database.variant_info(
        chrom=chrom, pos=pos,
        ref=ref, alt=alt,
        phenotype=_expand_tokens(phenotype), sex=sex,
        tech=_expand_tokens(tech),
        min_pass=min_pass, min_observed=min_observed,
        min_quality_evidence=min_quality_evidence,
    )

    _variant_key = (chrom, pos, ref or ".", alt or ".")
    _print_carriers(carriers, _variant_key, fmt)


@cli.command()
@click.option("--db",      required=True, help="Path to database directory.")
@click.option("--input",   "input_vcf",  required=True, help="Input VCF file (plain or .gz).")
@click.option("--output",  "output_vcf", required=True, help="Output annotated VCF file.")
@click.option("--phenotype",   multiple=True, help="Phenotype code to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--sex",     default="both", type=click.Choice(["male", "female", "both"]), help="Restrict to specified sex. Options: male, female, both. (default: both)")
@click.option("--tech", multiple=True, help="Technology filter to include. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--threads", default=None, type=int,
              help="Number of worker threads for parallel annotation. (default: all available CPU cores)")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output with per-item progress. (default: false)")
@click.option("--min-pass", default=0, type=int, help="Min PASS carriers per WES tech for hom-ref to be assumed (default: 0).")
@click.option("--min-observed", default=0, type=int, help="Min any-VCF entries per WES tech (default: 0).")
@click.option("--min-quality-evidence", default=0, type=int, help="Min quality-passing carriers per WES tech (Phase 2 DBs only) (default: 0).")
@click.option("--no-warn", is_flag=True, default=False, help="Suppress AfqueryWarning messages.")
def annotate(db, input_vcf, output_vcf, phenotype, sex, tech, threads, verbose,
             min_pass, min_observed, min_quality_evidence, no_warn):
    """Annotate a VCF with allele frequency INFO fields.

    The following INFO fields are added to each variant:

    \b
      AFQUERY_AC          allele count (per ALT)
      AFQUERY_AN          allele number (per site)
      AFQUERY_AF          allele frequency (per ALT)
      AFQUERY_N_HET       heterozygous carrier count (per ALT)
      AFQUERY_N_HOM_ALT   homozygous alt count (per ALT)
      AFQUERY_N_HOM_REF   homozygous ref count (per ALT)
      AFQUERY_N_FAIL      samples with FILTER!=PASS (per site)
      AFQUERY_N_NO_COVERAGE  eligible samples lacking coverage evidence (per ALT)
    """
    if no_warn:
        import warnings
        from .models import AfqueryWarning
        warnings.filterwarnings("ignore", category=AfqueryWarning)
    _configure_logging(verbose)
    database = Database(db)
    stats = database.annotate_vcf(
        input_vcf, output_vcf,
        phenotype=_expand_tokens(phenotype), sex=sex,
        tech=_expand_tokens(tech), n_workers=threads,
        min_pass=min_pass, min_observed=min_observed,
        min_quality_evidence=min_quality_evidence,
    )
    click.echo(
        f"Annotated {stats['n_annotated']} variants "
        f"({stats['n_uncovered']} uncovered, {stats['n_variants']} total)."
    )


@cli.command()
@click.option("--db",       required=True, help="Path to database directory.")
@click.option("--output", "-o", default=None, help="Output CSV file path. (default: stdout)")
@click.option("--chrom",   default=None, help="Restrict export to this chromosome.")
@click.option("--start",   default=None, type=int, help="1-based start position (inclusive). Requires --chrom.")
@click.option("--end",     default=None, type=int, help="1-based end position (inclusive). Requires --chrom.")
@click.option("--phenotype", multiple=True, help="Phenotype filter. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--sex",     default="both", type=click.Choice(["male", "female", "both"]), help="Restrict to specified sex. Options: male, female, both. (default: both)")
@click.option("--tech",    multiple=True, help="Technology filter. Repeatable; comma-separated or multiple flags. Use ^ prefix to exclude. (default: include all samples)")
@click.option("--by-sex",  is_flag=True, help="Disaggregate output by sex (adds AC_male/AC_female columns). (default: false)")
@click.option("--by-tech", is_flag=True, help="Disaggregate output by technology (adds AC_<tech> columns). (default: false)")
@click.option("--by-phenotype", multiple=True, help="Disaggregate by specific phenotype codes. Repeatable; comma-separated or multiple flags.")
@click.option("--all-groups", is_flag=True, help="Disaggregate by ALL sexes x technologies x phenotypes (Cartesian product). WARNING: can generate a very large number of columns. (default: false)")
@click.option("--threads", default=None, type=int, help="Number of worker threads for parallel export. (default: all available CPU cores)")
@click.option("--all-variants", is_flag=True, help="Include variants with AC=0 (covered but not observed). WARNING: may produce very large output on whole-genome databases. (default: false)")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output with per-item progress. (default: false)")
@click.option("--min-pass", default=0, type=int, help="Min PASS carriers per WES tech for hom-ref to be assumed (default: 0).")
@click.option("--min-observed", default=0, type=int, help="Min any-VCF entries per WES tech (default: 0).")
@click.option("--min-quality-evidence", default=0, type=int, help="Min quality-passing carriers per WES tech (Phase 2 DBs only) (default: 0).")
def dump(db, output, chrom, start, end, phenotype, sex, tech,
         by_sex, by_tech, by_phenotype, all_groups, threads, all_variants, verbose,
         min_pass, min_observed, min_quality_evidence):
    """Export allele frequencies to CSV.

    By default only variants with AC > 0 are exported. Use --all-variants to
    include positions that are covered but have no observed carriers (AC=0).

    Supports disaggregation by sex (--by-sex), technology (--by-tech),
    phenotype (--by-phenotype), or all combinations (--all-groups).
    """
    _configure_logging(verbose)

    if (start is not None or end is not None) and chrom is None:
        raise click.UsageError("--start/--end require --chrom.")
    if start is not None and end is not None and start > end:
        raise click.UsageError("--start must be <= --end.")

    database = Database(db)
    stats = database.dump(
        output=output,
        phenotype=_expand_tokens(phenotype),
        sex=sex,
        tech=_expand_tokens(tech),
        by_sex=by_sex,
        by_tech=by_tech,
        by_phenotype=_expand_tokens(by_phenotype),
        all_groups=all_groups,
        chrom=chrom,
        start=start,
        end=end,
        n_workers=threads,
        include_ac_zero=all_variants,
        min_pass=min_pass,
        min_observed=min_observed,
        min_quality_evidence=min_quality_evidence,
    )
    click.echo(
        f"{stats['n_rows']} row(s) exported from {stats['n_buckets']} bucket(s)"
        f" across {stats['n_chroms']} chrom(s).",
        err=True,
    )


@cli.command()
@click.option("--db",       required=True, help="Path to database directory.")
@click.option("--samples",  is_flag=True,  help="List all samples with metadata. (default: false)")
@click.option("--changelog", "show_changelog", is_flag=True,
              help="Show full changelog history. (default: false)")
@click.option("--format", "fmt", default="table",
              type=click.Choice(["table", "tsv", "json"]),
              help="Output format. Options: table, tsv, json. (default: table)")
def info(db, samples, show_changelog, fmt):
    """Show database metadata and introspection."""
    database = Database(db)

    if samples:
        sample_list = database.list_samples()
        if fmt == "json":
            click.echo(json.dumps(sample_list, indent=2))
        elif fmt == "tsv":
            click.echo("sample_id\tsample_name\tsex\ttech\tphenotypes\tvcf_path\tingested_at")
            for s in sample_list:
                click.echo(
                    f"{s['sample_id']}\t{s['sample_name']}\t{s['sex']}\t{s['tech']}\t"
                    f"{','.join(s['phenotypes'])}\t{s['vcf_path'] or ''}\t{s['ingested_at'] or ''}"
                )
        else:
            click.echo(f"{'ID':<6} {'Name':<30} {'Sex':<8} {'Tech':<15} Phenotypes")
            click.echo("-" * 80)
            for s in sample_list:
                click.echo(
                    f"{s['sample_id']:<6} {s['sample_name']:<30} {s['sex']:<8} "
                    f"{s['tech']:<15} {','.join(s['phenotypes'])}"
                )
        return

    if show_changelog:
        cl = database.changelog()
        if fmt == "json":
            click.echo(json.dumps(cl, indent=2))
        elif fmt == "tsv":
            click.echo("event_id\tevent_type\tevent_time\tsample_names\tnotes")
            for e in cl:
                names = ','.join(e['sample_names']) if e['sample_names'] else ''
                click.echo(
                    f"{e['event_id']}\t{e['event_type']}\t{e['event_time']}\t"
                    f"{names}\t{e['notes'] or ''}"
                )
        else:
            for e in cl:
                n_str = f" ({len(e['sample_names'])} samples)" if e['sample_names'] else ""
                click.echo(f"  [{e['event_time']}] {e['event_type']:<20} — {e['notes'] or ''}{n_str}")
        return

    data = database.info()
    if fmt == "json":
        click.echo(json.dumps(data, indent=2))
    elif fmt == "tsv":
        for k, v in data.items():
            click.echo(f"{k}\t{v}")
    else:
        click.echo(f"Database:     {data['db_path']}")
        click.echo(f"Version:      {data['db_version']}")
        click.echo(f"Genome build: {data['genome_build']}   Schema: {data['schema_version']}")
        click.echo(f"Created: {data.get('created_at') or 'N/A'}    Updated: {data.get('updated_at') or 'N/A'}")
        click.echo()
        click.echo(f"Samples: {data['sample_count']} total")
        if data.get("by_sex"):
            sex_str = "  ".join(f"{k}={v}" for k, v in sorted(data["by_sex"].items()))
            click.echo(f"  By sex:       {sex_str}")
        if data.get("by_tech"):
            tech_str = "  ".join(f"{k}={v}" for k, v in sorted(data["by_tech"].items()))
            click.echo(f"  By tech:      {tech_str}")
        if data.get("by_phenotype"):
            pheno_str = "  ".join(f"{k}={v}" for k, v in sorted(data["by_phenotype"].items()))
            click.echo(f"  By phenotype: {pheno_str}")
        recent = data.get("changelog_recent", [])
        if recent:
            click.echo()
            click.echo("Recent changes:")
            for e in recent:
                click.echo(f"  [{e['event_time']}] {e['event_type']:<20} — {e['notes'] or ''}")


@cli.group()
def version():
    """Manage the database version label."""


@version.command("show")
@click.option("--db", required=True, help="Path to database directory.")
def version_show(db):
    """Show the current database version."""
    database = Database(db)
    click.echo(database.info().get("db_version", "unknown"))


@version.command("set")
@click.option("--db", required=True, help="Path to database directory.")
@click.argument("new_version")
def version_set(db, new_version):
    """Set the database version label."""
    database = Database(db)
    database.set_db_version(new_version)
    click.echo(f"Database version set to: {new_version}")


@cli.command("create-db")
@click.option("--manifest",     required=True, help="Path to TSV manifest file.")
@click.option("--output-dir",   required=True, help="Path to output database directory.")
@click.option("--genome-build", required=True, type=click.Choice(["GRCh37", "GRCh38"]), help="Reference genome build. Options: GRCh37, GRCh38.")
@click.option("--threads",      default=None, type=int, help="Number of worker threads for parallel processing. (default: all available CPU cores)")
@click.option("--build-threads", default=None, type=int, help="Max parallel workers for the build phase. (default: min(cpu_count, n_buckets))")
@click.option("--build-memory", default="2GB", help="DuckDB memory limit per build worker. Increase for WGS or dense regions. (default: 2GB)")
@click.option("--tmp-dir",      default=None, help="Temporary directory for intermediate files. (default: {output_dir}/.tmp_preprocess)")
@click.option("--bed-dir",      default=None, help="Directory containing BED files for WES technologies.")
@click.option("--force", is_flag=True, default=False, help="Delete any partial results and restart from scratch. (default: False)")
@click.option("--db-version", "db_version", default="1.0", help="Version label for this database. (default: 1.0)")
@click.option("--min-dp",      default=0,   type=int,   help="Phase 2: minimum FORMAT/DP for a carrier to count as quality evidence. (default: 0 = disabled)")
@click.option("--min-gq",      default=0,   type=int,   help="Phase 2: minimum FORMAT/GQ for a carrier to count as quality evidence. (default: 0 = disabled)")
@click.option("--min-qual",    default=0.0, type=float, help="Phase 2: minimum VCF QUAL for a carrier to count as quality evidence. (default: 0 = disabled)")
@click.option("--min-covered", default=0,   type=int,   help="Phase 2: minimum quality-passing carriers per WES tech for hom-ref to be assumed. Triggers Phase 2 storage when > 0. (default: 0 = disabled)")
@click.option("--verbose", "-v", is_flag=True, help="Verbose output with per-item progress. (default: false)")
def create_db_command(manifest, output_dir, genome_build, threads: int | None, build_threads: int | None, build_memory: str, tmp_dir, bed_dir, force, db_version,
                      min_dp, min_gq, min_qual, min_covered, verbose):
    """Create a new query database from a VCF manifest."""
    _configure_logging(verbose)
    from .preprocess import run_preprocess
    from .preprocess.manifest import ManifestError
    from .preprocess.ingest import IngestError
    try:
        run_preprocess(
            manifest_path=manifest,
            output_dir=output_dir,
            genome_build=genome_build,
            bed_dir=bed_dir,
            threads=threads,
            build_threads=build_threads,
            build_memory=build_memory,
            tmp_dir=tmp_dir,
            force=force,
            db_version=db_version,
            min_dp=min_dp,
            min_gq=min_gq,
            min_qual=min_qual,
            min_covered=min_covered,
        )
        click.echo(f"Database written to {output_dir}")
    except (ManifestError, IngestError) as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1)


@cli.command("update-db")
@click.option("--db",             required=True, help="Path to database directory.")
@click.option("--remove-samples", multiple=True, help="Sample name(s) to remove. Repeatable; comma-separated or multiple flags.")
@click.option("--add-samples",    multiple=True, type=click.Path(exists=True), help="Manifest TSV of new samples to add. Repeatable for multiple manifests.")
@click.option("--compact",        is_flag=True, help="Remove dead bits from removed samples to reclaim disk space.")
@click.option("--update-sample",  default=None, help="Sample name to update (single-sample metadata mode).")
@click.option("--set-sex",        default=None, type=click.Choice(["male", "female"]), help="New sex for --update-sample. Options: male, female.")
@click.option("--set-phenotype",  default=None, help="New phenotype codes (comma-separated) for --update-sample. Replaces all current codes.")
@click.option("--update-samples-file", default=None, type=click.Path(exists=True), help="TSV file for batch metadata update. Header: sample_name<TAB>field<TAB>new_value.")
@click.option("--operator-note",  default=None, help="Free-text note appended to each changelog entry for this update.")
@click.option("--threads",        default=None, type=int, help="Number of worker threads for parallel processing. (default: all available CPU cores)")
@click.option("--tmp-dir",        default=None, help="Temporary directory for intermediate files. (default: system temp)")
@click.option("--bed-dir",        default=None, help="Directory containing BED files for WES technologies.")
@click.option("--db-version", "db_version", default=None, help="New version label after update. (default: auto-increment current version)")
@click.option("--verbose", "-v",  is_flag=True, help="Verbose output with per-item progress. (default: false)")
def update_db_command(
    db, remove_samples, add_samples, compact,
    update_sample, set_sex, set_phenotype, update_samples_file, operator_note,
    threads, tmp_dir, bed_dir, db_version, verbose,
):
    """Update an existing database: remove samples, add samples, update metadata, and/or compact.

    Operations are applied in order: remove → update-metadata → add → compact.

    \b
    Update metadata for a single sample:
      afquery update-db --db ./db/ --update-sample S1 --set-sex female
      afquery update-db --db ./db/ --update-sample S1 --set-phenotype "E11.9,I10"

    \b
    Batch update from a TSV file (header: sample_name<TAB>field<TAB>new_value):
      afquery update-db --db ./db/ --update-samples-file changes.tsv
    """
    has_metadata_update = bool(update_sample or update_samples_file)
    if not any([remove_samples, add_samples, compact, has_metadata_update]):
        raise click.UsageError(
            "At least one of --remove-samples, --add-samples, --compact, "
            "--update-sample, or --update-samples-file is required."
        )

    # Validate metadata-update option combinations
    if (set_sex or set_phenotype) and not update_sample:
        raise click.UsageError("--set-sex and --set-phenotype require --update-sample.")
    if update_sample and not (set_sex or set_phenotype):
        raise click.UsageError("--update-sample requires --set-sex and/or --set-phenotype.")
    if update_samples_file and update_sample:
        raise click.UsageError("--update-samples-file and --update-sample are mutually exclusive.")

    _configure_logging(verbose)
    from .preprocess.update import (
        add_samples as _add_samples,
        remove_samples as _remove_samples,
        update_sample_metadata as _update_metadata,
        parse_updates_tsv,
        UpdateError,
    )
    from .preprocess.manifest import ManifestError
    from .preprocess.ingest import IngestError
    from .preprocess.compact import compact_database
    from pathlib import Path

    try:
        if remove_samples:
            names = _expand_tokens(remove_samples)
            result = _remove_samples(db, names)
            click.echo(f"Removed {len(result['removed'])} sample(s): {', '.join(result['removed'])}")

        if has_metadata_update:
            if update_samples_file:
                meta_updates = parse_updates_tsv(update_samples_file)
            else:
                meta_updates = []
                if set_sex:
                    meta_updates.append(
                        {"sample_name": update_sample, "field": "sex", "new_value": set_sex}
                    )
                if set_phenotype:
                    meta_updates.append(
                        {"sample_name": update_sample, "field": "phenotype_codes", "new_value": set_phenotype}
                    )
            entries = _update_metadata(db, meta_updates, operator_note=operator_note)
            n_samples = len({e["sample"] for e in entries})
            click.echo(f"Updated {len(entries)} field(s) for {n_samples} sample(s).")
            for e in entries:
                click.echo(f"  {e['sample']}: {e['field']}  {e['old']!r} → {e['new']!r}")

        for manifest in add_samples:
            result = _add_samples(
                db, manifest, threads=threads, tmp_dir=tmp_dir,
                bed_dir=bed_dir, db_version=db_version,
            )
            click.echo(
                f"Added {result['new_samples']} sample(s): "
                f"{result['new_variants']} new variants, "
                f"{result['updated_variants']} updated variants."
            )

        if compact:
            stats = compact_database(Path(db))
            click.echo(
                f"Compact complete: {stats['files_rewritten']} file(s) rewritten, "
                f"{stats['rows_removed']} empty row(s) removed, "
                f"{stats['rows_processed']} row(s) processed. "
                f"Size: {stats['size_before']} → {stats['size_after']} bytes."
            )
    except (UpdateError, ManifestError, IngestError) as e:
        click.echo(f"Error: {e}", err=True)
        raise SystemExit(1)


@cli.command("check")
@click.option("--db", required=True, help="Path to database directory.")
def check_cmd(db):
    """Validate database integrity. Exits with code 1 if any errors are found."""
    from .preprocess.update import check_database
    results = check_database(db)
    has_error = False
    for r in results:
        if r.severity == "error":
            prefix = "[ERROR]"
            has_error = True
        elif r.severity == "warning":
            prefix = "[WARN ]"
        else:
            prefix = "[INFO ]"
        click.echo(f"{prefix} {r.message}")
    if has_error:
        raise SystemExit(1)


@cli.command()
@click.option("--n-samples",  default=1000, type=int, help="Number of synthetic samples to generate. (default: 1000)")
@click.option("--n-variants", default=10_000, type=int, help="Number of variants per chromosome. (default: 10000)")
@click.option("--output",     default="benchmark_report.json", help="Output path for JSON benchmark report. (default: benchmark_report.json)")
@click.option("--db",         default=None, help="Use existing database instead of generating synthetic data.")
def benchmark(n_samples, n_variants, output, db):
    """Run performance benchmark suite and output a timing report."""
    from .benchmark import run_benchmark, run_benchmark_with_synth
    from pathlib import Path

    if db is not None:
        click.echo(f"Benchmarking existing DB at {db} ...")
        results = run_benchmark(Path(db))
        import json as _json
        with open(output, "w") as f:
            _json.dump(results, f, indent=2)
    else:
        click.echo(
            f"Generating synthetic DB ({n_samples} samples, {n_variants} variants/chrom) ..."
        )
        results = run_benchmark_with_synth(
            n_samples=n_samples,
            n_variants=n_variants,
            output_report=output,
        )

    if "error" in results:
        click.echo(f"Benchmark error: {results['error']}", err=True)
        raise SystemExit(1)

    click.echo(f"Report written to {output}")
    click.echo(f"  point cold:  {results['point_query_cold_ms']:.1f} ms")
    click.echo(f"  point warm:  {results['point_query_warm_ms']:.1f} ms")
    click.echo(f"  batch-100:   {results['batch_100_ms']:.1f} ms")
    click.echo(f"  batch-1000:  {results['batch_1000_ms']:.1f} ms")
    targets = results.get("targets", {})
    ok = all(targets.values())
    click.echo("Targets: " + ("ALL MET" if ok else "SOME MISSED"))

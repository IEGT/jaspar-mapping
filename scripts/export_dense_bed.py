#!/usr/bin/env python3

"""Export thresholded BED rows from a packaged dense-score DuckDB database."""

from __future__ import annotations

import argparse
import gzip
import math
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import BinaryIO


class ExportError(RuntimeError):
    pass


def finite_float(value: str) -> float:
    try:
        number = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"expected a number, got {value!r}") from error
    if not math.isfinite(number):
        raise argparse.ArgumentTypeError("score values must be finite")
    return number


def nonnegative_float(value: str) -> float:
    number = finite_float(value)
    if number < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return number


def nonnegative_integer(value: str) -> int:
    try:
        number = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            f"expected a non-negative integer, got {value!r}"
        ) from error
    if number < 0:
        raise argparse.ArgumentTypeError("coordinates must be non-negative")
    return number


def positive_integer(value: str) -> int:
    number = nonnegative_integer(value)
    if number == 0:
        raise argparse.ArgumentTypeError("value must be greater than zero")
    return number


def decimal_places(value: str) -> int:
    number = nonnegative_integer(value)
    if number > 15:
        raise argparse.ArgumentTypeError("decimal places must be between 0 and 15")
    return number


def sql_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def sql_number(value: float) -> str:
    return format(value, ".17g")


def parse_chromosomes(values: list[str]) -> list[str]:
    chromosomes: list[str] = []
    for value in values:
        for chromosome in value.split(","):
            chromosome = chromosome.strip()
            if not chromosome:
                raise ExportError("--chrom contains an empty chromosome name")
            if chromosome not in chromosomes:
                chromosomes.append(chromosome)
    return chromosomes


def resolve_database(package: Path, requested: Path | None) -> Path:
    if not package.is_dir():
        raise ExportError(f"package directory not found: {package}")

    if requested is not None:
        database = requested if requested.is_absolute() else package / requested
        if not database.is_file():
            raise ExportError(f"DuckDB database not found: {database}")
        return database.resolve()

    databases = sorted(package.glob("*.duckdb"))
    if not databases:
        raise ExportError(f"no .duckdb database found in {package}")
    if len(databases) > 1:
        names = ", ".join(database.name for database in databases)
        raise ExportError(
            f"multiple DuckDB databases found ({names}); select one with --database"
        )
    return databases[0].resolve()


def duckdb_command(database: Path, sql: str, header: bool) -> list[str]:
    return [
        "duckdb",
        "-readonly",
        "-batch",
        "-no-init",
        "-list",
        "-separator",
        "\t",
        "-header" if header else "-noheader",
        str(database),
        "-c",
        sql,
    ]


def run_capture(package: Path, database: Path, sql: str,
                header: bool = False) -> str:
    result = subprocess.run(
        duckdb_command(database, sql, header),
        cwd=package,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if result.returncode != 0:
        detail = result.stderr.strip() or "DuckDB query failed"
        raise ExportError(detail)
    return result.stdout


def inventory_sql() -> str:
    return """
SELECT
    i.motif_id AS MotifID,
    m.motif_name AS Name,
    i.score_mode AS ScoreMode,
    i.pseudocount AS Pseudocount,
    i.chrom AS Chromosome,
    i.strand AS Strand,
    i.alignment_start_begin AS Start,
    i.alignment_start_end AS End,
    i.n_valid_windows AS ValidScores,
    i.n_skipped_windows AS SkippedScores
FROM dense_run_inventory i
JOIN motif_metadata m USING (motif_id)
ORDER BY i.motif_id, i.score_mode, i.pseudocount,
         CASE
             WHEN try_cast(i.chrom AS BIGINT) IS NOT NULL THEN 0
             WHEN upper(i.chrom) = 'X' THEN 1
             WHEN upper(i.chrom) = 'Y' THEN 2
             WHEN upper(i.chrom) IN ('M', 'MT') THEN 3
             ELSE 4
         END,
         try_cast(i.chrom AS BIGINT), i.chrom, i.strand;
"""


def requested_strands(strand: str) -> list[str]:
    return ["+", "-"] if strand == "both" else [strand]


def configuration_filter(args: argparse.Namespace,
                         chromosomes: list[str]) -> list[str]:
    strands = requested_strands(args.strand)
    filters = [
        f"b.motif_id = {sql_string(args.motif)}",
        f"b.score_mode = {sql_string(args.score_mode)}",
        f"b.pseudocount = {sql_number(args.pseudocount)}",
        "b.strand IN (" + ", ".join(sql_string(value) for value in strands) + ")",
    ]
    if chromosomes:
        filters.append(
            "b.chrom IN ("
            + ", ".join(sql_string(chromosome) for chromosome in chromosomes)
            + ")"
        )
    return filters


def validation_sql(args: argparse.Namespace, chromosomes: list[str]) -> str:
    filters = configuration_filter(args, chromosomes)
    inventory_filters = [condition.replace("b.", "i.") for condition in filters]
    return (
        "SELECT DISTINCT i.chrom, i.strand FROM dense_run_inventory i WHERE "
        + " AND ".join(inventory_filters)
        + " ORDER BY i.chrom, i.strand;"
    )


def validate_configuration(inventory: str, chromosomes: list[str],
                           strand: str) -> None:
    available: set[tuple[str, str]] = set()
    for line in inventory.splitlines():
        fields = line.split("\t")
        if len(fields) != 2:
            raise ExportError(f"unexpected DuckDB validation row: {line!r}")
        available.add((fields[0], fields[1]))

    if not available:
        raise ExportError(
            "no stored configuration matches the requested motif, mode, "
            "pseudocount, chromosome, and strand; use --list-configs"
        )

    selected_chromosomes = chromosomes or sorted(
        {chromosome for chromosome, _ in available}
    )
    missing = [
        f"{chromosome}:{orientation}"
        for chromosome in selected_chromosomes
        for orientation in requested_strands(strand)
        if (chromosome, orientation) not in available
    ]
    if missing:
        raise ExportError(
            "requested chromosome/strand partitions are absent: "
            + ", ".join(missing)
            + "; use --list-configs"
        )


def export_cte(args: argparse.Namespace, chromosomes: list[str]) -> str:
    filters = configuration_filter(args, chromosomes)
    if args.start is not None:
        filters.append(f"b.block_start + len(b.scores) > {args.start}")
    if args.end is not None:
        filters.append(f"b.block_start < {args.end}")

    row_filters = ["e.score IS NOT NULL"]
    if args.start is not None:
        row_filters.append(f"e.start >= {args.start}")
    if args.end is not None:
        row_filters.append(f"e.start < {args.end}")
    if not args.all_scores:
        minimum = 0.0 if args.min_score is None else args.min_score
        row_filters.append(f"e.score >= {sql_number(minimum)}")
    if args.max_score is not None:
        row_filters.append(f"e.score <= {sql_number(args.max_score)}")

    return f"""
WITH selected_blocks AS (
    SELECT b.*, m.motif_name, m.motif_length
    FROM motif_score_dense_block b
    JOIN motif_metadata m USING (motif_id)
    WHERE {' AND '.join(filters)}
),
expanded AS (
    SELECT
        b.chrom,
        b.block_start + CAST(u.offset_one - 1 AS BIGINT) AS start,
        b.block_start + CAST(u.offset_one - 1 AS BIGINT)
            + b.motif_length AS "end",
        b.motif_id,
        b.motif_name,
        b.strand,
        b.score_mode,
        b.pseudocount,
        u.score
    FROM selected_blocks b
    CROSS JOIN UNNEST(b.scores) WITH ORDINALITY AS u(score, offset_one)
),
filtered AS (
    SELECT *
    FROM expanded e
    WHERE {' AND '.join(row_filters)}
)
"""


def chromosome_order_sql(alias: str) -> str:
    return f"""
CASE
    WHEN try_cast({alias}.chrom AS BIGINT) IS NOT NULL THEN 0
    WHEN upper({alias}.chrom) = 'X' THEN 1
    WHEN upper({alias}.chrom) = 'Y' THEN 2
    WHEN upper({alias}.chrom) IN ('M', 'MT') THEN 3
    ELSE 4
END,
try_cast({alias}.chrom AS BIGINT), {alias}.chrom, {alias}.start,
CASE {alias}.strand WHEN '+' THEN 0 ELSE 1 END
"""


def bed_sql(args: argparse.Namespace, chromosomes: list[str]) -> str:
    score = "CAST(f.score AS VARCHAR)"
    if args.score_decimals is not None:
        score = f"printf('%.{args.score_decimals}f', f.score)"

    columns = [
        'f.chrom AS "Chromosome"',
        'f.start AS "From"',
        'f."end" AS "To"',
        'f.motif_name AS "Name"',
        f'{score} AS "Score"',
        'f.strand AS "Strand"',
    ]
    if args.columns == "provenance":
        columns.extend(
            [
                'f.motif_id AS "MotifID"',
                'f.score_mode AS "ScoreMode"',
                'f.pseudocount AS "Pseudocount"',
                "'bed' AS \"CoordinateMode\"",
            ]
        )

    limit = f"LIMIT {args.limit}" if args.limit is not None else ""
    return (
        export_cte(args, chromosomes)
        + "SELECT\n    "
        + ",\n    ".join(columns)
        + "\nFROM filtered f\nORDER BY "
        + chromosome_order_sql("f")
        + limit
        + ";"
    )


def count_sql(args: argparse.Namespace, chromosomes: list[str]) -> str:
    return export_cte(args, chromosomes) + 'SELECT COUNT(*) AS "Matches" FROM filtered;'


def run_to_stream(package: Path, database: Path, sql: str, header: bool,
                  output: BinaryIO) -> None:
    process = subprocess.Popen(
        duckdb_command(database, sql, header),
        cwd=package,
        stdout=subprocess.PIPE,
    )
    if process.stdout is None:
        process.terminate()
        process.wait()
        raise ExportError("could not read DuckDB export output")
    try:
        shutil.copyfileobj(process.stdout, output, length=1024 * 1024)
    except BaseException:
        process.stdout.close()
        process.terminate()
        process.wait()
        raise
    process.stdout.close()
    process.wait()
    if process.returncode != 0:
        raise ExportError(f"DuckDB export failed with exit code {process.returncode}")


def output_permissions() -> int:
    current_umask = os.umask(0)
    os.umask(current_umask)
    return 0o666 & ~current_umask


def write_output(package: Path, database: Path, sql: str,
                 args: argparse.Namespace) -> None:
    if args.output == "-":
        run_to_stream(package, database, sql, args.header, sys.stdout.buffer)
        return

    output = Path(args.output).expanduser()
    if output.exists() and not args.force:
        raise ExportError(f"output already exists (use --force to replace it): {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.name}.", suffix=".tmp", dir=output.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        if output.suffix == ".gz":
            with temporary.open("wb") as raw_output:
                with gzip.GzipFile(
                    filename="", mode="wb", fileobj=raw_output, mtime=0
                ) as compressed_output:
                    run_to_stream(
                        package, database, sql, args.header, compressed_output
                    )
        else:
            with temporary.open("wb") as plain_output:
                run_to_stream(package, database, sql, args.header, plain_output)

        if output.exists() and not args.force:
            raise ExportError(f"output appeared during export and was not replaced: {output}")
        os.chmod(temporary, output_permissions())
        os.replace(temporary, output)
        print(f"I: Wrote {output}", file=sys.stderr)
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Regenerate thresholded BED rows from dense motif-score Parquet "
            "through its packaged DuckDB query index."
        ),
        epilog=(
            "Ranges are BED-style, 0-based, half-open ranges of motif alignment "
            "starts. Scores at skipped N/sentinel positions are NULL and are "
            "never exported."
        ),
    )
    parser.add_argument("--package", required=True, type=Path,
                        help="package directory containing Parquet and DuckDB")
    parser.add_argument("--database", type=Path,
                        help="database path or package-relative name (auto-detected by default)")
    parser.add_argument("--list-configs", action="store_true",
                        help="list stored motifs/configurations and exit")
    parser.add_argument("--motif", help="motif ID to export, for example MA0861.2")
    parser.add_argument(
        "--score-mode", choices=("log2_relative_risk", "log_odds"),
        help="stored scoring transform; modes are never combined",
    )
    parser.add_argument("--pseudocount", type=nonnegative_float, default=1.0,
                        help="stored pseudocount configuration (default: 1)")
    parser.add_argument(
        "--chrom", action="append", default=[], metavar="CHROM[,CHROM...]",
        help="chromosome selection; repeat or use commas (default: all stored)",
    )
    parser.add_argument("--strand", choices=("+", "-", "both"), default="both",
                        help="orientation to export (default: both)")
    parser.add_argument("--from", "--start", dest="start", type=nonnegative_integer,
                        help="first motif alignment start to include")
    parser.add_argument("--to", "--end", dest="end", type=nonnegative_integer,
                        help="exclusive motif alignment-start bound")
    minimum = parser.add_mutually_exclusive_group()
    minimum.add_argument("--min-score", "--threshold", type=finite_float,
                         help="inclusive lower score bound (default: 0)")
    minimum.add_argument("--all-scores", action="store_true",
                         help="remove the default lower score bound")
    parser.add_argument("--max-score", type=finite_float,
                        help="inclusive upper score bound")
    parser.add_argument(
        "--columns", choices=("bed6", "provenance"), default="bed6",
        help="BED-compatible first six columns, optionally with run provenance",
    )
    parser.add_argument(
        "--score-decimals", type=decimal_places,
        help="fixed score decimals; omit to preserve the stored FLOAT representation",
    )
    parser.add_argument("--limit", type=positive_integer,
                        help="emit at most this many sorted rows (for previews)")
    parser.add_argument("--count-only", action="store_true",
                        help="print the number of matching rows without exporting BED")
    parser.add_argument("--output", "-o", default="-",
                        help="output path or - for stdout; .gz is compressed directly")
    parser.add_argument("--force", action="store_true",
                        help="replace an existing output file atomically")
    header = parser.add_mutually_exclusive_group()
    header.add_argument("--header", dest="header", action="store_true",
                        help="write column names (default)")
    header.add_argument("--no-header", dest="header", action="store_false",
                        help="omit column names for strict BED consumers")
    parser.set_defaults(header=True)
    return parser


def main() -> int:
    parser = argument_parser()
    args = parser.parse_args()

    try:
        if shutil.which("duckdb") is None:
            raise ExportError("duckdb CLI is required in PATH")

        package = args.package.expanduser().resolve()
        database = resolve_database(package, args.database)

        if args.list_configs:
            if args.output != "-":
                raise ExportError("--list-configs writes to stdout; omit --output")
            sys.stdout.write(run_capture(package, database, inventory_sql(), header=True))
            return 0

        if not args.motif:
            parser.error("--motif is required unless --list-configs is used")
        if not args.score_mode:
            parser.error("--score-mode is required unless --list-configs is used")
        if args.start is not None and args.end is not None and args.end <= args.start:
            parser.error("--to/--end must be greater than --from/--start")

        minimum = None if args.all_scores else (
            0.0 if args.min_score is None else args.min_score
        )
        if (minimum is not None and args.max_score is not None
                and args.max_score < minimum):
            parser.error("--max-score must be greater than or equal to --min-score")

        chromosomes = parse_chromosomes(args.chrom)
        available = run_capture(
            package, database, validation_sql(args, chromosomes)
        )
        validate_configuration(available, chromosomes, args.strand)

        if args.count_only:
            if args.output != "-":
                raise ExportError("--count-only writes to stdout; omit --output")
            sys.stdout.write(
                run_capture(package, database, count_sql(args, chromosomes), header=True)
            )
            return 0

        write_output(package, database, bed_sql(args, chromosomes), args)
        return 0
    except (ExportError, OSError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

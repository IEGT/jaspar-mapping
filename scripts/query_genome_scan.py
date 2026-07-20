#!/usr/bin/env python3

"""Query a finalized genome-scan package through exact inventory paths."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable


class QueryError(RuntimeError):
    pass


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_list(values: Iterable[str | Path]) -> str:
    return "[" + ",".join(sql_string(value) for value in values) + "]"


def package_paths(arguments: argparse.Namespace) -> tuple[Path, Path]:
    package = arguments.package.expanduser().resolve()
    manifest_path = package / "manifest.json"
    if not manifest_path.is_file():
        raise QueryError(f"package manifest not found: {manifest_path}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise QueryError(f"cannot read package manifest: {error}") from error
    database_name = arguments.database or manifest.get("database")
    if not database_name:
        raise QueryError("manifest does not identify a DuckDB catalog")
    database = Path(database_name)
    if not database.is_absolute():
        database = package / database
    database = database.resolve()
    if not database.is_file():
        raise QueryError(f"DuckDB catalog not found: {database}")
    return package, database


def run_duckdb(duckdb: str, database: Path, query: str, output_format: str,
               output: Path | None = None) -> str:
    command = [duckdb, "-readonly"]
    if output_format == "json":
        command.append("-json")
    elif output_format == "csv":
        command.append("-csv")
    command.extend([str(database), "-c", query])
    process = subprocess.run(command, text=True, capture_output=True, check=False)
    if process.returncode != 0:
        raise QueryError(process.stderr.strip() or "DuckDB query failed")
    if output is not None:
        output = output.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("x", encoding="utf-8") as stream:
            stream.write(process.stdout)
    return process.stdout


def query_json(duckdb: str, database: Path, query: str) -> list[dict[str, object]]:
    raw = run_duckdb(duckdb, database, query, "json")
    try:
        rows = json.loads(raw or "[]")
    except json.JSONDecodeError as error:
        raise QueryError("DuckDB returned invalid JSON") from error
    if not isinstance(rows, list):
        raise QueryError("DuckDB JSON result is not a row array")
    return rows


def inventory_where(arguments: argparse.Namespace) -> str:
    conditions: list[str] = []
    if getattr(arguments, "motif", None):
        conditions.append(f"motif_id = {sql_string(arguments.motif)}")
    if getattr(arguments, "chrom", None):
        conditions.append(f"CAST(chrom AS VARCHAR) = {sql_string(arguments.chrom)}")
    if getattr(arguments, "strand", None):
        conditions.append(f"strand = {sql_string(arguments.strand)}")
    return " AND ".join(conditions) if conditions else "TRUE"


def selected_files(arguments: argparse.Namespace, package: Path,
                   database: Path) -> list[dict[str, object]]:
    query = f"""
SELECT task_id, output_relative_path, motif_id, CAST(chrom AS VARCHAR) AS chrom,
       strand, emitted_hits, bytes, sha256
FROM scan_file_inventory
WHERE {inventory_where(arguments)}
ORDER BY motif_id, chrom, strand;
"""
    rows = query_json(arguments.duckdb, database, query)
    for row in rows:
        path = (
            package / "task_data" / f"task_id={row['task_id']}"
            / str(row["output_relative_path"])
        ).resolve()
        if not path.is_file():
            raise QueryError(f"inventory payload is missing: {path}")
        row["absolute_path"] = str(path)
    return rows


def print_python_rows(rows: list[dict[str, object]], output_format: str) -> None:
    if output_format == "json":
        print(json.dumps(rows, indent=2, sort_keys=True))
        return
    if not rows:
        return
    columns = list(rows[0])
    separator = "," if output_format == "csv" else "\t"
    print(separator.join(columns))
    for row in rows:
        print(separator.join(str(row.get(column, "")) for column in columns))


def summary(arguments: argparse.Namespace) -> None:
    _, database = package_paths(arguments)
    query = """
SELECT r.run_id, r.genome_id, r.motif_set_id, r.score_mode, r.pseudocount,
       r.task_count, r.file_count, r.emitted_hit_count, s.bytes AS parquet_bytes,
       s.n_motifs, s.n_sequence_regions
FROM scan_run r
JOIN scan_inventory_summary s USING (run_id, genome_id, motif_set_id);
"""
    rendered = run_duckdb(
        arguments.duckdb, database, query, arguments.format, arguments.output
    )
    if arguments.output is None:
        sys.stdout.write(rendered)


def files(arguments: argparse.Namespace) -> None:
    package, database = package_paths(arguments)
    rows = selected_files(arguments, package, database)
    if arguments.output is not None:
        output = arguments.output.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("x", encoding="utf-8") as stream:
            if arguments.format == "json":
                json.dump(rows, stream, indent=2, sort_keys=True)
                stream.write("\n")
            elif rows:
                columns = list(rows[0])
                separator = "," if arguments.format == "csv" else "\t"
                stream.write(separator.join(columns) + "\n")
                for row in rows:
                    stream.write(
                        separator.join(str(row.get(column, "")) for column in columns)
                        + "\n"
                    )
        return
    print_python_rows(rows, arguments.format)


def hits(arguments: argparse.Namespace) -> None:
    package, database = package_paths(arguments)
    rows = selected_files(arguments, package, database)
    if not rows:
        print("I: No inventory files match the requested motif/chromosome/strand.",
              file=sys.stderr)
        return
    paths = [str(row["absolute_path"]) for row in rows]
    conditions = ["TRUE"]
    if arguments.start is not None:
        conditions.append(f'"end" > {arguments.start}')
    if arguments.end is not None:
        conditions.append(f"start < {arguments.end}")
    if arguments.minimum_score is not None:
        conditions.append(f"score >= {arguments.minimum_score:.17g}")
    limit = f"LIMIT {arguments.limit}" if arguments.limit is not None else ""
    query = f"""
SELECT run_id, task_id, genome_id, motif_set_id, chrom, start, "end",
       motif_id, motif_name, strand, score, pwm_relative_score, score_mode,
       pseudocount, minimum_score
FROM motif_hit_files({sql_list(paths)})
WHERE {' AND '.join(conditions)}
ORDER BY chrom, start, "end", motif_id, strand
{limit};
"""
    rendered = run_duckdb(
        arguments.duckdb, database, query, arguments.format, arguments.output
    )
    if arguments.output is None:
        sys.stdout.write(rendered)


def metadata_sql(arguments: argparse.Namespace) -> None:
    _, database = package_paths(arguments)
    if (arguments.query is None) == (arguments.file is None):
        raise QueryError("sql requires exactly one of --query or --file")
    query = (
        arguments.query
        if arguments.query is not None
        else arguments.file.expanduser().read_text(encoding="utf-8")
    )
    rendered = run_duckdb(
        arguments.duckdb, database, query, arguments.format, arguments.output
    )
    if arguments.output is None:
        sys.stdout.write(rendered)


def add_common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--package", required=True, type=Path)
    parser.add_argument("--database", type=Path)
    parser.add_argument("--duckdb", default="duckdb")
    parser.add_argument("--format", choices=("table", "json", "csv"), default="table")
    parser.add_argument("--output", type=Path)


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect scan metadata or query exact motif/chromosome Parquet files "
            "without binding the package-wide payload tree."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    summary_parser = subparsers.add_parser("summary", help="show run and inventory totals")
    add_common(summary_parser)
    summary_parser.set_defaults(function=summary)

    files_parser = subparsers.add_parser("files", help="list exact inventory payload paths")
    add_common(files_parser)
    files_parser.add_argument("--motif")
    files_parser.add_argument("--chrom")
    files_parser.add_argument("--strand", choices=("+", "-"))
    files_parser.set_defaults(function=files)

    hits_parser = subparsers.add_parser("hits", help="query selected motif-hit payloads")
    add_common(hits_parser)
    hits_parser.add_argument("--motif", required=True)
    hits_parser.add_argument("--chrom")
    hits_parser.add_argument("--strand", choices=("+", "-"))
    hits_parser.add_argument("--start", type=int)
    hits_parser.add_argument("--end", type=int)
    hits_parser.add_argument("--minimum-score", type=float)
    hits_parser.add_argument("--limit", type=int)
    hits_parser.set_defaults(function=hits)

    sql_parser = subparsers.add_parser("sql", help="run SQL against catalog metadata")
    add_common(sql_parser)
    sql_parser.add_argument("--query")
    sql_parser.add_argument("--file", type=Path)
    sql_parser.set_defaults(function=metadata_sql)
    return parser


def main() -> int:
    arguments = argument_parser().parse_args()
    try:
        if shutil.which(arguments.duckdb) is None:
            raise QueryError(f"DuckDB executable not found: {arguments.duckdb}")
        arguments.function(arguments)
        return 0
    except (OSError, QueryError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

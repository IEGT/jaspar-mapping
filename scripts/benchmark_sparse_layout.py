#!/usr/bin/env python3

"""Compare split-strand and combined-strand sparse Parquet layouts."""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
import time
from pathlib import Path

from query_genome_scan import QueryError, package_paths, query_json, sql_list, sql_string


class BenchmarkError(RuntimeError):
    pass


def timed_query(duckdb: str, query: str, repetitions: int) -> list[float]:
    elapsed: list[float] = []
    for _ in range(repetitions):
        started = time.perf_counter()
        process = subprocess.run(
            [duckdb, ":memory:", "-bail", "-c", query],
            text=True,
            capture_output=True,
            check=False,
        )
        if process.returncode != 0:
            raise BenchmarkError(process.stderr.strip() or "benchmark query failed")
        elapsed.append(time.perf_counter() - started)
    return elapsed


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Materialize one temporary motif/chromosome file containing both "
            "strands and compare size and simple read latency with the current "
            "two-file layout. The source package is never modified."
        )
    )
    parser.add_argument("--package", required=True, type=Path)
    parser.add_argument("--database", type=Path)
    parser.add_argument("--motif", required=True)
    parser.add_argument("--chrom", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--duckdb", default="duckdb")
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--row-group-size", type=int, default=122880)
    arguments = parser.parse_args()

    try:
        if arguments.repetitions <= 0 or arguments.row_group_size <= 0:
            raise BenchmarkError("repetitions and row-group size must be positive")
        package, database = package_paths(arguments)
        inventory_query = f"""
SELECT task_id, output_relative_path, strand, emitted_hits, bytes
FROM scan_file_inventory
WHERE motif_id = {sql_string(arguments.motif)}
  AND CAST(chrom AS VARCHAR) = {sql_string(arguments.chrom)}
ORDER BY strand;
"""
        rows = query_json(arguments.duckdb, database, inventory_query)
        if len(rows) != 2 or {row["strand"] for row in rows} != {"+", "-"}:
            raise BenchmarkError("expected exactly one plus and one minus inventory file")
        paths = [
            (
                package / "task_data" / f"task_id={row['task_id']}"
                / str(row["output_relative_path"])
            ).resolve()
            for row in rows
        ]
        if not all(path.is_file() for path in paths):
            raise BenchmarkError("one or more inventory payloads are missing")
        output_dir = arguments.output_dir.expanduser().resolve()
        output_dir.mkdir(parents=True, exist_ok=False)
        combined = output_dir / "combined-strands.parquet"
        copy_query = f"""
COPY (
    SELECT start, "end", score, pwm_relative_score, matched_seq,
           CASE strand WHEN 'plus' THEN '+' WHEN 'minus' THEN '-' ELSE strand END
               AS strand
    FROM read_parquet({sql_list(paths)}, hive_partitioning = true)
    ORDER BY start, "end", strand
) TO {sql_string(combined)} (
    FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE {arguments.row_group_size}
);
"""
        process = subprocess.run(
            [arguments.duckdb, ":memory:", "-bail", "-c", copy_query],
            text=True,
            capture_output=True,
            check=False,
        )
        if process.returncode != 0:
            raise BenchmarkError(process.stderr.strip() or "layout conversion failed")

        split_query = (
            f"SELECT COUNT(*), SUM(score) FROM read_parquet({sql_list(paths)}, "
            "hive_partitioning=true) WHERE score >= 0;"
        )
        combined_query = (
            f"SELECT COUNT(*), SUM(score) FROM read_parquet({sql_string(combined)}) "
            "WHERE score >= 0;"
        )
        split_times = timed_query(arguments.duckdb, split_query, arguments.repetitions)
        combined_times = timed_query(
            arguments.duckdb, combined_query, arguments.repetitions
        )
        split_bytes = sum(path.stat().st_size for path in paths)
        combined_bytes = combined.stat().st_size
        result = {
            "motif_id": arguments.motif,
            "chrom": arguments.chrom,
            "rows": sum(int(row["emitted_hits"]) for row in rows),
            "split_file_count": 2,
            "split_bytes": split_bytes,
            "combined_file_count": 1,
            "combined_bytes": combined_bytes,
            "combined_to_split_size_ratio": combined_bytes / split_bytes,
            "repetitions": arguments.repetitions,
            "split_query_seconds_median": statistics.median(split_times),
            "combined_query_seconds_median": statistics.median(combined_times),
            "row_group_size": arguments.row_group_size,
        }
        result_path = output_dir / "benchmark.json"
        result_path.write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0
    except (BenchmarkError, OSError, QueryError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

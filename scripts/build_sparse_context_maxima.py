#!/usr/bin/env python3

"""Build per-anchor cofactor maxima from sparse motif-hit Parquet files."""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


class ContextMaximaError(RuntimeError):
    pass


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def parse_cofactors(
    values: list[list[str]],
) -> list[tuple[str, Path, Path]]:
    cofactors: list[tuple[str, Path, Path]] = []
    observed: set[str] = set()
    for motif_id, plus_value, minus_value in values:
        if motif_id in observed:
            raise ContextMaximaError(f"duplicate cofactor motif: {motif_id}")
        observed.add(motif_id)
        plus = Path(plus_value).expanduser().resolve()
        minus = Path(minus_value).expanduser().resolve()
        for strand, path in (("plus", plus), ("minus", minus)):
            if not path.is_file():
                raise ContextMaximaError(
                    f"{strand} Parquet file for {motif_id} not found: {path}"
                )
        cofactors.append((motif_id, plus, minus))
    return cofactors


def hit_union_sql(cofactors: list[tuple[str, Path, Path]]) -> str:
    selects: list[str] = []
    for motif_id, plus, minus in cofactors:
        for strand, path in (("+", plus), ("-", minus)):
            selects.append(
                "SELECT "
                f"{sql_string(motif_id)}::VARCHAR AS motif_id, "
                "start::BIGINT AS start, \"end\"::BIGINT AS hit_end, "
                "score::FLOAT AS score, "
                f"{sql_string(strand)}::VARCHAR AS strand "
                f"FROM read_parquet({sql_string(path)})"
            )
    return "\nUNION ALL\n".join(selects)


def motif_values_sql(cofactors: list[tuple[str, Path, Path]]) -> str:
    return ",".join(f"({sql_string(motif_id)})" for motif_id, _, _ in cofactors)


def build_sql(arguments: argparse.Namespace, cofactors: list[tuple[str, Path, Path]],
              staging_output: Path, temp_directory: Path) -> str:
    return f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;
SET temp_directory={sql_string(temp_directory)};
SET max_temp_directory_size={sql_string(arguments.max_temp_size)};

COPY (
WITH anchor_input AS (
    SELECT
        CAST(chrom AS VARCHAR) AS chrom,
        anchor_start::BIGINT AS anchor_start,
        anchor_end::BIGINT AS anchor_end
    FROM read_parquet({sql_string(arguments.anchor_parquet)})
),
hit_input AS (
    {hit_union_sql(cofactors)}
),
span_parameter AS (
    SELECT
        MAX(anchor_end - anchor_start)::BIGINT AS max_anchor_span_bp,
        (SELECT MAX(hit_end - start)::BIGINT FROM hit_input)
            AS max_context_span_bp
    FROM anchor_input
),
parameter AS (
    SELECT
        {arguments.flank}::BIGINT AS context_flank_bp,
        GREATEST(
            1,
            {arguments.flank} + max_anchor_span_bp + max_context_span_bp
        )::BIGINT AS capture_prefilter_center_bp,
        max_anchor_span_bp,
        max_context_span_bp
    FROM span_parameter
),
anchors AS (
    SELECT
        a.*,
        FLOOR(
            ((a.anchor_start + a.anchor_end) / 2.0)
            / p.capture_prefilter_center_bp
        )::BIGINT AS capture_bin
    FROM anchor_input a
    CROSS JOIN parameter p
),
candidate_hit AS (
    SELECT
        h.*,
        FLOOR(
            ((h.start + h.hit_end) / 2.0)
            / p.capture_prefilter_center_bp
        )::BIGINT + offsets.bin_offset AS anchor_capture_bin
    FROM hit_input h
    CROSS JOIN parameter p
    CROSS JOIN (VALUES (-1), (0), (1)) AS offsets(bin_offset)
),
maxima AS (
    SELECT
        a.chrom,
        a.anchor_start,
        a.anchor_end,
        h.motif_id,
        MAX(h.score)::FLOAT AS context_score
    FROM anchors a
    JOIN candidate_hit h
      ON a.capture_bin = h.anchor_capture_bin
    CROSS JOIN parameter p
    WHERE GREATEST(a.anchor_start, h.start)
          - LEAST(a.anchor_end, h.hit_end) <= p.context_flank_bp
    GROUP BY a.chrom, a.anchor_start, a.anchor_end, h.motif_id
),
motif(motif_id) AS (
    VALUES {motif_values_sql(cofactors)}
)
SELECT
    a.chrom,
    a.anchor_start,
    a.anchor_end,
    m.motif_id,
    x.context_score,
    {arguments.source_score_floor}::DOUBLE AS source_score_floor,
    p.context_flank_bp,
    p.capture_prefilter_center_bp,
    p.max_anchor_span_bp AS observed_max_anchor_span_bp,
    p.max_context_span_bp AS observed_max_context_span_bp,
    'signed_interval_edge_distance'::VARCHAR AS context_distance_metric
FROM anchors a
CROSS JOIN motif m
CROSS JOIN parameter p
LEFT JOIN maxima x
  ON a.chrom = x.chrom
 AND a.anchor_start = x.anchor_start
 AND a.anchor_end = x.anchor_end
 AND m.motif_id = x.motif_id
) TO {sql_string(staging_output)}
  (FORMAT PARQUET, COMPRESSION ZSTD);
"""


def run_duckdb(duckdb: str, sql: str) -> None:
    process = subprocess.run(
        [duckdb, "-light-mode", ":memory:", "-c", sql],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ContextMaximaError(
            process.stderr.strip() or "DuckDB context-maximum build failed"
        )


def validate_output(duckdb: str, output: Path, expected_motifs: int) -> dict[str, int]:
    query = f"""
WITH output AS (SELECT * FROM read_parquet({sql_string(output)}))
SELECT
    COUNT(*)::BIGINT AS rows,
    COUNT(DISTINCT (chrom, anchor_start, anchor_end))::BIGINT AS anchors,
    COUNT(DISTINCT chrom)::BIGINT AS chromosomes,
    COUNT(DISTINCT motif_id)::BIGINT AS motifs,
    COUNT(*) - COUNT(DISTINCT (chrom, anchor_start, anchor_end, motif_id))
        AS duplicate_keys,
    COUNT(*) FILTER (
        WHERE context_score IS NOT NULL AND context_score < source_score_floor
    ) AS scores_below_source_floor,
    MIN(capture_prefilter_center_bp)::BIGINT AS capture_prefilter_center_bp,
    MIN(observed_max_anchor_span_bp)::BIGINT AS observed_max_anchor_span_bp,
    MIN(observed_max_context_span_bp)::BIGINT AS observed_max_context_span_bp,
    COUNT(DISTINCT capture_prefilter_center_bp)::BIGINT AS prefilter_values,
    COUNT(DISTINCT observed_max_anchor_span_bp)::BIGINT AS anchor_span_values,
    COUNT(DISTINCT observed_max_context_span_bp)::BIGINT AS context_span_values
FROM output;
"""
    process = subprocess.run(
        [duckdb, "-light-mode", "-json", ":memory:", "-c", query],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ContextMaximaError(
            process.stderr.strip() or "DuckDB output validation failed"
        )
    try:
        rows = json.loads(process.stdout)
    except json.JSONDecodeError as error:
        raise ContextMaximaError("DuckDB returned invalid validation JSON") from error
    if len(rows) != 1:
        raise ContextMaximaError("DuckDB returned an invalid validation row count")
    if rows[0].get("rows") in (None, 0):
        raise ContextMaximaError("anchor input produced no context-maximum rows")
    if any(rows[0].get(key) is None for key in (
        "capture_prefilter_center_bp", "observed_max_anchor_span_bp",
        "observed_max_context_span_bp",
    )):
        raise ContextMaximaError("anchor or motif-hit inputs lack interval spans")
    result = {key: int(value) for key, value in rows[0].items()}
    if result["chromosomes"] != 1:
        raise ContextMaximaError(
            "anchor input must contain exactly one chromosome because sparse "
            "motif-hit partitions omit chromosome"
        )
    if result["motifs"] != expected_motifs:
        raise ContextMaximaError("output does not contain every requested motif")
    if result["rows"] != result["anchors"] * expected_motifs:
        raise ContextMaximaError("output is not rectangular by anchor and motif")
    if result["duplicate_keys"] != 0 or result["scores_below_source_floor"] != 0:
        raise ContextMaximaError("output failed key or source-score-floor validation")
    if any(result[key] != 1 for key in (
        "prefilter_values", "anchor_span_values", "context_span_values"
    )):
        raise ContextMaximaError("output contains inconsistent geometry provenance")
    return result


def write_run_config(
    path: Path,
    arguments: argparse.Namespace,
    cofactors: list[tuple[str, Path, Path]],
    result: dict[str, int],
) -> None:
    fieldnames = (
        "output", "anchor_parquet", "motif_id", "plus_parquet",
        "minus_parquet", "context_flank_bp", "source_score_floor",
        "capture_prefilter_center_bp", "observed_max_anchor_span_bp",
        "observed_max_context_span_bp", "rows", "anchors",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for motif_id, plus, minus in cofactors:
            writer.writerow({
                "output": arguments.output,
                "anchor_parquet": arguments.anchor_parquet,
                "motif_id": motif_id,
                "plus_parquet": plus,
                "minus_parquet": minus,
                "context_flank_bp": arguments.flank,
                "source_score_floor": arguments.source_score_floor,
                "capture_prefilter_center_bp":
                    result["capture_prefilter_center_bp"],
                "observed_max_anchor_span_bp":
                    result["observed_max_anchor_span_bp"],
                "observed_max_context_span_bp":
                    result["observed_max_context_span_bp"],
                "rows": result["rows"],
                "anchors": result["anchors"],
            })


def build(arguments: argparse.Namespace) -> None:
    if shutil.which(arguments.duckdb) is None:
        raise ContextMaximaError(f"DuckDB executable not found: {arguments.duckdb}")
    arguments.anchor_parquet = arguments.anchor_parquet.expanduser().resolve()
    arguments.output = arguments.output.expanduser().resolve()
    if not arguments.anchor_parquet.is_file():
        raise ContextMaximaError(
            f"anchor Parquet file not found: {arguments.anchor_parquet}"
        )
    if arguments.output.exists():
        raise ContextMaximaError(f"output already exists: {arguments.output}")
    run_config = Path(f"{arguments.output}.run_config.tsv")
    if run_config.exists():
        raise ContextMaximaError(f"run config already exists: {run_config}")
    cofactors = parse_cofactors(arguments.cofactor)

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    temporary_parent = None
    if arguments.temp_directory is not None:
        temporary_parent = arguments.temp_directory.expanduser().resolve()
        temporary_parent.mkdir(parents=True, exist_ok=True)
    working_directory = Path(tempfile.mkdtemp(
        prefix="jaspar-context-maxima-", dir=temporary_parent
    ))
    descriptor, staging_name = tempfile.mkstemp(
        prefix=f".{arguments.output.name}.",
        suffix=".tmp",
        dir=arguments.output.parent,
    )
    os.close(descriptor)
    staging_output = Path(staging_name)
    staging_output.unlink()
    staging_run_config = Path(f"{staging_output}.run_config.tsv")
    try:
        sql = build_sql(arguments, cofactors, staging_output, working_directory)
        run_duckdb(arguments.duckdb, sql)
        result = validate_output(
            arguments.duckdb, staging_output, len(cofactors)
        )
        write_run_config(staging_run_config, arguments, cofactors, result)
        os.replace(staging_output, arguments.output)
        os.replace(staging_run_config, run_config)
    finally:
        if staging_output.exists():
            staging_output.unlink()
        if staging_run_config.exists():
            staging_run_config.unlink()
        shutil.rmtree(working_directory, ignore_errors=True)

    print(
        "I: Wrote "
        f"{result['rows']} rows for {result['motifs']} motifs and "
        f"{result['anchors']} anchors: {arguments.output}",
        file=sys.stderr,
    )


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Derive each cofactor's best sparse motif score within an interval "
            "edge-distance flank of every anchor."
        )
    )
    parser.add_argument(
        "--anchor-parquet", required=True, type=Path,
        help="single-chromosome TP73 anchor/evidence Parquet",
    )
    parser.add_argument(
        "--cofactor",
        required=True,
        action="append",
        nargs=3,
        metavar=("MOTIF_ID", "PLUS_PARQUET", "MINUS_PARQUET"),
        help=(
            "stable motif accession and its exact, single-chromosome plus/minus "
            "sparse Parquet files; repeat"
        ),
    )
    parser.add_argument(
        "--output", required=True, type=Path,
        help="new rectangular long-form Parquet; existing output is refused",
    )
    parser.add_argument(
        "--flank", type=int, default=150,
        help="maximum signed interval-edge distance in bp (default: 150)",
    )
    parser.add_argument(
        "--source-score-floor", type=float, default=-1.0,
        help="inclusive storage floor shared by source motif scans (default: -1)",
    )
    parser.add_argument("--duckdb", default="duckdb", help="DuckDB CLI command")
    parser.add_argument("--threads", type=int, default=4, help="DuckDB threads")
    parser.add_argument(
        "--memory-limit", default="4GB", help="DuckDB memory limit"
    )
    parser.add_argument(
        "--max-temp-size", default="100GB",
        help="maximum DuckDB spill size (default: 100GB)",
    )
    parser.add_argument(
        "--temp-directory", type=Path,
        help="parent for a private temporary DuckDB spill directory",
    )
    return parser


def main() -> int:
    parser = argument_parser()
    arguments = parser.parse_args()
    if arguments.flank < 0:
        parser.error("--flank must be non-negative")
    if arguments.threads <= 0:
        parser.error("--threads must be positive")
    if not (-float("inf") < arguments.source_score_floor < float("inf")):
        parser.error("--source-score-floor must be finite")
    try:
        build(arguments)
        return 0
    except (ContextMaximaError, OSError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3

"""Build per-anchor cofactor maxima and optional thresholded locus counts."""

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


def motif_id_list_sql(cofactors: list[tuple[str, Path, Path]]) -> str:
    return ",".join(sql_string(motif_id) for motif_id, _, _ in cofactors)


def threshold_source_sql(
    arguments: argparse.Namespace,
    cofactors: list[tuple[str, Path, Path]],
) -> str:
    if arguments.threshold_parquet is None:
        return f"""
SELECT
    motif_id,
    motif_id AS motif_name,
    NULL::VARCHAR AS threshold_set_id,
    NULL::VARCHAR AS genome_id,
    NULL::VARCHAR AS motif_set_id,
    NULL::VARCHAR AS threshold_role,
    NULL::VARCHAR AS target_motif_id,
    NULL::VARCHAR AS score_mode,
    NULL::DOUBLE AS pseudocount,
    NULL::VARCHAR AS background_model_id,
    NULL::VARCHAR AS pseudocount_scheme,
    NULL::VARCHAR AS calibration_stratum_id,
    NULL::DOUBLE AS recommended_threshold,
    NULL::VARCHAR AS calibration_status,
    NULL::BIGINT AS context_min_interval_distance_bp,
    NULL::BIGINT AS context_max_interval_distance_bp,
    NULL::VARCHAR AS context_relation_filter
FROM (VALUES {motif_values_sql(cofactors)}) AS motif(motif_id)
"""
    return f"""
SELECT
    motif_id,
    motif_name,
    threshold_set_id,
    genome_id,
    motif_set_id,
    threshold_role,
    target_motif_id,
    score_mode,
    pseudocount::DOUBLE AS pseudocount,
    background_model_id,
    pseudocount_scheme,
    calibration_stratum_id,
    recommended_threshold::DOUBLE AS recommended_threshold,
    calibration_status,
    context_min_interval_distance_bp::BIGINT
        AS context_min_interval_distance_bp,
    context_max_interval_distance_bp::BIGINT
        AS context_max_interval_distance_bp,
    context_relation_filter
FROM read_parquet({sql_string(arguments.threshold_parquet)})
WHERE threshold_set_id = {sql_string(arguments.threshold_set_id)}
  AND threshold_role = {sql_string(arguments.threshold_role)}
  AND calibration_stratum_id = {sql_string(arguments.calibration_stratum_id)}
  AND motif_id IN ({motif_id_list_sql(cofactors)})
  AND recommended_threshold IS NOT NULL
"""


def build_sql(arguments: argparse.Namespace, cofactors: list[tuple[str, Path, Path]],
              staging_output: Path, temp_directory: Path) -> str:
    return f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;
SET temp_directory={sql_string(temp_directory)};
SET max_temp_directory_size={sql_string(arguments.max_temp_size)};

COPY (
WITH threshold AS (
    {threshold_source_sql(arguments, cofactors)}
),
anchor_input AS (
    SELECT DISTINCT
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
candidate_geometry AS (
    SELECT
        a.chrom,
        a.anchor_start,
        a.anchor_end,
        h.motif_id,
        h.start AS hit_start,
        h.hit_end,
        h.score,
        GREATEST(a.anchor_start, h.start)
            - LEAST(a.anchor_end, h.hit_end)
            AS interval_distance_bp,
        CASE
            WHEN a.anchor_start = h.start AND a.anchor_end = h.hit_end
                THEN 'identical_span'
            WHEN a.anchor_start <= h.start AND a.anchor_end >= h.hit_end
                THEN 'anchor_contains_neighbor'
            WHEN h.start <= a.anchor_start AND h.hit_end >= a.anchor_end
                THEN 'neighbor_contains_anchor'
            WHEN GREATEST(a.anchor_start, h.start)
                 - LEAST(a.anchor_end, h.hit_end) < 0
                THEN 'partial_overlap'
            WHEN GREATEST(a.anchor_start, h.start)
                 - LEAST(a.anchor_end, h.hit_end) = 0
                THEN 'abutting'
            ELSE 'disjoint'
        END AS interval_relation
    FROM anchors a
    JOIN candidate_hit h
      ON a.capture_bin = h.anchor_capture_bin
    CROSS JOIN parameter p
    WHERE GREATEST(a.anchor_start, h.start)
          - LEAST(a.anchor_end, h.hit_end) <= p.context_flank_bp
),
context_locus AS (
    SELECT
        chrom,
        anchor_start,
        anchor_end,
        motif_id,
        hit_start,
        hit_end,
        MAX(score)::FLOAT AS locus_score,
        MAX(interval_distance_bp)::BIGINT AS interval_distance_bp,
        MAX(interval_relation) AS interval_relation
    FROM candidate_geometry
    GROUP BY chrom, anchor_start, anchor_end, motif_id, hit_start, hit_end
),
features AS (
    SELECT
        a.chrom,
        a.anchor_start,
        a.anchor_end,
        t.*,
        MAX(l.locus_score)::FLOAT AS context_score,
        COUNT(*) FILTER (
            WHERE l.hit_start IS NOT NULL
              AND l.locus_score >= t.recommended_threshold
              AND (t.context_min_interval_distance_bp IS NULL
                   OR l.interval_distance_bp >=
                      t.context_min_interval_distance_bp)
              AND (t.context_max_interval_distance_bp IS NULL
                   OR l.interval_distance_bp <=
                      t.context_max_interval_distance_bp)
              AND (t.context_relation_filter = 'any'
                   OR l.interval_relation = t.context_relation_filter)
        )::BIGINT AS qualifying_locus_count
    FROM anchors a
    CROSS JOIN threshold t
    LEFT JOIN context_locus l
      ON a.chrom = l.chrom
     AND a.anchor_start = l.anchor_start
     AND a.anchor_end = l.anchor_end
     AND t.motif_id = l.motif_id
    GROUP BY ALL
)
SELECT
    1::INTEGER AS schema_version,
    f.chrom,
    CASE WHEN f.target_motif_id IS NULL THEN NULL
         ELSE md5(concat_ws('|', f.genome_id, f.motif_set_id, f.chrom,
                            f.anchor_start::VARCHAR, f.anchor_end::VARCHAR,
                            f.target_motif_id, f.score_mode,
                            printf('%.17g', f.pseudocount)))
    END AS anchor_locus_id,
    f.anchor_start,
    f.anchor_end,
    f.motif_id,
    f.motif_name,
    f.context_score,
    CASE WHEN f.recommended_threshold IS NULL THEN NULL
         ELSE f.qualifying_locus_count
    END AS n_neighbor_loci_above_threshold,
    CASE WHEN f.recommended_threshold IS NULL THEN NULL
         ELSE f.qualifying_locus_count > 0
    END AS has_neighbor_locus_above_threshold,
    f.threshold_set_id,
    f.genome_id,
    f.motif_set_id,
    f.threshold_role,
    f.target_motif_id,
    f.score_mode,
    f.pseudocount,
    f.background_model_id,
    f.pseudocount_scheme,
    f.calibration_stratum_id,
    f.recommended_threshold,
    f.calibration_status,
    f.context_min_interval_distance_bp,
    f.context_max_interval_distance_bp,
    f.context_relation_filter,
    {arguments.source_score_floor}::DOUBLE AS source_score_floor,
    p.context_flank_bp,
    p.capture_prefilter_center_bp,
    p.max_anchor_span_bp AS observed_max_anchor_span_bp,
    p.max_context_span_bp AS observed_max_context_span_bp,
    'signed_interval_edge_distance'::VARCHAR AS context_distance_metric
FROM features f
CROSS JOIN parameter p
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


def validate_threshold_input(
    arguments: argparse.Namespace,
    cofactors: list[tuple[str, Path, Path]],
) -> None:
    if arguments.threshold_parquet is None:
        return
    query = f"""
WITH expected(motif_id) AS (
    VALUES {motif_values_sql(cofactors)}
),
selected AS (
    SELECT *
    FROM read_parquet({sql_string(arguments.threshold_parquet)})
    WHERE threshold_set_id = {sql_string(arguments.threshold_set_id)}
      AND threshold_role = {sql_string(arguments.threshold_role)}
      AND calibration_stratum_id =
          {sql_string(arguments.calibration_stratum_id)}
      AND motif_id IN ({motif_id_list_sql(cofactors)})
      AND recommended_threshold IS NOT NULL
)
SELECT
    (SELECT COUNT(*) FROM selected)::BIGINT AS selected_rows,
    (SELECT COUNT(DISTINCT motif_id) FROM selected)::BIGINT
        AS selected_motifs,
    (SELECT COUNT(*) FROM expected e
     LEFT JOIN selected s USING (motif_id)
     WHERE s.motif_id IS NULL)::BIGINT AS missing_motifs,
    COUNT(*) FILTER (
        WHERE motif_name IS NULL
           OR motif_id = target_motif_id
           OR genome_id IS NULL
           OR motif_set_id IS NULL
           OR target_motif_id IS NULL
           OR score_mode IS NULL
           OR pseudocount IS NULL
           OR background_model_id IS NULL
           OR pseudocount_scheme IS NULL
           OR calibration_status IS NULL
           OR threshold_inclusive IS DISTINCT FROM TRUE
           OR recommended_threshold IS NULL
           OR NOT isfinite(recommended_threshold)
           OR context_distance_metric <>
              'signed_interval_edge_distance'
           OR context_max_interval_distance_bp IS NULL
           OR context_max_interval_distance_bp > {arguments.flank}
           OR (context_min_interval_distance_bp IS NOT NULL
               AND context_min_interval_distance_bp >
                   context_max_interval_distance_bp)
           OR context_relation_filter IS NULL
           OR context_relation_filter NOT IN (
               'any', 'identical_span', 'anchor_contains_neighbor',
               'neighbor_contains_anchor', 'partial_overlap', 'abutting',
               'disjoint'
           )
           OR source_minimum_score IS NULL
           OR ABS(source_minimum_score -
                  {arguments.source_score_floor}) > 1e-9
           OR recommended_threshold < source_minimum_score
    )::BIGINT AS invalid_rows,
    COUNT(DISTINCT genome_id)::BIGINT AS genomes,
    COUNT(DISTINCT motif_set_id)::BIGINT AS motif_sets,
    COUNT(DISTINCT target_motif_id)::BIGINT AS target_motifs,
    COUNT(DISTINCT score_mode)::BIGINT AS score_modes,
    COUNT(DISTINCT pseudocount)::BIGINT AS pseudocounts,
    COUNT(DISTINCT background_model_id)::BIGINT AS backgrounds,
    COUNT(DISTINCT pseudocount_scheme)::BIGINT AS pseudocount_schemes
FROM selected;
"""
    process = subprocess.run(
        [arguments.duckdb, "-light-mode", "-json", ":memory:", "-c", query],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ContextMaximaError(
            process.stderr.strip() or "DuckDB threshold validation failed"
        )
    try:
        rows = json.loads(process.stdout)
    except json.JSONDecodeError as error:
        raise ContextMaximaError(
            "DuckDB returned invalid threshold-validation JSON"
        ) from error
    if len(rows) != 1:
        raise ContextMaximaError("DuckDB returned an invalid threshold row count")
    result = {key: int(value) for key, value in rows[0].items()}
    expected = len(cofactors)
    if (result["selected_rows"] != expected
            or result["selected_motifs"] != expected
            or result["missing_motifs"] != 0):
        raise ContextMaximaError(
            "threshold selection must contain exactly one populated row for "
            "every requested motif"
        )
    if result["invalid_rows"] != 0:
        raise ContextMaximaError(
            "threshold rows are incompatible with the source floor, interval "
            "geometry, or inclusive-count contract"
        )
    identity_counts = (
        "genomes", "motif_sets", "target_motifs", "score_modes",
        "pseudocounts", "backgrounds", "pseudocount_schemes",
    )
    if any(result[key] != 1 for key in identity_counts):
        raise ContextMaximaError(
            "threshold rows do not share one target and scoring configuration"
        )


def validate_output(
    duckdb: str,
    output: Path,
    expected_motifs: int,
    threshold_configured: bool,
) -> dict[str, int]:
    query = f"""
WITH output AS (SELECT * FROM read_parquet({sql_string(output)}))
SELECT
    COUNT(*)::BIGINT AS rows,
    COUNT(DISTINCT (chrom, anchor_start, anchor_end))::BIGINT AS anchors,
    COUNT(DISTINCT chrom)::BIGINT AS chromosomes,
    COUNT(DISTINCT motif_id)::BIGINT AS motifs,
    COUNT(DISTINCT schema_version)::BIGINT AS schema_versions,
    MIN(schema_version)::BIGINT AS schema_version,
    COUNT(*) - COUNT(DISTINCT (chrom, anchor_start, anchor_end, motif_id))
        AS duplicate_keys,
    COUNT(*) FILTER (
        WHERE context_score IS NOT NULL AND context_score < source_score_floor
    ) AS scores_below_source_floor,
    COALESCE(SUM(COALESCE(n_neighbor_loci_above_threshold, 0)), 0)::BIGINT
        AS thresholded_loci,
    COUNT(*) FILTER (WHERE n_neighbor_loci_above_threshold IS NOT NULL)::BIGINT
        AS populated_count_rows,
    COUNT(*) FILTER (
        WHERE n_neighbor_loci_above_threshold > 0
          AND (context_score IS NULL
               OR context_score < recommended_threshold)
    )::BIGINT AS invalid_positive_counts,
    COUNT(*) FILTER (
        WHERE threshold_set_id IS NULL OR anchor_locus_id IS NULL
    )::BIGINT AS missing_threshold_identity,
    COUNT(DISTINCT threshold_set_id)::BIGINT AS threshold_sets,
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
    if result["schema_versions"] != 1 or result["schema_version"] != 1:
        raise ContextMaximaError("output has an unsupported feature schema version")
    if result["rows"] != result["anchors"] * expected_motifs:
        raise ContextMaximaError("output is not rectangular by anchor and motif")
    if result["duplicate_keys"] != 0 or result["scores_below_source_floor"] != 0:
        raise ContextMaximaError("output failed key or source-score-floor validation")
    if any(result[key] != 1 for key in (
        "prefilter_values", "anchor_span_values", "context_span_values"
    )):
        raise ContextMaximaError("output contains inconsistent geometry provenance")
    if threshold_configured:
        if (result["populated_count_rows"] != result["rows"]
                or result["invalid_positive_counts"] != 0
                or result["missing_threshold_identity"] != 0
                or result["threshold_sets"] != 1):
            raise ContextMaximaError(
                "thresholded output failed count or identity validation"
            )
    elif result["populated_count_rows"] != 0 or result["threshold_sets"] != 0:
        raise ContextMaximaError(
            "maximum-only output unexpectedly contains thresholded counts"
        )
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
        "threshold_parquet", "threshold_set_id", "threshold_role",
        "calibration_stratum_id", "thresholded_loci",
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
                "threshold_parquet": arguments.threshold_parquet or "",
                "threshold_set_id": arguments.threshold_set_id or "",
                "threshold_role": (
                    arguments.threshold_role
                    if arguments.threshold_parquet is not None else ""
                ),
                "calibration_stratum_id": (
                    arguments.calibration_stratum_id
                    if arguments.threshold_parquet is not None else ""
                ),
                "thresholded_loci": result["thresholded_loci"],
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
    if arguments.threshold_parquet is not None:
        arguments.threshold_parquet = (
            arguments.threshold_parquet.expanduser().resolve()
        )
    if not arguments.anchor_parquet.is_file():
        raise ContextMaximaError(
            f"anchor Parquet file not found: {arguments.anchor_parquet}"
        )
    if (arguments.threshold_parquet is not None
            and not arguments.threshold_parquet.is_file()):
        raise ContextMaximaError(
            f"threshold Parquet file not found: {arguments.threshold_parquet}"
        )
    if arguments.output.exists():
        raise ContextMaximaError(f"output already exists: {arguments.output}")
    run_config = Path(f"{arguments.output}.run_config.tsv")
    if run_config.exists():
        raise ContextMaximaError(f"run config already exists: {run_config}")
    cofactors = parse_cofactors(arguments.cofactor)
    validate_threshold_input(arguments, cofactors)

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
            arguments.duckdb,
            staging_output,
            len(cofactors),
            arguments.threshold_parquet is not None,
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
        f"{result['anchors']} anchors"
        + (
            f", counting {result['thresholded_loci']} threshold-qualified loci"
            if arguments.threshold_parquet is not None else ""
        )
        + f": {arguments.output}",
        file=sys.stderr,
    )


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Derive each cofactor's best sparse motif score within an interval "
            "edge-distance flank of every anchor. With a threshold registry, "
            "also count distinct physical motif loci meeting each motif's "
            "inclusive convenient threshold."
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
    parser.add_argument(
        "--threshold-parquet", type=Path,
        help=(
            "optional motif_score_threshold Parquet; requires "
            "--threshold-set-id and emits zero-complete per-anchor locus counts"
        ),
    )
    parser.add_argument(
        "--threshold-set-id",
        help="threshold_set_id selected from --threshold-parquet",
    )
    parser.add_argument(
        "--threshold-role", default="tp73_context_binary_feature",
        help=(
            "threshold role selected from the registry "
            "(default: tp73_context_binary_feature)"
        ),
    )
    parser.add_argument(
        "--calibration-stratum-id", default="all_tp73_anchors",
        help=(
            "calibration stratum selected from the registry "
            "(default: all_tp73_anchors)"
        ),
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
    if ((arguments.threshold_parquet is None)
            != (arguments.threshold_set_id is None)):
        parser.error(
            "--threshold-parquet and --threshold-set-id must be supplied together"
        )
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

#!/usr/bin/env python3
"""Aggregate one motif/chromosome into TP73 distance-band test components."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


DISTANCE_BANDS = (
    ("overlap", 0),
    ("adjacent_0_5", 1),
    ("gap_6_20", 2),
    ("gap_21_50", 3),
    ("gap_51_100", 4),
    ("gap_101_150", 5),
)


class DistanceCountError(RuntimeError):
    """Raised when the compact distance-count build cannot be validated."""


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_number(value: float) -> str:
    return format(value, ".17g")


def band_values_sql() -> str:
    return ", ".join(
        f"({sql_string(name)}, {order})" for name, order in DISTANCE_BANDS
    )


def build_sql(
    arguments: argparse.Namespace,
    block_output: Path,
    class_output: Path,
    temp_directory: Path,
) -> str:
    threshold = sql_number(arguments.positive_threshold)
    flank = arguments.flank
    return f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;
SET temp_directory={sql_string(temp_directory)};
SET max_temp_directory_size={sql_string(arguments.max_temp_size)};

CREATE TEMP TABLE anchor AS
SELECT
    CAST(chrom AS VARCHAR) AS chrom,
    anchor_start::BIGINT AS anchor_start,
    anchor_end::BIGINT AS anchor_end,
    anchor_score::DOUBLE AS anchor_score,
    FLOOR(anchor_start / {arguments.block_size})::BIGINT AS block_index,
    CASE
      WHEN anchor_score < -1 THEN '[-5,-1)'
      WHEN anchor_score < 0 THEN '[-1,0)'
      WHEN anchor_score < 1 THEN '[0,1)'
      WHEN anchor_score < 2 THEN '[1,2)'
      WHEN anchor_score < 5 THEN '[2,5)'
      WHEN anchor_score < 10 THEN '[5,10)'
      WHEN anchor_score < 15 THEN '[10,15)'
      ELSE '[15,+Inf)'
    END::VARCHAR AS tp73_score_stratum,
    supported_tp73_saos2_TA::BOOLEAN AS anti_saos2_TA,
    supported_negative_control_saos2_TA::BOOLEAN AS control_saos2_TA,
    supported_tp73_saos2_DN::BOOLEAN AS anti_saos2_DN,
    supported_negative_control_saos2_DN::BOOLEAN AS control_saos2_DN,
    supported_tp73_skmel29_2_TA::BOOLEAN AS anti_skmel29_2_TA,
    supported_negative_control_skmel29_2_TA::BOOLEAN
        AS control_skmel29_2_TA,
    supported_tp73_skmel29_2_DN::BOOLEAN AS anti_skmel29_2_DN,
    supported_negative_control_skmel29_2_DN::BOOLEAN
        AS control_skmel29_2_DN
FROM read_parquet({sql_string(arguments.anchors)})
WHERE CAST(chrom AS VARCHAR) = {sql_string(arguments.chrom)};

SELECT CASE WHEN (SELECT count(*) FROM anchor) = 0
    THEN error('anchor chromosome is empty') END;
SELECT CASE WHEN (SELECT count(*) - count(DISTINCT (anchor_start, anchor_end))
                  FROM anchor) <> 0
    THEN error('anchor chromosome has duplicate spans') END;

CREATE TEMP TABLE physical_hit AS
SELECT start::BIGINT AS start, "end"::BIGINT AS hit_end,
       MAX(score::DOUBLE) AS score
FROM (
    SELECT start, "end", score
    FROM read_parquet({sql_string(arguments.plus_hits)})
    UNION ALL
    SELECT start, "end", score
    FROM read_parquet({sql_string(arguments.minus_hits)})
) hit
GROUP BY start, "end";

CREATE TEMP TABLE band_feature AS
WITH parameter AS (
    SELECT GREATEST(
               1,
               {flank} +
               (SELECT MAX(anchor_end - anchor_start) FROM anchor) +
               COALESCE((SELECT MAX(hit_end - start) FROM physical_hit), 0)
           )::BIGINT AS capture_width
), anchor_bin AS (
    SELECT a.*,
           FLOOR(((anchor_start + anchor_end) / 2.0) / p.capture_width)::BIGINT
               AS capture_bin
    FROM anchor a CROSS JOIN parameter p
), hit_bin AS (
    SELECT h.*,
           FLOOR(((start + hit_end) / 2.0) / p.capture_width)::BIGINT
               + offset_value AS anchor_capture_bin
    FROM physical_hit h
    CROSS JOIN parameter p
    CROSS JOIN (VALUES (-1), (0), (1)) offset_table(offset_value)
), geometry AS (
    SELECT
        a.anchor_start,
        a.anchor_end,
        h.start AS hit_start,
        h.hit_end,
        h.score,
        GREATEST(a.anchor_start, h.start)
            - LEAST(a.anchor_end, h.hit_end) AS interval_distance_bp
    FROM anchor_bin a
    JOIN hit_bin h ON h.anchor_capture_bin = a.capture_bin
    WHERE GREATEST(a.anchor_start, h.start)
          - LEAST(a.anchor_end, h.hit_end) <= {flank}
), classified AS (
    SELECT *,
        CASE
          WHEN interval_distance_bp < 0 THEN 'overlap'
          WHEN interval_distance_bp <= 5 THEN 'adjacent_0_5'
          WHEN interval_distance_bp <= 20 THEN 'gap_6_20'
          WHEN interval_distance_bp <= 50 THEN 'gap_21_50'
          WHEN interval_distance_bp <= 100 THEN 'gap_51_100'
          ELSE 'gap_101_150'
        END::VARCHAR AS distance_band
    FROM geometry
)
SELECT anchor_start, anchor_end, distance_band, MAX(score)::DOUBLE AS best_score
FROM classified
GROUP BY anchor_start, anchor_end, distance_band;

COPY (
WITH bands(distance_band, distance_band_order) AS (
    VALUES {band_values_sql()}
)
SELECT
    {sql_string(arguments.motif_id)}::VARCHAR AS motif_id,
    {sql_string(arguments.motif_name)}::VARCHAR AS motif_name,
    {sql_string(arguments.chrom)}::VARCHAR AS chrom,
    b.distance_band,
    b.distance_band_order,
    -1.0::DOUBLE AS source_score_floor,
    {threshold}::DOUBLE AS positive_threshold,
    count(*)::BIGINT AS anchors_total,
    count(*) FILTER (WHERE f.anchor_start IS NOT NULL)::BIGINT
        AS anchors_source_present,
    count(*) FILTER (WHERE f.best_score >= {threshold})::BIGINT
        AS anchors_positive,
    (count(*) FILTER (WHERE f.anchor_start IS NOT NULL)
       - count(*) FILTER (WHERE f.best_score >= {threshold}))::BIGINT
        AS anchors_intermediate,
    count(*) FILTER (WHERE f.anchor_start IS NULL)::BIGINT AS anchors_negative
FROM anchor a
CROSS JOIN bands b
LEFT JOIN band_feature f
  ON f.anchor_start = a.anchor_start
 AND f.anchor_end = a.anchor_end
 AND f.distance_band = b.distance_band
GROUP BY b.distance_band, b.distance_band_order
ORDER BY b.distance_band_order
) TO {sql_string(class_output)} (FORMAT PARQUET, COMPRESSION ZSTD);

COPY (
WITH bands(distance_band, distance_band_order) AS (
    VALUES {band_values_sql()}
), discordant AS (
    SELECT
        a.anchor_start,
        a.anchor_end,
        a.block_index,
        a.tp73_score_stratum,
        sample.series_id,
        sample.isoform,
        sample.anti_supported,
        sample.control_supported
    FROM anchor a
    CROSS JOIN LATERAL (VALUES
      ('saos2', 'TA', a.anti_saos2_TA, a.control_saos2_TA),
      ('saos2', 'DN', a.anti_saos2_DN, a.control_saos2_DN),
      ('skmel29_2', 'TA', a.anti_skmel29_2_TA, a.control_skmel29_2_TA),
      ('skmel29_2', 'DN', a.anti_skmel29_2_DN, a.control_skmel29_2_DN)
    ) sample(series_id, isoform, anti_supported, control_supported)
    WHERE sample.anti_supported <> sample.control_supported
), baseline AS (
    SELECT
        block_index,
        tp73_score_stratum,
        series_id,
        isoform,
        count(*) FILTER (WHERE anti_supported)::BIGINT AS total_anti,
        count(*) FILTER (WHERE control_supported)::BIGINT AS total_control
    FROM discordant
    GROUP BY ALL
), observed AS (
    SELECT
        d.block_index,
        d.tp73_score_stratum,
        d.series_id,
        d.isoform,
        f.distance_band,
        count(*) FILTER (WHERE d.anti_supported)::BIGINT AS source_anti,
        count(*) FILTER (WHERE d.control_supported)::BIGINT AS source_control,
        count(*) FILTER (
            WHERE d.anti_supported AND f.best_score >= {threshold}
        )::BIGINT AS positive_anti,
        count(*) FILTER (
            WHERE d.control_supported AND f.best_score >= {threshold}
        )::BIGINT AS positive_control
    FROM discordant d
    JOIN band_feature f USING (anchor_start, anchor_end)
    GROUP BY ALL
), cell AS (
    SELECT
        b.block_index,
        b.tp73_score_stratum,
        b.series_id,
        b.isoform,
        band.distance_band,
        band.distance_band_order,
        b.total_anti,
        b.total_control,
        COALESCE(o.source_anti, 0)::BIGINT AS source_anti,
        COALESCE(o.source_control, 0)::BIGINT AS source_control,
        COALESCE(o.positive_anti, 0)::BIGINT AS positive_anti,
        COALESCE(o.positive_control, 0)::BIGINT AS positive_control,
        (b.total_anti - COALESCE(o.source_anti, 0))::BIGINT AS negative_anti,
        (b.total_control - COALESCE(o.source_control, 0))::BIGINT
            AS negative_control
    FROM baseline b
    CROSS JOIN bands band
    LEFT JOIN observed o
      ON o.block_index = b.block_index
     AND o.tp73_score_stratum = b.tp73_score_stratum
     AND o.series_id = b.series_id
     AND o.isoform = b.isoform
     AND o.distance_band = band.distance_band
), component AS (
    SELECT *,
        (positive_anti + positive_control + negative_anti + negative_control)
            ::DOUBLE AS stratum_total,
        CASE WHEN stratum_total > 0
          THEN positive_anti * negative_control / stratum_total
          ELSE 0 END::DOUBLE AS mh_numerator,
        CASE WHEN stratum_total > 0
          THEN positive_control * negative_anti / stratum_total
          ELSE 0 END::DOUBLE AS mh_denominator
    FROM cell
)
SELECT
    {sql_string(arguments.motif_id)}::VARCHAR AS motif_id,
    {sql_string(arguments.motif_name)}::VARCHAR AS motif_name,
    {sql_string(arguments.chrom)}::VARCHAR AS chrom,
    block_index,
    series_id,
    isoform,
    distance_band,
    distance_band_order,
    -1.0::DOUBLE AS source_score_floor,
    {threshold}::DOUBLE AS positive_threshold,
    SUM(mh_numerator)::DOUBLE AS mh_numerator,
    SUM(mh_denominator)::DOUBLE AS mh_denominator,
    SUM(positive_anti)::BIGINT AS positive_anti_discordant,
    SUM(positive_control)::BIGINT AS positive_control_discordant,
    SUM(negative_anti)::BIGINT AS negative_anti_discordant,
    SUM(negative_control)::BIGINT AS negative_control_discordant,
    SUM(source_anti)::BIGINT AS source_anti_discordant,
    SUM(source_control)::BIGINT AS source_control_discordant,
    SUM(total_anti)::BIGINT AS total_anti_discordant,
    SUM(total_control)::BIGINT AS total_control_discordant,
    count(*)::BIGINT AS tp73_score_strata
FROM component
GROUP BY block_index, series_id, isoform,
         distance_band, distance_band_order
ORDER BY block_index, series_id, isoform, distance_band_order
) TO {sql_string(block_output)} (FORMAT PARQUET, COMPRESSION ZSTD);
"""


def run_duckdb(executable: str, sql: str) -> None:
    process = subprocess.run(
        [executable, "-light-mode", "-batch", ":memory:"],
        input=sql,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise DistanceCountError(
            process.stderr.strip() or "DuckDB distance-count build failed"
        )


def build(arguments: argparse.Namespace) -> None:
    inputs = (arguments.anchors, arguments.plus_hits, arguments.minus_hits)
    for path in inputs:
        if not path.is_file():
            raise DistanceCountError(f"input file is missing: {path}")
    if shutil.which(arguments.duckdb) is None:
        raise DistanceCountError(f"DuckDB executable not found: {arguments.duckdb}")
    for output in (arguments.block_output, arguments.class_output):
        if output.exists():
            raise DistanceCountError(f"refusing to replace output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)

    temporary_parent = arguments.temp_directory
    temporary_parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix="tp73-distance-", dir=temporary_parent))
    try:
        staged_block = staging / "block_components.parquet"
        staged_class = staging / "class_counts.parquet"
        sql = build_sql(
            arguments, staged_block, staged_class, staging / "duckdb-spill"
        )
        run_duckdb(arguments.duckdb, sql)
        os.replace(staged_block, arguments.block_output)
        os.replace(staged_class, arguments.class_output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--anchors", type=Path, required=True)
    result.add_argument("--plus-hits", type=Path, required=True)
    result.add_argument("--minus-hits", type=Path, required=True)
    result.add_argument("--motif-id", required=True)
    result.add_argument("--motif-name", required=True)
    result.add_argument("--chrom", required=True)
    result.add_argument("--positive-threshold", type=float, required=True)
    result.add_argument("--block-output", type=Path, required=True)
    result.add_argument("--class-output", type=Path, required=True)
    result.add_argument("--block-size", type=int, default=5_000_000)
    result.add_argument("--flank", type=int, default=150)
    result.add_argument("--threads", type=int, default=2)
    result.add_argument("--memory-limit", default="16GB")
    result.add_argument("--max-temp-size", default="100GB")
    result.add_argument(
        "--temp-directory", type=Path, default=Path(tempfile.gettempdir())
    )
    result.add_argument("--duckdb", default="duckdb")
    return result


def main() -> int:
    arguments = parser().parse_args()
    if arguments.block_size <= 2 * arguments.flank:
        print("E: --block-size must exceed twice --flank", file=sys.stderr)
        return 2
    if arguments.flank != 150:
        print("E: this schema requires --flank 150", file=sys.stderr)
        return 2
    if arguments.threads <= 0:
        print("E: --threads must be positive", file=sys.stderr)
        return 2
    if not (arguments.positive_threshold >= -1):
        print("E: --positive-threshold must be at least -1", file=sys.stderr)
        return 2
    try:
        build(arguments)
    except (DistanceCountError, OSError, ValueError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

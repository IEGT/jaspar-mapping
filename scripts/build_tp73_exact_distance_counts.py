#!/usr/bin/env python3
"""Build exact TP73-to-motif distance counts for one chromosome and motif.

The output is deliberately compact: it contains a zero-complete conditional
distance inventory and 5 Mb block sufficient statistics, not anchor-level
rows.  Motif centers are represented in doubled-base coordinates so that
half-base offsets caused by odd/even motif-length combinations remain exact.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


class ExactDistanceError(RuntimeError):
    """Raised when exact-distance inputs or outputs violate their contract."""


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_number(value: float) -> str:
    return format(value, ".17g")


def build_sql(
    arguments: argparse.Namespace,
    inventory_output: Path,
    block_output: Path,
    temp_directory: Path,
) -> str:
    threshold = sql_number(arguments.positive_threshold)
    source_floor = sql_number(arguments.source_score_floor)
    exclude_identical = (
        "true" if arguments.motif_id == arguments.anchor_motif_id else "false"
    )
    return f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;
SET temp_directory={sql_string(temp_directory)};
SET max_temp_directory_size={sql_string(arguments.max_temp_size)};

CREATE TEMP TABLE anchor_evidence AS
SELECT
    CAST(chrom AS VARCHAR) AS chrom,
    anchor_start::BIGINT AS anchor_start,
    anchor_end::BIGINT AS anchor_end,
    anchor_score::DOUBLE AS anchor_score,
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
FROM read_parquet({sql_string(arguments.anchors)}, hive_partitioning=false)
WHERE CAST(chrom AS VARCHAR) = {sql_string(arguments.chrom)};

SELECT CASE WHEN (SELECT count(*) FROM anchor_evidence) = 0
    THEN error('anchor chromosome is empty') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM anchor_evidence
    WHERE anchor_start < 0 OR anchor_end <= anchor_start
) THEN error('anchor evidence contains an invalid interval') END;
SELECT CASE WHEN (SELECT count(*) - count(DISTINCT (anchor_start, anchor_end))
                  FROM anchor_evidence) <> 0
    THEN error('anchor chromosome has duplicate spans') END;

CREATE TEMP TABLE anchor_orientation AS
WITH orientation_record AS (
    SELECT start::BIGINT AS start, "end"::BIGINT AS hit_end,
           score::DOUBLE AS score, '+'::VARCHAR AS strand
    FROM read_parquet({sql_string(arguments.anchor_plus)},
                      hive_partitioning=false)
    UNION ALL
    SELECT start::BIGINT, "end"::BIGINT, score::DOUBLE, '-'::VARCHAR
    FROM read_parquet({sql_string(arguments.anchor_minus)},
                      hive_partitioning=false)
)
SELECT start, hit_end,
       max(score) FILTER (WHERE strand = '+')::DOUBLE AS plus_score,
       max(score) FILTER (WHERE strand = '-')::DOUBLE AS minus_score
FROM orientation_record
GROUP BY start, hit_end;

CREATE TEMP TABLE anchor AS
SELECT
    e.*,
    FLOOR(e.anchor_start / {arguments.block_size})::BIGINT AS block_index,
    CASE
      WHEN e.anchor_score < -1 THEN '[-5,-1)'
      WHEN e.anchor_score < 0 THEN '[-1,0)'
      WHEN e.anchor_score < 1 THEN '[0,1)'
      WHEN e.anchor_score < 2 THEN '[1,2)'
      WHEN e.anchor_score < 5 THEN '[2,5)'
      WHEN e.anchor_score < 10 THEN '[5,10)'
      WHEN e.anchor_score < 15 THEN '[10,15)'
      ELSE '[15,+Inf)'
    END::VARCHAR AS tp73_score_stratum,
    o.plus_score AS anchor_plus_score,
    o.minus_score AS anchor_minus_score,
    CASE
      WHEN o.plus_score IS NOT NULL AND o.minus_score IS NOT NULL
        THEN 'ambiguous'
      WHEN o.plus_score IS NOT NULL THEN 'plus_only'
      WHEN o.minus_score IS NOT NULL THEN 'minus_only'
      ELSE 'missing'
    END::VARCHAR AS anchor_orientation_state,
    CASE
      WHEN o.plus_score IS NOT NULL AND o.minus_score IS NULL THEN '+'
      WHEN o.minus_score IS NOT NULL AND o.plus_score IS NULL THEN '-'
      WHEN o.plus_score > o.minus_score THEN '+'
      WHEN o.minus_score > o.plus_score THEN '-'
      ELSE NULL
    END::VARCHAR AS anchor_best_strand,
    (e.anchor_start + e.anchor_end)::BIGINT AS anchor_center_twice
FROM anchor_evidence e
LEFT JOIN anchor_orientation o
  ON o.start = e.anchor_start AND o.hit_end = e.anchor_end;

SELECT CASE WHEN EXISTS (
    SELECT 1 FROM anchor WHERE anchor_orientation_state = 'missing'
) THEN error('an anchor span is absent from the orientation-specific scans') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM anchor
    WHERE abs(anchor_score - greatest(
        coalesce(anchor_plus_score, '-Infinity'::DOUBLE),
        coalesce(anchor_minus_score, '-Infinity'::DOUBLE)
    )) > 1e-5
) THEN error('anchor evidence score differs from orientation scans') END;

CREATE TEMP TABLE physical_hit AS
WITH orientation_record AS (
    SELECT start::BIGINT AS start, "end"::BIGINT AS hit_end,
           score::DOUBLE AS score, '+'::VARCHAR AS strand
    FROM read_parquet({sql_string(arguments.cofactor_plus)},
                      hive_partitioning=false)
    WHERE score >= {source_floor}
    UNION ALL
    SELECT start::BIGINT, "end"::BIGINT, score::DOUBLE, '-'::VARCHAR
    FROM read_parquet({sql_string(arguments.cofactor_minus)},
                      hive_partitioning=false)
    WHERE score >= {source_floor}
), collapsed AS (
    SELECT start, hit_end,
           max(score)::DOUBLE AS score,
           max(score) FILTER (WHERE strand = '+')::DOUBLE AS plus_score,
           max(score) FILTER (WHERE strand = '-')::DOUBLE AS minus_score
    FROM orientation_record
    GROUP BY start, hit_end
)
SELECT *,
       CASE
         WHEN plus_score IS NOT NULL AND minus_score IS NOT NULL
           THEN 'ambiguous'
         WHEN plus_score IS NOT NULL THEN 'plus_only'
         ELSE 'minus_only'
       END::VARCHAR AS orientation_state,
       CASE
         WHEN plus_score IS NOT NULL AND minus_score IS NULL THEN '+'
         WHEN minus_score IS NOT NULL AND plus_score IS NULL THEN '-'
         WHEN plus_score > minus_score THEN '+'
         WHEN minus_score > plus_score THEN '-'
         ELSE NULL
       END::VARCHAR AS best_strand,
       (start + hit_end)::BIGINT AS center_twice
FROM collapsed;

SELECT CASE WHEN EXISTS (
    SELECT 1 FROM physical_hit WHERE start < 0 OR hit_end <= start
) THEN error('cofactor scan contains an invalid interval') END;

CREATE TEMP TABLE capture_parameter AS
SELECT
    max(anchor_end - anchor_start)::BIGINT AS max_anchor_span,
    coalesce((SELECT max(hit_end - start) FROM physical_hit), 0)::BIGINT
        AS max_neighbor_span,
    greatest(
        1,
        {arguments.flank} + max(anchor_end - anchor_start)
          + coalesce((SELECT max(hit_end - start) FROM physical_hit), 0)
    )::BIGINT AS capture_width,
    (2 * {arguments.flank} + max(anchor_end - anchor_start)
       + coalesce((SELECT max(hit_end - start) FROM physical_hit), 0))::BIGINT
        AS max_center_offset_twice
FROM anchor;

CREATE TEMP TABLE locus_pair AS
WITH anchor_bin AS (
    SELECT a.*,
           floor((anchor_center_twice / 2.0) / p.capture_width)::BIGINT
               AS capture_bin
    FROM anchor a CROSS JOIN capture_parameter p
), hit_bin AS (
    SELECT h.*,
           floor((center_twice / 2.0) / p.capture_width)::BIGINT
               + offset_value AS anchor_capture_bin
    FROM physical_hit h
    CROSS JOIN capture_parameter p
    CROSS JOIN (VALUES (-1), (0), (1)) offsets(offset_value)
), geometry AS (
    SELECT
        a.*,
        h.start AS hit_start,
        h.hit_end,
        h.score,
        h.orientation_state AS neighbor_orientation_state,
        h.best_strand AS neighbor_best_strand,
        (h.center_twice - a.anchor_center_twice)::BIGINT
            AS genomic_center_offset_twice,
        CASE
          WHEN a.anchor_best_strand = '+'
            THEN h.center_twice - a.anchor_center_twice
          WHEN a.anchor_best_strand = '-'
            THEN a.anchor_center_twice - h.center_twice
          ELSE NULL
        END::BIGINT AS oriented_center_offset_twice,
        (greatest(a.anchor_start, h.start)
         - least(a.anchor_end, h.hit_end))::BIGINT AS interval_distance_bp,
        CASE
          WHEN a.anchor_best_strand IS NULL OR h.best_strand IS NULL
            THEN 'ambiguous'
          WHEN a.anchor_best_strand = h.best_strand THEN 'same'
          ELSE 'opposite'
        END::VARCHAR AS relative_orientation
    FROM anchor_bin a
    JOIN hit_bin h ON h.anchor_capture_bin = a.capture_bin
    WHERE greatest(a.anchor_start, h.start)
          - least(a.anchor_end, h.hit_end) <= {arguments.flank}
      AND (NOT {exclude_identical}
           OR h.start <> a.anchor_start OR h.hit_end <> a.anchor_end)
)
SELECT * FROM geometry;

CREATE TEMP TABLE exact_feature AS
WITH framed AS (
    SELECT anchor_start, anchor_end, block_index, tp73_score_stratum,
           'genomic'::VARCHAR AS distance_frame,
           'all'::VARCHAR AS relative_orientation,
           genomic_center_offset_twice AS center_offset_twice,
           score, interval_distance_bp
    FROM locus_pair
    UNION ALL
    SELECT anchor_start, anchor_end, block_index, tp73_score_stratum,
           'tp73_oriented', 'all', oriented_center_offset_twice,
           score, interval_distance_bp
    FROM locus_pair
    WHERE oriented_center_offset_twice IS NOT NULL
    UNION ALL
    SELECT anchor_start, anchor_end, block_index, tp73_score_stratum,
           'tp73_oriented', relative_orientation,
           oriented_center_offset_twice, score, interval_distance_bp
    FROM locus_pair
    WHERE oriented_center_offset_twice IS NOT NULL
)
SELECT anchor_start, anchor_end, block_index, tp73_score_stratum,
       distance_frame, relative_orientation, center_offset_twice,
       max(score)::DOUBLE AS best_score,
       count(*)::BIGINT AS physical_locus_pairs,
       min(interval_distance_bp)::BIGINT AS minimum_interval_distance_bp,
       max(interval_distance_bp)::BIGINT AS maximum_interval_distance_bp
FROM framed
GROUP BY anchor_start, anchor_end, block_index, tp73_score_stratum,
         distance_frame, relative_orientation, center_offset_twice;

CREATE TEMP TABLE frame_orientation AS
SELECT * FROM (VALUES
    ('genomic', 'all'),
    ('tp73_oriented', 'all'),
    ('tp73_oriented', 'same'),
    ('tp73_oriented', 'opposite'),
    ('tp73_oriented', 'ambiguous')
) AS value(distance_frame, relative_orientation);

COPY (
WITH offset_grid AS (
    SELECT offset_twice::BIGINT AS center_offset_twice
    FROM range(
        -(SELECT max_center_offset_twice FROM capture_parameter),
        (SELECT max_center_offset_twice FROM capture_parameter) + 1
    ) offsets(offset_twice)
), frame_total AS (
    SELECT 'genomic'::VARCHAR AS distance_frame,
           count(*)::BIGINT AS anchors_total
    FROM anchor
    UNION ALL
    SELECT 'tp73_oriented', count(*)::BIGINT
    FROM anchor WHERE anchor_best_strand IS NOT NULL
), observed AS (
    SELECT
        distance_frame,
        relative_orientation,
        center_offset_twice,
        count(*)::BIGINT AS anchors_source_present,
        count(*) FILTER (WHERE best_score >= {threshold})::BIGINT
            AS anchors_positive,
        sum(physical_locus_pairs)::BIGINT AS physical_locus_pair_count,
        min(minimum_interval_distance_bp)::BIGINT
            AS minimum_interval_distance_bp,
        max(maximum_interval_distance_bp)::BIGINT
            AS maximum_interval_distance_bp
    FROM exact_feature
    GROUP BY distance_frame, relative_orientation, center_offset_twice
), cell AS (
    SELECT
        f.distance_frame,
        f.relative_orientation,
        o.center_offset_twice,
        t.anchors_total,
        coalesce(x.anchors_source_present, 0)::BIGINT
            AS anchors_source_present,
        coalesce(x.anchors_positive, 0)::BIGINT AS anchors_positive,
        coalesce(x.physical_locus_pair_count, 0)::BIGINT
            AS physical_locus_pair_count,
        x.minimum_interval_distance_bp,
        x.maximum_interval_distance_bp
    FROM frame_orientation f
    CROSS JOIN offset_grid o
    JOIN frame_total t ON t.distance_frame = f.distance_frame
    LEFT JOIN observed x
      ON x.distance_frame = f.distance_frame
     AND x.relative_orientation = f.relative_orientation
     AND x.center_offset_twice = o.center_offset_twice
)
SELECT
    {sql_string(arguments.motif_id)}::VARCHAR AS motif_id,
    {sql_string(arguments.motif_name)}::VARCHAR AS motif_name,
    {sql_string(arguments.chrom)}::VARCHAR AS chrom,
    {sql_string(arguments.anchor_motif_id)}::VARCHAR AS anchor_motif_id,
    'tp73_anchor'::VARCHAR AS transaction_scope,
    'conditional_on_tp73_anchor'::VARCHAR AS conditioning_scope,
    'physical_interval_orientation_collapsed'::VARCHAR AS locus_identity_rule,
    'motif_center_doubled_base_exact'::VARCHAR AS distance_metric,
    'included_and_labelled_by_interval_distance'::VARCHAR AS overlap_policy,
    {arguments.flank}::BIGINT AS context_flank_bp,
    {arguments.block_size}::BIGINT AS block_size_bp,
    (SELECT max_anchor_span FROM capture_parameter)::BIGINT
        AS maximum_anchor_span_bp,
    (SELECT max_neighbor_span FROM capture_parameter)::BIGINT
        AS maximum_neighbor_span_bp,
    distance_frame,
    relative_orientation,
    center_offset_twice,
    center_offset_twice / 2.0 AS center_offset_bp,
    {source_floor}::DOUBLE AS source_score_floor,
    {threshold}::DOUBLE AS positive_threshold,
    anchors_total AS transaction_count,
    anchors_source_present,
    anchors_positive AS itemset_count,
    (anchors_source_present - anchors_positive)::BIGINT AS anchors_intermediate,
    (anchors_total - anchors_source_present)::BIGINT AS anchors_negative,
    physical_locus_pair_count,
    minimum_interval_distance_bp,
    maximum_interval_distance_bp,
    anchors_source_present / nullif(anchors_total, 0)::DOUBLE AS source_support,
    anchors_positive / nullif(anchors_total, 0)::DOUBLE AS support,
    anchors_positive / nullif(anchors_total, 0)::DOUBLE
        AS confidence_anchor_to_neighbor
FROM cell
ORDER BY distance_frame, relative_orientation, center_offset_twice
) TO {sql_string(inventory_output)} (FORMAT PARQUET, COMPRESSION ZSTD);

COPY (
WITH discordant AS (
    SELECT
        a.anchor_start,
        a.anchor_end,
        a.block_index,
        a.tp73_score_stratum,
        a.anchor_best_strand,
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
), framed_discordant AS (
    SELECT d.*, 'genomic'::VARCHAR AS distance_frame
    FROM discordant d
    UNION ALL
    SELECT d.*, 'tp73_oriented'::VARCHAR
    FROM discordant d WHERE anchor_best_strand IS NOT NULL
), baseline AS (
    SELECT block_index, tp73_score_stratum, series_id, isoform, distance_frame,
           count(*) FILTER (WHERE anti_supported)::BIGINT AS total_anti,
           count(*) FILTER (WHERE control_supported)::BIGINT AS total_control
    FROM framed_discordant
    GROUP BY ALL
), observed AS (
    SELECT
        d.block_index,
        d.tp73_score_stratum,
        d.series_id,
        d.isoform,
        f.distance_frame,
        f.relative_orientation,
        f.center_offset_twice,
        count(*) FILTER (WHERE d.anti_supported)::BIGINT AS source_anti,
        count(*) FILTER (WHERE d.control_supported)::BIGINT AS source_control,
        count(*) FILTER (
            WHERE d.anti_supported AND f.best_score >= {threshold}
        )::BIGINT AS positive_anti,
        count(*) FILTER (
            WHERE d.control_supported AND f.best_score >= {threshold}
        )::BIGINT AS positive_control
    FROM framed_discordant d
    JOIN exact_feature f
      ON f.anchor_start = d.anchor_start
     AND f.anchor_end = d.anchor_end
     AND f.distance_frame = d.distance_frame
    GROUP BY ALL
), component AS (
    SELECT
        'baseline'::VARCHAR AS component_type,
        block_index,
        tp73_score_stratum,
        series_id,
        isoform,
        distance_frame,
        NULL::VARCHAR AS relative_orientation,
        NULL::BIGINT AS center_offset_twice,
        total_anti,
        total_control,
        0::BIGINT AS source_anti,
        0::BIGINT AS source_control,
        0::BIGINT AS positive_anti,
        0::BIGINT AS positive_control
    FROM baseline
    UNION ALL
    SELECT
        'observed', block_index, tp73_score_stratum, series_id, isoform,
        distance_frame, relative_orientation, center_offset_twice,
        0::BIGINT, 0::BIGINT,
        source_anti, source_control, positive_anti, positive_control
    FROM observed
)
SELECT
    {sql_string(arguments.motif_id)}::VARCHAR AS motif_id,
    {sql_string(arguments.motif_name)}::VARCHAR AS motif_name,
    {sql_string(arguments.chrom)}::VARCHAR AS chrom,
    {sql_string(arguments.anchor_motif_id)}::VARCHAR AS anchor_motif_id,
    'conditional_on_tp73_anchor'::VARCHAR AS conditioning_scope,
    {arguments.flank}::BIGINT AS context_flank_bp,
    {arguments.block_size}::BIGINT AS block_size_bp,
    component_type,
    distance_frame,
    relative_orientation,
    center_offset_twice,
    center_offset_twice / 2.0 AS center_offset_bp,
    block_index,
    tp73_score_stratum,
    series_id,
    isoform,
    {source_floor}::DOUBLE AS source_score_floor,
    {threshold}::DOUBLE AS positive_threshold,
    sum(positive_anti)::BIGINT AS positive_anti_discordant,
    sum(positive_control)::BIGINT AS positive_control_discordant,
    sum(source_anti)::BIGINT AS source_anti_discordant,
    sum(source_control)::BIGINT AS source_control_discordant,
    sum(total_anti)::BIGINT AS total_anti_discordant,
    sum(total_control)::BIGINT AS total_control_discordant,
    count(*)::BIGINT AS tp73_score_strata
FROM component
GROUP BY component_type, distance_frame, relative_orientation, center_offset_twice,
         block_index, tp73_score_stratum, series_id, isoform
ORDER BY component_type, distance_frame, relative_orientation, center_offset_twice,
         block_index, tp73_score_stratum, series_id, isoform
) TO {sql_string(block_output)} (FORMAT PARQUET, COMPRESSION ZSTD);
"""


def run_duckdb(executable: str, sql: str) -> None:
    process = subprocess.run(
        [executable, "-batch", "-bail", ":memory:"],
        input=sql,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ExactDistanceError(
            process.stderr.strip() or "DuckDB exact-distance build failed"
        )


def build(arguments: argparse.Namespace) -> None:
    inputs = (
        arguments.anchors,
        arguments.anchor_plus,
        arguments.anchor_minus,
        arguments.cofactor_plus,
        arguments.cofactor_minus,
    )
    for path in inputs:
        if not path.is_file():
            raise ExactDistanceError(f"input file is missing: {path}")
    if shutil.which(arguments.duckdb) is None:
        raise ExactDistanceError(f"DuckDB executable not found: {arguments.duckdb}")
    for output in (arguments.inventory_output, arguments.block_output):
        if output.exists():
            raise ExactDistanceError(f"refusing to replace output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)

    arguments.temp_directory.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix="tp73-exact-distance-", dir=arguments.temp_directory)
    )
    try:
        inventory = staging / "distance_inventory.parquet"
        blocks = staging / "block_components.parquet"
        sql = build_sql(arguments, inventory, blocks, staging / "duckdb-spill")
        run_duckdb(arguments.duckdb, sql)
        os.replace(inventory, arguments.inventory_output)
        os.replace(blocks, arguments.block_output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--anchors", type=Path, required=True)
    result.add_argument("--anchor-plus", type=Path, required=True)
    result.add_argument("--anchor-minus", type=Path, required=True)
    result.add_argument("--cofactor-plus", type=Path, required=True)
    result.add_argument("--cofactor-minus", type=Path, required=True)
    result.add_argument("--motif-id", required=True)
    result.add_argument("--motif-name", required=True)
    result.add_argument("--anchor-motif-id", default="MA0861.2")
    result.add_argument("--chrom", required=True)
    result.add_argument("--source-score-floor", type=float, required=True)
    result.add_argument("--positive-threshold", type=float, required=True)
    result.add_argument("--inventory-output", type=Path, required=True)
    result.add_argument("--block-output", type=Path, required=True)
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
    if arguments.flank <= 0:
        print("E: --flank must be positive", file=sys.stderr)
        return 2
    if arguments.threads <= 0:
        print("E: --threads must be positive", file=sys.stderr)
        return 2
    if arguments.positive_threshold < arguments.source_score_floor:
        print(
            "E: --positive-threshold must be at least --source-score-floor",
            file=sys.stderr,
        )
        return 2
    try:
        build(arguments)
    except (ExactDistanceError, OSError, ValueError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

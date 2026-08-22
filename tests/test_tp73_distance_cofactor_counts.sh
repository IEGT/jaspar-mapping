#!/usr/bin/env bash

set -euo pipefail

if ! command -v duckdb >/dev/null 2>&1; then
    echo "E: duckdb CLI is required for this test." >&2
    exit 1
fi

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/tp73-distance-counts.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb -light-mode -batch :memory: >/dev/null <<SQL
COPY (
  WITH base AS (
    SELECT i, (100 + 200 * i)::BIGINT AS anchor_start,
           (116 + 200 * i)::BIGINT AS anchor_end,
           CASE WHEN i < 4 THEN 'positive'
                WHEN i = 4 THEN 'intermediate' ELSE 'negative' END AS class
    FROM range(8) r(i)
  )
  SELECT '1'::VARCHAR AS chrom, anchor_start, anchor_end,
         3.0::FLOAT AS anchor_score,
         (class = 'positive' AND i <> 3 OR class = 'negative' AND i = 5)
           AS supported_tp73_saos2_TA,
         (class = 'positive' AND i = 3 OR class = 'negative' AND i <> 5)
           AS supported_negative_control_saos2_TA,
         (class = 'positive' AND i = 3 OR class = 'negative' AND i <> 7)
           AS supported_tp73_saos2_DN,
         (class = 'positive' AND i <> 3 OR class = 'negative' AND i = 7)
           AS supported_negative_control_saos2_DN,
         (class = 'positive' AND i <> 3 OR class = 'negative' AND i = 5)
           AS supported_tp73_skmel29_2_TA,
         (class = 'positive' AND i = 3 OR class = 'negative' AND i <> 5)
           AS supported_negative_control_skmel29_2_TA,
         (class = 'positive' AND i = 3 OR class = 'negative' AND i <> 7)
           AS supported_tp73_skmel29_2_DN,
         (class = 'positive' AND i <> 3 OR class = 'negative' AND i = 7)
           AS supported_negative_control_skmel29_2_DN
  FROM base
) TO '$temporary/anchors.parquet' (FORMAT PARQUET);

COPY (
  SELECT (116 + 200 * i)::BIGINT AS start,
         (125 + 200 * i)::BIGINT AS "end",
         CASE WHEN i < 4 THEN 3.0 ELSE 0.0 END::FLOAT AS score
  FROM range(5) r(i)
) TO '$temporary/plus.parquet' (FORMAT PARQUET);

COPY (
  SELECT 116::BIGINT AS start, 125::BIGINT AS "end", 2.5::FLOAT AS score
) TO '$temporary/minus.parquet' (FORMAT PARQUET);
SQL

"$repository_root/scripts/build_tp73_distance_cofactor_counts.py" \
    --anchors "$temporary/anchors.parquet" \
    --plus-hits "$temporary/plus.parquet" \
    --minus-hits "$temporary/minus.parquet" \
    --motif-id MA9999.1 --motif-name SYNTHETIC --chrom 1 \
    --positive-threshold 2 \
    --block-output "$temporary/block.parquet" \
    --class-output "$temporary/class.parquet" \
    --block-size 10000 --memory-limit 1GB --max-temp-size 1GB \
    --temp-directory "$temporary/spill"

duckdb -light-mode -batch :memory: >/dev/null <<SQL
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM read_parquet('$temporary/class.parquet')
  WHERE distance_band = 'adjacent_0_5'
    AND anchors_total = 8 AND anchors_source_present = 5
    AND anchors_positive = 4 AND anchors_intermediate = 1
    AND anchors_negative = 3
) THEN error('adjacent class counts are incorrect') END;
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM read_parquet('$temporary/class.parquet')
  WHERE distance_band <> 'adjacent_0_5'
    AND (anchors_source_present <> 0 OR anchors_negative <> 8)
) THEN error('empty distance bands are not zero-complete') END;
SELECT CASE WHEN ABS((
  SELECT SUM(mh_numerator) / SUM(mh_denominator)
  FROM read_parquet('$temporary/block.parquet')
  WHERE isoform = 'TA'
) - 6.0) > 1e-9 THEN error('TA common odds ratio is incorrect') END;
SELECT CASE WHEN ABS((
  SELECT SUM(mh_numerator) / SUM(mh_denominator)
  FROM read_parquet('$temporary/block.parquet')
  WHERE isoform = 'DN'
) - (1.0 / 6.0)) > 1e-9
THEN error('DN common odds ratio is incorrect') END;
SELECT CASE WHEN (SELECT count(DISTINCT series_id)
                  FROM read_parquet('$temporary/block.parquet')) <> 2
THEN error('the two valid cell-line series were not retained') END;
SQL

echo "TP73 distance cofactor count tests passed."

#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-cofactor-thresholds.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping cofactor-threshold test." >&2
    exit 0
}
command -v Rscript >/dev/null 2>&1 || {
    echo "I: Rscript unavailable; skipping cofactor-threshold test." >&2
    exit 0
}
Rscript -e 'library(data.table); library(ggplot2)' >/dev/null 2>&1 || {
    echo "I: Required R packages unavailable; skipping cofactor-threshold test." >&2
    exit 0
}

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT
        '1'::VARCHAR AS chrom,
        (i * 1000)::BIGINT AS anchor_start,
        (i * 1000 + 16)::BIGINT AS anchor_end,
        (((i * 7) % 13) / 2.0)::FLOAT AS anchor_score,
        (i % 2 = 0) AS supported_anti_saos2_TA,
        (i % 2 = 0) AS supported_anti_saos2_DN,
        (i % 2 = 0) AS supported_anti_skmel29_1_TA,
        (i % 2 = 0) AS supported_anti_skmel29_1_DN,
        (i % 2 = 0) AS supported_anti_skmel29_2_TA,
        (i % 2 = 0) AS supported_anti_skmel29_2_DN,
        (i % 2 = 1) AS supported_control_saos2_TA,
        (i % 2 = 1) AS supported_control_saos2_DN,
        (i % 2 = 1) AS supported_control_skmel29_1_TA,
        (i % 2 = 1) AS supported_control_skmel29_1_DN,
        (i % 2 = 1) AS supported_control_skmel29_2_TA,
        (i % 2 = 1) AS supported_control_skmel29_2_DN
    FROM range(100) AS r(i)
) TO '$temporary/anchors.parquet' (FORMAT PARQUET);

COPY (
    SELECT
        '1'::VARCHAR AS chrom,
        (i * 1000)::BIGINT AS anchor_start,
        (i * 1000 + 16)::BIGINT AS anchor_end,
        motif_id,
        CASE
            WHEN motif_id = 'M_INHIB' AND i % 10 IN (1, 2, 3, 5, 6, 7)
                THEN 5.0::FLOAT
            WHEN motif_id = 'M_NOISE' AND i % 5 IN (0, 1, 2)
                THEN 5.0::FLOAT
            WHEN motif_id = 'M_NOISE' AND i % 10 = 3
                THEN 2.0::FLOAT
            WHEN motif_id = 'M_NOISE' AND i % 10 = 8
                THEN -1.0::FLOAT
            ELSE NULL
        END AS context_score,
        -1.0::DOUBLE AS source_score_floor,
        150::BIGINT AS context_flank_bp,
        180::BIGINT AS capture_prefilter_center_bp,
        16::BIGINT AS observed_max_anchor_span_bp,
        14::BIGINT AS observed_max_context_span_bp,
        'signed_interval_edge_distance'::VARCHAR AS context_distance_metric
    FROM range(100) AS r(i)
    CROSS JOIN (VALUES ('M_INHIB'), ('M_NOISE')) AS m(motif_id)
) TO '$temporary/maxima.parquet' (FORMAT PARQUET);
SQL

"$repository_root/scripts/evaluate_tp73_cofactor_thresholds.R" \
    --anchor-evidence "$temporary/anchors.parquet" \
    --cofactor-maxima "$temporary/maxima.parquet" \
    --output-prefix "$temporary/result" \
    --thresholds 0,4 --folds 5 --chrom-size 100000 \
    --spline-df 3 --duckdb "$duckdb"

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW metrics AS
SELECT * FROM read_csv_auto('$temporary/result_threshold_metrics.tsv',
                            delim='\t', header=true);
CREATE VIEW folds AS
SELECT * FROM read_csv_auto('$temporary/result_fold_manifest.tsv',
                            delim='\t', header=true);

SELECT CASE WHEN (SELECT COUNT(*) FROM metrics) <> 4
    THEN error('threshold evaluator emitted the wrong number of rows') END;
SELECT CASE WHEN (SELECT COUNT(*) FROM folds) <> 5
    THEN error('threshold evaluator did not preserve contiguous folds') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM metrics
    WHERE motif_id = 'M_INHIB' AND threshold = 0
      AND retained_fraction > 0 AND retained_fraction < 1
      AND delta_macro_roc_auc > 0
      AND median_adjusted_odds_ratio < 1
) THEN error('known inhibitory cofactor was not recovered') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM metrics
    WHERE NOT isfinite(augmented_macro_roc_auc)
       OR NOT isfinite(median_adjusted_odds_ratio)
) THEN error('threshold evaluator emitted non-finite metrics') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1
    FROM metrics low
    JOIN metrics high USING (motif_id)
    WHERE low.motif_id = 'M_NOISE'
      AND low.threshold = 0 AND high.threshold = 4
      AND low.comparison_mode = 'fixed-negative-reference'
      AND high.comparison_mode = 'fixed-negative-reference'
      AND low.anchors_negative_reference = high.anchors_negative_reference
      AND low.anchors_negative_reference = 20
      AND low.anchors_intermediate = 10
      AND high.anchors_intermediate = 20
      AND high.negative_reference_threshold = -1
      AND NOT high.negative_reference_inclusive
      AND high.intermediate_excluded
) THEN error('fixed negative reference changed with the positive threshold') END;
SQL

[[ -s $temporary/result_threshold_sweep.png ]] || {
    echo "E: Threshold evaluator did not create its plot." >&2
    exit 1
}

"$repository_root/scripts/evaluate_tp73_cofactor_thresholds.R" \
    --anchor-evidence "$temporary/anchors.parquet" \
    --cofactor-maxima "$temporary/maxima.parquet" \
    --output-prefix "$temporary/result_zero_reference" \
    --thresholds 0,4 --folds 5 --chrom-size 100000 \
    --negative-reference-threshold 0 --spline-df 3 --compact-output \
    --duckdb "$duckdb" >/dev/null

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW zero_reference AS
SELECT * FROM read_csv_auto(
    '$temporary/result_zero_reference_threshold_metrics.tsv',
    delim='\t', header=true, nullstr='NA');
SELECT CASE WHEN NOT EXISTS (
    SELECT 1
    FROM zero_reference low
    JOIN zero_reference high USING (motif_id)
    WHERE low.motif_id = 'M_NOISE'
      AND low.threshold = 0 AND high.threshold = 4
      AND low.anchors_negative_reference = 30
      AND high.anchors_negative_reference = 30
      AND low.anchors_intermediate = 0
      AND high.anchors_intermediate = 10
      AND low.negative_reference_threshold = 0
      AND NOT low.negative_reference_inclusive
) THEN error('score-zero compatibility contrast is incorrect') END;
SQL

"$repository_root/scripts/evaluate_tp73_cofactor_thresholds.R" \
    --anchor-evidence "$temporary/anchors.parquet" \
    --cofactor-maxima "$temporary/maxima.parquet" \
    --output-prefix "$temporary/result_auto" \
    --thresholds auto --folds 5 --chrom-size 100000 \
    --spline-df 3 --duckdb "$duckdb" >/dev/null

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW auto_metrics AS
SELECT * FROM read_csv_auto('$temporary/result_auto_threshold_metrics.tsv',
                            delim='\t', header=true);
CREATE VIEW auto_config AS
SELECT * FROM read_csv_auto('$temporary/result_auto_run_config.tsv',
                            delim='\t', header=true);

SELECT CASE WHEN (SELECT MIN(threshold) FROM auto_metrics) <> 0
                  OR (SELECT MAX(threshold) FROM auto_metrics) <> 5
    THEN error('automatic threshold grid did not follow observed scores') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM auto_metrics
    WHERE samples_total <> 6 OR sample_fold_cells <> 30
) THEN error('threshold evidence counts are missing') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM auto_config
    WHERE threshold_mode = 'observed_integer_grid'
      AND threshold_specification = 'auto'
      AND comparison_mode = 'fixed-negative-reference'
      AND negative_reference_threshold = -1
      AND NOT negative_reference_inclusive
      AND intermediate_handling = 'excluded_from_positive_vs_negative_fit'
) THEN error('automatic threshold provenance is incomplete') END;
SQL

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT chrom, anchor_start, anchor_end, anchor_score,
           supported_anti_saos2_TA AS supported_tp73_saos2_TA,
           supported_anti_saos2_DN AS supported_tp73_saos2_DN,
           supported_anti_skmel29_2_TA AS supported_tp73_skmel29_2_TA,
           supported_anti_skmel29_2_DN AS supported_tp73_skmel29_2_DN,
           supported_control_saos2_TA AS supported_negative_control_saos2_TA,
           supported_control_saos2_DN AS supported_negative_control_saos2_DN,
           supported_control_skmel29_2_TA
               AS supported_negative_control_skmel29_2_TA,
           supported_control_skmel29_2_DN
               AS supported_negative_control_skmel29_2_DN
    FROM read_parquet('$temporary/anchors.parquet')
) TO '$temporary/manifest-anchors.parquet' (FORMAT PARQUET);
SQL
"$repository_root/scripts/evaluate_tp73_cofactor_thresholds.R" \
    --anchor-evidence "$temporary/manifest-anchors.parquet" \
    --cofactor-maxima "$temporary/maxima.parquet" \
    --output-prefix "$temporary/result_manifest_samples" \
    --thresholds 4 --folds 5 --chrom-size 100000 \
    --spline-df 3 --compact-output --duckdb "$duckdb" >/dev/null
"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW manifest_metrics AS SELECT * FROM read_csv_auto(
    '$temporary/result_manifest_samples_threshold_metrics.tsv',
    delim='\t', header=true, nullstr='NA');
CREATE VIEW manifest_config AS SELECT * FROM read_csv_auto(
    '$temporary/result_manifest_samples_run_config.tsv',
    delim='\t', header=true, nullstr='NA');
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM manifest_metrics
    WHERE samples_total <> 4 OR sample_fold_cells <> 20
) THEN error('manifest-selected evidence did not discover four samples') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM manifest_config
    WHERE evidence_column_scheme = 'supported_tp73_and_negative_control'
      AND sample_count = 4
      AND sample_ids NOT LIKE '%skmel29_1%'
) THEN error('manifest-selected evidence provenance is incomplete') END;
SQL

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT
        '1'::VARCHAR AS chrom,
        (i * 1000)::BIGINT AS anchor_start,
        (i * 1000 + 16)::BIGINT AS anchor_end,
        'M_NEGATIVE'::VARCHAR AS motif_id,
        -0.5::FLOAT AS context_score,
        -1.0::DOUBLE AS source_score_floor,
        150::BIGINT AS context_flank_bp,
        180::BIGINT AS capture_prefilter_center_bp,
        16::BIGINT AS observed_max_anchor_span_bp,
        14::BIGINT AS observed_max_context_span_bp,
        'signed_interval_edge_distance'::VARCHAR AS context_distance_metric
    FROM range(100) AS r(i)
) TO '$temporary/negative-maxima.parquet' (FORMAT PARQUET);
SQL

"$repository_root/scripts/evaluate_tp73_cofactor_thresholds.R" \
    --anchor-evidence "$temporary/anchors.parquet" \
    --cofactor-maxima "$temporary/negative-maxima.parquet" \
    --output-prefix "$temporary/result_negative" \
    --thresholds auto --folds 5 --chrom-size 100000 \
    --spline-df 3 --minimum-class-fraction 0.01 --compact-output \
    --duckdb "$duckdb" >/dev/null

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW negative_metrics AS
SELECT * FROM read_csv_auto('$temporary/result_negative_threshold_metrics.tsv',
                            delim='\t', header=true);
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM negative_metrics
    WHERE motif_id = 'M_NEGATIVE' AND threshold = 0
      AND anchors_retained = 0
      AND anchors_negative_reference = 0
      AND anchors_intermediate = 100
      AND evaluation_status = 'empty_comparison_class'
      AND delta_macro_roc_auc IS NULL
) THEN error('intermediate-only motif did not receive an explicit unevaluable row') END;
SQL

"$repository_root/scripts/evaluate_tp73_cofactor_thresholds.R" \
    --anchor-evidence "$temporary/anchors.parquet" \
    --cofactor-maxima "$temporary/negative-maxima.parquet" \
    --output-prefix "$temporary/result_source_floor" \
    --thresholds auto-source-floor --folds 5 --chrom-size 100000 \
    --comparison-mode threshold-complement \
    --spline-df 3 --minimum-class-fraction 0.01 --compact-output \
    --duckdb "$duckdb" >/dev/null

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW source_floor_metrics AS
SELECT * FROM read_csv_auto(
    '$temporary/result_source_floor_threshold_metrics.tsv',
    delim='\t', header=true, nullstr='NA');
CREATE VIEW source_floor_config AS
SELECT * FROM read_csv_auto(
    '$temporary/result_source_floor_run_config.tsv',
    delim='\t', header=true, nullstr='NA');

SELECT CASE WHEN (SELECT list(threshold ORDER BY threshold)
                  FROM source_floor_metrics) <> [-1.0, 0.0]
    THEN error('source-floor grid did not include every integer threshold') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM source_floor_metrics
    WHERE evaluation_status <> 'empty_comparison_class'
) THEN error('source-floor complement rows did not retain their status') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM source_floor_config
    WHERE threshold_mode = 'source_floor_integer_grid'
      AND threshold_specification = 'auto-source-floor'
      AND source_score_floor = -1
      AND comparison_mode = 'threshold-complement'
) THEN error('source-floor threshold provenance is incomplete') END;
SQL

[[ ! -e $temporary/result_negative_sample_fold_metrics.tsv &&
   ! -e $temporary/result_negative_threshold_sweep.png ]] || {
    echo "E: Compact threshold evaluation wrote detailed artifacts." >&2
    exit 1
}

echo "TP73 cofactor-threshold tests passed."

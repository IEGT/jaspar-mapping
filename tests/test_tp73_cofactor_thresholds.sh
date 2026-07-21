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
SQL

[[ -s $temporary/result_threshold_sweep.png ]] || {
    echo "E: Threshold evaluator did not create its plot." >&2
    exit 1
}

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
) THEN error('automatic threshold provenance is incomplete') END;
SQL

echo "TP73 cofactor-threshold tests passed."

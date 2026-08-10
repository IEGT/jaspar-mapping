#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-cofactor-enrichment.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping cofactor-enrichment test." >&2
    exit 0
}
command -v Rscript >/dev/null 2>&1 || {
    echo "I: Rscript unavailable; skipping cofactor-enrichment test." >&2
    exit 0
}
Rscript -e 'library(data.table)' >/dev/null 2>&1 || {
    echo "I: data.table unavailable; skipping cofactor-enrichment test." >&2
    exit 0
}

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT
        '1'::VARCHAR AS chrom,
        (i * 1000)::BIGINT AS anchor_start,
        (i * 1000 + 16)::BIGINT AS anchor_end,
        (((i * 7) % 17) / 2.0)::FLOAT AS anchor_score,
        ((i % 8 IN (0, 1) AND floor(i / 8) % 5 <> 0)
          OR (i % 8 IN (4, 5) AND floor(i / 8) % 5 = 0)
          OR i % 8 = 6) AS supported_tp73_s1,
        ((i % 8 IN (4, 5) AND floor(i / 8) % 5 <> 0)
          OR (i % 8 IN (0, 1) AND floor(i / 8) % 5 = 0)
          OR i % 8 = 2) AS supported_negative_control_s1,
        CASE WHEN ((i % 8 IN (0, 1) AND floor(i / 8) % 5 <> 0)
                        OR (i % 8 IN (4, 5) AND floor(i / 8) % 5 = 0)
                        OR i % 8 = 6)
             THEN 1 + (i % 10) ELSE 0 END::FLOAT
            AS depth_tp73_s1,
        CASE WHEN ((i % 8 IN (4, 5) AND floor(i / 8) % 5 <> 0)
                        OR (i % 8 IN (0, 1) AND floor(i / 8) % 5 = 0)
                        OR i % 8 = 2)
             THEN 1 + (i % 7) ELSE 0 END::FLOAT
            AS depth_negative_control_s1,
        ((i % 8 IN (0, 1) AND floor(i / 8) % 5 <> 1)
          OR (i % 8 IN (4, 5) AND floor(i / 8) % 5 = 1)
          OR i % 8 = 7) AS supported_tp73_s2,
        ((i % 8 IN (4, 5) AND floor(i / 8) % 5 <> 1)
          OR (i % 8 IN (0, 1) AND floor(i / 8) % 5 = 1)
          OR i % 8 = 3) AS supported_negative_control_s2,
        CASE WHEN ((i % 8 IN (0, 1) AND floor(i / 8) % 5 <> 1)
                        OR (i % 8 IN (4, 5) AND floor(i / 8) % 5 = 1)
                        OR i % 8 = 7)
             THEN 1 ELSE 0 END::FLOAT
            AS depth_tp73_s2,
        CASE WHEN ((i % 8 IN (4, 5) AND floor(i / 8) % 5 <> 1)
                        OR (i % 8 IN (0, 1) AND floor(i / 8) % 5 = 1)
                        OR i % 8 = 3)
             THEN 1 ELSE 0 END::FLOAT
            AS depth_negative_control_s2
    FROM range(400) AS r(i)
) TO '$temporary/anchors.parquet' (FORMAT PARQUET);

COPY (
    SELECT
        '1'::VARCHAR AS chrom,
        (i * 1000)::BIGINT AS anchor_start,
        (i * 1000 + 16)::BIGINT AS anchor_end,
        motif_id,
        CASE
            WHEN motif_id = 'M_ENRICHED' AND i % 8 IN (0, 1) THEN 5.0
            WHEN motif_id = 'M_ENRICHED' AND i % 8 = 2 THEN -0.5
            WHEN motif_id = 'M_ENRICHED' AND i % 8 = 3 THEN -1.0
            WHEN motif_id = 'M_DEPLETED' AND i % 8 IN (4, 5) THEN 5.0
            WHEN motif_id = 'M_DEPLETED' AND i % 8 = 6 THEN -0.5
            WHEN motif_id = 'M_DEPLETED' AND i % 8 = 7 THEN -1.0
            ELSE NULL
        END::FLOAT AS context_score,
        -1.0::DOUBLE AS source_score_floor,
        150::BIGINT AS context_flank_bp,
        'signed_interval_edge_distance'::VARCHAR AS context_distance_metric
    FROM range(400) AS r(i)
    CROSS JOIN (VALUES ('M_DEPLETED'), ('M_ENRICHED')) AS m(motif_id)
) TO '$temporary/maxima.parquet' (FORMAT PARQUET);
SQL

printf 'motif_id\tpositive_threshold\tfactor_name\tpositive_threshold_source\tselection_semantics\nM_DEPLETED\t4\tSynthetic depleted\ttest_fixture\tprespecified\nM_ENRICHED\t4\tSynthetic enriched\ttest_fixture\tprespecified\n' \
    > "$temporary/thresholds.tsv"

Rscript "$repository_root/scripts/analyze_tp73_cofactor_enrichment.R" \
    --anchor-evidence "$temporary/anchors.parquet" \
    --cofactor-maxima "$temporary/maxima.parquet" \
    --thresholds "$temporary/thresholds.tsv" \
    --output-prefix "$temporary/result" \
    --block-size 50000 --spline-df 3

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW classes AS SELECT * FROM read_csv_auto(
    '$temporary/result_class_counts.tsv', delim='\t', header=true);
CREATE VIEW descriptive AS SELECT * FROM read_csv_auto(
    '$temporary/result_descriptive.tsv', delim='\t', header=true);
CREATE VIEW macro AS SELECT * FROM read_csv_auto(
    '$temporary/result_macro_summary.tsv', delim='\t', header=true);
CREATE VIEW primary_test AS SELECT * FROM read_csv_auto(
    '$temporary/result_primary_occupancy.tsv', delim='\t', header=true,
    nullstr='NA');
CREATE VIEW depth_manifest AS SELECT * FROM read_csv_auto(
    '$temporary/result_depth_tier_manifest.tsv', delim='\t', header=true);
CREATE VIEW run_config AS SELECT * FROM read_csv_auto(
    '$temporary/result_run_config.tsv', delim='\t', header=true);

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM classes
    WHERE motif_id = 'M_ENRICHED' AND tp73_score_stratum = 'all'
      AND negative_reference_threshold = -1
      AND anchors_positive = 100 AND anchors_negative_reference = 200
      AND anchors_intermediate = 100
) THEN error('strict <-1 class partition is wrong') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM classes
    WHERE motif_id = 'M_ENRICHED' AND tp73_score_stratum = 'all'
      AND negative_reference_threshold = 0
      AND anchors_positive = 100 AND anchors_negative_reference = 300
      AND anchors_intermediate = 0
) THEN error('cumulative <0 class partition is wrong') END;
SELECT CASE WHEN (SELECT count(*) FROM depth_manifest) <> 24
    THEN error('per-track depth tier manifest is incomplete') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM depth_manifest
    WHERE depth_tier <> 'strict_immersion'
      AND (positive_depth_quantile IS NULL OR effective_depth_cutoff <= 0)
) THEN error('positive-depth quantile tier is invalid') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM depth_manifest
    WHERE sample_id = 's2' AND depth_tier = 'positive_depth_q50'
      AND duplicates_previous_tier
) THEN error('duplicate depth event tiers were not identified') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM macro
    WHERE motif_id = 'M_ENRICHED' AND tp73_score_stratum = 'all'
      AND negative_reference_threshold = -1
      AND depth_tier = 'strict_immersion'
      AND mean_log2_anti_control_specificity_ratio_jeffreys > 0
) THEN error('known anti-p73 enrichment was not recovered') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM macro
    WHERE motif_id = 'M_DEPLETED' AND tp73_score_stratum = 'all'
      AND negative_reference_threshold = -1
      AND depth_tier = 'strict_immersion'
      AND mean_log2_anti_control_specificity_ratio_jeffreys < 0
) THEN error('known anti-p73 depletion was not recovered') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM primary_test
    WHERE motif_id = 'M_ENRICHED' AND evaluation_status = 'ok'
      AND factor_name = 'Synthetic enriched'
      AND comparison_label = 'score>=4 versus score<-1'
      AND negative_reference_threshold = -1
      AND adjusted_odds_ratio > 1 AND p_value IS NOT NULL
) THEN error('primary matched test missed enrichment') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM primary_test
    WHERE motif_id = 'M_DEPLETED' AND evaluation_status = 'ok'
      AND adjusted_odds_ratio < 1 AND p_value IS NOT NULL
) THEN error('primary matched test missed depletion') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM run_config
    WHERE negative_reference_semantics = 'strict context_score < N or absent'
      AND primary_negative_reference = -1
      AND context_geometry = 'signed_interval_edge_distance'
      AND evidence_column_scheme = 'supported_tp73_and_negative_control'
      AND sample_count = 2
) THEN error('run provenance omits comparison semantics') END;
SQL

echo "I: TP73 cofactor enrichment synthetic test passed."

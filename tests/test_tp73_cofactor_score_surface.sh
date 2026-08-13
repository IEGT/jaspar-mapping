#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-score-surface.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping TP73/cofactor score-surface test." >&2
    exit 0
}
command -v Rscript >/dev/null 2>&1 || {
    echo "I: Rscript unavailable; skipping TP73/cofactor score-surface test." >&2
    exit 0
}
Rscript -e 'library(data.table)' >/dev/null 2>&1 || {
    echo "I: data.table unavailable; skipping TP73/cofactor score-surface test." >&2
    exit 0
}

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    WITH base AS (
        SELECT
            i,
            CASE WHEN i < 4000 THEN '1' ELSE '2' END::VARCHAR AS chrom,
            (i * 10000)::BIGINT AS anchor_start,
            (-1.0 + 16.0 * (((i * 17) % 1000) / 999.0))::FLOAT
                AS anchor_score,
            CASE WHEN i % 5 = 0 THEN NULL
                 ELSE -1.0 + 11.0 * (((i * 37) % 1000) / 999.0)
            END::DOUBLE AS context_score,
            ((hash(i, 's1') % 100000)::DOUBLE / 100000.0) AS u1,
            ((hash(i, 's2') % 100000)::DOUBLE / 100000.0) AS u2,
            ((hash(i, 's3') % 100000)::DOUBLE / 100000.0) AS u3
        FROM range(8000) AS r(i)
    ), probability AS (
        SELECT *,
            1.0 / (1.0 + exp(-(
                -0.5 + 0.02 * anchor_score +
                0.06 * coalesce(context_score + 1.0, 0.0) +
                0.01 * anchor_score * coalesce(context_score + 1.0, 0.0)
            ))) AS p1,
            1.0 / (1.0 + exp(-(
                -0.35 + 0.02 * anchor_score +
                0.06 * coalesce(context_score + 1.0, 0.0) +
                0.01 * anchor_score * coalesce(context_score + 1.0, 0.0)
            ))) AS p2,
            1.0 / (1.0 + exp(-(
                -0.2 + 0.02 * anchor_score +
                0.06 * coalesce(context_score + 1.0, 0.0) +
                0.01 * anchor_score * coalesce(context_score + 1.0, 0.0)
            ))) AS p3
        FROM base
    )
    SELECT
        chrom,
        anchor_start,
        anchor_start + 16 AS anchor_end,
        anchor_score,
        u1 < p1 AS supported_tp73_s1,
        NOT (u1 < p1) AS supported_negative_control_s1,
        u2 < p2 AS supported_tp73_s2,
        NOT (u2 < p2) AS supported_negative_control_s2,
        u3 < p3 AS supported_tp73_s3,
        NOT (u3 < p3) AS supported_negative_control_s3
    FROM probability
) TO '$temporary/anchors.parquet' (FORMAT PARQUET);

COPY (
    WITH base AS (
        SELECT
            i,
            CASE WHEN i < 4000 THEN '1' ELSE '2' END::VARCHAR AS chrom,
            (i * 10000)::BIGINT AS anchor_start,
            CASE WHEN i % 5 = 0 THEN NULL
                 ELSE -1.0 + 11.0 * (((i * 37) % 1000) / 999.0)
            END::FLOAT AS context_score
        FROM range(8000) AS r(i)
    )
    SELECT
        chrom,
        anchor_start,
        anchor_start + 16 AS anchor_end,
        motif_id,
        context_score,
        CASE WHEN motif_id = 'M_CENSORED' THEN 2.0 ELSE -2.0 END::DOUBLE
            AS source_score_floor,
        150::BIGINT AS context_flank_bp,
        'signed_interval_edge_distance'::VARCHAR AS context_distance_metric
    FROM base
    CROSS JOIN (VALUES ('M_CENSORED'), ('M_RESPONSE')) AS m(motif_id)
) TO '$temporary/maxima.parquet' (FORMAT PARQUET);
SQL

printf 'motif_id\tpositive_threshold\tfactor_name\nM_CENSORED\t4\tSynthetic censored\nM_RESPONSE\t4\tSynthetic response\n' \
    > "$temporary/thresholds.tsv"
printf 'motif_id\tpositive_threshold\tfactor_name\nM_CENSORED\t4\tSynthetic censored\n' \
    > "$temporary/censored-thresholds.tsv"

Rscript "$repository_root/scripts/analyze_tp73_cofactor_score_surface.R" \
    --anchor-evidence "$temporary/anchors.parquet" \
    --cofactor-maxima "$temporary/maxima.parquet" \
    --thresholds "$temporary/thresholds.tsv" \
    --output-prefix "$temporary/result" \
    --exclude-samples s3 \
    --source-dirty true \
    --block-size 5000000 --spline-df 3 --minimum-band-count 50 \
    --inference-status synthetic_known_tp73_by_cofactor_score_response

Rscript "$repository_root/scripts/analyze_tp73_cofactor_score_surface.R" \
    --anchor-evidence "$temporary/anchors.parquet" \
    --cofactor-maxima "$temporary/maxima.parquet" \
    --thresholds "$temporary/censored-thresholds.tsv" \
    --output-prefix "$temporary/censored" \
    --exclude-samples s3 --block-size 5000000 --spline-df 3 \
    --inference-status synthetic_all_censored_contract
for suffix in motif_status score_band score_surface tp73_score_contrast run_config; do
    head -n 1 "$temporary/censored_${suffix}.tsv" | grep -q '^motif_id\|^anchor_evidence' || {
        echo "E: all-censored output lacks a header: $suffix" >&2
        exit 1
    }
done

if [[ ${VERBOSE_TEST_OUTPUT:-0} == 1 ]]; then
    "$duckdb" -box -c "SELECT * FROM read_csv_auto(
        '$temporary/result_motif_status.tsv', delim='\\t', header=true,
        nullstr='NA');"
    "$duckdb" -box -c "SELECT * FROM read_csv_auto(
        '$temporary/result_tp73_score_contrast.tsv', delim='\\t', header=true,
        nullstr='NA') WHERE motif_id='M_RESPONSE';"
fi

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW status AS SELECT * FROM read_csv_auto(
    '$temporary/result_motif_status.tsv', delim='\t', header=true, nullstr='NA');
CREATE VIEW band AS SELECT * FROM read_csv_auto(
    '$temporary/result_score_band.tsv', delim='\t', header=true, nullstr='NA');
CREATE VIEW surface AS SELECT * FROM read_csv_auto(
    '$temporary/result_score_surface.tsv', delim='\t', header=true, nullstr='NA');
CREATE VIEW tp73_contrast AS SELECT * FROM read_csv_auto(
    '$temporary/result_tp73_score_contrast.tsv', delim='\t', header=true,
    nullstr='NA');
CREATE VIEW config AS SELECT * FROM read_csv_auto(
    '$temporary/result_run_config.tsv', delim='\t', header=true, nullstr='NA');

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM status
    WHERE motif_id = 'M_CENSORED'
      AND evaluation_status = 'negative_reference_below_source_floor'
      AND NOT negative_reference_observable
      AND anchors_negative_reference = 0
) THEN error('censored source floor was treated as an observed negative') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM status
    WHERE motif_id = 'M_RESPONSE' AND evaluation_status = 'ok'
      AND negative_reference_threshold = -1
      AND negative_reference_observable
      AND try_cast(tp73_by_cofactor_interaction_p_value AS DOUBLE) < 0.05
) THEN error('known TP73-by-cofactor score interaction was not recovered') END;

SELECT CASE WHEN (SELECT count(*) FROM band WHERE motif_id = 'M_RESPONSE') < 4
    THEN error('empirical cofactor score bands were not emitted') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM band
    WHERE motif_id = 'M_RESPONSE' AND score_band <> 'negative'
      AND anchors < 50
) THEN error('an emitted positive score band violates minimum support') END;

SELECT CASE WHEN NOT EXISTS (
    WITH highest AS (
        SELECT * FROM surface
        WHERE motif_id = 'M_RESPONSE'
          AND score_band_order = (
              SELECT max(score_band_order) FROM surface
              WHERE motif_id = 'M_RESPONSE'
          )
    ), ends AS (
        SELECT
            arg_min(
                try_cast(adjusted_odds_ratio_vs_negative AS DOUBLE),
                try_cast(tp73_score AS DOUBLE)
            ) AS low_score_or,
            arg_max(
                try_cast(adjusted_odds_ratio_vs_negative AS DOUBLE),
                try_cast(tp73_score AS DOUBLE)
            ) AS high_score_or
        FROM highest
    )
    SELECT 1 FROM ends
    WHERE high_score_or > low_score_or AND high_score_or > 1
) THEN error('the fitted high cofactor band does not strengthen with TP73 score') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM surface
    WHERE motif_id = 'M_RESPONSE'
      AND try_cast(adjusted_odds_ratio_vs_negative AS DOUBLE) > 1
      AND try_cast(adjusted_odds_ratio_vs_previous_band AS DOUBLE) IS NOT NULL
      AND try_cast(confidence_interval_95_lower_vs_negative AS DOUBLE) IS NOT NULL
) THEN error('score surface lacks reference and adjacent-band contrasts') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM surface
    WHERE motif_id = 'M_RESPONSE'
      AND score_band_order = (
          SELECT max(score_band_order) FROM surface
          WHERE motif_id = 'M_RESPONSE'
      )
      AND try_cast(adjusted_odds_ratio_vs_weakest_positive AS DOUBLE) > 1
      AND try_cast(p_value_vs_weakest_positive AS DOUBLE) IS NOT NULL
) THEN error('strongest-versus-weakest cofactor contrast is absent') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73_contrast
    WHERE motif_id = 'M_RESPONSE'
      AND score_band_order = (
          SELECT max(score_band_order) FROM tp73_contrast
          WHERE motif_id = 'M_RESPONSE'
      )
      AND try_cast(adjusted_ratio_of_odds_ratios AS DOUBLE) > 1
      AND try_cast(p_value AS DOUBLE) IS NOT NULL
) THEN error('high-versus-low TP73 score contrast missed the known increase') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM config
    WHERE negative_reference_semantics = 'strict context_score < N or absent'
      AND score_band_semantics LIKE '%not cumulative thresholds%'
      AND model_formula LIKE '%ns(anchor_score%* score_band%'
      AND minimum_class_fraction = 0.01
      AND sample_ids = 's1,s2' AND excluded_sample_ids = 's3'
      AND source_dirty
      AND context_geometry = 'signed_interval_edge_distance'
      AND inference_status = 'synthetic_known_tp73_by_cofactor_score_response'
) THEN error('score-surface provenance is incomplete') END;
SQL

echo "I: TP73/cofactor score-surface synthetic test passed."

#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-motif-thresholds.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping motif-score-threshold test." >&2
    exit 0
}

cat > "$temporary/jaspar.txt" <<'EOF'
>MA0001.1	INHIB
A [ 1 ]
C [ 1 ]
G [ 1 ]
T [ 1 ]
>MA0002.1	NO_GAIN
A [ 1 ]
C [ 1 ]
G [ 1 ]
T [ 1 ]
>MA0003.1	PENDING
A [ 1 ]
C [ 1 ]
G [ 1 ]
T [ 1 ]
>MA0004.1	TOO_RARE
A [ 1 ]
C [ 1 ]
G [ 1 ]
T [ 1 ]
>MA0861.2	TP73
A [ 1 ]
C [ 1 ]
G [ 1 ]
T [ 1 ]
EOF

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT * FROM (VALUES
        ('MA0001.1', 0.0, 1000, 900, 0.90, 2000, 0.50, 0.51, 0.01,
         0.60, 0.61, 0.01, 0.70, 0.69, -0.01, 0.20, 0.19, -0.01,
         0.80, 0.75, 0.85, 0.70, 6, 6, 30, 25, 24),
        ('MA0001.1', 1.0, 1000, 600, 0.60, 2000, 0.50, 0.55, 0.05,
         0.60, 0.64, 0.04, 0.70, 0.66, -0.04, 0.20, 0.18, -0.02,
         0.70, 0.65, 0.75, 0.60, 6, 6, 30, 29, 26),
        ('MA0001.1', 2.0, 1000, 400, 0.40, 2000, 0.50, 0.548, 0.048,
         0.60, 0.638, 0.038, 0.70, 0.662, -0.038, 0.20, 0.181, -0.019,
         0.72, 0.67, 0.77, 0.62, 6, 6, 30, 28, 25),
        ('MA0001.1', 3.0, 1000, 200, 0.20, 2000, 0.50, 0.52, 0.02,
         0.60, 0.615, 0.015, 0.70, 0.685, -0.015, 0.20, 0.195, -0.005,
         0.78, 0.72, 0.84, 0.68, 6, 6, 30, 20, 18),
        ('MA0002.1', 0.0, 1000, 500, 0.50, 2000, 0.50, 0.49, -0.01,
         0.60, 0.59, -0.01, 0.70, 0.71, 0.01, 0.20, 0.21, 0.01,
         1.10, 1.05, 1.15, 1.08, 0, 6, 30, 5, 4),
        ('MA0002.1', 1.0, 1000, 300, 0.30, 2000, 0.50, 0.495, -0.005,
         0.60, 0.595, -0.005, 0.70, 0.705, 0.005, 0.20, 0.205, 0.005,
         1.05, 1.01, 1.10, 1.04, 0, 6, 30, 8, 7),
        ('MA0004.1', 8.0, 10000, 5, 0.0005, 2000, 0.50, 0.60, 0.10,
         0.60, 0.70, 0.10, 0.70, 0.60, -0.10, 0.20, 0.15, -0.05,
         0.20, 0.10, 0.30, 0.15, 6, 6, 30, 30, 30)
    ) AS t(
        motif_id, threshold, anchors_total, anchors_retained,
        retained_fraction, discordant_observations,
        baseline_macro_roc_auc, augmented_macro_roc_auc, delta_macro_roc_auc,
        baseline_macro_average_precision, augmented_macro_average_precision,
        delta_macro_average_precision, baseline_macro_log_loss,
        augmented_macro_log_loss, delta_macro_log_loss,
        baseline_macro_brier_score, augmented_macro_brier_score,
        delta_macro_brier_score, median_adjusted_odds_ratio,
        minimum_adjusted_odds_ratio, maximum_adjusted_odds_ratio,
        median_raw_sample_odds_ratio, samples_with_raw_odds_ratio_below_one,
        samples_total, sample_fold_cells,
        sample_fold_cells_with_roc_auc_gain,
        sample_fold_cells_with_log_loss_gain
    )
) TO '$temporary/metrics.parquet' (FORMAT PARQUET);
SQL

"$repository_root/scripts/build_motif_score_thresholds.py" \
    --metrics "$temporary/metrics.parquet" \
    --jaspar "$temporary/jaspar.txt" \
    --motif MA0001.1 --motif MA0002.1 --motif MA0003.1 --motif MA0004.1 \
    --output "$temporary/thresholds.parquet" \
    --threshold-set-id synthetic_v1 \
    --calibration-run-id synthetic_run \
    --genome-id genome1 --motif-set-id motifs1 \
    --calibration-scope chromosome_1 \
    --evidence-dataset-id synthetic_cutrun \
    --source-commit 0000000000000000000000000000000000000000 \
    --source-dirty --duckdb "$duckdb"

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW t AS SELECT * FROM read_parquet('$temporary/thresholds.parquet');

SELECT CASE WHEN (SELECT COUNT(*) FROM t) <> 4
    THEN error('threshold registry omitted an expected motif') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM t
    WHERE motif_id = 'MA0001.1'
      AND recommended_threshold = 1
      AND useful_threshold_min = 1
      AND useful_threshold_max = 2
      AND selected_metric_gain = 0.05
      AND association_direction = 'inhibitory_association'
      AND calibration_status = 'exploratory_positive_gain'
) THEN error('positive-gain motif threshold was selected incorrectly') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM t
    WHERE motif_id = 'MA0002.1'
      AND recommended_threshold IS NULL
      AND calibration_status = 'no_positive_gain'
) THEN error('no-gain motif did not retain explicit null semantics') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM t
    WHERE motif_id = 'MA0003.1'
      AND recommended_threshold IS NULL
      AND calibration_status = 'pending'
) THEN error('missing calibration did not become pending') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM t
    WHERE motif_id = 'MA0004.1'
      AND recommended_threshold IS NULL
      AND eligible_threshold_count = 0
      AND calibration_status = 'insufficient_class_support'
) THEN error('rare candidate bypassed the minimum class-support rule') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM t
    WHERE threshold_inclusive IS DISTINCT FROM TRUE
       OR schema_version <> 1
       OR calibration_stratum_id <> 'all_tp73_anchors'
       OR context_max_interval_distance_bp <> 150
       OR source_metrics_sha256 IS NULL
       OR jaspar_sha256 IS NULL
) THEN error('threshold provenance is incomplete') END;
SQL

mkdir -p "$temporary/package/tables/jaspar2026/motif_score_threshold"
cp "$temporary/thresholds.parquet" \
    "$temporary/package/tables/jaspar2026/motif_score_threshold/part-000000.parquet"
awk '
    /^CREATE OR REPLACE VIEW motif_score_threshold AS/ { capture = 1 }
    /^CREATE OR REPLACE VIEW scan_task AS/ { capture = 0 }
    capture { print }
' "$repository_root/sql/schema.sql" > "$temporary/threshold-schema.sql"
(
    cd "$temporary/package"
    "$duckdb" -batch threshold.duckdb < "$temporary/threshold-schema.sql"
    "$duckdb" -batch threshold.duckdb >/dev/null <<SQL
SELECT CASE WHEN (SELECT COUNT(*) FROM motif_score_threshold) <> 4
    THEN error('schema threshold view omitted rows') END;
SELECT CASE WHEN (SELECT COUNT(*) FROM motif_convenient_threshold) <> 1
    THEN error('convenient-threshold view did not enforce status/null semantics') END;
SQL

    {
        cat <<'SQL'
CREATE TABLE motif_context_pair (
    anchor_hit_id VARCHAR,
    genome_id VARCHAR,
    motif_set_id VARCHAR,
    chrom VARCHAR,
    anchor_start BIGINT,
    anchor_end BIGINT,
    anchor_strand VARCHAR,
    neighbor_hit_id VARCHAR,
    neighbor_start BIGINT,
    neighbor_end BIGINT,
    neighbor_motif_id VARCHAR,
    neighbor_motif_name VARCHAR,
    neighbor_strand VARCHAR,
    neighbor_score DOUBLE,
    anchor_neighbor_interval_distance_bp BIGINT,
    interval_relation VARCHAR,
    interval_distance_band VARCHAR,
    relative_orientation VARCHAR,
    anchor_oriented_side VARCHAR,
    score_mode VARCHAR,
    pseudocount DOUBLE,
    background_model_id VARCHAR,
    pseudocount_scheme VARCHAR
);
INSERT INTO motif_context_pair VALUES
    ('A1', 'genome1', 'motifs1', '1', 100, 116, '+',
     'N1', 120, 130, 'MA0001.1', 'INHIB', '-', 2.0,
     4, 'disjoint', 'adjacent_0_5', 'opposite', 'downstream',
     'log2_relative_risk', 1.0, 'uniform_acgt_v1', 'additive_per_base'),
    ('A2', 'genome1', 'motifs1', '1', 1000, 1016, '+',
     'N2', 1167, 1177, 'MA0001.1', 'INHIB', '-', 2.0,
     151, 'disjoint', 'out_of_range', 'opposite', 'downstream',
     'log2_relative_risk', 1.0, 'uniform_acgt_v1', 'additive_per_base'),
    ('A3', 'genome1', 'motifs1', '1', 2000, 2016, '+',
     'N3', 2020, 2030, 'MA0001.1', 'INHIB', '-', 0.0,
     4, 'disjoint', 'adjacent_0_5', 'opposite', 'downstream',
     'log2_relative_risk', 1.0, 'uniform_acgt_v1', 'additive_per_base');

PREPARE q16 AS
SQL
        awk '
            /^-- Q16\./ { capture = 1 }
            /^-- Q17\./ { capture = 0 }
            capture { print }
        ' "$repository_root/sql/queries.sql"
        cat <<'SQL'
EXECUTE q16(
    threshold_set_id := 'synthetic_v1', genome_id := 'genome1',
    motif_set_id := 'motifs1', score_mode := 'log2_relative_risk',
    pseudocount := 1.0, threshold_role := 'tp73_context_binary_feature',
    target_motif_id := 'MA0861.2',
    calibration_stratum_id := 'all_tp73_anchors'
);

PREPARE q17 AS
SQL
        awk '
            /^-- Q17\./ { capture = 1 }
            /^-- Q18\./ { capture = 0 }
            capture { print }
        ' "$repository_root/sql/queries.sql"
        cat <<'SQL'
EXECUTE q17(
    threshold_set_id := 'synthetic_v1', genome_id := 'genome1',
    motif_set_id := 'motifs1', chrom := '1',
    threshold_role := 'tp73_context_binary_feature',
    target_motif_id := 'MA0861.2', score_mode := 'log2_relative_risk',
    pseudocount := 1.0, calibration_stratum_id := 'all_tp73_anchors'
);
SQL
    } | "$duckdb" -csv threshold.duckdb > "$temporary/query-output.csv"
)

grep -Fq 'N1' "$temporary/query-output.csv" || {
    echo "E: In-range motif above its convenient threshold was omitted." >&2
    exit 1
}
if grep -Eq 'N2|N3' "$temporary/query-output.csv"; then
    echo "E: Convenient-threshold query ignored distance or score bounds." >&2
    exit 1
fi

if "$repository_root/scripts/build_motif_score_thresholds.py" \
    --metrics "$temporary/metrics.parquet" --jaspar "$temporary/jaspar.txt" \
    --motif MA0001.1 --output "$temporary/thresholds.parquet" \
    --threshold-set-id reuse --calibration-run-id reuse \
    --genome-id genome1 --motif-set-id motifs1 \
    --calibration-scope chromosome_1 --evidence-dataset-id synthetic \
    --source-commit 0000000000000000000000000000000000000000 \
    >"$temporary/reuse.out" 2>"$temporary/reuse.err"
then
    echo "E: Existing motif threshold output was overwritten." >&2
    exit 1
fi
grep -Fq 'output already exists' "$temporary/reuse.err" || {
    echo "E: Existing-output failure was not explained." >&2
    exit 1
}

echo "Motif score-threshold tests passed."

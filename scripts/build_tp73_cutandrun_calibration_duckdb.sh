#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: build_tp73_cutandrun_calibration_duckdb.sh --run-root DIR
       [--output FILE.duckdb]

Materialize the compact anti-p73, IgG, and threshold-recommendation TSV files
from a TP73 chromosome calibration into one read-only-friendly DuckDB file.
The input directory must contain risk_p0, risk_p1, log_odds_p1 and matching
igg_* result directories. The default output is DIR/tp73_calibration.duckdb.
An existing output is never replaced.

Options:
  --run-root DIR   Calibration result root
  --output FILE    New DuckDB path
  -h, --help       Show this help
EOF
}

run_root=
output=
while (($#)); do
    case "$1" in
        --run-root)
            [[ $# -ge 2 ]] || { echo "E: --run-root requires a value" >&2; exit 2; }
            run_root=$2
            shift 2
            ;;
        --output)
            [[ $# -ge 2 ]] || { echo "E: --output requires a value" >&2; exit 2; }
            output=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "E: unknown argument: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

[[ -n "$run_root" ]] || { echo "E: --run-root is required" >&2; exit 2; }
[[ -d "$run_root" ]] || { echo "E: run root not found: $run_root" >&2; exit 1; }
command -v duckdb >/dev/null 2>&1 || {
    echo "E: DuckDB CLI is required." >&2
    exit 1
}

run_root=$(cd "$run_root" && pwd)
if [[ -z "$output" ]]; then
    output="$run_root/tp73_calibration.duckdb"
fi
[[ ! -e "$output" ]] || {
    echo "E: refusing to replace existing output: $output" >&2
    exit 1
}

for run in risk_p0 risk_p1 log_odds_p1 \
           igg_risk_p0 igg_risk_p1 igg_log_odds_p1; do
    [[ -f "$run_root/$run/threshold_curve.tsv" ]] || {
        echo "E: missing calibration run: $run_root/$run" >&2
        exit 1
    }
done
for suffix in aggregate recommendations; do
    summary_path="$run_root/tp73_threshold_$suffix.tsv"
    [[ -f "$summary_path" ]] || {
        echo "E: missing threshold summary: $summary_path" >&2
        exit 1
    }
done

sql_root=$(printf '%s' "$run_root" | sed "s/'/''/g")

duckdb "$output" -bail -c "
CREATE TABLE threshold_curve AS
WITH source AS (
    SELECT regexp_extract(filename, '/([^/]+)/threshold_curve.tsv', 1) AS run_name,
           * EXCLUDE(filename)
    FROM read_csv('$sql_root/*/threshold_curve.tsv',
                  delim='\t', header=true, filename=true)
)
SELECT CASE WHEN starts_with(run_name, 'igg_') THEN 'igg'
            ELSE 'anti_p73' END AS evidence_type,
       replace(run_name, 'igg_', '') AS run_name,
       * EXCLUDE(run_name)
FROM source;

CREATE TABLE score_histogram AS
WITH source AS (
    SELECT regexp_extract(filename, '/([^/]+)/score_histogram.tsv', 1) AS run_name,
           * EXCLUDE(filename)
    FROM read_csv('$sql_root/*/score_histogram.tsv',
                  delim='\t', header=true, filename=true)
)
SELECT CASE WHEN starts_with(run_name, 'igg_') THEN 'igg'
            ELSE 'anti_p73' END AS evidence_type,
       replace(run_name, 'igg_', '') AS run_name,
       * EXCLUDE(run_name)
FROM source;

CREATE TABLE calibration_summary AS
WITH source AS (
    SELECT regexp_extract(filename, '/([^/]+)/summary.tsv', 1) AS run_name,
           * EXCLUDE(filename)
    FROM read_csv('$sql_root/*/summary.tsv',
                  delim='\t', header=true, filename=true)
)
SELECT CASE WHEN starts_with(run_name, 'igg_') THEN 'igg'
            ELSE 'anti_p73' END AS evidence_type,
       replace(run_name, 'igg_', '') AS run_name,
       * EXCLUDE(run_name)
FROM source;

CREATE TABLE run_config AS
WITH source AS (
    SELECT regexp_extract(filename, '/([^/]+)/run_config.tsv', 1) AS run_name,
           * EXCLUDE(filename)
    FROM read_csv('$sql_root/*/run_config.tsv',
                  delim='\t', header=true, filename=true)
)
SELECT CASE WHEN starts_with(run_name, 'igg_') THEN 'igg'
            ELSE 'anti_p73' END AS evidence_type,
       replace(run_name, 'igg_', '') AS run_name,
       * EXCLUDE(run_name)
FROM source;

CREATE TABLE threshold_aggregate AS
SELECT * FROM read_csv('$sql_root/tp73_threshold_aggregate.tsv',
                       delim='\t', header=true);

CREATE TABLE threshold_recommendation AS
SELECT * FROM read_csv('$sql_root/tp73_threshold_recommendations.tsv',
                       delim='\t', header=true);

CREATE INDEX threshold_curve_lookup
    ON threshold_curve(evidence_type, run_name, sample_id, threshold);
CREATE INDEX score_histogram_lookup
    ON score_histogram(evidence_type, run_name, sample_id, threshold);
CREATE INDEX summary_lookup
    ON calibration_summary(evidence_type, run_name, sample_id);
CHECKPOINT;"

duckdb -readonly "$output" -c "
SELECT recommendation_kind, mode, threshold, selected_motifs,
       round(minimum_support_ratio, 3) AS min_support_ratio,
       round(median_support_ratio, 3) AS median_support_ratio,
       round(minimum_depth_ratio, 3) AS min_depth_ratio,
       round(median_depth_ratio, 3) AS median_depth_ratio
FROM threshold_recommendation
ORDER BY recommendation_kind, mode;"

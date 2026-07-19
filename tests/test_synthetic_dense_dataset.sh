#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
fixture_dir="$repository_root/test_files/synthetic_dense"
scanner="$repository_root/pssm_scan_parquet"
temporary=""

if [[ $# -gt 1 ]]; then
    echo "Usage: test_synthetic_dense_dataset.sh [OUTPUT_DIR]" >&2
    exit 2
fi

if [[ $# -eq 1 ]]; then
    output=$1
    [[ $output = /* ]] || output="$repository_root/$output"
    mkdir -p "$output"
else
    temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-synthetic-dense.XXXXXX")
    output=$temporary
    trap 'rm -rf "$temporary"' EXIT HUP INT TERM
fi
output=$(cd "$output" && pwd -P)

command -v duckdb >/dev/null 2>&1 || {
    echo "E: duckdb CLI is required." >&2
    exit 1
}
[[ -x $scanner ]] || {
    echo "E: Build the Arrow scanner with 'make pssm_scan_parquet'." >&2
    exit 1
}

genome="$fixture_dir/genome.fna"
pssm="$fixture_dir/JASPAR2026_synthetic.jaspar"
expected="$fixture_dir/expected_scores.tsv"
fasta_index="$output/genome.fna.fai"

python3 "$repository_root/scripts/build_fasta_index.py" \
    "$genome" --output "$fasta_index" >/dev/null 2>&1

for mode in log2_relative_risk log_odds; do
    "$scanner" \
        --dense-scores \
        --dense-block-size 3 \
        --outdir "$output" \
        --genome "$genome" \
        --fasta-index "$fasta_index" \
        --pssm "$pssm" \
        --motif MA9999.1 \
        --motif-set-id synthetic_jaspar2026 \
        --genome-id synthetic_acgtn_v1 \
        --chr 1 \
        --strand both \
        --coordinate-mode bed \
        --score-mode "$mode" \
        --pseudocount 0 \
        --skip-N >/dev/null
done

metadata_dir="$output/tables/jaspar2026/motif_metadata"
metadata_file="$metadata_dir/part-000000.parquet"
mkdir -p "$metadata_dir"
if [[ ! -f $metadata_file ]]; then
    (
        cd "$output"
        duckdb :memory: -bail -c "COPY (
            SELECT 'MA9999.1'::VARCHAR AS motif_id,
                   'synthetic_jaspar2026'::VARCHAR AS motif_set_id,
                   'A_RICH'::VARCHAR AS motif_name,
                   1::INTEGER AS motif_length,
                   2026::INTEGER AS jaspar_version,
                   'synthetic-fixture'::VARCHAR AS source_sha256
        ) TO 'tables/jaspar2026/motif_metadata/part-000000.parquet'
          (FORMAT PARQUET, COMPRESSION ZSTD);" >/dev/null
    )
fi

expected_sql=${expected//\'/\'\'}
if source_commit=$(git -C "$repository_root" rev-parse HEAD 2>/dev/null); then
    :
elif [[ -s $repository_root/SOURCE_COMMIT ]]; then
    source_commit=$(<"$repository_root/SOURCE_COMMIT")
else
    echo "E: Cannot determine the source commit from Git or SOURCE_COMMIT." >&2
    exit 1
fi
database=synthetic_dense.duckdb

(
    cd "$output"
    duckdb "$database" -bail -f \
        "$repository_root/sql/chr1_dense_dry_run_schema.sql" >/dev/null
    duckdb "$database" -bail -c "
        CREATE OR REPLACE TABLE run_manifest AS
        SELECT 'synthetic_ACGTN'::VARCHAR AS run_id,
               '$source_commit'::VARCHAR AS source_commit,
               'ACGTN'::VARCHAR AS reference_sequence,
               0.0::DOUBLE AS pseudocount;

        CREATE OR REPLACE TABLE synthetic_expected AS
        SELECT * FROM read_csv(
            '$expected_sql', delim='\\t', header=true, nullstr='NULL'
        );

        CREATE OR REPLACE VIEW synthetic_actual AS
        SELECT
            start,
            MIN(\"end\") AS \"end\",
            MAX(score) FILTER (
                WHERE score_mode = 'log2_relative_risk' AND strand = '+'
            ) AS risk_plus,
            MAX(score) FILTER (
                WHERE score_mode = 'log2_relative_risk' AND strand = '-'
            ) AS risk_minus,
            MAX(score) FILTER (
                WHERE score_mode = 'log_odds' AND strand = '+'
            ) AS odds_plus,
            MAX(score) FILTER (
                WHERE score_mode = 'log_odds' AND strand = '-'
            ) AS odds_minus
        FROM motif_score_dense
        WHERE motif_id = 'MA9999.1'
        GROUP BY start;

        CREATE OR REPLACE VIEW synthetic_comparison AS
        SELECT
            e.start,
            e.reference_base,
            e.risk_plus AS expected_risk_plus,
            a.risk_plus AS actual_risk_plus,
            e.risk_minus AS expected_risk_minus,
            a.risk_minus AS actual_risk_minus,
            e.odds_plus AS expected_odds_plus,
            a.odds_plus AS actual_odds_plus,
            e.odds_minus AS expected_odds_minus,
            a.odds_minus AS actual_odds_minus,
            (
                (ABS(e.risk_plus - a.risk_plus) < 1e-5 OR
                    (e.risk_plus IS NULL AND a.risk_plus IS NULL)) AND
                (ABS(e.risk_minus - a.risk_minus) < 1e-5 OR
                    (e.risk_minus IS NULL AND a.risk_minus IS NULL)) AND
                (ABS(e.odds_plus - a.odds_plus) < 1e-5 OR
                    (e.odds_plus IS NULL AND a.odds_plus IS NULL)) AND
                (ABS(e.odds_minus - a.odds_minus) < 1e-5 OR
                    (e.odds_minus IS NULL AND a.odds_minus IS NULL))
            ) AS matches
        FROM synthetic_expected e
        FULL OUTER JOIN synthetic_actual a USING (start);

        SELECT CASE
            WHEN (SELECT COUNT(*) FROM synthetic_comparison) <> 5
              OR NOT (SELECT BOOL_AND(matches) FROM synthetic_comparison)
            THEN error('synthetic scores differ from the inspectable fixture')
        END;

        SELECT CASE
            WHEN (SELECT COUNT(*) FROM dense_run_inventory) <> 4
              OR EXISTS (
                  SELECT 1 FROM dense_run_inventory
                  WHERE n_windows <> 5
                     OR n_valid_windows <> 4
                     OR n_skipped_windows <> 1
                     OR blocks <> 2
              )
            THEN error('synthetic dense inventory is incorrect')
        END;

        SELECT CASE
            WHEN EXISTS (SELECT 1 FROM synthetic_actual WHERE \"end\" <> start + 1)
            THEN error('one-base motif produced a non-unit interval')
        END;" >/dev/null

    duckdb -readonly -box "$database" -c "
        SELECT start, reference_base,
               ROUND(actual_risk_plus, 6) AS risk_plus,
               ROUND(actual_risk_minus, 6) AS risk_minus,
               ROUND(actual_odds_plus, 6) AS odds_plus,
               ROUND(actual_odds_minus, 6) AS odds_minus,
               matches
        FROM synthetic_comparison
        ORDER BY start;"
)

echo "Synthetic dense-score test passed: $output"

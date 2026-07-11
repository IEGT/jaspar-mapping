#!/usr/bin/env bash

set -euo pipefail

command -v duckdb >/dev/null 2>&1 || {
    echo "E: duckdb CLI is required for this test." >&2
    exit 1
}

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-dense-duckdb.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

metadata_dir="$temporary/tables/jaspar2026/motif_metadata"
dense_dir="$temporary/tables/jaspar2026/motif_score_dense/motif_id=MA0001.1/score_mode=log_odds/pseudocount=1/chrom=1/strand=plus"
mkdir -p "$metadata_dir" "$dense_dir"

(
    cd "$temporary"
    duckdb :memory: -bail -c "
        COPY (
            SELECT 'MA0001.1'::VARCHAR AS motif_id,
                   'TEST'::VARCHAR AS motif_name,
                   4::INTEGER AS motif_length,
                   2026::INTEGER AS jaspar_version,
                   'fixture'::VARCHAR AS source_sha256
        ) TO 'tables/jaspar2026/motif_metadata/part-000000.parquet'
          (FORMAT PARQUET, COMPRESSION ZSTD);
        COPY (
            SELECT 100::BIGINT AS block_start,
                   [1.0::FLOAT, NULL::FLOAT, -1.0::FLOAT] AS scores
        ) TO 'tables/jaspar2026/motif_score_dense/motif_id=MA0001.1/score_mode=log_odds/pseudocount=1/chrom=1/strand=plus/part-from=100-to=107-n_policy=skip-000000.parquet'
          (FORMAT PARQUET, COMPRESSION ZSTD);" >/dev/null

    duckdb dense.duckdb -bail -f \
        "$repository_root/sql/chr1_dense_dry_run_schema.sql" >/dev/null
    duckdb dense.duckdb -bail -c \
        "CREATE TABLE run_manifest AS SELECT 'fixture'::VARCHAR AS run_id;" >/dev/null

    duckdb dense.duckdb -bail -c "
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM dense_run_inventory
            WHERE n_windows = 3
              AND n_valid_windows = 2
              AND n_skipped_windows = 1
              AND alignment_start_begin = 100
              AND alignment_start_end = 103
        ) THEN error('dense inventory is incorrect') END;

        SELECT CASE WHEN (SELECT COUNT(*) FROM dense_scores_region(
            'MA0001.1', 'log_odds', 1.0, '1', '+', 100, 103
        )) <> 3 THEN error('region macro returned the wrong row count') END;

        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM dense_scores_region(
                'MA0001.1', 'log_odds', 1.0, '1', '+', 100, 103
            ) WHERE start = 101 AND score IS NULL AND \"end\" = 105
        ) THEN error('region macro did not preserve a skipped score') END;

        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM dense_score_summary(
                'MA0001.1', 'log_odds', 1.0, '1', '+', 100, 103
            ) WHERE n_windows = 3 AND n_valid_windows = 2 AND n_skipped_windows = 1
        ) THEN error('summary macro is incorrect') END;

        SELECT CASE WHEN (SELECT SUM(n_windows) FROM dense_score_histogram(
            'MA0001.1', 'log_odds', 1.0, '1', '+', 100, 103, 1.0
        )) <> 2 THEN error('histogram macro counted skipped scores') END;

        SELECT CASE WHEN (SELECT SUM(n_windows) FROM dense_score_calibration_bins(
            'MA0001.1', 'log_odds', 1.0, '1', '+', 100, 103
        )) <> 2 THEN error('calibration bins did not cover valid scores') END;

        SELECT CASE WHEN (SELECT COUNT(*) FROM dense_score_calibration_bins(
            'MA0001.1', 'log_odds', 1.0, '1', '+', 100, 103
        )) <> 15 THEN error('calibration bins omitted empty intervals') END;" >/dev/null
)

cp "$temporary/dense.duckdb" \
    "$temporary/jaspar2026_chr1_patz1_tp73.duckdb"
bash "$repository_root/scripts/inspect_chr1_dense_dry_run.sh" \
    --package "$temporary" summary MA0001.1 log_odds + 100 103 >/dev/null

echo "Dense chr1 DuckDB inspection tests passed."

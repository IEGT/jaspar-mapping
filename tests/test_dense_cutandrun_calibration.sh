#!/usr/bin/env bash

set -euo pipefail

command -v duckdb >/dev/null 2>&1 || {
    echo "E: duckdb CLI is required for this test." >&2
    exit 1
}

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-cutandrun-dense.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

for strand in plus minus; do
    mkdir -p "$temporary/package/tables/jaspar2026/motif_score_dense/"\
"motif_id=TEST.1/score_mode=log2_relative_risk/pseudocount=0/chrom=1/strand=$strand"
done

duckdb :memory: -bail -c "
    COPY (
        SELECT 0::BIGINT AS block_start,
               [5.0::FLOAT, '-Infinity'::FLOAT, NULL::FLOAT, 1.0::FLOAT,
                '-Infinity'::FLOAT] AS scores
    ) TO '$temporary/package/tables/jaspar2026/motif_score_dense/motif_id=TEST.1/score_mode=log2_relative_risk/pseudocount=0/chrom=1/strand=plus/part-from=0-to=7-n_policy=skip-000000.parquet'
      (FORMAT PARQUET, COMPRESSION ZSTD);
    COPY (
        SELECT 0::BIGINT AS block_start,
               [4.0::FLOAT, '-Infinity'::FLOAT, NULL::FLOAT, 2.0::FLOAT,
                0.0::FLOAT] AS scores
    ) TO '$temporary/package/tables/jaspar2026/motif_score_dense/motif_id=TEST.1/score_mode=log2_relative_risk/pseudocount=0/chrom=1/strand=minus/part-from=0-to=7-n_policy=skip-000000.parquet'
      (FORMAT PARQUET, COMPRESSION ZSTD);" >/dev/null

printf 'chr1\t0\t4\t7.5\n' > "$temporary/coverage.bedGraph"

"$repository_root/scripts/analyze_dense_cutandrun_coverage.py" \
    --package "$temporary/package" \
    --coverage-bed "$temporary/coverage.bedGraph" \
    --output-dir "$temporary/result" \
    --sample-id synthetic_negative_infinity \
    --motif-id TEST.1 \
    --motif-length 2 \
    --chrom 1 \
    --pseudocount 0 \
    --score-mode log2_relative_risk \
    --bin-width 0.1 \
    --threads 1 \
    --memory-limit 256MB \
    --max-temp-size 256MB >/dev/null

duckdb -readonly "$temporary/result/calibration.duckdb" -bail -c "
    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM calibration_summary
        WHERE score_mode = 'log2_relative_risk'
          AND n_valid_motifs = 4
          AND n_supported_motifs = 1
          AND n_score_bins = 4
          AND roc_auc_binned = 0.0
          AND average_precision = 0.25
    ) THEN error('negative-infinity calibration summary is incorrect') END;

    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM threshold_curve
        WHERE isinf(threshold) AND threshold < 0
          AND selected_motifs = 4
          AND supported_selected_motifs = 1
          AND motif_recall = 1.0
          AND motif_false_positive_rate = 1.0
    ) THEN error('negative-infinity threshold endpoint is missing') END;

    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM immersed_motif
        WHERE start = 1 AND isinf(score) AND score < 0
          AND effective_max_depth = 7.5
    ) THEN error('bedGraph depth was not retained for an immersed score') END;

    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM run_config
        WHERE coverage_format = 'bedgraph'
          AND coverage_depth_semantics = 'column_4_signal'
    ) THEN error('bedGraph depth provenance is missing') END;

    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM calibration_summary
        WHERE coverage_format = 'bedgraph'
          AND coverage_depth_semantics = 'column_4_signal'
    ) THEN error('bedGraph provenance is missing from the summary') END;" >/dev/null

echo "Dense CUT&RUN calibration tests passed."

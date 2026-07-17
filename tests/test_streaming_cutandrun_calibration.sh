#!/usr/bin/env bash

set -euo pipefail

command -v duckdb >/dev/null 2>&1 || {
    echo "E: duckdb CLI is required for this test." >&2
    exit 1
}

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "/tmp/jaspar-cutandrun-streaming.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb :memory: -bail -c "
    COPY (
        SELECT 0::BIGINT AS block_start,
               [0.0::FLOAT, 1.0::FLOAT, 2.0::FLOAT, 3.0::FLOAT,
                4.0::FLOAT, 5.0::FLOAT, 6.0::FLOAT,
                '-Infinity'::FLOAT] AS scores
    ) TO '$temporary/plus.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
    COPY (
        SELECT 0::BIGINT AS block_start,
               [0.0::FLOAT, 0.0::FLOAT, 0.0::FLOAT, 0.0::FLOAT,
                0.0::FLOAT, 0.0::FLOAT, 0.0::FLOAT,
                '-Infinity'::FLOAT] AS scores
    ) TO '$temporary/minus.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
    COPY (
        SELECT 0::BIGINT AS block_start,
               [0.0::FLOAT, 10.0::FLOAT, 1.0::FLOAT, 2.0::FLOAT,
                3.0::FLOAT, 4.0::FLOAT, 5.0::FLOAT, 6.0::FLOAT,
                7.0::FLOAT] AS scores
    ) TO '$temporary/context-plus.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
    COPY (
        SELECT 0::BIGINT AS block_start,
               [0.0::FLOAT, 0.0::FLOAT, 0.0::FLOAT, 0.0::FLOAT,
                0.0::FLOAT, 0.0::FLOAT, 0.0::FLOAT, 0.0::FLOAT,
                0.0::FLOAT] AS scores
    ) TO '$temporary/context-minus.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);" >/dev/null

cat > "$temporary/coverage.bedGraph" <<'EOF'
chr1	0	2	1
chr1	2	4	3
chr1	5	9	2
EOF

"$repository_root/cutandrun_score_calibration" \
    --plus-parquet "$temporary/plus.parquet" \
    --minus-parquet "$temporary/minus.parquet" \
    --coverage test="$temporary/coverage.bedGraph" \
    --output-dir "$temporary/result" \
    --motif-id TEST.1 \
    --motif-length 2 \
    --chrom 1 \
    --score-mode log2_relative_risk \
    --pseudocount 0 \
    --bin-width 1 \
    --progress-every 0 >/dev/null

"$repository_root/cutandrun_score_calibration" \
    --plus-parquet "$temporary/plus.parquet" \
    --minus-parquet "$temporary/minus.parquet" \
    --context-plus-parquet "$temporary/context-plus.parquet" \
    --context-minus-parquet "$temporary/context-minus.parquet" \
    --coverage test="$temporary/coverage.bedGraph" \
    --output-dir "$temporary/joint-result" \
    --feature-parquet "$temporary/joint-result/anchor_features.parquet" \
    --motif-id TEST.1 --motif-length 2 \
    --context-motif-id CONTEXT.1 --context-motif-length 1 \
    --context-flank 1 --minimum-anchor-score -100 \
    --chrom 1 --score-mode log2_relative_risk \
    --pseudocount 0 --context-pseudocount 0 \
    --bin-width 1 --progress-every 0 >/dev/null

duckdb :memory: -bail -c "
    CREATE VIEW histogram AS
    SELECT * FROM read_csv('$temporary/result/score_histogram.tsv',
                           delim='\t', header=true);
    CREATE VIEW summary AS
    SELECT * FROM read_csv('$temporary/result/summary.tsv',
                           delim='\t', header=true);
    CREATE VIEW curve AS
    SELECT * FROM read_csv('$temporary/result/threshold_curve.tsv',
                           delim='\t', header=true);

    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM summary
        WHERE sample_id = 'test'
          AND n_alignment_starts = 8
          AND n_valid_motifs = 8
          AND n_supported_motifs = 2
          AND n_coverage_intervals = 3
          AND n_coverage_components = 2
          AND covered_bases = 8
          AND possible_immersed_starts = 2
    ) THEN error('streaming summary is incorrect') END;

    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM histogram
        WHERE threshold = 1
          AND n_supported = 1
          AND effective_depth_sum = 3
    ) THEN error('adjacent-run immersion or maximum depth is incorrect') END;

    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM histogram
        WHERE threshold = 6
          AND n_supported = 1
          AND effective_depth_sum = 2
    ) THEN error('second component support is incorrect') END;

    SELECT CASE WHEN EXISTS (
        SELECT 1 FROM histogram
        WHERE threshold IN (0, 2) AND n_supported <> 0
    ) THEN error('strict boundary rule was not enforced') END;

    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM curve
        WHERE threshold = '-Inf'
          AND selected_motifs = 8
          AND supported_selected_motifs = 2
          AND motif_recall = 1
          AND motif_false_positive_rate = 1
    ) THEN error('negative-infinity endpoint is incorrect') END;" >/dev/null

duckdb :memory: -bail -c "
    CREATE VIEW joint AS
    SELECT * FROM read_csv('$temporary/joint-result/joint_score_histogram.tsv',
                           delim='\t', header=true);
    CREATE VIEW config AS
    SELECT * FROM read_csv('$temporary/joint-result/joint_run_config.tsv',
                           delim='\t', header=true);
    CREATE VIEW features AS
    SELECT * FROM read_parquet(
        '$temporary/joint-result/anchor_features.parquet');

    SELECT CASE WHEN (SELECT SUM(n_total) FROM joint) <> 7
        THEN error('joint candidate count is incorrect') END;
    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM joint
        WHERE anchor_threshold = 1 AND context_threshold = 10
          AND n_total = 1 AND n_supported = 1
          AND n_same_orientation = 1 AND n_supported_same_orientation = 1
          AND mean_signed_center_distance_bp = -0.5
    ) THEN error('sliding context maximum or geometry is incorrect') END;
    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM joint
        WHERE anchor_threshold = 6 AND context_threshold = 6
          AND n_total = 1 AND n_supported = 1
          AND effective_depth_sum = 2
    ) THEN error('joint coverage evidence is incorrect') END;
    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM config
        WHERE anchor_motif_id = 'TEST.1' AND context_motif_id = 'CONTEXT.1'
          AND context_flank_bp = 1 AND minimum_anchor_score = -100
          AND n_joint_candidates = 7
          AND feature_parquet = '$temporary/joint-result/anchor_features.parquet'
    ) THEN error('joint run provenance is incomplete') END;" >/dev/null

duckdb :memory: -bail -c "
    CREATE VIEW features AS
    SELECT * FROM read_parquet(
        '$temporary/joint-result/anchor_features.parquet');

    SELECT CASE WHEN (SELECT COUNT(*) FROM features) <> 7
        THEN error('exact feature row count is incorrect') END;
    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM features
        WHERE chrom = '1'
          AND anchor_start = 1 AND anchor_end = 3
          AND anchor_score = 1 AND anchor_strand = 1
          AND context_start = 1 AND context_end = 2
          AND context_score = 10 AND context_strand = 1
          AND context_center_distance_bp = -0.5
          AND same_orientation
          AND supported_test AND depth_test = 3
    ) THEN error('exact feature score, geometry, or coverage is incorrect') END;
    SELECT CASE WHEN NOT EXISTS (
        SELECT 1 FROM features
        WHERE anchor_start = 0 AND NOT supported_test AND depth_test = 0
    ) THEN error('unsupported feature row is incorrect') END;" >/dev/null

echo "Streaming CUT&RUN calibration tests passed."

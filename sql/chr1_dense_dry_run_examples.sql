-- Run from a generated package directory:
--   duckdb -readonly jaspar2026_chr1_patz1_tp73.duckdb \
--     -f /path/to/jaspar-mapping/sql/chr1_dense_dry_run_examples.sql

SELECT * FROM run_manifest;
SELECT * FROM motif_metadata ORDER BY motif_id;
SELECT * FROM dense_run_inventory ORDER BY motif_id, score_mode, strand;

SELECT *
FROM dense_score_summary(
    'MA0861.2', 'log2_relative_risk', 1.0, '1', '+', 3651800, 3652600
);

SELECT *
FROM dense_scores_region(
    'MA1961.2', 'log_odds', 1.0, '1', '+', 3651800, 3651820
)
ORDER BY start;

SELECT *
FROM dense_score_histogram(
    'MA0861.2', 'log_odds', 1.0, '1', '+', 3600000, 3700000, 0.5
);

SELECT *
FROM dense_score_calibration_bins(
    'MA0861.2', 'log2_relative_risk', 1.0, '1', '+', 3600000, 3700000
);

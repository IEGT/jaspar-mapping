-- Query contract for a finalized sparse whole-genome scan package.
-- Build from the package root. Metadata is copied into the small DuckDB catalog
-- so it remains queryable from any working directory. Coordinates are BED
-- 0-based half-open.

CREATE OR REPLACE TABLE motif_set AS
SELECT * FROM read_parquet('tables/jaspar2026/motif_set/part-000000.parquet');

CREATE OR REPLACE TABLE genome AS
SELECT * FROM read_parquet('tables/jaspar2026/genome/part-000000.parquet');

CREATE OR REPLACE TABLE sequence_region AS
SELECT * FROM read_parquet('tables/jaspar2026/sequence_region/part-000000.parquet');

CREATE OR REPLACE TABLE motif_metadata AS
SELECT * FROM read_parquet('tables/jaspar2026/motif_metadata/part-000000.parquet');

CREATE OR REPLACE TABLE scan_run AS
SELECT * FROM read_parquet('tables/jaspar2026/scan_run/part-000000.parquet');

CREATE OR REPLACE TABLE scan_threshold_policy AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/scan_threshold_policy/part-000000.parquet'
);

CREATE OR REPLACE TABLE scan_task AS
SELECT * FROM read_parquet('tables/jaspar2026/scan_task/part-000000.parquet');

CREATE OR REPLACE TABLE scan_file_inventory AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/scan_file_inventory/part-000000.parquet'
);

-- Backward-compatible fallback for packages finalized before detailed
-- per-motif threshold provenance was stored. The catalog builder replaces this
-- table from scan_motif_threshold Parquet when that file is present.
CREATE OR REPLACE TABLE scan_motif_threshold AS
SELECT
    run_id,
    genome_id,
    motif_set_id,
    ('legacy_' || run_id)::VARCHAR AS threshold_set_id,
    NULL::VARCHAR AS threshold_registry_sha256,
    motif_id,
    NULL::DOUBLE AS informative_threshold,
    'legacy_scan_inventory'::VARCHAR AS informative_source,
    NULL::DOUBLE AS default_minimum_score,
    MIN(CAST(minimum_score AS DOUBLE)) AS candidate_minimum_score,
    NULL::DOUBLE AS density_minimum_spacing_bp,
    NULL::BIGINT AS density_maximum_loci,
    NULL::DOUBLE AS density_threshold,
    MIN(CAST(minimum_score AS DOUBLE)) AS final_minimum_score,
    NULL::BOOLEAN AS density_limited,
    NULL::VARCHAR AS density_chrom,
    NULL::BIGINT AS valid_locus_starts,
    NULL::BIGINT AS skipped_locus_starts,
    NULL::BIGINT AS loci_at_candidate_threshold,
    NULL::BIGINT AS loci_at_final_threshold,
    NULL::DOUBLE AS mean_spacing_bp_at_final_threshold,
    NULL::VARCHAR AS distribution_sha256
FROM scan_file_inventory
GROUP BY run_id, genome_id, motif_set_id, motif_id;

-- Hit payloads are intentionally not bound here. Opening a package-wide glob
-- made DuckDB inspect every Parquet footer (131,650 files in the GRCh38 run)
-- merely to create the catalog. Call motif_hit_files() with the exact paths
-- selected from scan_file_inventory, normally through query_genome_scan.py.
CREATE OR REPLACE MACRO motif_hit_files(file_paths) AS TABLE
SELECT
    t.run_id,
    h.task_id,
    h.genome_id,
    h.motif_set_id,
    CAST(h.chrom AS VARCHAR) AS chrom,
    h.start,
    h."end",
    h.motif_id,
    m.motif_name,
    CASE h.strand WHEN 'plus' THEN '+' WHEN 'minus' THEN '-'
         ELSE h.strand END AS strand,
    CAST(h.score AS DOUBLE) AS score,
    h.score_mode,
    CAST(h.pseudocount AS DOUBLE) AS pseudocount,
    h.background_model_id,
    h.pseudocount_scheme,
    TRY_CAST(h.minimum_score AS DOUBLE) AS minimum_score,
    TRY_CAST(h.minimum_pwm_relative_score AS DOUBLE)
        AS minimum_pwm_relative_score,
    TRY_CAST(h.maximum_pwm_relative_score AS DOUBLE)
        AS maximum_pwm_relative_score,
    h.n_policy,
    h.matched_sequence AS matched_sequence_policy,
    'bed'::VARCHAR AS coordinate_mode,
    CAST(h.pwm_relative_score AS DOUBLE) AS pwm_relative_score,
    h.matched_seq
FROM read_parquet(file_paths, hive_partitioning = true) h
JOIN motif_metadata m USING (motif_set_id, motif_id)
JOIN scan_task t ON t.task_id = h.task_id;

CREATE OR REPLACE VIEW scan_inventory_summary AS
SELECT
    run_id,
    genome_id,
    motif_set_id,
    state,
    COUNT(*) AS n_files,
    COUNT(DISTINCT task_id) AS n_tasks,
    COUNT(DISTINCT motif_id) AS n_motifs,
    COUNT(DISTINCT chrom) AS n_sequence_regions,
    SUM(expected_windows) AS expected_windows,
    SUM(skipped_n_windows) AS skipped_n_windows,
    SUM(sentinel_score_windows) AS sentinel_score_windows,
    SUM(below_minimum_score_windows) AS below_minimum_score_windows,
    SUM(pwm_filtered_windows) AS pwm_filtered_windows,
    SUM(emitted_hits) AS emitted_hits,
    SUM(bytes) AS bytes
FROM scan_file_inventory
GROUP BY run_id, genome_id, motif_set_id, state;

CREATE OR REPLACE VIEW scan_task_completeness AS
SELECT
    t.run_id,
    t.task_id,
    t.chrom,
    t.policy_id,
    t.minimum_score,
    t.motif_count,
    t.state,
    COUNT(i.output_relative_path) AS inventory_files,
    2 * t.motif_count AS expected_files,
    COUNT(i.output_relative_path) = 2 * t.motif_count AS complete
FROM scan_task t
LEFT JOIN scan_file_inventory i USING (run_id, task_id)
GROUP BY t.run_id, t.task_id, t.chrom, t.policy_id, t.minimum_score,
         t.motif_count, t.state;

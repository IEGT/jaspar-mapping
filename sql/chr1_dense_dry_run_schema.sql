-- Dense chr1 PATZ1/TP73 dry-run query surface.
-- Run this from the package directory containing tables/jaspar2026/.

CREATE OR REPLACE VIEW motif_metadata AS
SELECT motif_id, motif_name, motif_length, jaspar_version, source_sha256
FROM read_parquet('tables/jaspar2026/motif_metadata/*.parquet');

CREATE OR REPLACE VIEW motif_score_dense_block AS
SELECT
    motif_id,
    score_mode,
    CAST(pseudocount AS DOUBLE) AS pseudocount,
    CAST(chrom AS VARCHAR) AS chrom,
    CASE
        WHEN strand IN ('plus', '+') THEN '+'
        WHEN strand IN ('minus', '-') THEN '-'
        ELSE strand
    END AS strand,
    block_start,
    scores,
    filename AS source_file
FROM read_parquet(
    'tables/jaspar2026/motif_score_dense/*/*/*/*/*/*.parquet',
    hive_partitioning = true,
    filename = true,
    union_by_name = true
);

-- This inventory never expands score lists, so it remains fast for full chr1.
CREATE OR REPLACE VIEW dense_run_inventory AS
SELECT
    motif_id,
    score_mode,
    pseudocount,
    chrom,
    strand,
    COUNT(DISTINCT source_file) AS part_files,
    COUNT(*) AS blocks,
    MIN(block_start) AS alignment_start_begin,
    MAX(block_start + len(scores)) AS alignment_start_end,
    SUM(len(scores)) AS n_windows,
    SUM(list_count(scores)) AS n_valid_windows,
    SUM(len(scores) - list_count(scores)) AS n_skipped_windows
FROM motif_score_dense_block
GROUP BY motif_id, score_mode, pseudocount, chrom, strand;

-- Convenient logical row view. Filter by all partition columns and a small
-- region; do not SELECT the complete view for a full chromosome.
CREATE OR REPLACE VIEW motif_score_dense AS
SELECT
    b.chrom,
    b.block_start + CAST(u.offset_one - 1 AS BIGINT) AS start,
    b.block_start + CAST(u.offset_one - 1 AS BIGINT) + m.motif_length AS "end",
    b.motif_id,
    m.motif_name,
    b.strand,
    b.score_mode,
    b.pseudocount,
    u.score
FROM motif_score_dense_block b
JOIN motif_metadata m USING (motif_id)
CROSS JOIN UNNEST(b.scores) WITH ORDINALITY AS u(score, offset_one);

CREATE OR REPLACE TABLE score_calibration_bin AS
SELECT *
FROM (
    VALUES
        (1,  '(-Inf,-10000)', CAST(NULL AS DOUBLE), -10000.0),
        (2,  '[-10000,-200)', -10000.0, -200.0),
        (3,  '[-200,-50)', -200.0, -50.0),
        (4,  '[-50,-10)', -50.0, -10.0),
        (5,  '[-10,-5)', -10.0, -5.0),
        (6,  '[-5,-2)', -5.0, -2.0),
        (7,  '[-2,-1)', -2.0, -1.0),
        (8,  '[-1,-0.5)', -1.0, -0.5),
        (9,  '[-0.5,0)', -0.5, 0.0),
        (10, '[0,0.5)', 0.0, 0.5),
        (11, '[0.5,1)', 0.5, 1.0),
        (12, '[1,2)', 1.0, 2.0),
        (13, '[2,5)', 2.0, 5.0),
        (14, '[5,10)', 5.0, 10.0),
        (15, '[10,+Inf)', 10.0, CAST(NULL AS DOUBLE))
) AS t(bin_order, bin_label, lower_bound, upper_bound);

CREATE OR REPLACE MACRO dense_scores_region(
    p_motif_id,
    p_score_mode,
    p_pseudocount,
    p_chrom,
    p_strand,
    p_start,
    p_end
) AS TABLE
WITH selected_blocks AS (
    SELECT b.*, m.motif_name, m.motif_length
    FROM motif_score_dense_block b
    JOIN motif_metadata m USING (motif_id)
    WHERE b.motif_id = p_motif_id
      AND b.score_mode = p_score_mode
      AND b.pseudocount = p_pseudocount
      AND b.chrom = p_chrom
      AND b.strand = p_strand
      AND b.block_start < p_end
      AND b.block_start + len(b.scores) > p_start
)
SELECT
    b.chrom,
    b.block_start + CAST(u.offset_one - 1 AS BIGINT) AS start,
    b.block_start + CAST(u.offset_one - 1 AS BIGINT) + b.motif_length AS "end",
    b.motif_id,
    b.motif_name,
    b.strand,
    b.score_mode,
    b.pseudocount,
    u.score
FROM selected_blocks b
CROSS JOIN UNNEST(b.scores) WITH ORDINALITY AS u(score, offset_one)
WHERE b.block_start + CAST(u.offset_one - 1 AS BIGINT) >= p_start
  AND b.block_start + CAST(u.offset_one - 1 AS BIGINT) < p_end;

CREATE OR REPLACE MACRO dense_score_summary(
    p_motif_id,
    p_score_mode,
    p_pseudocount,
    p_chrom,
    p_strand,
    p_start,
    p_end
) AS TABLE
SELECT
    COUNT(*) AS n_windows,
    COUNT(score) AS n_valid_windows,
    COUNT(*) - COUNT(score) AS n_skipped_windows,
    COUNT(*) FILTER (WHERE isinf(score) AND score < 0) AS n_negative_infinity_windows,
    COUNT(*) FILTER (WHERE isfinite(score)) AS n_finite_windows,
    MIN(score) AS min_score,
    AVG(score) AS mean_score,
    MAX(score) AS max_score,
    approx_quantile(score, [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]) AS quantiles
FROM dense_scores_region(
    p_motif_id, p_score_mode, p_pseudocount,
    p_chrom, p_strand, p_start, p_end
);

CREATE OR REPLACE MACRO dense_score_histogram(
    p_motif_id,
    p_score_mode,
    p_pseudocount,
    p_chrom,
    p_strand,
    p_start,
    p_end,
    p_bin_width
) AS TABLE
WITH checked_width AS (
    SELECT CASE
        WHEN p_bin_width > 0 THEN CAST(p_bin_width AS DOUBLE)
        ELSE error('bin width must be greater than zero')
    END AS bin_width
),
valid_scores AS (
    SELECT score
    FROM dense_scores_region(
        p_motif_id, p_score_mode, p_pseudocount,
        p_chrom, p_strand, p_start, p_end
    )
    WHERE score IS NOT NULL
)
SELECT
    CASE WHEN isinf(score) AND score < 0
         THEN '-Infinity'::DOUBLE
         ELSE ROUND(FLOOR(score / bin_width) * bin_width, 12)
    END AS bin_start,
    CASE WHEN isinf(score) AND score < 0
         THEN '-Infinity'::DOUBLE
         ELSE ROUND((FLOOR(score / bin_width) + 1) * bin_width, 12)
    END AS bin_end,
    COUNT(*) AS n_windows,
    COUNT(*)::DOUBLE / SUM(COUNT(*)) OVER () AS fraction
FROM valid_scores
CROSS JOIN checked_width
GROUP BY bin_start, bin_end
ORDER BY bin_start;

CREATE OR REPLACE MACRO dense_score_calibration_bins(
    p_motif_id,
    p_score_mode,
    p_pseudocount,
    p_chrom,
    p_strand,
    p_start,
    p_end
) AS TABLE
WITH valid_scores AS (
    SELECT score
    FROM dense_scores_region(
        p_motif_id, p_score_mode, p_pseudocount,
        p_chrom, p_strand, p_start, p_end
    )
    WHERE score IS NOT NULL
),
bin_counts AS (
    SELECT b.bin_order, COUNT(*) AS n_windows
    FROM valid_scores s
    JOIN score_calibration_bin b
      ON (b.lower_bound IS NULL OR s.score >= b.lower_bound)
     AND (b.upper_bound IS NULL OR s.score < b.upper_bound)
    GROUP BY b.bin_order
),
score_total AS (
    SELECT COUNT(*) AS n_windows
    FROM valid_scores
)
SELECT
    b.bin_order,
    b.bin_label,
    b.lower_bound,
    b.upper_bound,
    COALESCE(c.n_windows, 0) AS n_windows,
    CASE
        WHEN t.n_windows = 0 THEN 0.0
        ELSE COALESCE(c.n_windows, 0)::DOUBLE / t.n_windows
    END AS fraction
FROM score_calibration_bin b
LEFT JOIN bin_counts c USING (bin_order)
CROSS JOIN score_total t
ORDER BY b.bin_order;

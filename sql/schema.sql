-- ============================================================================
-- schema.sql — DuckDB contract over the versioned JASPAR-2026 table package
-- ----------------------------------------------------------------------------
-- DRAFT. Targets the layout described in docs/promoter_architecture_ml_schema.md.
-- Nothing here materializes data that does not yet exist; it defines the query
-- surface that the exporter must fill.
--
-- Conventions:
--   * All coordinates are BED 0-based half-open (run pssm_scan --coordinate-mode bed).
--   * Parquet is partitioned by chrom (hive style: .../chrom=1/part-*.parquet),
--     sorted within partition by start, ZSTD-compressed.
--   * VIEWs are zero-copy over Parquet. TABLEs are materialized because they
--     carry an interval join or an aggregation we do not want to repeat per query.
--   * Run duckdb from the exported package root (the directory containing
--     tables/jaspar2026/). The Parquet paths below are intentionally literal:
--     DuckDB session variables are not persisted inside views.
-- ============================================================================

-- ---------------------------------------------------------------------------
-- Layer 1: raw motif mapping (long form, one row per hit)
--   Produced directly by pssm_scan. Deterministic from (genome, JASPAR, params).
-- ---------------------------------------------------------------------------
CREATE OR REPLACE VIEW motif_metadata AS
SELECT
    motif_id,
    motif_name,
    motif_length,
    jaspar_version,
    source_sha256
FROM read_parquet('tables/jaspar2026/motif_metadata/*.parquet');

CREATE OR REPLACE VIEW motif_hit AS
SELECT
    CAST(chrom AS VARCHAR) AS chrom, -- e.g. '1', 'X'  (partition key)
    start,                -- 0-based half-open
    "end",
    motif_id,             -- JASPAR accession, e.g. 'MA0861.2'
    motif_name,           -- e.g. 'TP73'
    strand,               -- '+' / '-'
    score,                -- PSSM score in the run's score_mode
    score_mode,           -- log2_relative_risk | log_odds | raw_counts
    pseudocount,          -- count added to each A/C/G/T motif entry before normalization
    pwm_relative_score,   -- 0..1, comparable across motifs
    matched_seq           -- nullable; present only if scanned with --show-sequence
FROM read_parquet('tables/jaspar2026/motif_hit/*/*.parquet', hive_partitioning = 1);

-- Dense per-base score blocks for calibration runs. Physical Parquet stores
-- only block_start and a FLOAT[] score vector; identity lives in hive
-- partitions: motif_id, score_mode, pseudocount, chrom, strand.
CREATE OR REPLACE VIEW motif_score_dense_block AS
SELECT
    motif_id,
    score_mode,
    CAST(pseudocount AS DOUBLE) AS pseudocount,
    CAST(chrom AS VARCHAR) AS chrom,
    CASE WHEN strand IN ('plus', '+') THEN '+'
         WHEN strand IN ('minus', '-') THEN '-'
         ELSE strand END AS strand,
    block_start,
    scores
FROM read_parquet('tables/jaspar2026/motif_score_dense/*/*/*/*/*/*.parquet', hive_partitioning = 1);

-- Logical row view over dense score blocks. Use this for region inspection;
-- aggregate directly from blocks for whole-chromosome calibration jobs.
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

-- Shared score bins for the dense chr1 CUT&RUN calibration. Intervals are
-- lower-inclusive and upper-exclusive, except the first open lower tail and
-- the final open upper tail.
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

-- Precomputed dense-score × CUT&RUN summaries. These are the calibration
-- products used to choose genome-wide storage thresholds.
CREATE OR REPLACE VIEW motif_cutandrun_score_bin_stats AS
SELECT
    motif_id,
    motif_name,
    CAST(chrom AS VARCHAR) AS chrom,
    strand,
    score_mode,
    CAST(pseudocount AS DOUBLE) AS pseudocount,
    bin_order,
    bin_label,
    lower_bound,
    upper_bound,
    cell_line,
    isoform,
    antibody,
    replicate,
    n_windows,
    n_covered_windows,
    overlap_fraction,
    baseline_fraction,
    enrichment_ratio,
    log2_enrichment,
    mean_signal,
    max_signal
FROM read_parquet('tables/jaspar2026/motif_cutandrun_score_bin_stats/*/*/*/*/*.parquet', hive_partitioning = 1);

-- ---------------------------------------------------------------------------
-- Layer 2a: motif-level architecture dimension (small, one row per motif_id)
--   Curated for the TP53 family; 'unknown'/NULL elsewhere. See schema notes.
-- ---------------------------------------------------------------------------
CREATE OR REPLACE VIEW motif_architecture AS
SELECT
    motif_id,
    motif_name,
    tf_family,                 -- e.g. 'TP53_family'
    binding_unit_model,        -- monomer|homodimer|heterodimer|dimer_of_dimers|tetramer|unknown
    oligomeric_state,          -- integer where known (e.g. 4)
    dna_binding_unit_count,
    half_site_count,
    half_site_pattern,         -- e.g. 'RRRCWWGYYY'
    spacer_min_bp,
    spacer_max_bp,
    architecture_source,       -- citation / PMCID
    architecture_confidence    -- curated|database_inferred|unknown
FROM read_parquet('tables/jaspar2026/motif_architecture/*.parquet');

-- ---------------------------------------------------------------------------
-- Layer 2b: per-hit architecture (optional; keyed to motif_hit)
--   Derived from the matched sequence, never written back into layer 1.
-- ---------------------------------------------------------------------------
CREATE OR REPLACE VIEW hit_architecture AS
SELECT
    CAST(chrom AS VARCHAR) AS chrom, start, motif_id, strand,
    spacer_bp,
    half_site_1_score,
    half_site_2_score,
    architecture_class,        -- full_site|half_site|degenerate_full_site
    oligomer_compatible,       -- boolean
    nucleosome_context_model   -- not_assessed|edge_accessible|dyad_centered|orientation_sensitive
FROM read_parquet('tables/jaspar2026/hit_architecture/*/*.parquet', hive_partitioning = 1);

-- ---------------------------------------------------------------------------
-- Layer 3: experimental evidence (CUT&RUN), long form
--   One row per (locus, sample). Sample decomposed into its facets so TA/DN
--   contrasts are a GROUP BY, not a column-name parse.
-- ---------------------------------------------------------------------------
CREATE OR REPLACE VIEW cutandrun AS
SELECT
    CAST(chrom AS VARCHAR) AS chrom,
    start, "end",             -- the reference locus the signal was mapped onto
    cell_line,                 -- 'saos2' | 'skmel29'
    isoform,                   -- 'TA' | 'DN' | 'GFP'
    antibody,                  -- e.g. 'p73' | 'pos_ctrl'
    replicate,                 -- 'R1', ...
    signal                     -- bedtools map value (null -> 0 upstream)
FROM read_parquet('tables/jaspar2026/cutandrun/*/*.parquet', hive_partitioning = 1);

-- ---------------------------------------------------------------------------
-- Layer 4: regulatory / gene context
-- ---------------------------------------------------------------------------
CREATE OR REPLACE VIEW gene AS
SELECT gene_id, gene_name, chrom, strand, gene_start, gene_end, tss, biotype
FROM read_parquet('tables/jaspar2026/gene/*.parquet');

CREATE OR REPLACE VIEW promoter AS
SELECT gene_id, transcript_id, chrom, strand, promoter_start, promoter_end, tss
FROM read_parquet('tables/jaspar2026/promoter/*.parquet');

CREATE OR REPLACE VIEW gene_set AS
SELECT gene_id, set_name           -- many-to-many (HALLMARK_*, KEGG_*, ...)
FROM read_parquet('tables/jaspar2026/gene_set/*.parquet');

-- Expression label source. `expression` is tidy per (gene, cell_line, isoform);
-- `expression_differential` is the precomputed TA-vs-DN contrast = the ML label.
CREATE OR REPLACE VIEW expression AS
SELECT gene_id, cell_line, isoform, value, replicate
FROM read_parquet('tables/jaspar2026/expression/*.parquet');

CREATE OR REPLACE VIEW expression_differential AS
SELECT
    gene_id, cell_line,
    log2fc_ta_vs_dn,           -- + means higher under TAp73 than DNp73
    lfc_se, pvalue, padj, base_mean
FROM read_parquet('tables/jaspar2026/expression_differential/*.parquet');

-- ===========================================================================
-- Materialized joins — the two things we never want to recompute per query.
-- ===========================================================================

-- Bridge: which motif hits fall inside which promoter, with strand-aware TSS
-- distance. This is the one interval join in the whole system; do it once.
CREATE OR REPLACE TABLE promoter_motif_hit AS
SELECT
    p.gene_id,
    p.transcript_id,
    h.motif_id,
    h.strand,
    h.score,
    h.score_mode,
    h.pseudocount,
    h.pwm_relative_score,
    h.start,
    h."end",
    -- signed distance from TSS, positive = downstream of TSS in gene direction
    CASE WHEN p.strand = '+' THEN h.start - p.tss
         ELSE p.tss - h."end" END AS tss_distance
FROM promoter p
JOIN motif_hit h
  ON  h.chrom = p.chrom
  AND h.start >= p.promoter_start
  AND h."end" <= p.promoter_end;

-- Per-(gene, motif) architecture summary in LONG form (no ~8000-col matrix).
-- Pivot to wide only at the ML boundary (see queries.sql / the Python builder).
CREATE OR REPLACE TABLE promoter_arch_feature AS
SELECT
    b.gene_id,
    b.motif_id,
    b.score_mode,
    b.pseudocount,
    a.tf_family,
    a.binding_unit_model,
    COUNT(*)                                    AS n_hits,
    MAX(b.pwm_relative_score)                   AS max_rel_score,
    SUM(CASE WHEN b.pwm_relative_score >= 0.8 THEN 1 ELSE 0 END) AS n_strong_hits,
    MIN(ABS(b.tss_distance))                    AS min_abs_tss_distance,
    BOOL_OR(a.binding_unit_model = 'dimer_of_dimers') AS has_dimer_of_dimers
FROM promoter_motif_hit b
LEFT JOIN motif_architecture a USING (motif_id)
GROUP BY b.gene_id, b.motif_id, b.score_mode, b.pseudocount, a.tf_family, a.binding_unit_model;

-- The ML surface: features (long) + the TA-vs-DN label. Keep long; pivot at read.
CREATE OR REPLACE VIEW ml_ta_vs_dn AS
SELECT
    f.gene_id,
    d.cell_line,
    f.motif_id,
    f.score_mode,
    f.pseudocount,
    f.tf_family,
    f.n_hits,
    f.max_rel_score,
    f.n_strong_hits,
    f.min_abs_tss_distance,
    f.has_dimer_of_dimers,
    d.log2fc_ta_vs_dn,
    d.padj,
    CASE WHEN d.log2fc_ta_vs_dn > 0 THEN 1 ELSE 0 END AS label_ta_up   -- classification
FROM promoter_arch_feature f
JOIN expression_differential d USING (gene_id);

-- Hot agent lookup: one compact, low-latency promoter summary per gene.
-- This is the first stop for OpenClaw / Claude / Codex before deeper joins.
CREATE OR REPLACE TABLE promoter_card AS
WITH promoter_span AS (
    SELECT
        gene_id,
        MIN(promoter_start) AS promoter_start,
        MAX(promoter_end) AS promoter_end,
        COUNT(DISTINCT transcript_id) AS n_promoter_transcripts
    FROM promoter
    GROUP BY gene_id
),
feature_summary AS (
    SELECT
        gene_id,
        MIN(score_mode) AS score_mode,
        MIN(pseudocount) AS pseudocount,
        COUNT(DISTINCT motif_id) AS n_motifs_with_hits,
        SUM(n_hits) AS n_motif_hits,
        MAX(max_rel_score) AS max_pwm_relative_score,
        MIN(min_abs_tss_distance) AS min_abs_tss_distance,
        BOOL_OR(COALESCE(has_dimer_of_dimers, false)) AS has_dimer_of_dimers,
        COUNT(DISTINCT tf_family) FILTER (WHERE tf_family IS NOT NULL) AS n_tf_families
    FROM promoter_arch_feature
    GROUP BY gene_id
)
SELECT
    g.gene_id,
    g.gene_name,
    g.chrom,
    g.strand,
    g.tss,
    ps.promoter_start,
    ps.promoter_end,
    fs.score_mode,
    fs.pseudocount,
    COALESCE(ps.n_promoter_transcripts, 0) AS n_promoter_transcripts,
    COALESCE(fs.n_motifs_with_hits, 0) AS n_motifs_with_hits,
    COALESCE(fs.n_motif_hits, 0) AS n_motif_hits,
    fs.max_pwm_relative_score,
    fs.min_abs_tss_distance,
    COALESCE(fs.has_dimer_of_dimers, false) AS has_dimer_of_dimers,
    COALESCE(fs.n_tf_families, 0) AS n_tf_families
FROM gene g
LEFT JOIN promoter_span ps USING (gene_id)
LEFT JOIN feature_summary fs USING (gene_id);

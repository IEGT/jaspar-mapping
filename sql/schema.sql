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
    motif_set_id,
    motif_id,
    motif_name,
    motif_length,
    jaspar_version,
    source_sha256
FROM read_parquet('tables/jaspar2026/motif_metadata/*.parquet');

CREATE OR REPLACE VIEW motif_hit AS
SELECT
    h.task_id,
    h.genome_id,
    h.motif_set_id,
    CAST(h.chrom AS VARCHAR) AS chrom, -- e.g. '1', 'X'  (partition key)
    h.start,                -- 0-based half-open
    h."end",
    h.motif_id,             -- JASPAR accession, e.g. 'MA0861.2' (partition key)
    m.motif_name,           -- e.g. 'TP73' (small dimension, not repeated per hit)
    CASE WHEN h.strand IN ('plus', '+') THEN '+'
         WHEN h.strand IN ('minus', '-') THEN '-'
         ELSE h.strand END AS strand,
    CAST(h.score AS DOUBLE) AS score,
    h.score_mode,           -- log2_relative_risk | log_odds | raw_counts
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
    h.matched_seq           -- nullable; NULL unless scanned with --show-sequence
FROM read_parquet(
    'task_data/task_id=*/tables/jaspar2026/motif_hit/**/*.parquet',
    hive_partitioning = 1,
    union_by_name = 1
) h
JOIN motif_metadata m USING (motif_set_id, motif_id);

-- Immutable run and file inventory. scan_file_inventory is the completeness
-- authority: one row per promoted Parquet file, including zero-hit files.
CREATE OR REPLACE VIEW motif_set AS
SELECT * FROM read_parquet('tables/jaspar2026/motif_set/*.parquet');

CREATE OR REPLACE VIEW genome AS
SELECT * FROM read_parquet('tables/jaspar2026/genome/*.parquet');

CREATE OR REPLACE VIEW sequence_region AS
SELECT * FROM read_parquet('tables/jaspar2026/sequence_region/*.parquet');

CREATE OR REPLACE VIEW scan_run AS
SELECT * FROM read_parquet('tables/jaspar2026/scan_run/*.parquet');

CREATE OR REPLACE VIEW scan_threshold_policy AS
SELECT * FROM read_parquet('tables/jaspar2026/scan_threshold_policy/*.parquet');

CREATE OR REPLACE VIEW scan_task AS
SELECT * FROM read_parquet('tables/jaspar2026/scan_task/*.parquet');

CREATE OR REPLACE VIEW scan_file_inventory AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/scan_file_inventory/**/*.parquet',
    hive_partitioning = 1,
    union_by_name = 1
);

-- Dense score blocks for calibration runs. Each score belongs to one PSSM
-- alignment start, not to an asserted TF footprint. Physical Parquet stores
-- only block_start and a FLOAT[] score vector; skipped/sentinel alignments are
-- NULL list elements. Identity lives in hive partitions: motif_id, score_mode,
-- pseudocount, chrom, strand.
CREATE OR REPLACE VIEW motif_score_dense_block AS
SELECT
    genome_id,
    motif_set_id,
    motif_id,
    score_mode,
    CAST(pseudocount AS DOUBLE) AS pseudocount,
    background_model_id,
    pseudocount_scheme,
    CAST(chrom AS VARCHAR) AS chrom,
    CASE WHEN strand IN ('plus', '+') THEN '+'
         WHEN strand IN ('minus', '-') THEN '-'
         ELSE strand END AS strand,
    block_start,
    scores
FROM read_parquet(
    'tables/jaspar2026/motif_score_dense/**/*.parquet',
    hive_partitioning = 1,
    union_by_name = 1
);

-- Logical row view over dense score blocks. The start/end interval is the DNA
-- span scored by the PSSM alignment. Use this for region inspection; aggregate
-- directly from blocks for whole-chromosome calibration jobs.
CREATE OR REPLACE VIEW motif_score_dense AS
SELECT
    b.genome_id,
    b.motif_set_id,
    b.chrom,
    b.block_start + CAST(u.offset_one - 1 AS BIGINT) AS start,
    b.block_start + CAST(u.offset_one - 1 AS BIGINT) + m.motif_length AS "end",
    b.motif_id,
    m.motif_name,
    b.strand,
    b.score_mode,
    b.pseudocount,
    b.background_model_id,
    b.pseudocount_scheme,
    u.score
FROM motif_score_dense_block b
JOIN motif_metadata m USING (motif_set_id, motif_id)
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
    genome_id,
    motif_set_id,
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
FROM read_parquet(
    'tables/jaspar2026/motif_cutandrun_score_bin_stats/**/*.parquet',
    hive_partitioning = 1,
    union_by_name = 1
);

-- ---------------------------------------------------------------------------
-- Layer 2a: motif-level architecture dimension (small, one row per motif_id)
--   Curated for the TP53 family; 'unknown'/NULL elsewhere. See schema notes.
-- ---------------------------------------------------------------------------
CREATE OR REPLACE VIEW motif_architecture AS
SELECT
    motif_set_id,
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
    genome_id, motif_set_id,
    CAST(chrom AS VARCHAR) AS chrom, start, motif_id, strand,
    spacer_bp,
    half_site_1_score,
    half_site_2_score,
    architecture_class,        -- full_site|half_site|degenerate_full_site
    oligomer_compatible,       -- boolean
    nucleosome_context_model   -- not_assessed|edge_accessible|dyad_centered|orientation_sensitive
FROM read_parquet('tables/jaspar2026/hit_architecture/*/*.parquet', hive_partitioning = 1);

-- ---------------------------------------------------------------------------
-- Layer 2c: local motif syntax around TP73 anchors
--   Schema v4 uses signed interval distance for capture and bands. Center
--   offsets remain directional fields and are not used to decide proximity.
-- ---------------------------------------------------------------------------
CREATE OR REPLACE VIEW motif_context_run_config AS
SELECT
    schema_version,
    genome_id,
    motif_set_id,
    anchor_motif_id,
    score_mode,
    pseudocount,
    background_model_id,
    pseudocount_scheme,
    anchor_minimum_score,
    partner_minimum_score,
    anchor_selection_mode,
    anchor_local_peak_flank_bp,
    anchor_local_peak_rule,
    tandem_score_rule,
    capture_flank_bp,
    context_flank_bp,
    tandem_flank_bp,
    cofactor_pair_flank_bp,
    output_tier,
    capture_geometry,
    distance_metric,
    center_distance_metric,
    capture_prefilter_rule,
    capture_prefilter_center_bp,
    cofactor_pair_prefilter_center_bp,
    observed_max_anchor_span_bp,
    observed_max_neighbor_span_bp,
    distance_band_rule,
    tandem_distance_metric,
    oriented_distance_rule,
    tandem_identity_rule,
    partner_locus_identity_rule,
    pair_class_rule,
    motif_hit_source,
    gtf_source
FROM read_parquet('tables/jaspar2026/context_run_config.parquet');

-- Backward-compatible short name used by the standalone context packages.
CREATE OR REPLACE VIEW context_run_config AS
SELECT * FROM motif_context_run_config;

CREATE OR REPLACE VIEW motif_context_pair AS
SELECT
    anchor_hit_id,
    genome_id,
    motif_set_id,
    CAST(chrom AS VARCHAR) AS chrom,
    anchor_start,
    anchor_end,
    anchor_center_bp,
    anchor_strand,
    anchor_score,
    anchor_pwm_relative_score,
    neighbor_hit_id,
    neighbor_start,
    neighbor_end,
    neighbor_center_bp,
    neighbor_motif_id,
    neighbor_motif_name,
    neighbor_strand,
    neighbor_score,
    neighbor_pwm_relative_score,
    genomic_center_distance_bp,
    anchor_oriented_center_distance_bp,
    absolute_center_distance_bp,
    genomic_side,               -- left | right | coincident_center
    relative_orientation,       -- same | opposite
    anchor_oriented_side,       -- upstream | downstream | coincident_center
    same_alignment_span,
    same_anchor_motif_span,     -- orientation alternative, not another locus
    anchor_neighbor_interval_distance_bp,
    interval_overlap_bp,
    inter_motif_gap_bp,
    interval_relation,          -- containment | overlap | abutting | disjoint
    interval_distance_band,
    within_5,
    within_20,
    within_50,
    within_100,
    within_150,
    within_context_flank,
    is_tandem_tp73,
    score_mode,
    pseudocount,
    background_model_id,
    pseudocount_scheme,
    anchor_minimum_score,
    partner_minimum_score,
    capture_flank_bp,
    context_flank_bp,
    tandem_flank_bp,
    cofactor_pair_flank_bp
FROM read_parquet(
    'tables/jaspar2026/motif_context_pair/**/*.parquet',
    hive_partitioning = 1
);

-- One row per TP73 alignment span. Orientation-specific records remain in
-- motif_context_pair; this locus grain avoids duplicating physical labels.
CREATE OR REPLACE VIEW tp73_anchor_locus AS
SELECT * REPLACE (
    CAST(genome_id AS VARCHAR) AS genome_id,
    CAST(chrom AS VARCHAR) AS chrom
)
FROM read_parquet(
    'tables/jaspar2026/tp73_anchor_locus/**/*.parquet',
    hive_partitioning = 1
);

-- Nonempty long-form groups for all-motif screening. These partitions are the
-- scalable surface; raw relationship rows are the selected-motif tier.
CREATE OR REPLACE VIEW tp73_motif_context_summary AS
SELECT * REPLACE (
    CAST(genome_id AS VARCHAR) AS genome_id,
    CAST(chrom AS VARCHAR) AS chrom,
    CAST(neighbor_motif_id AS VARCHAR) AS neighbor_motif_id
)
FROM read_parquet(
    'tables/jaspar2026/tp73_motif_context_summary/**/*.parquet',
    hive_partitioning = 1
);

-- TP73-context-relevant same-motif cofactor architecture. Identical-span
-- strand alternatives are collapsed before canonical pairs are formed. Loci
-- include context members plus outside partners; pairs require at least one
-- context member, and locus features are emitted only for context members.
CREATE OR REPLACE VIEW cofactor_motif_locus AS
SELECT * FROM read_parquet('tables/jaspar2026/cofactor_motif_locus.parquet');

CREATE OR REPLACE VIEW cofactor_motif_pair AS
SELECT * FROM read_parquet('tables/jaspar2026/cofactor_motif_pair.parquet');

CREATE OR REPLACE VIEW cofactor_locus_pair_feature AS
SELECT *
FROM read_parquet('tables/jaspar2026/cofactor_locus_pair_feature.parquet');

CREATE OR REPLACE VIEW tp73_cofactor_pair_context AS
SELECT *
FROM read_parquet('tables/jaspar2026/tp73_cofactor_pair_context.parquet');

CREATE OR REPLACE VIEW tp73_cofactor_pair_summary AS
SELECT *
FROM read_parquet('tables/jaspar2026/tp73_cofactor_pair_summary.parquet');

-- One row per TP73 anchor. Distinct strand records at the same neighboring
-- alignment span are collapsed into one partner locus; dual-strand support is
-- represented as orientation-ambiguous. These fields describe sequence
-- architecture compatible with a complex, not an observed quaternary state.
CREATE OR REPLACE VIEW tp73_pair_feature AS
SELECT
    anchor_hit_id,
    genome_id,
    motif_set_id,
    CAST(chrom AS VARCHAR) AS chrom,
    start,
    "end",
    center_bp,
    motif_id,
    motif_name,
    strand,
    score,
    anchor_locus_id,
    anchor_selection_class,
    anchor_locus_best_score,
    best_other_anchor_locus_score,
    anchor_locus_score_prominence,
    anchor_locus_is_local_peak,
    score_mode,
    pseudocount,
    background_model_id,
    pseudocount_scheme,
    anchor_minimum_score,
    partner_minimum_score,
    anchor_local_peak_flank_bp,
    pwm_relative_score,
    pair_class,
    n_tandem_tp73_partner_loci,
    n_same_orientation_partner_loci,
    n_opposite_orientation_partner_loci,
    n_ambiguous_orientation_partner_loci,
    has_multiple_tandem_partner_loci,
    nearest_tandem_inter_motif_gap_bp,
    nearest_tandem_absolute_center_distance_bp,
    nearest_same_orientation_gap_bp,
    nearest_opposite_orientation_gap_bp,
    nearest_ambiguous_orientation_gap_bp,
    best_partner_score,
    best_partner_pwm_relative_score,
    best_same_orientation_partner_score,
    best_opposite_orientation_partner_score,
    best_ambiguous_orientation_partner_score,
    best_pair_min_score,
    best_pair_sum_score,
    best_pair_min_pwm_relative_score,
    tandem_flank_bp
FROM read_parquet(
    'tables/jaspar2026/tp73_pair_feature/**/*.parquet',
    hive_partitioning = 1
);

-- Raw neighboring occurrences decorated with the anchor's pair class. This is
-- the per-TP73-site feature surface used by CUT&RUN models; unlike the promoter
-- summaries below, it covers every retained anchor on the requested genome.
CREATE OR REPLACE VIEW tp73_context_pair_feature AS
SELECT
    p.*,
    f.motif_id AS anchor_motif_id,
    f.pair_class AS anchor_pair_class,
    f.n_tandem_tp73_partner_loci,
    f.n_same_orientation_partner_loci,
    f.n_opposite_orientation_partner_loci,
    f.n_ambiguous_orientation_partner_loci,
    CASE
        WHEN p.is_tandem_tp73 THEN 'same_motif_tandem'
        WHEN p.neighbor_motif_id = f.motif_id AND p.interval_overlap_bp > 0
            THEN 'same_motif_overlapping_alignment'
        WHEN p.neighbor_motif_id = f.motif_id THEN 'same_motif_context'
        ELSE 'heterotypic_context'
    END AS pair_relation,
    FLOOR(p.anchor_oriented_center_distance_bp / 5.0) * 5.0
        AS oriented_distance_bin_start_bp
FROM motif_context_pair p
JOIN tp73_pair_feature f USING (anchor_hit_id);

-- Published-work compatibility only. The historical context used motif start
-- coordinates and a 100 bp cutoff; new analyses use interval geometry.
CREATE OR REPLACE VIEW legacy_tp73_context_100 AS
SELECT
    p.*,
    p.neighbor_start - p.anchor_start AS legacy_genomic_start_distance_bp,
    CASE WHEN p.anchor_strand = '-'
         THEN p.anchor_start - p.neighbor_start
         ELSE p.neighbor_start - p.anchor_start END
        AS legacy_anchor_oriented_start_distance_bp
FROM motif_context_pair p
WHERE ABS(p.neighbor_start - p.anchor_start) <= 100;

CREATE OR REPLACE VIEW tp73_context_anchor AS
SELECT
    anchor_hit_id,
    genome_id,
    motif_set_id,
    CAST(chrom AS VARCHAR) AS chrom,
    start,
    "end",
    center_bp,
    motif_id,
    motif_name,
    strand,
    score,
    anchor_locus_id,
    anchor_selection_class,
    anchor_locus_best_score,
    best_other_anchor_locus_score,
    anchor_locus_score_prominence,
    anchor_locus_is_local_peak,
    score_mode,
    pseudocount,
    background_model_id,
    pseudocount_scheme,
    anchor_minimum_score,
    partner_minimum_score,
    anchor_local_peak_flank_bp,
    pwm_relative_score,
    n_context_neighbor_hits,
    n_context_neighbor_loci,
    n_context_motifs,
    has_tandem_tp73,
    n_tandem_tp73_partners,
    pair_class,
    n_tandem_tp73_partner_loci,
    n_same_orientation_partner_loci,
    n_opposite_orientation_partner_loci,
    n_ambiguous_orientation_partner_loci,
    has_multiple_tandem_partner_loci,
    nearest_tandem_hit_id,
    nearest_tandem_oriented_distance_bp,
    nearest_tandem_genomic_distance_bp,
    nearest_tandem_absolute_distance_bp,
    nearest_tandem_relative_orientation,
    nearest_tandem_score,
    gene_annotation_available,
    nearest_gene_id,
    nearest_gene_name,
    nearest_transcript_id,
    nearest_tss_distance_bp,
    in_any_transcript,
    in_any_intron,
    overlaps_any_intron_boundary,
    primary_transcript_region,
    capture_flank_bp,
    context_flank_bp,
    tandem_flank_bp
FROM read_parquet(
    'tables/jaspar2026/tp73_context_anchor/**/*.parquet',
    hive_partitioning = 1
);

-- ---------------------------------------------------------------------------
-- Layer 3: experimental evidence (CUT&RUN), long form
--   One row per (locus, sample). Sample decomposed into its facets so TA/DN
--   contrasts are a GROUP BY, not a column-name parse.
-- ---------------------------------------------------------------------------
CREATE OR REPLACE VIEW cutandrun AS
SELECT
    genome_id,
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
SELECT genome_id, gene_id, gene_name, chrom, strand, gene_start, gene_end, tss, biotype
FROM read_parquet('tables/jaspar2026/gene/*.parquet');

CREATE OR REPLACE VIEW promoter AS
SELECT genome_id, gene_id, transcript_id, chrom, strand, promoter_start, promoter_end, tss
FROM read_parquet('tables/jaspar2026/promoter/*.parquet');

-- Transcript and intron intervals are regenerated directly from the pinned
-- GTF. A motif can be intronic for one transcript and non-intronic for another;
-- keep that ambiguity in the long per-transcript bridge rather than flattening
-- it into one premature gene-level label.
CREATE OR REPLACE VIEW transcript AS
SELECT genome_id, gene_id, transcript_id, chrom, strand, transcript_start, transcript_end,
       tss, gene_name, biotype
FROM read_parquet('tables/jaspar2026/transcript.parquet');

CREATE OR REPLACE VIEW intron AS
SELECT genome_id, gene_id, transcript_id, chrom, strand, start, "end", intron_number
FROM read_parquet('tables/jaspar2026/intron.parquet');

CREATE OR REPLACE VIEW motif_transcript_context AS
SELECT
    anchor_hit_id,
    genome_id,
    motif_set_id,
    chrom,
    anchor_start,
    anchor_end,
    gene_id,
    gene_name,
    transcript_id,
    transcript_strand,
    tss,
    signed_tss_distance_bp,
    fully_within_intron,
    overlaps_intron,
    transcript_region
FROM read_parquet('tables/jaspar2026/motif_transcript_context.parquet');

-- Transcript-oriented TP73/cofactor direction is derived here rather than
-- multiplying the chromosome-wide feature table once per transcript.
CREATE OR REPLACE VIEW motif_transcript_context_pair AS
SELECT
    p.*,
    t.gene_id,
    t.gene_name,
    t.transcript_id,
    t.transcript_strand,
    t.tss,
    t.signed_tss_distance_bp AS anchor_signed_tss_distance_bp,
    CASE WHEN t.transcript_strand = '+' THEN p.neighbor_center_bp - t.tss
         ELSE t.tss - p.neighbor_center_bp END
        AS neighbor_signed_tss_distance_bp,
    CASE WHEN t.transcript_strand = '+' THEN p.genomic_center_distance_bp
         ELSE -p.genomic_center_distance_bp END
        AS transcript_oriented_center_distance_bp,
    CASE
        WHEN (CASE WHEN t.transcript_strand = '+'
                        THEN p.genomic_center_distance_bp
                   ELSE -p.genomic_center_distance_bp END) < 0 THEN 'upstream'
        WHEN (CASE WHEN t.transcript_strand = '+'
                        THEN p.genomic_center_distance_bp
                   ELSE -p.genomic_center_distance_bp END) > 0 THEN 'downstream'
        ELSE 'coincident_center'
    END AS transcript_oriented_side,
    p.anchor_start <= t.tss AND p.anchor_end > t.tss AS anchor_spans_tss,
    p.neighbor_start <= t.tss AND p.neighbor_end > t.tss AS neighbor_spans_tss
FROM motif_context_pair p
JOIN motif_transcript_context t USING (anchor_hit_id);

CREATE OR REPLACE VIEW gene_set AS
SELECT genome_id, gene_id, set_name -- many-to-many (HALLMARK_*, KEGG_*, ...)
FROM read_parquet('tables/jaspar2026/gene_set/*.parquet');

-- Expression label source. `expression` is tidy per (gene, cell_line, isoform);
-- `expression_differential` is the precomputed TA-vs-DN contrast = the ML label.
CREATE OR REPLACE VIEW expression AS
SELECT genome_id, gene_id, cell_line, isoform, value, replicate
FROM read_parquet('tables/jaspar2026/expression/*.parquet');

CREATE OR REPLACE VIEW expression_differential AS
SELECT
    genome_id, gene_id, cell_line,
    log2fc_ta_vs_dn,           -- + means higher under TAp73 than DNp73
    lfc_se, pvalue, padj, base_mean
FROM read_parquet('tables/jaspar2026/expression_differential/*.parquet');

-- ===========================================================================
-- Materialized joins and summaries that we never want to recompute per query.
-- ===========================================================================

-- Bridge: which motif hits fall inside which gene's promoter set, with
-- strand-aware TSS distance. A hit can overlap several transcript promoters
-- for one gene, but appears here only once per gene and scoring configuration.
-- The closest transcript supplies the signed distance; the overlap count keeps
-- the discarded transcript multiplicity visible without multiplying features.
-- This is the one interval join in the whole system; do it once.
CREATE OR REPLACE TABLE promoter_motif_hit AS
WITH transcript_overlap AS (
    SELECT
        p.genome_id,
        h.motif_set_id,
        p.gene_id,
        p.transcript_id,
        h.chrom,
        h.motif_id,
        h.strand,
        h.score,
        h.score_mode,
        h.pseudocount,
        h.pwm_relative_score,
        h.start AS hit_start,
        h."end" AS hit_end,
        -- signed distance from TSS, positive = downstream in gene direction
        CASE WHEN p.strand = '+' THEN h.start - p.tss
             ELSE p.tss - h."end" END AS tss_distance
    FROM promoter p
    JOIN motif_hit h
      ON  h.genome_id = p.genome_id
      AND h.chrom = p.chrom
      AND h.start >= p.promoter_start
      AND h."end" <= p.promoter_end
),
ranked_overlap AS (
    SELECT
        *,
        COUNT(DISTINCT transcript_id) OVER hit_identity
            AS n_overlapping_transcripts,
        ROW_NUMBER() OVER (
            hit_identity
            ORDER BY ABS(tss_distance), transcript_id
        ) AS transcript_distance_rank
    FROM transcript_overlap
    WINDOW hit_identity AS (
        PARTITION BY genome_id, motif_set_id, gene_id, chrom, hit_start,
                     hit_end, motif_id, strand,
                     score_mode, pseudocount, score, pwm_relative_score
    )
)
SELECT
    genome_id,
    motif_set_id,
    gene_id,
    transcript_id AS closest_transcript_id,
    n_overlapping_transcripts,
    chrom,
    motif_id,
    strand,
    score,
    score_mode,
    pseudocount,
    pwm_relative_score,
    hit_start AS start,
    hit_end AS "end",
    tss_distance
FROM ranked_overlap
WHERE transcript_distance_rank = 1;

-- Per-(gene, motif) architecture summary in LONG form (no ~8000-col matrix).
-- Pivot to wide only at the ML boundary (see queries.sql / the Python builder).
CREATE OR REPLACE TABLE promoter_arch_feature AS
SELECT
    b.genome_id,
    b.motif_set_id,
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
LEFT JOIN motif_architecture a USING (motif_set_id, motif_id)
GROUP BY b.genome_id, b.motif_set_id, b.gene_id, b.motif_id,
         b.score_mode, b.pseudocount, a.tf_family, a.binding_unit_model;

-- Pair state of promoter-associated TP73 anchors. Every TP73 occurrence is
-- represented, including singletons. Keeping pair_class in the table grain
-- lets an ML exporter compare singleton, same-, opposite-, mixed-, and
-- orientation-ambiguous sequence architectures without duplicating labels.
CREATE OR REPLACE TABLE promoter_pair_feature AS
SELECT
    b.genome_id,
    b.motif_set_id,
    b.gene_id,
    f.motif_id AS anchor_motif_id,
    f.score_mode,
    f.pseudocount,
    f.pair_class,
    COUNT(*)::BIGINT AS n_anchor_orientation_hits,
    COUNT(DISTINCT concat_ws('|', f.chrom, f.start::VARCHAR, f."end"::VARCHAR))::BIGINT
        AS n_anchor_loci,
    SUM(CASE WHEN f.pwm_relative_score >= 0.8 THEN 1 ELSE 0 END)::BIGINT
        AS n_strong_anchor_hits,
    MAX(f.score) AS max_anchor_score,
    MAX(f.pwm_relative_score) AS max_anchor_pwm_relative_score,
    SUM(f.n_tandem_tp73_partner_loci)::BIGINT AS n_tandem_partner_loci,
    MAX(f.n_tandem_tp73_partner_loci)::BIGINT AS max_partners_per_anchor,
    MIN(f.nearest_tandem_inter_motif_gap_bp)::BIGINT
        AS nearest_tandem_inter_motif_gap_bp,
    MAX(f.best_partner_score) AS best_partner_score,
    MAX(f.best_partner_pwm_relative_score) AS best_partner_pwm_relative_score,
    MAX(f.best_pair_min_score) AS best_pair_min_score,
    MAX(f.best_pair_sum_score) AS best_pair_sum_score,
    MAX(f.best_pair_min_pwm_relative_score) AS best_pair_min_pwm_relative_score
FROM promoter_motif_hit b
JOIN tp73_pair_feature f
  ON  f.genome_id = b.genome_id
  AND f.motif_set_id = b.motif_set_id
  AND f.chrom = b.chrom
  AND f.start = b.start
  AND f."end" = b."end"
  AND f.motif_id = b.motif_id
  AND f.strand = b.strand
  AND f.score IS NOT DISTINCT FROM b.score
  AND f.score_mode = b.score_mode
  AND f.pseudocount IS NOT DISTINCT FROM b.pseudocount
GROUP BY b.genome_id, b.motif_set_id, b.gene_id, f.motif_id,
         f.score_mode, f.pseudocount, f.pair_class;

-- Directed TP73-anchor/neighbor features for promoter models. The raw pair
-- geometry remains in motif_context_pair; this table only aggregates the
-- provisional context radius into stable 5-bp oriented-distance bins. Pair
-- scores from different motif models remain separate; only PWM-relative
-- scores are combined through a cross-motif minimum.
CREATE OR REPLACE TABLE promoter_motif_pair_feature AS
WITH promoter_context_pair AS (
    SELECT
        b.genome_id,
        b.motif_set_id,
        b.gene_id,
        f.motif_id AS anchor_motif_id,
        p.neighbor_motif_id,
        p.neighbor_motif_name,
        f.pair_class AS anchor_pair_class,
        p.relative_orientation,
        p.anchor_oriented_side,
        CASE
            WHEN p.is_tandem_tp73 THEN 'same_motif_tandem'
            WHEN p.neighbor_motif_id = f.motif_id AND p.interval_overlap_bp > 0
                THEN 'same_motif_overlapping_alignment'
            WHEN p.neighbor_motif_id = f.motif_id THEN 'same_motif_context'
            ELSE 'heterotypic_context'
        END AS pair_relation,
        FLOOR(p.anchor_oriented_center_distance_bp / 5.0) * 5.0
            AS oriented_distance_bin_start_bp,
        f.anchor_hit_id,
        p.neighbor_hit_id,
        p.neighbor_start,
        p.neighbor_end,
        p.anchor_score,
        p.anchor_pwm_relative_score,
        p.neighbor_score,
        p.neighbor_pwm_relative_score,
        p.absolute_center_distance_bp,
        p.score_mode,
        p.pseudocount
    FROM promoter_motif_hit b
    JOIN tp73_pair_feature f
      ON  f.genome_id = b.genome_id
      AND f.motif_set_id = b.motif_set_id
      AND f.chrom = b.chrom
      AND f.start = b.start
      AND f."end" = b."end"
      AND f.motif_id = b.motif_id
      AND f.strand = b.strand
      AND f.score IS NOT DISTINCT FROM b.score
      AND f.score_mode = b.score_mode
      AND f.pseudocount IS NOT DISTINCT FROM b.pseudocount
    JOIN motif_context_pair p
      ON p.anchor_hit_id = f.anchor_hit_id
     AND p.genome_id = f.genome_id
     AND p.motif_set_id = f.motif_set_id
    WHERE p.within_context_flank
      AND NOT p.same_anchor_motif_span
)
SELECT
    genome_id,
    motif_set_id,
    gene_id,
    anchor_motif_id,
    neighbor_motif_id,
    MAX(neighbor_motif_name) AS neighbor_motif_name,
    score_mode,
    pseudocount,
    anchor_pair_class,
    pair_relation,
    relative_orientation,
    anchor_oriented_side,
    oriented_distance_bin_start_bp,
    oriented_distance_bin_start_bp + 5.0 AS oriented_distance_bin_end_bp,
    5::INTEGER AS distance_bin_width_bp,
    COUNT(*)::BIGINT AS n_directed_orientation_pairs,
    COUNT(DISTINCT concat_ws(
        '|', anchor_hit_id, neighbor_motif_id,
        neighbor_start::VARCHAR, neighbor_end::VARCHAR
    ))::BIGINT AS n_anchor_neighbor_loci,
    COUNT(DISTINCT anchor_hit_id)::BIGINT AS n_anchor_hits,
    MAX(anchor_score) AS max_anchor_score,
    MAX(anchor_pwm_relative_score) AS max_anchor_pwm_relative_score,
    MAX(neighbor_score) AS max_neighbor_score,
    MAX(neighbor_pwm_relative_score) AS max_neighbor_pwm_relative_score,
    MAX(LEAST(anchor_pwm_relative_score, neighbor_pwm_relative_score))
        AS max_pair_min_pwm_relative_score,
    MIN(absolute_center_distance_bp) AS min_absolute_center_distance_bp
FROM promoter_context_pair
GROUP BY genome_id, motif_set_id, gene_id, anchor_motif_id,
         neighbor_motif_id, score_mode, pseudocount,
         anchor_pair_class, pair_relation, relative_orientation,
         anchor_oriented_side, oriented_distance_bin_start_bp;

-- The ML surface: features (long) + the TA-vs-DN label. Keep long; pivot at read.
CREATE OR REPLACE VIEW ml_ta_vs_dn AS
SELECT
    f.genome_id,
    f.motif_set_id,
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
JOIN expression_differential d USING (genome_id, gene_id);

-- Pair-aware sequence-only ML surfaces. They remain separate from the ordinary
-- per-motif frame so callers opt into pair stratification deliberately.
CREATE OR REPLACE VIEW ml_ta_vs_dn_pair AS
SELECT
    f.*,
    d.cell_line,
    d.log2fc_ta_vs_dn,
    d.padj,
    d.log2fc_ta_vs_dn > 0 AS label_ta_up
FROM promoter_pair_feature f
JOIN expression_differential d USING (genome_id, gene_id);

CREATE OR REPLACE VIEW ml_ta_vs_dn_motif_pair AS
SELECT
    f.*,
    d.cell_line,
    d.log2fc_ta_vs_dn,
    d.padj,
    d.log2fc_ta_vs_dn > 0 AS label_ta_up
FROM promoter_motif_pair_feature f
JOIN expression_differential d USING (genome_id, gene_id);

-- Hot agent lookup: one compact, low-latency promoter summary per gene and
-- scoring configuration. Keeping score_mode and pseudocount in the table grain
-- prevents totals from unrelated scoring runs being combined under MIN labels.
-- This is the first stop for OpenClaw / Claude / Codex before deeper joins.
CREATE OR REPLACE TABLE promoter_card AS
WITH promoter_span AS (
    SELECT
        genome_id,
        gene_id,
        MIN(promoter_start) AS promoter_start,
        MAX(promoter_end) AS promoter_end,
        COUNT(DISTINCT transcript_id) AS n_promoter_transcripts
    FROM promoter
    GROUP BY genome_id, gene_id
),
score_configuration AS (
    SELECT DISTINCT genome_id, motif_set_id, score_mode, pseudocount
    FROM motif_hit
),
feature_summary AS (
    SELECT
        genome_id,
        motif_set_id,
        gene_id,
        score_mode,
        pseudocount,
        COUNT(DISTINCT motif_id) AS n_motifs_with_hits,
        SUM(n_hits) AS n_motif_hits,
        MAX(max_rel_score) AS max_pwm_relative_score,
        MIN(min_abs_tss_distance) AS min_abs_tss_distance,
        BOOL_OR(COALESCE(has_dimer_of_dimers, false)) AS has_dimer_of_dimers,
        COUNT(DISTINCT tf_family) FILTER (WHERE tf_family IS NOT NULL) AS n_tf_families
    FROM promoter_arch_feature
    GROUP BY genome_id, motif_set_id, gene_id, score_mode, pseudocount
),
pair_summary AS (
    SELECT
        genome_id,
        motif_set_id,
        gene_id,
        score_mode,
        pseudocount,
        SUM(n_anchor_orientation_hits)::BIGINT AS n_tp73_anchor_hits,
        SUM(n_anchor_orientation_hits) FILTER (WHERE pair_class = 'singleton')::BIGINT
            AS n_singleton_tp73_anchor_hits,
        SUM(n_anchor_orientation_hits) FILTER (WHERE pair_class <> 'singleton')::BIGINT
            AS n_tandem_tp73_anchor_hits,
        COUNT(DISTINCT pair_class)::BIGINT AS n_tp73_pair_classes
    FROM promoter_pair_feature
    GROUP BY genome_id, motif_set_id, gene_id, score_mode, pseudocount
)
SELECT
    g.genome_id,
    sc.motif_set_id,
    g.gene_id,
    g.gene_name,
    g.chrom,
    g.strand,
    g.tss,
    ps.promoter_start,
    ps.promoter_end,
    sc.score_mode,
    sc.pseudocount,
    COALESCE(ps.n_promoter_transcripts, 0) AS n_promoter_transcripts,
    COALESCE(fs.n_motifs_with_hits, 0) AS n_motifs_with_hits,
    COALESCE(fs.n_motif_hits, 0) AS n_motif_hits,
    fs.max_pwm_relative_score,
    fs.min_abs_tss_distance,
    COALESCE(fs.has_dimer_of_dimers, false) AS has_dimer_of_dimers,
    COALESCE(fs.n_tf_families, 0) AS n_tf_families,
    COALESCE(pa.n_tp73_anchor_hits, 0) AS n_tp73_anchor_hits,
    COALESCE(pa.n_singleton_tp73_anchor_hits, 0) AS n_singleton_tp73_anchor_hits,
    COALESCE(pa.n_tandem_tp73_anchor_hits, 0) AS n_tandem_tp73_anchor_hits,
    COALESCE(pa.n_tp73_pair_classes, 0) AS n_tp73_pair_classes
FROM gene g
JOIN score_configuration sc USING (genome_id)
LEFT JOIN promoter_span ps USING (genome_id, gene_id)
LEFT JOIN feature_summary fs
  ON  fs.genome_id = g.genome_id
  AND fs.motif_set_id = sc.motif_set_id
  AND fs.gene_id = g.gene_id
  AND fs.score_mode = sc.score_mode
  AND fs.pseudocount IS NOT DISTINCT FROM sc.pseudocount
LEFT JOIN pair_summary pa
  ON  pa.genome_id = g.genome_id
  AND pa.motif_set_id = sc.motif_set_id
  AND pa.gene_id = g.gene_id
  AND pa.score_mode = sc.score_mode
  AND pa.pseudocount IS NOT DISTINCT FROM sc.pseudocount;

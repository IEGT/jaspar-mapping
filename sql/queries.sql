-- ============================================================================
-- queries.sql — canned, low-latency queries over schema.sql
-- ----------------------------------------------------------------------------
-- These are the queries agents (OpenClaw / Claude / Codex) and humans reach for.
-- Each is a single statement, parameterized with $placeholders so it can be
-- prepared once and reused. Run `.read sql/schema.sql` first.
-- ============================================================================

-- Q0. Hot promoter card — first lookup for an agent answering "what do we know
--     about this gene's promoter architecture?"
--     Params: $genome_id, $motif_set_id, $gene_name, $score_mode, $pseudocount
SELECT *
FROM promoter_card
WHERE genome_id = $genome_id
  AND motif_set_id = $motif_set_id
  AND gene_name = $gene_name
  AND score_mode = $score_mode
  AND pseudocount = $pseudocount;

-- Q1. Region browse — every motif hit in a window. Prunes by chrom partition
--     + start zonemap, so latency is ~row-group scan, not full table.
--     Params: $genome_id, $motif_set_id, $chrom, $start, $end
SELECT genome_id, motif_set_id, chrom, start, "end", motif_id, motif_name,
       strand, score_mode, pseudocount, minimum_score, pwm_relative_score
FROM motif_hit
WHERE genome_id = $genome_id AND motif_set_id = $motif_set_id
  AND chrom = $chrom AND start >= $start AND "end" <= $end
ORDER BY start;

-- Q2. Promoter architecture for one gene — what the ML sees, human-readable.
--     Params: $genome_id, $motif_set_id, $gene_name, $score_mode, $pseudocount
SELECT f.motif_id, f.score_mode, f.pseudocount, f.tf_family, f.binding_unit_model,
       f.n_hits, f.n_strong_hits, f.max_rel_score, f.min_abs_tss_distance
FROM promoter_arch_feature f
JOIN gene g ON g.genome_id = f.genome_id AND g.gene_id = f.gene_id
WHERE f.genome_id = $genome_id
  AND f.motif_set_id = $motif_set_id
  AND g.gene_name = $gene_name
  AND f.score_mode = $score_mode
  AND f.pseudocount = $pseudocount
ORDER BY f.n_strong_hits DESC, f.max_rel_score DESC;

-- Q3. Prepared-safe ML pull in LONG family-feature form. DuckDB must discover
--     data-driven PIVOT columns while binding a statement, so a dynamic PIVOT
--     cannot contain parameters in its source. Bind this query first, then
--     pivot its result in the ML client/exporter.
--     Params: $genome_id, $motif_set_id, $cell_line, $score_mode, $pseudocount
SELECT
    m.genome_id,
    m.motif_set_id,
    m.gene_id,
    m.cell_line,
    m.score_mode,
    m.pseudocount,
    m.tf_family,
    SUM(m.n_hits) AS n_hits,
    SUM(m.n_strong_hits) AS n_strong_hits,
    MAX(m.max_rel_score) AS max_rel_score,
    MIN(m.min_abs_tss_distance) AS min_abs_tss_distance,
    BOOL_OR(COALESCE(m.has_dimer_of_dimers, false)) AS has_dimer_of_dimers,
    m.log2fc_ta_vs_dn,
    m.padj,
    m.label_ta_up
FROM ml_ta_vs_dn m
WHERE m.genome_id = $genome_id
  AND m.motif_set_id = $motif_set_id
  AND m.cell_line = $cell_line
  AND m.score_mode = $score_mode
  AND m.pseudocount = $pseudocount
  AND m.tf_family IS NOT NULL
GROUP BY m.genome_id, m.motif_set_id, m.gene_id, m.cell_line,
         m.score_mode, m.pseudocount, m.tf_family,
         m.log2fc_ta_vs_dn, m.padj, m.label_ta_up
ORDER BY m.gene_id, m.tf_family;

-- Q4. Experimental cross-check — TA-vs-DN differential CUT&RUN signal at a
--     gene's promoter loci, kept separate by antibody. Keep this OUT of the
--     sequence-only feature set (leakage); use it to validate predictions.
--     Params: $genome_id, $gene_name, $cell_line
SELECT c.genome_id, c.chrom, c.start, c."end", c.cell_line, c.antibody,
       AVG(c.signal) FILTER (WHERE c.isoform = 'TA') AS ta_signal,
       AVG(c.signal) FILTER (WHERE c.isoform = 'DN') AS dn_signal,
       AVG(c.signal) FILTER (WHERE c.isoform = 'TA')
     - AVG(c.signal) FILTER (WHERE c.isoform = 'DN') AS ta_minus_dn
FROM cutandrun c
JOIN promoter p
  ON c.genome_id = p.genome_id AND c.chrom = p.chrom
 AND c.start < p.promoter_end AND c."end" > p.promoter_start
JOIN gene g ON g.genome_id = p.genome_id AND g.gene_id = p.gene_id
WHERE c.genome_id = $genome_id
  AND g.gene_name = $gene_name AND c.cell_line = $cell_line
GROUP BY c.genome_id, c.chrom, c.start, c."end", c.cell_line, c.antibody
ORDER BY ABS(ta_minus_dn) DESC;

-- Q5. Architecture enrichment — do dimer-of-dimers TP53-family sites associate
--     with TA-up genes? A quick 2x2 signal before any modelling.
--     Params: $genome_id, $motif_set_id, $cell_line, $score_mode, $pseudocount
SELECT has_arch,
       AVG(label_ta_up)          AS frac_ta_up,
       COUNT(*)                  AS n_genes
FROM (
    SELECT m.gene_id,
           MAX(CASE WHEN m.has_dimer_of_dimers THEN 1 ELSE 0 END) AS has_arch,
           MAX(m.label_ta_up)                                     AS label_ta_up
    FROM ml_ta_vs_dn m
    WHERE m.genome_id = $genome_id
      AND m.motif_set_id = $motif_set_id
      AND m.cell_line = $cell_line
      AND m.score_mode = $score_mode
      AND m.pseudocount = $pseudocount
    GROUP BY m.gene_id
)
GROUP BY has_arch;

-- Q6. Gene-set slice — restrict the ML frame to a pathway (e.g. an EMT set) so a
--     model can be trained per program.
--     Params: $genome_id, $motif_set_id, $set_name, $cell_line, $score_mode, $pseudocount
SELECT m.*
FROM ml_ta_vs_dn m
JOIN gene_set s ON s.genome_id = m.genome_id AND s.gene_id = m.gene_id
WHERE m.genome_id = $genome_id
  AND m.motif_set_id = $motif_set_id
  AND s.set_name = $set_name
  AND m.cell_line = $cell_line
  AND m.score_mode = $score_mode
  AND m.pseudocount = $pseudocount;

-- Q7. Dense score slice — reconstruct PSSM alignment-score rows from dense
--     score blocks for a small genomic interval. Params: $genome_id,
--     $motif_set_id, $motif_id, $chrom, $strand,
--     $score_mode, $pseudocount, $start, $end
SELECT chrom, start, "end", motif_id, motif_name, strand,
       score_mode, pseudocount, score
FROM motif_score_dense
WHERE genome_id = $genome_id
  AND motif_set_id = $motif_set_id
  AND motif_id = $motif_id
  AND chrom = $chrom
  AND strand = $strand
  AND score_mode = $score_mode
  AND pseudocount = $pseudocount
  AND start >= $start
  AND start < $end
ORDER BY start;

-- Q8. Dense calibration result — compare CUT&RUN support by score bin for one
--     motif/scoring run. Params: $genome_id, $motif_set_id, $motif_id, $chrom,
--     $score_mode, $pseudocount
SELECT motif_id, motif_name, chrom, strand, score_mode, pseudocount,
       bin_order, bin_label, lower_bound, upper_bound,
       cell_line, isoform, antibody, replicate,
       n_windows, n_covered_windows, overlap_fraction, baseline_fraction,
       enrichment_ratio, log2_enrichment, mean_signal, max_signal
FROM motif_cutandrun_score_bin_stats
WHERE genome_id = $genome_id
  AND motif_set_id = $motif_set_id
  AND motif_id = $motif_id
  AND chrom = $chrom
  AND score_mode = $score_mode
  AND pseudocount = $pseudocount
ORDER BY cell_line, antibody, isoform, replicate, strand, bin_order;

-- Q9. Calibration-bin ladder — the common bins used for log2-relative-risk and
--     log-odds dense chr1 comparisons.
SELECT *
FROM score_calibration_bin
ORDER BY bin_order;

-- Q10. TP73 tandem partners for one anchor. An exact same-span orientation
--      alternative is retained in motif_context_pair but is not a tandem site.
--      Params: $genome_id, $motif_set_id, $chrom, $start, $strand,
--      $score_mode, $pseudocount
SELECT
    anchor_hit_id,
    neighbor_hit_id,
    neighbor_start,
    neighbor_end,
    neighbor_strand,
    neighbor_score,
    anchor_neighbor_interval_distance_bp,
    interval_relation,
    interval_distance_band,
    genomic_center_distance_bp,
    anchor_oriented_center_distance_bp,
    relative_orientation,
    anchor_oriented_side
FROM motif_context_pair
WHERE genome_id = $genome_id
  AND motif_set_id = $motif_set_id
  AND chrom = $chrom
  AND anchor_start = $start
  AND anchor_strand = $strand
  AND score_mode = $score_mode
  AND pseudocount = $pseudocount
  AND is_tandem_tp73
ORDER BY absolute_center_distance_bp, neighbor_start, neighbor_strand;

-- Q11. Compact TP73 context for anchors nearest a named gene. This is the
--      sequence/transcript feature surface; CUT&RUN remains a separate layer.
--      Params: $genome_id, $motif_set_id, $gene_name, $score_mode, $pseudocount
SELECT
    chrom,
    start,
    "end",
    strand,
    score,
    pwm_relative_score,
    anchor_selection_class,
    anchor_locus_best_score,
    best_other_anchor_locus_score,
    anchor_locus_score_prominence,
    anchor_locus_is_local_peak,
    nearest_tss_distance_bp,
    primary_transcript_region,
    in_any_intron,
    has_tandem_tp73,
    n_tandem_tp73_partners,
    pair_class,
    n_tandem_tp73_partner_loci,
    n_same_orientation_partner_loci,
    n_opposite_orientation_partner_loci,
    n_ambiguous_orientation_partner_loci,
    nearest_tandem_oriented_distance_bp,
    nearest_tandem_relative_orientation,
    n_context_neighbor_loci,
    n_context_motifs
FROM tp73_context_anchor
WHERE genome_id = $genome_id
  AND motif_set_id = $motif_set_id
  AND nearest_gene_name = $gene_name
  AND score_mode = $score_mode
  AND pseudocount = $pseudocount
ORDER BY ABS(nearest_tss_distance_bp), start, strand;

-- Q12. Pair-stratified TP73 promoter features for the expression model. This
--      includes singleton anchors; pair classes describe sequence-compatible
--      architecture, not an observed protein complex.
--      Params: $genome_id, $motif_set_id, $gene_name, $score_mode, $pseudocount
SELECT
    f.anchor_motif_id,
    f.pair_class,
    f.n_anchor_orientation_hits,
    f.n_anchor_loci,
    f.n_strong_anchor_hits,
    f.max_anchor_score,
    f.max_anchor_pwm_relative_score,
    f.n_tandem_partner_loci,
    f.max_partners_per_anchor,
    f.nearest_tandem_inter_motif_gap_bp,
    f.best_partner_score,
    f.best_partner_pwm_relative_score,
    f.best_pair_min_score,
    f.best_pair_sum_score,
    f.best_pair_min_pwm_relative_score
FROM promoter_pair_feature f
JOIN gene g USING (genome_id, gene_id)
WHERE f.genome_id = $genome_id
  AND f.motif_set_id = $motif_set_id
  AND g.gene_name = $gene_name
  AND f.score_mode = $score_mode
  AND f.pseudocount = $pseudocount
ORDER BY f.pair_class;

-- Q13. Neighbor-motif predictors stratified by the TP73 anchor's pair class,
--      relative orientation, side, and 5-bp oriented-distance bin. This is the
--      long feature frame for interactions such as PATZ1 x TP73-pair-class.
--      Params: $genome_id, $motif_set_id, $gene_name, $neighbor_motif_id,
--      $score_mode, $pseudocount
SELECT
    f.anchor_motif_id,
    f.neighbor_motif_id,
    f.neighbor_motif_name,
    f.anchor_pair_class,
    f.pair_relation,
    f.relative_orientation,
    f.anchor_oriented_side,
    f.oriented_distance_bin_start_bp,
    f.oriented_distance_bin_end_bp,
    f.n_directed_orientation_pairs,
    f.n_anchor_neighbor_loci,
    f.n_anchor_hits,
    f.max_anchor_pwm_relative_score,
    f.max_neighbor_pwm_relative_score,
    f.max_pair_min_pwm_relative_score,
    f.min_absolute_center_distance_bp
FROM promoter_motif_pair_feature f
JOIN gene g USING (genome_id, gene_id)
WHERE f.genome_id = $genome_id
  AND f.motif_set_id = $motif_set_id
  AND g.gene_name = $gene_name
  AND f.neighbor_motif_id = $neighbor_motif_id
  AND f.score_mode = $score_mode
  AND f.pseudocount = $pseudocount
ORDER BY f.anchor_pair_class, f.oriented_distance_bin_start_bp,
         f.relative_orientation;

-- Q14. Chromosome-wide, per-anchor neighboring motif records decorated with
--      TP73 pair class. Use this surface to stratify a CUT&RUN predictor before
--      any promoter restriction. Params: $genome_id, $motif_set_id, $chrom,
--      $start, $strand,
--      $neighbor_motif_id, $score_mode, $pseudocount
SELECT
    anchor_hit_id,
    anchor_motif_id,
    anchor_start,
    anchor_end,
    anchor_strand,
    anchor_score,
    anchor_pwm_relative_score,
    anchor_pair_class,
    n_tandem_tp73_partner_loci,
    neighbor_hit_id,
    neighbor_motif_id,
    neighbor_motif_name,
    neighbor_start,
    neighbor_end,
    neighbor_strand,
    neighbor_score,
    neighbor_pwm_relative_score,
    pair_relation,
    anchor_neighbor_interval_distance_bp,
    interval_relation,
    interval_distance_band,
    relative_orientation,
    anchor_oriented_side,
    anchor_oriented_center_distance_bp,
    oriented_distance_bin_start_bp
FROM tp73_context_pair_feature
WHERE genome_id = $genome_id
  AND motif_set_id = $motif_set_id
  AND chrom = $chrom
  AND anchor_start = $start
  AND anchor_strand = $strand
  AND neighbor_motif_id = $neighbor_motif_id
  AND score_mode = $score_mode
  AND pseudocount = $pseudocount
  AND within_context_flank
  AND NOT same_anchor_motif_span
ORDER BY neighbor_score DESC, absolute_center_distance_bp,
         neighbor_start, neighbor_strand;

-- Q15. Canonical same-motif cofactor pairs attributed to one TP73 anchor.
--      Each pair appears once even when both members are inside context.
--      Params: $genome_id, $motif_set_id, $chrom, $start, $strand,
--      $neighbor_motif_id, $score_mode, $pseudocount
SELECT
    p.cofactor_pair_id,
    p.cofactor_motif_id,
    p.left_locus_id,
    p.right_locus_id,
    p.pair_member_interval_distance_bp,
    p.pair_member_distance_band,
    p.pair_arrangement,
    p.n_pair_members_in_context,
    p.pair_fully_within_context,
    p.nearest_member_locus_id,
    p.nearest_member_anchor_neighbor_interval_distance_bp,
    p.nearest_member_anchor_distance_band,
    p.nearest_member_anchor_oriented_side,
    p.nearest_member_relative_orientation,
    p.pair_min_score,
    p.pair_sum_score
FROM tp73_cofactor_pair_context p
JOIN tp73_pair_feature a USING (anchor_hit_id)
WHERE p.genome_id = $genome_id
  AND p.motif_set_id = $motif_set_id
  AND p.chrom = $chrom
  AND a.start = $start
  AND a.strand = $strand
  AND p.cofactor_motif_id = $neighbor_motif_id
  AND p.score_mode = $score_mode
  AND p.pseudocount = $pseudocount
ORDER BY p.nearest_member_anchor_neighbor_interval_distance_bp,
         p.pair_member_interval_distance_bp, p.cofactor_pair_id;

-- Q16. Inspect one explicitly versioned convenient-threshold set. Pending and
--      no-gain motifs remain visible rather than silently falling back to zero.
--      Params: $threshold_set_id, $genome_id, $motif_set_id, $score_mode,
--      $pseudocount, $threshold_role, $target_motif_id,
--      $calibration_stratum_id
SELECT
    motif_id,
    motif_name,
    recommended_threshold,
    useful_threshold_min,
    useful_threshold_max,
    selected_retained_fraction,
    selected_metric_gain,
    selected_adjusted_odds_ratio,
    association_direction,
    calibration_status,
    calibration_scope,
    evidence_dataset_id
FROM motif_score_threshold
WHERE threshold_set_id = $threshold_set_id
  AND genome_id = $genome_id
  AND motif_set_id = $motif_set_id
  AND score_mode = $score_mode
  AND pseudocount = $pseudocount
  AND threshold_role = $threshold_role
  AND target_motif_id = $target_motif_id
  AND calibration_stratum_id = $calibration_stratum_id
ORDER BY motif_name, motif_id;

-- Q17. Apply one threshold set to raw TP73-neighbor relationships. This view of
--      the data remains derived: motif_context_pair and its continuous scores
--      are not discarded or rewritten.
--      Params: $threshold_set_id, $genome_id, $motif_set_id, $chrom,
--      $threshold_role, $target_motif_id, $score_mode, $pseudocount,
--      $calibration_stratum_id
SELECT
    p.anchor_hit_id,
    p.chrom,
    p.anchor_start,
    p.anchor_end,
    p.anchor_strand,
    p.neighbor_hit_id,
    p.neighbor_start,
    p.neighbor_end,
    p.neighbor_motif_id,
    p.neighbor_motif_name,
    p.neighbor_strand,
    p.neighbor_score,
    t.recommended_threshold,
    p.anchor_neighbor_interval_distance_bp,
    p.interval_relation,
    p.interval_distance_band,
    p.relative_orientation,
    p.anchor_oriented_side
FROM motif_context_pair p
JOIN motif_convenient_threshold t
  ON t.genome_id = p.genome_id
 AND t.motif_set_id = p.motif_set_id
 AND t.motif_id = p.neighbor_motif_id
 AND t.score_mode = p.score_mode
 AND t.pseudocount = p.pseudocount
 AND t.background_model_id = p.background_model_id
 AND t.pseudocount_scheme = p.pseudocount_scheme
WHERE t.threshold_set_id = $threshold_set_id
  AND t.genome_id = $genome_id
  AND t.motif_set_id = $motif_set_id
  AND p.chrom = $chrom
  AND t.threshold_role = $threshold_role
  AND t.target_motif_id = $target_motif_id
  AND t.calibration_stratum_id = $calibration_stratum_id
  AND t.score_mode = $score_mode
  AND t.pseudocount = $pseudocount
  AND t.threshold_inclusive
  AND t.context_distance_metric = 'signed_interval_edge_distance'
  AND (t.context_min_interval_distance_bp IS NULL
       OR p.anchor_neighbor_interval_distance_bp >=
          t.context_min_interval_distance_bp)
  AND (t.context_max_interval_distance_bp IS NULL
       OR p.anchor_neighbor_interval_distance_bp <=
          t.context_max_interval_distance_bp)
  AND (t.context_relation_filter = 'any'
       OR p.interval_relation = t.context_relation_filter)
  AND p.neighbor_score >= t.recommended_threshold
ORDER BY p.anchor_start, p.neighbor_motif_id, p.neighbor_score DESC;

-- Q18. One zero-complete count for every physical TP73 anchor and one selected
--      neighboring motif. Counts are distinct alignment spans meeting the
--      motif's inclusive convenient threshold and its recorded interval bounds;
--      strand alternatives at one span count once.
--      Params: $threshold_set_id, $genome_id, $motif_set_id, $chrom,
--      $threshold_role, $target_motif_id, $neighbor_motif_id, $score_mode,
--      $pseudocount, $calibration_stratum_id
SELECT
    anchor_locus_id,
    chrom,
    anchor_start,
    anchor_end,
    target_motif_id,
    motif_id AS neighbor_motif_id,
    motif_name AS neighbor_motif_name,
    context_score AS best_unthresholded_neighbor_score,
    recommended_threshold,
    n_neighbor_loci_above_threshold,
    has_neighbor_locus_above_threshold,
    context_min_interval_distance_bp,
    context_max_interval_distance_bp,
    context_relation_filter,
    calibration_status
FROM tp73_motif_threshold_count
WHERE threshold_set_id = $threshold_set_id
  AND genome_id = $genome_id
  AND motif_set_id = $motif_set_id
  AND chrom = $chrom
  AND threshold_role = $threshold_role
  AND target_motif_id = $target_motif_id
  AND motif_id = $neighbor_motif_id
  AND score_mode = $score_mode
  AND pseudocount = $pseudocount
  AND calibration_stratum_id = $calibration_stratum_id
ORDER BY anchor_start, anchor_end, anchor_locus_id;

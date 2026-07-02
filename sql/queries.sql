-- ============================================================================
-- queries.sql — canned, low-latency queries over schema.sql
-- ----------------------------------------------------------------------------
-- These are the queries agents (OpenClaw / Claude / Codex) and humans reach for.
-- Each is a single statement, parameterized with $placeholders so it can be
-- prepared once and reused. Run `.read sql/schema.sql` first.
-- ============================================================================

-- Q0. Hot promoter card — first lookup for an agent answering "what do we know
--     about this gene's promoter architecture?" Param: $gene_name
SELECT *
FROM promoter_card
WHERE gene_name = $gene_name;

-- Q1. Region browse — every motif hit in a window. Prunes by chrom partition
--     + start zonemap, so latency is ~row-group scan, not full table.
--     Params: $chrom, $start, $end
SELECT chrom, start, "end", motif_id, motif_name, strand, score_mode, pseudocount, pwm_relative_score
FROM motif_hit
WHERE chrom = $chrom AND start >= $start AND "end" <= $end
ORDER BY start;

-- Q2. Promoter architecture for one gene — what the ML sees, human-readable.
--     Param: $gene_name
SELECT f.motif_id, f.score_mode, f.pseudocount, f.tf_family, f.binding_unit_model,
       f.n_hits, f.n_strong_hits, f.max_rel_score, f.min_abs_tss_distance
FROM promoter_arch_feature f
JOIN gene g ON g.gene_id = f.gene_id
WHERE g.gene_name = $gene_name
ORDER BY f.n_strong_hits DESC, f.max_rel_score DESC;

-- Q3. The ML pull, WIDE — one row per gene, one column per (tf_family feature),
--     plus the TA-vs-DN label. This is the training frame; PIVOT happens here so
--     the stored tables stay long. Param: $cell_line
PIVOT (
    SELECT m.gene_id, m.tf_family, m.n_strong_hits, m.log2fc_ta_vs_dn, m.label_ta_up
    FROM ml_ta_vs_dn m
    WHERE m.cell_line = $cell_line AND m.tf_family IS NOT NULL
)
ON tf_family
USING SUM(n_strong_hits)
GROUP BY gene_id, log2fc_ta_vs_dn, label_ta_up;

-- Q4. Experimental cross-check — TA-vs-DN differential CUT&RUN signal at a
--     gene's promoter loci, kept separate by antibody. Keep this OUT of the
--     sequence-only feature set (leakage); use it to validate predictions.
--     Params: $gene_name, $cell_line
SELECT c.chrom, c.start, c."end", c.cell_line, c.antibody,
       AVG(c.signal) FILTER (WHERE c.isoform = 'TA') AS ta_signal,
       AVG(c.signal) FILTER (WHERE c.isoform = 'DN') AS dn_signal,
       AVG(c.signal) FILTER (WHERE c.isoform = 'TA')
     - AVG(c.signal) FILTER (WHERE c.isoform = 'DN') AS ta_minus_dn
FROM cutandrun c
JOIN promoter p
  ON c.chrom = p.chrom AND c.start < p.promoter_end AND c."end" > p.promoter_start
JOIN gene g ON g.gene_id = p.gene_id
WHERE g.gene_name = $gene_name AND c.cell_line = $cell_line
GROUP BY c.chrom, c.start, c."end", c.cell_line, c.antibody
ORDER BY ABS(ta_minus_dn) DESC;

-- Q5. Architecture enrichment — do dimer-of-dimers TP53-family sites associate
--     with TA-up genes? A quick 2x2 signal before any modelling. Param: $cell_line
SELECT has_arch,
       AVG(label_ta_up)          AS frac_ta_up,
       COUNT(*)                  AS n_genes
FROM (
    SELECT m.gene_id,
           MAX(CASE WHEN m.has_dimer_of_dimers THEN 1 ELSE 0 END) AS has_arch,
           MAX(m.label_ta_up)                                     AS label_ta_up
    FROM ml_ta_vs_dn m
    WHERE m.cell_line = $cell_line
    GROUP BY m.gene_id
)
GROUP BY has_arch;

-- Q6. Gene-set slice — restrict the ML frame to a pathway (e.g. an EMT set) so a
--     model can be trained per program. Params: $set_name, $cell_line
SELECT m.*
FROM ml_ta_vs_dn m
JOIN gene_set s ON s.gene_id = m.gene_id
WHERE s.set_name = $set_name AND m.cell_line = $cell_line;

-- Q7. Dense score slice — reconstruct per-base rows from dense score blocks
--     for a small genomic interval. Params: $motif_id, $chrom, $strand,
--     $score_mode, $pseudocount, $start, $end
SELECT chrom, start, "end", motif_id, motif_name, strand,
       score_mode, pseudocount, score
FROM motif_score_dense
WHERE motif_id = $motif_id
  AND chrom = $chrom
  AND strand = $strand
  AND score_mode = $score_mode
  AND pseudocount = $pseudocount
  AND start >= $start
  AND start < $end
ORDER BY start;

-- Q8. Dense calibration result — compare CUT&RUN support by score bin for one
--     motif/scoring run. Params: $motif_id, $chrom, $score_mode, $pseudocount
SELECT motif_id, motif_name, chrom, strand, score_mode, pseudocount,
       bin_order, bin_label, lower_bound, upper_bound,
       cell_line, isoform, antibody, replicate,
       n_windows, n_covered_windows, overlap_fraction, baseline_fraction,
       enrichment_ratio, log2_enrichment, mean_signal, max_signal
FROM motif_cutandrun_score_bin_stats
WHERE motif_id = $motif_id
  AND chrom = $chrom
  AND score_mode = $score_mode
  AND pseudocount = $pseudocount
ORDER BY cell_line, antibody, isoform, replicate, strand, bin_order;

-- Q9. Calibration-bin ladder — the common bins used for log2-relative-risk and
--     log-odds dense chr1 comparisons.
SELECT *
FROM score_calibration_bin
ORDER BY bin_order;

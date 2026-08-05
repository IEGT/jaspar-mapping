#!/usr/bin/env bash

set -euo pipefail

if ! command -v duckdb >/dev/null 2>&1; then
    echo "E: duckdb CLI is required for this test." >&2
    exit 1
fi

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-motif-context.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

cat > "$temporary/hits.sql" <<'SQL'
COPY (
    SELECT t.*,
           CASE WHEN motif_id = 'MTHRESH.1' THEN 3.0 ELSE 0.0 END
               AS minimum_score
    FROM (VALUES
        ('1', 140::BIGINT, 156::BIGINT, 'MA0861.2', 'TP73', '+', 10.0, 'log2_relative_risk', 1.0, 0.95),
        ('1', 140::BIGINT, 156::BIGINT, 'MA0861.2', 'TP73', '-', 9.5,  'log2_relative_risk', 1.0, 0.99),
        ('1', 145::BIGINT, 161::BIGINT, 'MA0861.2', 'TP73', '+', 6.5,  'log2_relative_risk', 1.0, 0.86),
        ('1', 158::BIGINT, 174::BIGINT, 'MA0861.2', 'TP73', '+', 8.0,  'log2_relative_risk', 1.0, 0.90),
        ('1', 200::BIGINT, 216::BIGINT, 'MA0861.2', 'TP73', '+', 7.0,  'log2_relative_risk', 1.0, 0.88),
        ('1', 220::BIGINT, 236::BIGINT, 'MA0861.2', 'TP73', '-', 9.0,  'log2_relative_risk', 1.0, 0.92),
        ('1', 400::BIGINT, 416::BIGINT, 'MA0861.2', 'TP73', '+', 5.0,  'log2_relative_risk', 1.0, 0.81),
        ('1', 600::BIGINT, 616::BIGINT, 'MA0861.2', 'TP73', '+', -2.0, 'log2_relative_risk', 1.0, 0.30),
        ('1', 620::BIGINT, 636::BIGINT, 'MA0861.2', 'TP73', '+', -0.5, 'log2_relative_risk', 1.0, 0.45),
        ('1', 1000::BIGINT, 1016::BIGINT, 'MA0861.2', 'TP73', '+', 4.0, 'log2_relative_risk', 1.0, 0.79),
        ('1', 130::BIGINT, 139::BIGINT, 'MA0079.5', 'SP1',  '+', 6.0,  'log2_relative_risk', 1.0, 0.85),
        ('1', 410::BIGINT, 419::BIGINT, 'MA0079.5', 'SP1',  '+', 4.0,  'log2_relative_risk', 1.0, 0.80),
        ('1', 499::BIGINT, 508::BIGINT, 'MA0079.5', 'SP1',  '+', 3.0,  'log2_relative_risk', 1.0, 0.70),
        ('1', 156::BIGINT, 165::BIGINT, 'MABUT.1', 'ABUT', '+', 2.0, 'log2_relative_risk', 1.0, 0.70),
        ('1', 144::BIGINT, 150::BIGINT, 'MNEST.1', 'NEST', '+', 2.0, 'log2_relative_risk', 1.0, 0.70),
        ('1', 130::BIGINT, 170::BIGINT, 'MCONTAIN.1', 'CONTAIN', '-', 2.0, 'log2_relative_risk', 1.0, 0.70),
        ('1', 1166::BIGINT, 1206::BIGINT, 'MWIDE.1', 'WIDE', '+', 2.0, 'log2_relative_risk', 1.0, 0.70),
        ('1', 300::BIGINT, 306::BIGINT, 'MPAIR.1', 'PAIR', '+', 3.0, 'log2_relative_risk', 1.0, 0.75),
        ('1', 300::BIGINT, 306::BIGINT, 'MPAIR.1', 'PAIR', '-', 2.5, 'log2_relative_risk', 1.0, 0.72),
        ('1', 306::BIGINT, 312::BIGINT, 'MPAIR.1', 'PAIR', '+', 4.0, 'log2_relative_risk', 1.0, 0.80),
        ('1', 318::BIGINT, 324::BIGINT, 'MPAIR.1', 'PAIR', '-', 5.0, 'log2_relative_risk', 1.0, 0.85),
        ('1', 1166::BIGINT, 1176::BIGINT, 'MBOUND.1', 'BOUND', '+', 3.0, 'log2_relative_risk', 1.0, 0.75),
        ('1', 1177::BIGINT, 1187::BIGINT, 'MBOUND.1', 'BOUND', '-', 4.0, 'log2_relative_risk', 1.0, 0.80),
        ('1', 422::BIGINT, 428::BIGINT, 'MBEST.1', 'BEST', '+', 2.0, 'log2_relative_risk', 1.0, 0.60),
        ('1', 430::BIGINT, 436::BIGINT, 'MBEST.1', 'BEST', '-', 5.0, 'log2_relative_risk', 1.0, 0.85),
        ('1', 436::BIGINT, 442::BIGINT, 'MBEST.1', 'BEST', '+', 1.0, 'log2_relative_risk', 1.0, 0.55),
        ('1', 440::BIGINT, 446::BIGINT, 'MTIE.1', 'TIE', '+', 3.0, 'log2_relative_risk', 1.0, 0.70),
        ('1', 440::BIGINT, 446::BIGINT, 'MTIE.1', 'TIE', '-', 3.0, 'log2_relative_risk', 1.0, 0.70),
        ('1', 464::BIGINT, 470::BIGINT, 'MTIE.1', 'TIE', '+', 3.0, 'log2_relative_risk', 1.0, 0.70),
        ('1', 448::BIGINT, 454::BIGINT, 'MTHRESH.1', 'THRESH', '+', 2.0, 'log2_relative_risk', 1.0, 0.60),
        ('1', 456::BIGINT, 462::BIGINT, 'MTHRESH.1', 'THRESH', '+', 4.0, 'log2_relative_risk', 1.0, 0.80),
        ('1', 456::BIGINT, 462::BIGINT, 'MTHRESH.1', 'THRESH', '-', 3.5, 'log2_relative_risk', 1.0, 0.99),
        ('1', 2000::BIGINT, 2006::BIGINT, 'MFAR.1', 'FAR', '+', 3.0, 'log2_relative_risk', 1.0, 0.75),
        ('1', 2006::BIGINT, 2012::BIGINT, 'MFAR.1', 'FAR', '-', 4.0, 'log2_relative_risk', 1.0, 0.80)
    ) AS t(chrom, start, "end", motif_id, motif_name, strand, score,
           score_mode, pseudocount, pwm_relative_score)
) TO 'motif_hit.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);

-- Match the direct sparse writer: no motif_name column and word-form strands.
COPY (
    SELECT * FROM (VALUES
        ('1', 100::BIGINT, 116::BIGINT, 'MA0861.2', 'plus', 5.0,
         'log2_relative_risk', 1::BIGINT, 0.80),
        ('1', 120::BIGINT, 136::BIGINT, 'MA0861.2', 'minus', 4.0,
         'log2_relative_risk', 1::BIGINT, 0.75)
    ) AS t(chrom, start, "end", motif_id, strand, score,
           score_mode, pseudocount, pwm_relative_score)
) TO 'direct_sparse.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL

(
    cd "$temporary"
    duckdb :memory: < hits.sql >/dev/null
)

cat > "$temporary/annotation.gtf" <<'EOF'
##gtf-version 3
1	test	transcript	51	300	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENE1"; transcript_biotype "protein_coding";
1	test	exon	51	120	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENE1"; exon_number "1";
1	test	exon	181	300	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENE1"; exon_number "2";
EOF

"$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/motif_hit.parquet" \
    --gtf "$temporary/annotation.gtf" \
    --output "$temporary/context_package" \
    --anchor-motif MA0861.2 \
    --motif-set-id synthetic_jaspar2026 --genome-id synthetic_grch38_v1 \
    --anchor-minimum-score 0 --partner-minimum-score 0 \
    --score-mode log2_relative_risk --pseudocount 1 \
    --chrom 1 --capture-flank 150 --context-flank 150 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB \
    --temp-directory "$temporary/external_duckdb_spill" >/dev/null

[[ -d "$temporary/external_duckdb_spill" ]] || {
    echo "E: external DuckDB spill directory was not created" >&2
    exit 1
}

"$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/direct_sparse.parquet" \
    --output "$temporary/direct_sparse_context" \
    --anchor-motif MA0861.2 \
    --motif-set-id synthetic_jaspar2026 --genome-id synthetic_grch38_v1 \
    --anchor-minimum-score 0 --partner-minimum-score 0 \
    --score-mode log2_relative_risk --pseudocount 1 \
    --chrom 1 --capture-flank 100 --context-flank 50 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null

(
    cd "$temporary/direct_sparse_context"
    duckdb context.duckdb -c \
        "SELECT CASE WHEN NOT EXISTS (
             SELECT 1 FROM tp73_pair_feature
             WHERE start = 100 AND strand = '+'
               AND motif_name = 'MA0861.2'
               AND pair_class = 'tandem_opposite_orientation'
               AND nearest_tandem_inter_motif_gap_bp = 4
         ) THEN error('direct sparse strand/name normalization failed') END;" \
        >/dev/null
)

"$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/motif_hit.parquet" \
    --output "$temporary/context_summary" --output-tier summary \
    --anchor-motif MA0861.2 \
    --motif-set-id synthetic_jaspar2026 --genome-id synthetic_grch38_v1 \
    --anchor-minimum-score 0 --partner-minimum-score 0 \
    --score-mode log2_relative_risk --pseudocount 1 \
    --chrom 1 --capture-flank 150 --context-flank 150 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null

(
    cd "$temporary/context_summary"
    duckdb context.duckdb -bail -c "
        SELECT CASE WHEN (SELECT output_tier FROM context_run_config) <> 'summary'
            THEN error('summary output tier was not recorded') END;
        SELECT CASE WHEN (SELECT COUNT(*) FROM motif_context_pair) <> 0
            THEN error('summary tier retained raw motif-context pairs') END;
        SELECT CASE WHEN (SELECT COUNT(*) FROM cofactor_motif_pair) <> 0
            THEN error('summary tier retained raw cofactor pairs') END;
        SELECT CASE WHEN (SELECT COUNT(*) FROM tp73_cofactor_pair_context) <> 0
            THEN error('summary tier retained raw pair attribution') END;
        SELECT CASE WHEN (SELECT COUNT(*) FROM tp73_motif_context_summary) = 0
            THEN error('summary tier dropped compact motif features') END;
        SELECT CASE WHEN (SELECT COUNT(*) FROM anchor_motif_band_feature) = 0
            THEN error('summary tier dropped strongest-hit band features') END;
        SELECT CASE WHEN (SELECT COUNT(*) FROM tp73_cofactor_pair_summary) = 0
            THEN error('summary tier dropped compact cofactor-pair features') END;" \
        >/dev/null
)

"$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/motif_hit.parquet" \
    --motif-hit-source-label durable/context/input/**/*.parquet \
    --output "$temporary/context_band" --output-tier band \
    --anchor-motif MA0861.2 \
    --motif-set-id synthetic_jaspar2026 --genome-id synthetic_grch38_v1 \
    --anchor-minimum-score 0 --partner-minimum-score 0 \
    --score-mode log2_relative_risk --pseudocount 1 \
    --chrom 1 --capture-flank 150 --context-flank 150 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null

(
    cd "$temporary/context_band"
    duckdb context.duckdb -bail -c "
        SELECT CASE WHEN (SELECT output_tier FROM context_run_config) <> 'band'
            THEN error('band output tier was not recorded') END;
        SELECT CASE WHEN (SELECT motif_hit_source FROM context_run_config)
                <> 'durable/context/input/**/*.parquet'
            THEN error('durable motif-hit provenance label was not recorded') END;
        SELECT CASE WHEN (SELECT COUNT(*) FROM anchor_motif_band_feature) = 0
            THEN error('band tier dropped per-band features') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM anchor_motif_band_feature
            WHERE neighbor_motif_id = 'MA0861.2'
        ) THEN error('band tier duplicated the TP73 anchor motif') END;
        SELECT CASE WHEN
               (SELECT COUNT(*) FROM motif_context_pair) <> 0
            OR (SELECT COUNT(*) FROM tp73_anchor_locus) <> 0
            OR (SELECT COUNT(*) FROM tp73_motif_context_summary) <> 0
            OR (SELECT COUNT(*) FROM tp73_cofactor_pair_summary) <> 0
            OR (SELECT COUNT(*) FROM tp73_pair_feature) <> 0
            OR (SELECT COUNT(*) FROM tp73_context_anchor) <> 0
            THEN error('band tier retained a shared or compatibility surface') END;" \
        >/dev/null
)

"$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/direct_sparse.parquet" \
    --output "$temporary/context_band_empty" --output-tier band \
    --anchor-motif MA0861.2 \
    --motif-set-id synthetic_jaspar2026 --genome-id synthetic_grch38_v1 \
    --anchor-minimum-score 0 --partner-minimum-score 0 \
    --score-mode log2_relative_risk --pseudocount 1 \
    --chrom 1 --capture-flank 150 --context-flank 150 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null
(
    cd "$temporary/context_band_empty"
    duckdb context.duckdb -bail -c "
        SELECT CASE WHEN (SELECT COUNT(*) FROM anchor_motif_band_feature) <> 0
            THEN error('empty band package retained a below-threshold row') END;" \
        >/dev/null
)

cat > "$temporary/assertions.sql" <<'SQL'
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 140 AND anchor_strand = '+' AND neighbor_start = 158
      AND genomic_center_distance_bp = 18
      AND anchor_oriented_center_distance_bp = 18
      AND relative_orientation = 'same' AND inter_motif_gap_bp = 2
      AND interval_overlap_bp = 0 AND is_tandem_tp73
) THEN error('same-orientation tandem pair was not retained') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 220 AND anchor_strand = '-' AND neighbor_start = 200
      AND genomic_center_distance_bp = -20
      AND anchor_oriented_center_distance_bp = 20
      AND relative_orientation = 'opposite' AND inter_motif_gap_bp = 4
      AND interval_overlap_bp = 0 AND is_tandem_tp73
) THEN error('reverse-anchor tandem orientation was not normalized') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 140 AND anchor_strand = '+' AND neighbor_start = 140
      AND neighbor_strand = '-' AND same_alignment_span
      AND same_anchor_motif_span AND NOT is_tandem_tp73
) THEN error('same-span opposite orientation was mistaken for a tandem site') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 140 AND anchor_strand = '+' AND neighbor_start = 145
      AND interval_overlap_bp = 11 AND inter_motif_gap_bp = 0
      AND NOT is_tandem_tp73
) THEN error('overlapping shifted TP73 match was mistaken for a tandem site') END;

-- This pair crosses a center-prefilter bin boundary. The +/-1-bin expansion
-- must retain it before exact interval filtering.
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 400 AND anchor_strand = '+' AND neighbor_start = 499
      AND absolute_center_distance_bp = 95.5 AND capture_flank_bp = 150
) THEN error('valid pair across a capture-bin boundary was dropped') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 140 AND anchor_strand = '+'
      AND neighbor_motif_id = 'MABUT.1'
      AND anchor_neighbor_interval_distance_bp = 0
      AND interval_relation = 'abutting'
      AND interval_distance_band = 'adjacent_0_5'
      AND interval_overlap_bp = 0 AND inter_motif_gap_bp = 0
      AND within_5 AND within_20 AND within_50 AND within_100 AND within_150
) THEN error('abutting interval geometry is incorrect') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 140 AND anchor_strand = '+'
      AND neighbor_motif_id = 'MA0079.5' AND neighbor_start = 130
      AND anchor_neighbor_interval_distance_bp = 1
      AND interval_relation = 'disjoint'
      AND interval_distance_band = 'adjacent_0_5'
) THEN error('one-base interval gap is incorrect') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 140 AND anchor_strand = '+'
      AND neighbor_start = 145 AND neighbor_motif_id = 'MA0861.2'
      AND anchor_neighbor_interval_distance_bp = -11
      AND interval_relation = 'partial_overlap'
      AND interval_distance_band = 'overlap'
) THEN error('partial-overlap interval geometry is incorrect') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 140 AND anchor_strand = '+'
      AND neighbor_motif_id = 'MNEST.1'
      AND anchor_neighbor_interval_distance_bp = -6
      AND interval_relation = 'anchor_contains_neighbor'
) THEN error('anchor containment was not distinguished from partial overlap') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 140 AND anchor_strand = '+'
      AND neighbor_motif_id = 'MCONTAIN.1'
      AND anchor_neighbor_interval_distance_bp = -16
      AND interval_relation = 'neighbor_contains_anchor'
) THEN error('neighbor containment was not distinguished from partial overlap') END;

-- Interval distance, not center distance, defines the 150 bp boundary. The
-- 40 bp neighbor is exactly 150 bp away despite a center distance of 178 bp.
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 1000 AND anchor_strand = '+'
      AND neighbor_motif_id = 'MWIDE.1'
      AND anchor_neighbor_interval_distance_bp = 150
      AND absolute_center_distance_bp = 178
      AND interval_distance_band = 'gap_101_150'
      AND within_150 AND NOT within_100
) THEN error('wide interval-near motif was lost by the center prefilter') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73_anchor_locus
    WHERE start = 140 AND "end" = 156
      AND orientation_state = 'ambiguous'
      AND n_orientation_records = 2
      AND plus_score = 10.0 AND minus_score = 9.5
) THEN error('TP73 orientation records were not collapsed at locus grain') END;

SELECT CASE WHEN (SELECT COUNT(*) FROM cofactor_motif_locus
                  WHERE motif_id = 'MPAIR.1') <> 3
    THEN error('same-span cofactor strand alternatives formed duplicate loci') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM cofactor_motif_locus
    WHERE motif_id = 'MPAIR.1' AND start = 300 AND "end" = 306
      AND orientation_state = 'ambiguous' AND n_orientation_records = 2
) THEN error('cofactor locus orientation ambiguity was lost') END;

SELECT CASE WHEN (SELECT COUNT(*) FROM cofactor_motif_pair
                  WHERE motif_id = 'MPAIR.1') <> 3
    THEN error('three cofactor loci did not produce three canonical pairs') END;

SELECT CASE WHEN EXISTS (
    SELECT 1 FROM cofactor_motif_pair WHERE left_locus_id = right_locus_id
) THEN error('a cofactor locus was paired with itself') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM cofactor_motif_pair
    WHERE motif_id = 'MPAIR.1' AND left_start = 300 AND right_start = 306
      AND pair_member_interval_distance_bp = 0
      AND pair_member_interval_relation = 'abutting'
      AND pair_arrangement = 'ambiguous'
) THEN error('ambiguous cofactor pair was classified incorrectly') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM cofactor_motif_pair
    WHERE motif_id = 'MPAIR.1' AND left_start = 306 AND right_start = 318
      AND pair_member_interval_distance_bp = 6
      AND pair_member_distance_band = 'gap_6_20'
      AND pair_arrangement = 'convergent'
) THEN error('convergent cofactor pair was classified incorrectly') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM cofactor_locus_pair_feature
    WHERE motif_id = 'MPAIR.1' AND start = 306
      AND n_same_motif_partner_loci = 2
      AND n_ambiguous_pairs = 1 AND n_convergent_pairs = 1
) THEN error('cofactor locus pair features are incorrect') END;

SELECT CASE WHEN (SELECT COUNT(*) FROM tp73_cofactor_pair_context
                  WHERE anchor_hit_id = (
                      SELECT anchor_hit_id FROM tp73_pair_feature
                      WHERE start = 400 AND strand = '+'
                  ) AND cofactor_motif_id = 'MPAIR.1'
                  AND n_pair_members_in_context = 2) <> 3
    THEN error('canonical cofactor pairs were multiplied during TP73 attribution') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73_cofactor_pair_context
    WHERE anchor_hit_id = (
              SELECT anchor_hit_id FROM tp73_pair_feature
              WHERE start = 1000 AND strand = '+'
          )
      AND cofactor_motif_id = 'MBOUND.1'
      AND n_pair_members_in_context = 1
      AND left_member_in_context AND NOT right_member_in_context
      AND nearest_member_anchor_neighbor_interval_distance_bp = 150
      AND nearest_member_anchor_distance_band = 'gap_101_150'
) THEN error('one-member cofactor-pair boundary attribution is incorrect') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM cofactor_motif_locus
    WHERE motif_id = 'MBOUND.1' AND start = 1166
      AND in_any_tp73_context
) OR NOT EXISTS (
    SELECT 1 FROM cofactor_motif_locus
    WHERE motif_id = 'MBOUND.1' AND start = 1177
      AND NOT in_any_tp73_context
) THEN error('context seed and outside pair partner were not distinguished') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM cofactor_locus_pair_feature
    WHERE motif_id = 'MBOUND.1' AND start = 1166
      AND in_any_tp73_context AND n_same_motif_partner_loci = 1
) OR EXISTS (
    SELECT 1 FROM cofactor_locus_pair_feature
    WHERE motif_id = 'MBOUND.1' AND start = 1177
) THEN error('cofactor locus features were not scoped to TP73 context loci') END;

SELECT CASE WHEN EXISTS (
    SELECT 1 FROM cofactor_motif_locus WHERE motif_id = 'MFAR.1'
) OR EXISTS (
    SELECT 1 FROM cofactor_motif_pair WHERE motif_id = 'MFAR.1'
) THEN error('cofactor pair unrelated to every TP73 context was retained') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73_motif_context_summary
    WHERE anchor_hit_id = (
              SELECT anchor_hit_id FROM tp73_pair_feature
              WHERE start = 400 AND strand = '+'
          )
      AND neighbor_motif_id = 'MPAIR.1'
      AND relative_orientation_state = 'ambiguous'
      AND n_neighbor_loci = 1 AND n_orientation_records = 2
) THEN error('compact motif-context summary did not collapse strand alternatives') END;

-- The strongest score and every descriptor below must come from one physical
-- neighboring locus. A nearer but weaker locus must not donate its distance.
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM anchor_motif_band_feature
    WHERE anchor_start = 400 AND anchor_end = 416
      AND neighbor_motif_id = 'MBEST.1'
      AND interval_distance_band = 'gap_6_20'
      AND interval_distance_band_order = 2
      AND qualifying_threshold = 0 AND threshold_inclusive
      AND n_neighbor_loci_above_threshold = 3
      AND n_best_score_ties = 1
      AND best_neighbor_start = 430 AND best_neighbor_end = 436
      AND best_neighbor_score = 5.0
      AND best_neighbor_pwm_relative_score = 0.85
      AND best_interval_distance_bp = 14
      AND best_interval_overlap_bp = 0
      AND best_inter_motif_gap_bp = 14
      AND best_interval_relation = 'disjoint'
      AND best_genomic_center_distance_bp = 25
      AND best_anchor_oriented_center_distance_bp = 25
      AND best_genomic_side = 'right'
      AND best_anchor_oriented_side = 'downstream'
      AND anchor_orientation_state = 'plus'
      AND anchor_best_orientation_state = 'plus'
      AND best_neighbor_orientation_state = 'minus'
      AND best_relative_orientation_state = 'opposite'
      AND best_hit_pair_architecture_assessed
      AND best_hit_has_same_motif_partner
      AND best_hit_n_same_motif_partner_loci = 2
      AND best_hit_n_convergent_pairs = 1
      AND best_hit_n_divergent_pairs = 1
      AND best_hit_nearest_pair_member_distance_bp = 0
      AND best_hit_nearest_pair_arrangement = 'divergent'
      AND best_hit_nearest_partner_score = 1.0
      AND best_hit_best_pair_min_score = 2.0
      AND best_hit_best_pair_sum_score = 7.0
) THEN error('strongest per-band locus fields were uncoupled or incomplete') END;

-- Same-span strand alternatives count once. A score tie between two physical
-- loci is reported and resolved deterministically by distance for the winner.
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM anchor_motif_band_feature
    WHERE anchor_start = 400 AND anchor_end = 416
      AND neighbor_motif_id = 'MTIE.1'
      AND interval_distance_band = 'gap_21_50'
      AND n_neighbor_loci_above_threshold = 2
      AND n_best_score_ties = 2
      AND best_neighbor_start = 440 AND best_neighbor_end = 446
      AND best_neighbor_score = 3.0
      AND best_neighbor_n_orientation_records = 2
      AND best_neighbor_orientation_state = 'ambiguous'
      AND best_relative_orientation_state = 'ambiguous'
) THEN error('physical-locus counting or strongest-score tie handling failed') END;

-- The input's per-motif minimum_score is inclusive. The score-2 locus is not
-- counted and cannot make the qualifying score-4 locus appear paired.
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM anchor_motif_band_feature
    WHERE anchor_start = 400 AND anchor_end = 416
      AND neighbor_motif_id = 'MTHRESH.1'
      AND interval_distance_band = 'gap_21_50'
      AND qualifying_threshold = 3.0
      AND n_neighbor_loci_above_threshold = 1
      AND best_neighbor_start = 456 AND best_neighbor_end = 462
      AND best_neighbor_score = 4.0
      AND best_neighbor_pwm_relative_score = 0.80
      AND best_neighbor_n_orientation_records = 2
      AND best_neighbor_plus_score = 4.0
      AND best_neighbor_minus_score = 3.5
      AND best_neighbor_orientation_state = 'plus'
      AND best_hit_pair_architecture_assessed
      AND NOT best_hit_has_same_motif_partner
      AND best_hit_n_same_motif_partner_loci = 0
) THEN error('per-motif qualifying threshold was not applied consistently') END;

-- TP73 neighbors remain useful context observations, but their tandem state is
-- represented by tp73_pair_feature rather than the generic cofactor pair layer.
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM anchor_motif_band_feature
    WHERE anchor_start = 140 AND anchor_end = 156
      AND neighbor_motif_id = 'MA0861.2'
      AND interval_distance_band = 'overlap'
      AND NOT best_hit_pair_architecture_assessed
      AND best_hit_has_same_motif_partner IS NULL
      AND best_hit_n_same_motif_partner_loci IS NULL
) THEN error('unassessed anchor-motif pair architecture was encoded as absent') END;

-- Physical TP73 anchors are not duplicated by orientation alternatives. Keep
-- both their raw ambiguity and the orientation selected by the stronger score.
SELECT CASE WHEN (SELECT COUNT(*) FROM anchor_motif_band_feature
                  WHERE anchor_start = 140 AND anchor_end = 156
                    AND neighbor_motif_id = 'MA0079.5'
                    AND interval_distance_band = 'adjacent_0_5') <> 1
    OR NOT EXISTS (
        SELECT 1 FROM anchor_motif_band_feature
        WHERE anchor_start = 140 AND anchor_end = 156
          AND neighbor_motif_id = 'MA0079.5'
          AND interval_distance_band = 'adjacent_0_5'
          AND anchor_orientation_state = 'ambiguous'
          AND anchor_best_orientation_state = 'plus'
          AND anchor_best_score = 10.0
          AND anchor_best_pwm_relative_score = 0.95
    )
    THEN error('physical anchor grain or anchor orientation state is incorrect') END;

SELECT CASE WHEN EXISTS (
    SELECT anchor_locus_id, genome_id, motif_set_id, neighbor_motif_id,
           interval_distance_band
    FROM anchor_motif_band_feature
    GROUP BY ALL
    HAVING COUNT(*) <> 1
) THEN error('anchor/motif/distance-band feature key is not unique') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM legacy_tp73_context_100
    WHERE anchor_start = 400 AND neighbor_start = 499
      AND legacy_genomic_start_distance_bp = 99
) THEN error('historical start-distance compatibility view is incorrect') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73_context_anchor
    WHERE start = 140 AND strand = '+'
      AND has_tandem_tp73 AND n_tandem_tp73_partners = 1
      AND pair_class = 'tandem_same_orientation'
      AND n_tandem_tp73_partner_loci = 1
      AND nearest_tandem_oriented_distance_bp = 18
      AND nearest_tandem_relative_orientation = 'same'
      AND nearest_gene_id = 'G1' AND nearest_tss_distance_bp = 98
      AND in_any_intron AND primary_transcript_region = 'intron'
) THEN error('anchor summary lost tandem or intron/TSS context') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73_pair_feature
    WHERE start = 140 AND strand = '+'
      AND pair_class = 'tandem_same_orientation'
      AND n_tandem_tp73_partner_loci = 1
      AND n_same_orientation_partner_loci = 1
      AND n_opposite_orientation_partner_loci = 0
      AND nearest_tandem_inter_motif_gap_bp = 2
      AND best_partner_score = 8.0
      AND best_pair_min_score = 8.0
      AND best_pair_sum_score = 18.0
) THEN error('same-orientation pair features are incorrect') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73_pair_feature
    WHERE start = 200 AND strand = '+'
      AND pair_class = 'tandem_opposite_orientation'
      AND n_opposite_orientation_partner_loci = 1
) THEN error('opposite-orientation pair class is missing') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73_pair_feature
    WHERE start = 158 AND strand = '+'
      AND pair_class = 'tandem_orientation_ambiguous'
      AND n_tandem_tp73_partner_loci = 1
      AND n_ambiguous_orientation_partner_loci = 1
) THEN error('dual-strand partner locus was not collapsed as ambiguous') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73_pair_feature
    WHERE start = 400 AND strand = '+'
      AND pair_class = 'singleton'
      AND n_tandem_tp73_partner_loci = 0
      AND best_partner_score IS NULL
      AND best_pair_min_score IS NULL
) THEN error('singleton pair features are incorrect') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_transcript_context
    WHERE anchor_start = 140 AND transcript_id = 'T1'
      AND fully_within_intron AND transcript_region = 'intron'
) THEN error('per-transcript intron classification is missing') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_transcript_context_pair
    WHERE anchor_start = 140 AND anchor_strand = '+'
      AND neighbor_motif_id = 'MA0079.5' AND neighbor_start = 130
      AND transcript_id = 'T1'
      AND anchor_signed_tss_distance_bp = 98
      AND neighbor_signed_tss_distance_bp = 84.5
      AND transcript_oriented_center_distance_bp = -13.5
      AND transcript_oriented_side = 'upstream'
) THEN error('transcript-oriented cofactor direction is incorrect') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_run_config
    WHERE schema_version = 5
      AND builder_source_commit = 'unknown'
      AND genome_id = 'synthetic_grch38_v1'
      AND motif_set_id = 'synthetic_jaspar2026'
      AND anchor_minimum_score = 0 AND partner_minimum_score = 0
      AND anchor_selection_mode = 'threshold'
      AND anchor_local_peak_flank_bp = 150
      AND anchor_local_peak_rule =
          'physical_locus_best_score_no_stronger_start_within_flank'
      AND tandem_score_rule =
          'both_orientation_specific_scores_at_or_above_tandem_minimum'
      AND capture_flank_bp = 150 AND context_flank_bp = 150
      AND tandem_flank_bp = 20 AND cofactor_pair_flank_bp = 150
      AND output_tier = 'selected'
      AND capture_geometry = 'interval'
      AND distance_metric = 'signed_interval_edge_distance'
      AND center_distance_metric = 'motif_center'
      AND capture_prefilter_center_bp = 206
      AND cofactor_pair_prefilter_center_bp = 230
      AND observed_max_anchor_span_bp = 16
      AND observed_max_neighbor_span_bp = 40
      AND tandem_distance_metric = 'nonoverlapping_edge_gap'
      AND partner_locus_identity_rule = 'same_alignment_span_collapses_orientation_records'
      AND cofactor_pair_scope = 'at_least_one_member_is_a_tp73_context_locus'
      AND cofactor_motif_locus_scope = 'tp73_context_loci_plus_their_pair_partners'
      AND cofactor_locus_pair_feature_scope = 'tp73_context_loci_only'
      AND default_neighbor_minimum_score = 0
      AND neighbor_qualifying_threshold_rule =
          'source_minimum_score_else_default_inclusive'
      AND anchor_motif_band_feature_grain =
          'physical_anchor_locus_neighbor_motif_interval_distance_band'
      AND anchor_motif_band_winner_rule =
          'highest_score_then_nearest_center_then_coordinates'
      AND cofactor_pair_score_rule =
          'both_locus_best_scores_at_or_above_neighbor_qualifying_threshold'
) THEN error('context run provenance is incomplete') END;

SELECT CASE WHEN
       (SELECT schema_version FROM context_run_config)
           <> (SELECT schema_version FROM motif_context_run_config)
    OR (SELECT COUNT(*) FROM context_run_config)
           <> (SELECT COUNT(*) FROM motif_context_run_config)
    THEN error('context run-config compatibility alias is inconsistent') END;

-- Hive infers an all-numeric chromosome partition as BIGINT unless told
-- otherwise. Lock the package contract to VARCHAR so chr1-only output can be
-- combined with X/Y and with packages from other species.
SELECT CASE WHEN
       (SELECT TYPEOF(chrom) FROM motif_context_pair LIMIT 1) <> 'VARCHAR'
    OR (SELECT TYPEOF(chrom) FROM tp73_anchor_locus LIMIT 1) <> 'VARCHAR'
    OR (SELECT TYPEOF(chrom) FROM tp73_motif_context_summary LIMIT 1) <> 'VARCHAR'
    OR (SELECT TYPEOF(chrom) FROM anchor_motif_band_feature LIMIT 1) <> 'VARCHAR'
    OR (SELECT TYPEOF(chrom) FROM tp73_pair_feature LIMIT 1) <> 'VARCHAR'
    OR (SELECT TYPEOF(chrom) FROM tp73_context_anchor LIMIT 1) <> 'VARCHAR'
    THEN error('partitioned context views did not preserve chromosome text') END;

SELECT CASE WHEN EXISTS (
    SELECT 1 FROM tp73_pair_feature
    WHERE genome_id <> 'synthetic_grch38_v1'
       OR motif_set_id <> 'synthetic_jaspar2026'
       OR anchor_minimum_score <> 0
       OR partner_minimum_score <> 0
) THEN error('pair feature lost genome identity or score floors') END;

-- The pair model joins on anchor_hit_id, so a context package must contain
-- exactly one feature row for each configuration-aware anchor identity.
SELECT CASE WHEN EXISTS (
    SELECT anchor_hit_id
    FROM tp73_pair_feature
    GROUP BY anchor_hit_id
    HAVING COUNT(*) <> 1
) THEN error('tp73_pair_feature anchor_hit_id is not unique') END;

-- Partner-locus orientation classes form an exhaustive, exclusive partition.
-- Counts are added here; subtracting them would invert the intended invariant.
SELECT CASE WHEN EXISTS (
    SELECT 1
    FROM tp73_pair_feature
    WHERE n_tandem_tp73_partner_loci < 0
       OR n_same_orientation_partner_loci < 0
       OR n_opposite_orientation_partner_loci < 0
       OR n_ambiguous_orientation_partner_loci < 0
       OR n_tandem_tp73_partner_loci <>
            n_same_orientation_partner_loci
            + n_opposite_orientation_partner_loci
            + n_ambiguous_orientation_partner_loci
       OR pair_class NOT IN (
            'singleton',
            'tandem_same_orientation',
            'tandem_opposite_orientation',
            'tandem_mixed_orientation',
            'tandem_orientation_ambiguous'
       )
       OR (pair_class = 'singleton') <>
            (n_tandem_tp73_partner_loci = 0)
       OR (pair_class = 'tandem_orientation_ambiguous') <>
            (n_tandem_tp73_partner_loci > 0
             AND n_ambiguous_orientation_partner_loci > 0)
       OR (pair_class = 'tandem_mixed_orientation') <>
            (n_ambiguous_orientation_partner_loci = 0
             AND n_same_orientation_partner_loci > 0
             AND n_opposite_orientation_partner_loci > 0)
       OR (pair_class = 'tandem_same_orientation') <>
            (n_ambiguous_orientation_partner_loci = 0
             AND n_same_orientation_partner_loci > 0
             AND n_opposite_orientation_partner_loci = 0)
       OR (pair_class = 'tandem_opposite_orientation') <>
            (n_ambiguous_orientation_partner_loci = 0
             AND n_same_orientation_partner_loci = 0
             AND n_opposite_orientation_partner_loci > 0)
       OR has_multiple_tandem_partner_loci <>
            (n_tandem_tp73_partner_loci > 1)
) THEN error('pair_class / partner-locus counts are inconsistent') END;

-- NULL means that no partner exists. Models and exports rely on zero partners
-- remaining distinguishable from a partner whose measured score is zero.
SELECT CASE WHEN EXISTS (
    SELECT 1
    FROM tp73_pair_feature
    WHERE (pair_class = 'singleton' AND (
              best_partner_score IS NOT NULL
           OR best_partner_pwm_relative_score IS NOT NULL
           OR best_pair_min_score IS NOT NULL
           OR best_pair_sum_score IS NOT NULL
           OR best_pair_min_pwm_relative_score IS NOT NULL
          ))
       OR (pair_class <> 'singleton' AND best_partner_score IS NULL)
) THEN error('singleton/non-singleton partner-score nullability is inconsistent') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73_context_pair_feature
    WHERE anchor_start = 140 AND anchor_strand = '+'
      AND neighbor_motif_id = 'MA0079.5'
      AND anchor_pair_class = 'tandem_same_orientation'
      AND pair_relation = 'heterotypic_context'
      AND oriented_distance_bin_start_bp = -15.0
) THEN error('chromosome-wide pair-stratified context view is incorrect') END;
SQL

{
    echo "PREPARE q10 AS"
    awk '
        /^-- Q10\./ { capture = 1; next }
        /^-- Q11\./ { capture = 0 }
        capture { print }
    ' "$repository_root/sql/queries.sql"
    cat <<'SQL'
EXECUTE q10(
    genome_id := 'synthetic_grch38_v1', motif_set_id := 'synthetic_jaspar2026',
    chrom := '1', start := 140, strand := '+',
    score_mode := 'log2_relative_risk', pseudocount := 1.0
);
PREPARE q11 AS
SQL
    awk '
        /^-- Q11\./ { capture = 1; next }
        /^-- Q12\./ { capture = 0 }
        capture { print }
    ' "$repository_root/sql/queries.sql"
    cat <<'SQL'
EXECUTE q11(
    genome_id := 'synthetic_grch38_v1', motif_set_id := 'synthetic_jaspar2026',
    gene_name := 'GENE1', score_mode := 'log2_relative_risk', pseudocount := 1.0
);
PREPARE q14 AS
SQL
    awk '
        /^-- Q14\./ { capture = 1; next }
        /^-- Q15\./ { capture = 0 }
        capture { print }
    ' "$repository_root/sql/queries.sql"
    cat <<'SQL'
EXECUTE q14(
    genome_id := 'synthetic_grch38_v1', motif_set_id := 'synthetic_jaspar2026',
    chrom := '1', start := 140, strand := '+',
    neighbor_motif_id := 'MA0079.5',
    score_mode := 'log2_relative_risk', pseudocount := 1.0
);
PREPARE q15 AS
SQL
    awk '
        /^-- Q15\./ { capture = 1; next }
        /^-- Q16\./ { capture = 0 }
        capture { print }
    ' "$repository_root/sql/queries.sql"
    cat <<'SQL'
EXECUTE q15(
    genome_id := 'synthetic_grch38_v1', motif_set_id := 'synthetic_jaspar2026',
    chrom := '1', start := 400, strand := '+',
    neighbor_motif_id := 'MPAIR.1',
    score_mode := 'log2_relative_risk', pseudocount := 1.0
);
PREPARE q19 AS
SQL
    awk '
        /^-- Q19\./ { capture = 1 }
        capture { print }
    ' "$repository_root/sql/queries.sql"
    cat <<'SQL'
EXECUTE q19(
    genome_id := 'synthetic_grch38_v1', motif_set_id := 'synthetic_jaspar2026',
    chrom := '1', anchor_start := 400, anchor_end := 416,
    score_mode := 'log2_relative_risk', pseudocount := 1.0
);
SQL
} >> "$temporary/assertions.sql"

(
    cd "$temporary/context_package"
    duckdb context.duckdb < "$temporary/assertions.sql" >/dev/null
)

for partner_floor in 0 -1; do
    package="$temporary/context_partner_floor_${partner_floor//-/_minus_}"
    "$repository_root/scripts/build_motif_context.py" \
        --motif-hits "$temporary/motif_hit.parquet" \
        --output "$package" --anchor-motif MA0861.2 \
        --motif-set-id synthetic_jaspar2026 --genome-id synthetic_grch38_v1 \
        --anchor-minimum-score -5 --tandem-minimum-score "$partner_floor" \
        --anchor-selection-mode threshold \
        --score-mode log2_relative_risk --pseudocount 1 --chrom 1 \
        --capture-flank 100 --context-flank 50 --tandem-flank 20 \
        --memory-limit 1GB --max-temp-size 1GB >/dev/null
done

(
    cd "$temporary/context_partner_floor_0"
    duckdb context.duckdb -bail -c "
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM tp73_pair_feature
            WHERE start = 600 AND strand = '+' AND pair_class = 'singleton'
        ) THEN error('tandem floor 0 did not retain the negative anchor as a singleton') END;" \
        >/dev/null
)
(
    cd "$temporary/context_partner_floor__minus_1"
    duckdb context.duckdb -bail -c "
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM tp73_pair_feature
            WHERE start = 600 AND strand = '+' AND pair_class = 'singleton'
        ) THEN error('tandem floor -1 ignored the anchor score below its floor') END;
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM tp73_pair_feature
            WHERE start = 620 AND strand = '+' AND pair_class = 'singleton'
        ) THEN error('tandem floor -1 ignored the partner score below its floor') END;" \
        >/dev/null
)

"$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/motif_hit.parquet" \
    --output "$temporary/context_local_peak" \
    --anchor-motif MA0861.2 \
    --motif-set-id synthetic_jaspar2026 --genome-id synthetic_grch38_v1 \
    --anchor-minimum-score -5 --tandem-minimum-score 0 \
    --anchor-selection-mode local_peak --anchor-local-peak-flank 150 \
    --score-mode log2_relative_risk --pseudocount 1 --chrom 1 \
    --capture-flank 150 --context-flank 150 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null

(
    cd "$temporary/context_local_peak"
    duckdb context.duckdb -bail -c "
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM tp73_pair_feature WHERE start = 600
        ) THEN error('negative non-peak survived anchor selection') END;
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM tp73_pair_feature
            WHERE start = 620 AND strand = '+'
              AND anchor_selection_class = 'local_peak'
              AND anchor_locus_best_score = -0.5
              AND best_other_anchor_locus_score = -2.0
              AND anchor_locus_score_prominence = 1.5
              AND anchor_locus_is_local_peak
              AND pair_class = 'singleton'
        ) THEN error('negative regional peak provenance is incorrect') END;
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM motif_context_pair
            WHERE anchor_start = 620 AND neighbor_start = 600
              AND neighbor_score = -2.0
        ) THEN error('rejected anchor was lost from retained-neighbor context') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM tp73_pair_feature WHERE start = 145
        ) THEN error('subordinate nonnegative locus survived local-peak selection') END;
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM tp73_pair_feature
            WHERE start = 140 AND strand = '+' AND score = 10.0
              AND anchor_selection_class = 'local_peak'
        ) THEN error('strongest nonnegative regional peak was not retained') END;" \
        >/dev/null
)

"$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/motif_hit.parquet" \
    --output "$temporary/context_without_gtf" \
    --anchor-motif MA0861.2 \
    --motif-set-id synthetic_jaspar2026 --genome-id synthetic_grch38_v1 \
    --anchor-minimum-score 0 --partner-minimum-score 0 \
    --score-mode log2_relative_risk --pseudocount 1 \
    --chrom 1 --capture-flank 100 --context-flank 50 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null

if "$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/motif_hit.parquet" \
    --output "$temporary/context_without_gtf" \
    --anchor-motif MA0861.2 \
    --motif-set-id synthetic_jaspar2026 --genome-id synthetic_grch38_v1 \
    --anchor-minimum-score 0 --partner-minimum-score 0 \
    --score-mode log2_relative_risk --pseudocount 1 \
    --chrom 1 --capture-flank 100 --context-flank 50 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null 2>&1; then
    echo "E: motif-context build replaced an output without --force." >&2
    exit 1
fi

"$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/motif_hit.parquet" \
    --output "$temporary/context_without_gtf" --force \
    --anchor-motif MA0861.2 \
    --motif-set-id synthetic_jaspar2026 --genome-id synthetic_grch38_v1 \
    --anchor-minimum-score 0 --partner-minimum-score 0 \
    --score-mode log2_relative_risk --pseudocount 1 \
    --chrom 1 --capture-flank 100 --context-flank 50 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null

(
    cd "$temporary/context_without_gtf"
    duckdb context.duckdb -c \
        "SELECT CASE WHEN EXISTS (
             SELECT 1 FROM tp73_context_anchor
             WHERE gene_annotation_available
                OR primary_transcript_region <> 'not_assessed'
         ) THEN error('missing GTF was represented as negative annotation') END;" \
        >/dev/null
)

if find "$temporary/context_package" -type f \( -name '*.tsv' -o -name '*.bed' \) -print -quit | grep -q .; then
    echo "E: motif-context build created an intermediate BED/TSV product." >&2
    exit 1
fi

"$repository_root/scripts/build_motif_context.py" --help | grep -q '^usage:'

echo "Motif-context Parquet tests passed."

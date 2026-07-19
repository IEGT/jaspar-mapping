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
    SELECT * FROM (VALUES
        ('1', 140::BIGINT, 156::BIGINT, 'MA0861.2', 'TP73', '+', 10.0, 'log2_relative_risk', 1.0, 0.95),
        ('1', 140::BIGINT, 156::BIGINT, 'MA0861.2', 'TP73', '-', 9.5,  'log2_relative_risk', 1.0, 0.94),
        ('1', 145::BIGINT, 161::BIGINT, 'MA0861.2', 'TP73', '+', 6.5,  'log2_relative_risk', 1.0, 0.86),
        ('1', 158::BIGINT, 174::BIGINT, 'MA0861.2', 'TP73', '+', 8.0,  'log2_relative_risk', 1.0, 0.90),
        ('1', 200::BIGINT, 216::BIGINT, 'MA0861.2', 'TP73', '+', 7.0,  'log2_relative_risk', 1.0, 0.88),
        ('1', 220::BIGINT, 236::BIGINT, 'MA0861.2', 'TP73', '-', 9.0,  'log2_relative_risk', 1.0, 0.92),
        ('1', 400::BIGINT, 416::BIGINT, 'MA0861.2', 'TP73', '+', 5.0,  'log2_relative_risk', 1.0, 0.81),
        ('1', 600::BIGINT, 616::BIGINT, 'MA0861.2', 'TP73', '+', -2.0, 'log2_relative_risk', 1.0, 0.30),
        ('1', 620::BIGINT, 636::BIGINT, 'MA0861.2', 'TP73', '+', -0.5, 'log2_relative_risk', 1.0, 0.45),
        ('1', 130::BIGINT, 139::BIGINT, 'MA0079.5', 'SP1',  '+', 6.0,  'log2_relative_risk', 1.0, 0.85),
        ('1', 410::BIGINT, 419::BIGINT, 'MA0079.5', 'SP1',  '+', 4.0,  'log2_relative_risk', 1.0, 0.80),
        ('1', 499::BIGINT, 508::BIGINT, 'MA0079.5', 'SP1',  '+', 3.0,  'log2_relative_risk', 1.0, 0.70)
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
    --chrom 1 --capture-flank 100 --context-flank 50 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null

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

-- Anchor center 408 is in capture bin 4; neighbor center 503.5 is in bin 5.
-- This guards the +/-1-bin expansion used to make the bounded join lossless.
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_context_pair
    WHERE anchor_start = 400 AND anchor_strand = '+' AND neighbor_start = 499
      AND absolute_center_distance_bp = 95.5 AND capture_flank_bp = 100
) THEN error('valid pair across a capture-bin boundary was dropped') END;

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
    SELECT 1 FROM context_run_config
    WHERE schema_version = 3
      AND genome_id = 'synthetic_grch38_v1'
      AND motif_set_id = 'synthetic_jaspar2026'
      AND anchor_minimum_score = 0 AND partner_minimum_score = 0
      AND capture_flank_bp = 100 AND context_flank_bp = 50
      AND tandem_flank_bp = 20 AND distance_metric = 'motif_center'
      AND tandem_distance_metric = 'nonoverlapping_edge_gap'
      AND partner_locus_identity_rule = 'same_alignment_span_collapses_orientation_records'
) THEN error('context run provenance is incomplete') END;

SELECT CASE WHEN EXISTS (
    SELECT 1 FROM tp73_pair_feature
    WHERE genome_id <> 'synthetic_grch38_v1'
       OR motif_set_id <> 'synthetic_jaspar2026'
       OR anchor_minimum_score <> 0
       OR partner_minimum_score <> 0
) THEN error('pair feature lost genome identity or score floors') END;

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
        capture { print }
    ' "$repository_root/sql/queries.sql"
    cat <<'SQL'
EXECUTE q14(
    genome_id := 'synthetic_grch38_v1', motif_set_id := 'synthetic_jaspar2026',
    chrom := '1', start := 140, strand := '+',
    neighbor_motif_id := 'MA0079.5',
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
        --anchor-minimum-score -5 --partner-minimum-score "$partner_floor" \
        --score-mode log2_relative_risk --pseudocount 1 --chrom 1 \
        --capture-flank 100 --context-flank 50 --tandem-flank 20 \
        --memory-limit 1GB --max-temp-size 1GB >/dev/null
done

(
    cd "$temporary/context_partner_floor_0"
    duckdb context.duckdb -bail -c "
        SELECT CASE WHEN (SELECT pair_class FROM tp73_pair_feature
                          WHERE start = 600 AND strand = '+') <> 'singleton'
            THEN error('partner floor 0 retained a -0.5 tandem partner') END;" \
        >/dev/null
)
(
    cd "$temporary/context_partner_floor__minus_1"
    duckdb context.duckdb -bail -c "
        SELECT CASE WHEN (SELECT pair_class FROM tp73_pair_feature
                          WHERE start = 600 AND strand = '+') <> 'tandem_same_orientation'
            THEN error('partner floor -1 dropped a -0.5 tandem partner') END;" \
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

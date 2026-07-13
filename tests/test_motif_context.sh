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
        ('1', 130::BIGINT, 139::BIGINT, 'MA0079.5', 'SP1',  '+', 6.0,  'log2_relative_risk', 1.0, 0.85)
    ) AS t(chrom, start, "end", motif_id, motif_name, strand, score,
           score_mode, pseudocount, pwm_relative_score)
) TO 'motif_hit.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
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
    --score-mode log2_relative_risk --pseudocount 1 \
    --chrom 1 --capture-flank 100 --context-flank 50 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null

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

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73_context_anchor
    WHERE start = 140 AND strand = '+'
      AND has_tandem_tp73 AND n_tandem_tp73_partners = 1
      AND nearest_tandem_oriented_distance_bp = 18
      AND nearest_tandem_relative_orientation = 'same'
      AND nearest_gene_id = 'G1' AND nearest_tss_distance_bp = 98
      AND in_any_intron AND primary_transcript_region = 'intron'
) THEN error('anchor summary lost tandem or intron/TSS context') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM motif_transcript_context
    WHERE anchor_start = 140 AND transcript_id = 'T1'
      AND fully_within_intron AND transcript_region = 'intron'
) THEN error('per-transcript intron classification is missing') END;

SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM context_run_config
    WHERE capture_flank_bp = 100 AND context_flank_bp = 50
      AND tandem_flank_bp = 20 AND distance_metric = 'motif_center'
      AND tandem_distance_metric = 'nonoverlapping_edge_gap'
) THEN error('context run provenance is incomplete') END;
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
    chrom := '1', start := 140, strand := '+',
    score_mode := 'log2_relative_risk', pseudocount := 1.0
);
PREPARE q11 AS
SQL
    awk '
        /^-- Q11\./ { capture = 1; next }
        capture { print }
    ' "$repository_root/sql/queries.sql"
    cat <<'SQL'
EXECUTE q11(
    gene_name := 'GENE1', score_mode := 'log2_relative_risk', pseudocount := 1.0
);
SQL
} >> "$temporary/assertions.sql"

(
    cd "$temporary/context_package"
    duckdb context.duckdb < "$temporary/assertions.sql" >/dev/null
)

"$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/motif_hit.parquet" \
    --output "$temporary/context_without_gtf" \
    --anchor-motif MA0861.2 \
    --score-mode log2_relative_risk --pseudocount 1 \
    --chrom 1 --capture-flank 100 --context-flank 50 --tandem-flank 20 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null

if "$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/motif_hit.parquet" \
    --output "$temporary/context_without_gtf" \
    --anchor-motif MA0861.2 \
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

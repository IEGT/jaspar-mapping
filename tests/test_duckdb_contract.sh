#!/usr/bin/env bash

set -euo pipefail

if ! command -v duckdb >/dev/null 2>&1; then
    echo "E: duckdb CLI is required for this test." >&2
    exit 1
fi

repository_root=$(cd "$(dirname "$0")/.." && pwd)
test_sql=$(mktemp "${TMPDIR:-/tmp}/jaspar-duckdb-contract.XXXXXX.sql")
trap 'rm -f "$test_sql"' EXIT HUP INT TERM

cat > "$test_sql" <<'SQL'
CREATE TABLE promoter (
    genome_id VARCHAR,
    gene_id VARCHAR,
    transcript_id VARCHAR,
    chrom VARCHAR,
    strand VARCHAR,
    promoter_start BIGINT,
    promoter_end BIGINT,
    tss BIGINT
);
INSERT INTO promoter VALUES
    ('genome1', 'G1', 'T1', '1', '+', 0, 100, 50),
    ('genome1', 'G1', 'T2', '1', '+', 0, 100, 60),
    ('genome1', 'G2', 'T3', '1', '+', 200, 300, 250);

CREATE TABLE motif_hit (
    genome_id VARCHAR,
    motif_set_id VARCHAR,
    chrom VARCHAR,
    start BIGINT,
    "end" BIGINT,
    motif_id VARCHAR,
    strand VARCHAR,
    score DOUBLE,
    score_mode VARCHAR,
    pseudocount DOUBLE,
    pwm_relative_score DOUBLE
);
INSERT INTO motif_hit VALUES
    ('genome1', 'motifs1', '1', 55, 60, 'M1', '+', 1.0, 'log_odds', 1.0, 0.90),
    ('genome1', 'motifs1', '1', 70, 75, 'M2', '+', 2.0, 'log_odds', 1.0, 0.85),
    ('genome1', 'motifs1', '1', 55, 60, 'M1', '+', 3.0, 'log2_relative_risk', 1.0, 0.88);

CREATE TABLE motif_architecture (
    motif_set_id VARCHAR,
    motif_id VARCHAR,
    tf_family VARCHAR,
    binding_unit_model VARCHAR
);
INSERT INTO motif_architecture VALUES
    ('motifs1', 'M1', 'FAMILY', 'dimer_of_dimers'),
    ('motifs1', 'M2', 'FAMILY', 'monomer');

CREATE TABLE tp73_pair_feature (
    anchor_hit_id VARCHAR,
    genome_id VARCHAR,
    motif_set_id VARCHAR,
    chrom VARCHAR,
    start BIGINT,
    "end" BIGINT,
    motif_id VARCHAR,
    strand VARCHAR,
    score DOUBLE,
    score_mode VARCHAR,
    pseudocount DOUBLE,
    pwm_relative_score DOUBLE,
    pair_class VARCHAR,
    n_tandem_tp73_partner_loci BIGINT,
    nearest_tandem_inter_motif_gap_bp BIGINT,
    best_partner_score DOUBLE,
    best_partner_pwm_relative_score DOUBLE,
    best_pair_min_score DOUBLE,
    best_pair_sum_score DOUBLE,
    best_pair_min_pwm_relative_score DOUBLE
);
INSERT INTO tp73_pair_feature VALUES
    ('H1_LO', 'genome1', 'motifs1', '1', 55, 60, 'M1', '+', 1.0, 'log_odds', 1.0, 0.90,
     'tandem_same_orientation', 1, 10, 2.0, 0.85, 1.0, 3.0, 0.85),
    ('H1_RR', 'genome1', 'motifs1', '1', 55, 60, 'M1', '+', 3.0, 'log2_relative_risk', 1.0, 0.88,
     'singleton', 0, NULL, NULL, NULL, NULL, NULL, NULL);

CREATE TABLE motif_context_pair (
    anchor_hit_id VARCHAR,
    genome_id VARCHAR,
    motif_set_id VARCHAR,
    neighbor_hit_id VARCHAR,
    neighbor_start BIGINT,
    neighbor_end BIGINT,
    neighbor_motif_id VARCHAR,
    neighbor_motif_name VARCHAR,
    anchor_score DOUBLE,
    anchor_pwm_relative_score DOUBLE,
    neighbor_score DOUBLE,
    neighbor_pwm_relative_score DOUBLE,
    anchor_oriented_center_distance_bp DOUBLE,
    absolute_center_distance_bp DOUBLE,
    relative_orientation VARCHAR,
    anchor_oriented_side VARCHAR,
    interval_overlap_bp BIGINT,
    within_context_flank BOOLEAN,
    same_anchor_motif_span BOOLEAN,
    is_tandem_tp73 BOOLEAN,
    score_mode VARCHAR,
    pseudocount DOUBLE
);
INSERT INTO motif_context_pair VALUES
    ('H1_LO', 'genome1', 'motifs1', 'H2_LO', 70, 75, 'M2', 'MOTIF2', 1.0, 0.90, 2.0, 0.85,
     15.0, 15.0, 'same', 'downstream', 0, true, false, false,
     'log_odds', 1.0);

CREATE TABLE expression_differential (
    genome_id VARCHAR,
    gene_id VARCHAR,
    cell_line VARCHAR,
    log2fc_ta_vs_dn DOUBLE,
    lfc_se DOUBLE,
    pvalue DOUBLE,
    padj DOUBLE,
    base_mean DOUBLE
);
INSERT INTO expression_differential VALUES
    ('genome1', 'G1', 'saos2', 1.5, 0.1, 0.001, 0.01, 100.0),
    ('genome1', 'G2', 'saos2', -0.5, 0.2, 0.05, 0.10, 50.0);

CREATE TABLE gene (
    genome_id VARCHAR,
    gene_id VARCHAR,
    gene_name VARCHAR,
    chrom VARCHAR,
    strand VARCHAR,
    tss BIGINT
);
INSERT INTO gene VALUES
    ('genome1', 'G1', 'GENE1', '1', '+', 50),
    ('genome1', 'G2', 'GENE2', '1', '+', 250);
SQL

# Exercise the actual materialized feature/card definitions rather than a copy,
# then extract Q3 through the Q4 marker with its placeholders unchanged.
{
    sed -n '/^CREATE OR REPLACE TABLE promoter_motif_hit AS/,$p' \
        "$repository_root/sql/schema.sql"

    cat <<'SQL'
SELECT CASE
    WHEN (SELECT COUNT(*) FROM promoter_motif_hit) <> 3
    THEN error('promoter_motif_hit multiplied hits by transcript count')
END;

SELECT CASE
    WHEN EXISTS (
        SELECT 1
        FROM promoter_motif_hit
        WHERE gene_id = 'G1' AND n_overlapping_transcripts <> 2
    )
    THEN error('transcript overlap provenance was not retained')
END;

SELECT CASE
    WHEN EXISTS (
        SELECT 1
        FROM promoter_motif_hit
        WHERE gene_id = 'G1'
          AND motif_id = 'M1'
          AND score_mode = 'log_odds'
          AND (closest_transcript_id <> 'T1' OR tss_distance <> 5)
    )
    THEN error('closest-transcript TSS distance is not deterministic')
END;

SELECT CASE
    WHEN (SELECT SUM(n_hits) FROM promoter_arch_feature WHERE gene_id = 'G1') <> 3
    THEN error('promoter architecture hit counts are inflated')
END;

SELECT CASE
    WHEN (SELECT COUNT(*) FROM promoter_card WHERE gene_id = 'G1') <> 2
    THEN error('promoter_card did not preserve scoring configurations')
END;

SELECT CASE
    WHEN (SELECT n_motif_hits FROM promoter_card
          WHERE gene_id = 'G1' AND score_mode = 'log_odds' AND pseudocount = 1.0) <> 2
    THEN error('promoter_card mixed or miscounted log-odds features')
END;

SELECT CASE
    WHEN EXISTS (
        SELECT 1
        FROM promoter_card
        WHERE gene_id = 'G2' AND n_motif_hits <> 0
    ) OR (SELECT COUNT(*) FROM promoter_card WHERE gene_id = 'G2') <> 2
    THEN error('zero-hit genes are missing configuration-specific cards')
END;

SELECT CASE
    WHEN NOT EXISTS (
        SELECT 1 FROM promoter_pair_feature
        WHERE gene_id = 'G1' AND score_mode = 'log_odds'
          AND pair_class = 'tandem_same_orientation'
          AND n_anchor_orientation_hits = 1
          AND n_tandem_partner_loci = 1
          AND best_pair_min_pwm_relative_score = 0.85
    ) OR NOT EXISTS (
        SELECT 1 FROM promoter_pair_feature
        WHERE gene_id = 'G1' AND score_mode = 'log2_relative_risk'
          AND pair_class = 'singleton'
          AND n_anchor_orientation_hits = 1
    )
    THEN error('promoter TP73 pair classes were not preserved')
END;

SELECT CASE
    WHEN NOT EXISTS (
        SELECT 1 FROM promoter_motif_pair_feature
        WHERE gene_id = 'G1' AND neighbor_motif_id = 'M2'
          AND anchor_pair_class = 'tandem_same_orientation'
          AND pair_relation = 'heterotypic_context'
          AND relative_orientation = 'same'
          AND oriented_distance_bin_start_bp = 15.0
          AND n_directed_orientation_pairs = 1
          AND max_pair_min_pwm_relative_score = 0.85
    )
    THEN error('pair-stratified neighboring-motif feature is incorrect')
END;

SELECT CASE
    WHEN (SELECT n_tandem_tp73_anchor_hits FROM promoter_card
          WHERE gene_id = 'G1' AND score_mode = 'log_odds') <> 1
      OR (SELECT n_singleton_tp73_anchor_hits FROM promoter_card
          WHERE gene_id = 'G1' AND score_mode = 'log2_relative_risk') <> 1
    THEN error('promoter card omitted TP73 pair-state counts')
END;

PREPARE q3 AS
SQL

    awk '
        /^-- Q3\./ { capture = 1 }
        /^-- Q4\./ { capture = 0 }
        capture { print }
    ' "$repository_root/sql/queries.sql"

    cat <<'SQL'
EXECUTE q3(
    genome_id := 'genome1',
    motif_set_id := 'motifs1',
    cell_line := 'saos2',
    score_mode := 'log_odds',
    pseudocount := 1.0
);
SQL

    echo "PREPARE q12 AS"
    awk '
        /^-- Q12\./ { capture = 1 }
        /^-- Q13\./ { capture = 0 }
        capture { print }
    ' "$repository_root/sql/queries.sql"
    cat <<'SQL'
EXECUTE q12(
    genome_id := 'genome1',
    motif_set_id := 'motifs1',
    gene_name := 'GENE1',
    score_mode := 'log_odds',
    pseudocount := 1.0
);
PREPARE q13 AS
SQL
    awk '
        /^-- Q13\./ { capture = 1 }
        /^-- Q14\./ { capture = 0 }
        capture { print }
    ' "$repository_root/sql/queries.sql"
    cat <<'SQL'
EXECUTE q13(
    genome_id := 'genome1',
    motif_set_id := 'motifs1',
    gene_name := 'GENE1',
    neighbor_motif_id := 'M2',
    score_mode := 'log_odds',
    pseudocount := 1.0
);
SQL
} >> "$test_sql"

duckdb :memory: < "$test_sql" >/dev/null
echo "DuckDB promoter/ML contract tests passed."

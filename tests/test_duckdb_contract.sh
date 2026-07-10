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
    gene_id VARCHAR,
    transcript_id VARCHAR,
    chrom VARCHAR,
    strand VARCHAR,
    promoter_start BIGINT,
    promoter_end BIGINT,
    tss BIGINT
);
INSERT INTO promoter VALUES
    ('G1', 'T1', '1', '+', 0, 100, 50),
    ('G1', 'T2', '1', '+', 0, 100, 60),
    ('G2', 'T3', '1', '+', 200, 300, 250);

CREATE TABLE motif_hit (
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
    ('1', 55, 60, 'M1', '+', 1.0, 'log_odds', 1.0, 0.90),
    ('1', 70, 75, 'M2', '+', 2.0, 'log_odds', 1.0, 0.85),
    ('1', 55, 60, 'M1', '+', 3.0, 'log2_relative_risk', 1.0, 0.88);

CREATE TABLE motif_architecture (
    motif_id VARCHAR,
    tf_family VARCHAR,
    binding_unit_model VARCHAR
);
INSERT INTO motif_architecture VALUES
    ('M1', 'FAMILY', 'dimer_of_dimers'),
    ('M2', 'FAMILY', 'monomer');

CREATE TABLE expression_differential (
    gene_id VARCHAR,
    cell_line VARCHAR,
    log2fc_ta_vs_dn DOUBLE,
    lfc_se DOUBLE,
    pvalue DOUBLE,
    padj DOUBLE,
    base_mean DOUBLE
);
INSERT INTO expression_differential VALUES
    ('G1', 'saos2', 1.5, 0.1, 0.001, 0.01, 100.0),
    ('G2', 'saos2', -0.5, 0.2, 0.05, 0.10, 50.0);

CREATE TABLE gene (
    gene_id VARCHAR,
    gene_name VARCHAR,
    chrom VARCHAR,
    strand VARCHAR,
    tss BIGINT
);
INSERT INTO gene VALUES
    ('G1', 'GENE1', '1', '+', 50),
    ('G2', 'GENE2', '1', '+', 250);
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

PREPARE q3 AS
SQL

    awk '
        /^-- Q3\./ { capture = 1 }
        /^-- Q4\./ { capture = 0 }
        capture { print }
    ' "$repository_root/sql/queries.sql"

    cat <<'SQL'
EXECUTE q3(
    cell_line := 'saos2',
    score_mode := 'log_odds',
    pseudocount := 1.0
);
SQL
} >> "$test_sql"

duckdb :memory: < "$test_sql" >/dev/null
echo "DuckDB promoter/ML contract tests passed."

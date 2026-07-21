#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-context-maxima.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping sparse context-maximum test." >&2
    exit 0
}

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT * FROM (VALUES
        ('1', 100::BIGINT, 110::BIGINT),
        ('1', 500::BIGINT, 510::BIGINT)
    ) AS t(chrom, anchor_start, anchor_end)
) TO '$temporary/anchors.parquet' (FORMAT PARQUET);

COPY (
    SELECT * FROM (VALUES
        (50::BIGINT, 60::BIGINT, 5.5::FLOAT),
        (260::BIGINT, 270::BIGINT, 5.0::FLOAT),
        (661::BIGINT, 671::BIGINT, 9.0::FLOAT)
    ) AS t(start, "end", score)
) TO '$temporary/m1-plus.parquet' (FORMAT PARQUET);

COPY (
    SELECT * FROM (VALUES
        (50::BIGINT, 60::BIGINT, 6.0::FLOAT),
        (340::BIGINT, 350::BIGINT, 4.0::FLOAT)
    ) AS t(start, "end", score)
) TO '$temporary/m1-minus.parquet' (FORMAT PARQUET);

COPY (
    SELECT * FROM (VALUES
        (260::BIGINT, 300::BIGINT, 7.0::FLOAT)
    ) AS t(start, "end", score)
) TO '$temporary/m2-plus.parquet' (FORMAT PARQUET);

COPY (
    SELECT * FROM (VALUES
        (650::BIGINT, 660::BIGINT, 8.0::FLOAT)
    ) AS t(start, "end", score)
) TO '$temporary/m2-minus.parquet' (FORMAT PARQUET);

COPY (
    SELECT * FROM (VALUES
        ('M1', 'MOTIF1', 'synthetic_v1', 'genome1', 'motifs1',
         'tp73_context_binary_feature', 'MA0861.2',
         'log2_relative_risk', 1.0::DOUBLE, 'uniform_acgt_v1',
         'additive_per_base', 'all_tp73_anchors', 5.0::DOUBLE,
         'exploratory_positive_gain', NULL::BIGINT, 50::BIGINT, 'any',
         true, 'signed_interval_edge_distance', -1.0::DOUBLE),
        ('M2', 'MOTIF2', 'synthetic_v1', 'genome1', 'motifs1',
         'tp73_context_binary_feature', 'MA0861.2',
         'log2_relative_risk', 1.0::DOUBLE, 'uniform_acgt_v1',
         'additive_per_base', 'all_tp73_anchors', 8.0::DOUBLE,
         'exploratory_positive_gain', NULL::BIGINT, 150::BIGINT, 'any',
         true, 'signed_interval_edge_distance', -1.0::DOUBLE)
    ) AS t(
        motif_id, motif_name, threshold_set_id, genome_id, motif_set_id,
        threshold_role, target_motif_id, score_mode, pseudocount,
        background_model_id, pseudocount_scheme, calibration_stratum_id,
        recommended_threshold, calibration_status,
        context_min_interval_distance_bp, context_max_interval_distance_bp,
        context_relation_filter, threshold_inclusive,
        context_distance_metric, source_minimum_score
    )
) TO '$temporary/thresholds.parquet' (FORMAT PARQUET);
SQL

"$repository_root/scripts/build_sparse_context_maxima.py" \
    --anchor-parquet "$temporary/anchors.parquet" \
    --cofactor M1 "$temporary/m1-plus.parquet" "$temporary/m1-minus.parquet" \
    --cofactor M2 "$temporary/m2-plus.parquet" "$temporary/m2-minus.parquet" \
    --output "$temporary/maxima.parquet" \
    --flank 150 --source-score-floor -1 \
    --duckdb "$duckdb" --threads 1 --memory-limit 256MB \
    --temp-directory "$temporary"

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW maxima AS SELECT * FROM read_parquet('$temporary/maxima.parquet');

SELECT CASE WHEN (SELECT COUNT(*) FROM maxima) <> 4
    THEN error('sparse context maxima are not rectangular') END;
SELECT CASE WHEN (SELECT COUNT(DISTINCT (chrom, anchor_start, motif_id))
                  FROM maxima) <> 4
    THEN error('sparse context maxima keys are not unique') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM maxima
    WHERE anchor_start = 100 AND motif_id = 'M1' AND context_score = 6
) THEN error('highest M1 score was not selected') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM maxima
    WHERE anchor_start = 500 AND motif_id = 'M1' AND context_score = 4
) THEN error('a distance-151 M1 hit was admitted') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM maxima
    WHERE anchor_start = 100 AND motif_id = 'M2' AND context_score = 7
) THEN error('wide interval-near M2 hit was lost by the center prefilter') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM maxima
    WHERE anchor_start = 500 AND motif_id = 'M2' AND context_score = 8
) THEN error('second M2 anchor maximum is incorrect') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM maxima
    WHERE capture_prefilter_center_bp <> 200
       OR observed_max_anchor_span_bp <> 10
       OR observed_max_context_span_bp <> 40
       OR context_distance_metric <> 'signed_interval_edge_distance'
) THEN error('derived interval-prefilter provenance is incorrect') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM maxima
    WHERE n_neighbor_loci_above_threshold IS NOT NULL
       OR has_neighbor_locus_above_threshold IS NOT NULL
) THEN error('maximum-only output fabricated thresholded counts') END;
SQL

"$repository_root/scripts/build_sparse_context_maxima.py" \
    --anchor-parquet "$temporary/anchors.parquet" \
    --cofactor M1 "$temporary/m1-plus.parquet" "$temporary/m1-minus.parquet" \
    --cofactor M2 "$temporary/m2-plus.parquet" "$temporary/m2-minus.parquet" \
    --threshold-parquet "$temporary/thresholds.parquet" \
    --threshold-set-id synthetic_v1 \
    --output "$temporary/threshold-counts.parquet" \
    --flank 150 --source-score-floor -1 \
    --duckdb "$duckdb" --threads 1 --memory-limit 256MB \
    --temp-directory "$temporary"

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW counts AS
SELECT * FROM read_parquet('$temporary/threshold-counts.parquet');

SELECT CASE WHEN (SELECT COUNT(*) FROM counts) <> 4
    THEN error('threshold counts are not zero-complete') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM counts
    WHERE anchor_start = 100 AND motif_id = 'M1'
      AND recommended_threshold = 5
      AND n_neighbor_loci_above_threshold = 1
      AND has_neighbor_locus_above_threshold
) THEN error('same-span strand alternatives inflated the M1 locus count') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM counts
    WHERE anchor_start = 500 AND motif_id = 'M1'
      AND n_neighbor_loci_above_threshold = 0
      AND NOT has_neighbor_locus_above_threshold
) THEN error('an explicit zero threshold count was omitted') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM counts
    WHERE anchor_start = 100 AND motif_id = 'M2'
      AND n_neighbor_loci_above_threshold = 0
) OR NOT EXISTS (
    SELECT 1 FROM counts
    WHERE anchor_start = 500 AND motif_id = 'M2'
      AND n_neighbor_loci_above_threshold = 1
) THEN error('per-motif score or interval thresholds were ignored') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM counts
    WHERE schema_version <> 1
       OR anchor_locus_id IS NULL
       OR threshold_set_id <> 'synthetic_v1'
       OR genome_id <> 'genome1'
       OR motif_set_id <> 'motifs1'
       OR target_motif_id <> 'MA0861.2'
       OR calibration_stratum_id <> 'all_tp73_anchors'
) THEN error('threshold count identity or provenance is incomplete') END;
SELECT CASE WHEN (SELECT COUNT(DISTINCT anchor_locus_id) FROM counts) <> 2
    THEN error('physical TP73 anchor IDs are unstable across motifs') END;
SQL

mkdir -p "$temporary/package/tables/jaspar2026/tp73_motif_threshold_count"
cp "$temporary/threshold-counts.parquet" \
    "$temporary/package/tables/jaspar2026/tp73_motif_threshold_count/part-000000.parquet"
awk '
    /^CREATE OR REPLACE VIEW tp73_motif_threshold_count AS/ { capture = 1 }
    /^-- TP73-context-relevant same-motif cofactor architecture/ { capture = 0 }
    capture { print }
' "$repository_root/sql/schema.sql" > "$temporary/threshold-count-schema.sql"
(
    cd "$temporary/package"
    "$duckdb" -batch threshold-count.duckdb \
        < "$temporary/threshold-count-schema.sql"
    "$duckdb" -batch threshold-count.duckdb >/dev/null <<SQL
SELECT CASE WHEN (SELECT COUNT(*) FROM tp73_motif_threshold_count) <> 4
    THEN error('schema view omitted threshold-count rows') END;
SELECT CASE WHEN (SELECT COUNT(*) FROM tp73_motif_threshold_count
                  WHERE n_neighbor_loci_above_threshold = 0) <> 2
    THEN error('schema view lost explicit zero counts') END;
SQL

    {
        echo 'PREPARE q18 AS'
        awk '/^-- Q18\./ { capture = 1 } capture { print }' \
            "$repository_root/sql/queries.sql"
        cat <<'SQL'
EXECUTE q18(
    threshold_set_id := 'synthetic_v1', genome_id := 'genome1',
    motif_set_id := 'motifs1', chrom := '1',
    threshold_role := 'tp73_context_binary_feature',
    target_motif_id := 'MA0861.2', neighbor_motif_id := 'M1',
    score_mode := 'log2_relative_risk', pseudocount := 1.0,
    calibration_stratum_id := 'all_tp73_anchors'
);
SQL
    } | "$duckdb" -csv threshold-count.duckdb \
        > "$temporary/q18-output.csv"
)

[[ $(grep -c '^' "$temporary/q18-output.csv") -eq 3 ]] || {
    echo "E: Q18 did not return one header plus every physical anchor." >&2
    exit 1
}

[[ -s $temporary/maxima.parquet.run_config.tsv ]] || {
    echo "E: Sparse context maxima lack their source-file run config." >&2
    exit 1
}
grep -Fq $'M1\t' "$temporary/maxima.parquet.run_config.tsv" || {
    echo "E: Sparse context run config lacks the motif accession." >&2
    exit 1
}

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT * FROM (VALUES
        ('1', 100::BIGINT, 110::BIGINT),
        ('2', 500::BIGINT, 510::BIGINT)
    ) AS t(chrom, anchor_start, anchor_end)
) TO '$temporary/multichrom-anchors.parquet' (FORMAT PARQUET);
SQL

if "$repository_root/scripts/build_sparse_context_maxima.py" \
    --anchor-parquet "$temporary/multichrom-anchors.parquet" \
    --cofactor M1 "$temporary/m1-plus.parquet" "$temporary/m1-minus.parquet" \
    --output "$temporary/multichrom.parquet" --duckdb "$duckdb" \
    >"$temporary/multichrom.out" 2>"$temporary/multichrom.err"
then
    echo "E: Multi-chromosome anchors were accepted for chromosome-free hits." >&2
    exit 1
fi
grep -Fq 'exactly one chromosome' "$temporary/multichrom.err" || {
    echo "E: Multi-chromosome failure was not explained." >&2
    exit 1
}

if "$repository_root/scripts/build_sparse_context_maxima.py" \
    --anchor-parquet "$temporary/anchors.parquet" \
    --cofactor M1 "$temporary/m1-plus.parquet" "$temporary/m1-minus.parquet" \
    --output "$temporary/maxima.parquet" --duckdb "$duckdb" \
    >"$temporary/reuse.out" 2>"$temporary/reuse.err"
then
    echo "E: Existing output was overwritten." >&2
    exit 1
fi
grep -Fq 'output already exists' "$temporary/reuse.err" || {
    echo "E: Existing-output failure was not explained." >&2
    exit 1
}

echo "Sparse context-maximum tests passed."

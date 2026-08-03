#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
scanner="$repository_root/pssm_scan_parquet"
fixture="$repository_root/test_files/synthetic_dense"
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-context-maxima.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

command -v duckdb >/dev/null 2>&1 || {
    echo "E: duckdb CLI is required." >&2
    exit 1
}
[[ -x $scanner ]] || {
    echo "E: Build the Arrow scanner with 'make pssm_scan_parquet'." >&2
    exit 1
}

python3 "$repository_root/scripts/build_fasta_index.py" \
    "$fixture/genome.fna" --output "$temporary/genome.fna.fai" >/dev/null
cat > "$temporary/anchors.tsv" <<'EOF'
Chromosome	From	To	Name
1	1	2	anchor_1
1	3	4	anchor_2
EOF

duckdb -batch :memory: >/dev/null <<SQL
COPY (
    SELECT * FROM (VALUES
        ('1'::VARCHAR, 1::BIGINT, 2::BIGINT),
        ('1'::VARCHAR, 3::BIGINT, 4::BIGINT)
    ) AS anchors(chrom, anchor_start, anchor_end)
) TO '$temporary/anchors.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL

common=(
    --genome "$fixture/genome.fna"
    --fasta-index "$temporary/genome.fna.fai"
    --pssm "$fixture/JASPAR2026_synthetic.jaspar"
    --motif MA9999.1
    --strand both
    --coordinate-mode bed
    --score-mode log2_relative_risk
    --pseudocount 0
    --threshold -20
    --skip-N
)

direct="$temporary/direct.parquet"
"$scanner" "${common[@]}" --regions "$temporary/anchors.tsv" \
    --context-maxima "$direct" --context-flank 0 >/dev/null

sparse_root="$temporary/sparse/task_data/task_id=synthetic"
"$scanner" "${common[@]}" --chr 1 --sparse-parquet \
    --motif-set-id synthetic_jaspar2026 --genome-id synthetic_acgtn_v1 \
    --outdir "$sparse_root" >/dev/null
plus=$(find "$sparse_root" -path '*/strand=plus/*.parquet' -print -quit)
minus=$(find "$sparse_root" -path '*/strand=minus/*.parquet' -print -quit)
[[ -n $plus && -n $minus ]] || {
    echo "E: Synthetic sparse inputs were not produced." >&2
    exit 1
}

established="$temporary/established.parquet"
"$repository_root/scripts/build_sparse_context_maxima.py" \
    --anchor-parquet "$temporary/anchors.parquet" \
    --cofactor MA9999.1 "$plus" "$minus" \
    --output "$established" --flank 0 --source-score-floor -20 \
    --duckdb duckdb --threads 1 --memory-limit 256MB \
    --max-temp-size 1GB --temp-directory "$temporary/duckdb-spill" >/dev/null

duckdb -batch :memory: >/dev/null <<SQL
CREATE VIEW direct AS SELECT * FROM read_parquet('$direct');
CREATE VIEW established AS SELECT * FROM read_parquet('$established');

SELECT CASE WHEN (SELECT count(*) FROM direct) <> 2
    THEN error('direct context output is not rectangular by anchor') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM direct
    WHERE schema_version <> 1
       OR motif_id <> 'MA9999.1'
       OR abs(context_score - 1.0) > 1e-6
       OR source_score_floor <> -20
       OR context_flank_bp <> 0
       OR capture_prefilter_center_bp <> 2
       OR observed_max_anchor_span_bp <> 1
       OR observed_max_context_span_bp <> 1
       OR context_distance_metric <> 'signed_interval_edge_distance'
) THEN error('direct context values or provenance are incorrect') END;

WITH direct_common AS (
    SELECT chrom, anchor_start, anchor_end, motif_id,
           round(context_score, 6) AS context_score,
           source_score_floor, context_flank_bp,
           capture_prefilter_center_bp,
           observed_max_anchor_span_bp, observed_max_context_span_bp,
           context_distance_metric
    FROM direct
), established_common AS (
    SELECT chrom, anchor_start, anchor_end, motif_id,
           round(context_score, 6) AS context_score,
           source_score_floor, context_flank_bp,
           capture_prefilter_center_bp,
           observed_max_anchor_span_bp, observed_max_context_span_bp,
           context_distance_metric
    FROM established
), differences AS (
    (SELECT * FROM direct_common EXCEPT ALL SELECT * FROM established_common)
    UNION ALL
    (SELECT * FROM established_common EXCEPT ALL SELECT * FROM direct_common)
)
SELECT CASE WHEN (SELECT count(*) FROM differences) <> 0
    THEN error('direct and sparse-derived context maxima disagree') END;
SQL

if "$scanner" "${common[@]}" --regions "$temporary/anchors.tsv" \
    --context-maxima "$direct" --context-flank 0 \
    >"$temporary/retry.log" 2>&1; then
    echo "E: Direct context writer replaced an existing result." >&2
    exit 1
fi
grep -Fq 'Refusing to replace existing context-maximum Parquet file' \
    "$temporary/retry.log"

echo "Direct context-maximum Parquet tests passed."

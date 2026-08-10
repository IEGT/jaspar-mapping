#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-context-finalize.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM
run_root="$temporary/run"
mkdir -p "$run_root/plan" "$run_root/staging"

{
    printf 'task_index\tchrom\tcofactor_motif_ids\toutput_tier\tbuilder_source_commit\tcontext_schema_version\tgtf_size_bytes\tgtf_sha256\tannotation_release\tpromoter_definition_id\tpromoter_upstream_bp\tpromoter_downstream_bp\n'
    printf '0\t1\tMA0001.1\tband\tabc123\t7\t0\tnone\tensembl_113\ttss_upstream_2000_downstream_500_v1\t2000\t500\n'
    printf '1\tX\tMA0002.1\tband\tabc123\t7\t0\tnone\tensembl_113\ttss_upstream_2000_downstream_500_v1\t2000\t500\n'
} > "$run_root/plan/context_tasks.tsv"

build_task() {
    local task_index=$1 chrom=$2 cofactor=$3
    local hits="$temporary/hits-$task_index.parquet"
    local input="$run_root/inputs/task-$task_index-chrom-$chrom"
    local package="$run_root/packages/chrom-$chrom/task-$task_index"
    duckdb :memory: -bail -c "
COPY (
    SELECT * FROM (VALUES
        ('$chrom', 100::BIGINT, 116::BIGINT, 'MA0861.2', 'TP73', '+',
         2.0::DOUBLE, 'log2_relative_risk', 1.0::DOUBLE, 0.8::DOUBLE,
         -5.0::DOUBLE,
         'homo_sapiens_grch38_ensembl113_primary',
         'jaspar2026_core_nonredundant', 'uniform_acgt_v1',
         'additive_per_base'),
        ('$chrom', 120::BIGINT, 130::BIGINT, '$cofactor', 'Synthetic', '+',
         1.0::DOUBLE, 'log2_relative_risk', 1.0::DOUBLE, 0.7::DOUBLE,
         0.0::DOUBLE,
         'homo_sapiens_grch38_ensembl113_primary',
         'jaspar2026_core_nonredundant', 'uniform_acgt_v1',
         'additive_per_base')
    ) AS t(chrom, start, \"end\", motif_id, motif_name, strand, score,
           score_mode, pseudocount, pwm_relative_score, minimum_score,
           genome_id, motif_set_id, background_model_id, pseudocount_scheme)
) TO '$hits' (FORMAT PARQUET);" >/dev/null

    "$repository_root/scripts/build_motif_context.py" \
        --motif-hits "$hits" \
        --motif-hit-source-label "$input/**/*.parquet" \
        --input-uniqueness validated_scan_inventory \
        --output "$package" --output-tier band \
        --anchor-motif MA0861.2 --source-commit abc123 \
        --motif-set-id jaspar2026_core_nonredundant \
        --genome-id homo_sapiens_grch38_ensembl113_primary \
        --annotation-release ensembl_113 \
        --promoter-definition-id tss_upstream_2000_downstream_500_v1 \
        --promoter-upstream-bp 2000 --promoter-downstream-bp 500 \
        --anchor-minimum-score -1 --partner-minimum-score 0 \
        --anchor-selection-mode local_peak --anchor-local-peak-flank 150 \
        --score-mode log2_relative_risk --pseudocount 1 \
        --chrom "$chrom" --capture-flank 150 --context-flank 150 \
        --tandem-flank 20 --cofactor-pair-flank 150 \
        --memory-limit 1GB --max-temp-size 1GB >/dev/null

    mkdir -p "$input"
    cat > "$input/input_manifest.json" <<EOF
{
  "format_version": 1,
  "source_package": "/synthetic/scan",
  "source_manifest_sha256": "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
  "motifs": ["MA0861.2", "$cofactor"],
  "chromosomes": ["$chrom"],
  "files": []
}
EOF
    cp "$input/input_manifest.json" "$package/input_manifest.json"
}

build_task 0 1 MA0001.1
if "$repository_root/scripts/finalize_motif_context_run.py" \
    --run-root "$run_root" --task-file "$run_root/plan/context_tasks.tsv" \
    --duckdb duckdb >"$temporary/incomplete.out" 2>"$temporary/incomplete.err"; then
    echo "E: Incomplete context run was finalized." >&2
    exit 1
fi
grep -Fq '1 of 2 tasks failed' "$temporary/incomplete.err"
[[ ! -e $run_root/final/manifest.json ]]

build_task 1 X MA0002.1
mkdir -p "$run_root/staging/task-0/orphan-attempt"
touch "$run_root/staging/task-0/orphan-attempt/marker"
"$repository_root/scripts/finalize_motif_context_run.py" \
    --run-root "$run_root" --task-file "$run_root/plan/context_tasks.tsv" \
    --duckdb duckdb >/dev/null
"$repository_root/scripts/finalize_motif_context_run.py" \
    --run-root "$run_root" --task-file "$run_root/plan/context_tasks.tsv" \
    --duckdb duckdb >/dev/null

[[ -f $run_root/final/manifest.json ]]
[[ -f $run_root/final/context.duckdb ]]
[[ -f $run_root/locks/finalizer.lock ]]
[[ -e $run_root/staging/task-0/orphan-attempt/marker ]]
grep -Fq 'staging/task-0/orphan-attempt' \
    "$run_root/final/orphan_staging_paths.json"

feature_0="$run_root/packages/chrom-1/task-0/tables/jaspar2026/anchor_motif_band_feature/data.parquet"
feature_1="$run_root/packages/chrom-X/task-1/tables/jaspar2026/anchor_motif_band_feature/data.parquet"
duckdb "$run_root/final/context.duckdb" -readonly -bail -c "
SELECT CASE WHEN (SELECT COUNT(*) FROM context_task_inventory) <> 2
    THEN error('final task inventory is incomplete') END;
SELECT CASE WHEN (SELECT COUNT(*) FROM context_file_inventory
                  WHERE dataset = 'anchor_motif_band_feature') <> 2
    THEN error('feature-file inventory is incomplete') END;
SELECT CASE WHEN (SELECT COUNT(*) FROM anchor_motif_band_feature_files(
                    ['$feature_0', '$feature_1'])) <> 2
    THEN error('exact-file feature macro did not read both task packages') END;" \
    >/dev/null

echo "Motif-context finalizer tests passed."

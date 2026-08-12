#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-context-finalize.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM
run_root="$temporary/run"
mkdir -p "$run_root/plan" "$run_root/staging"

{
    printf 'task_index\tchrom\tcofactor_motif_ids\toutput_tier\tbuilder_source_commit\tcontext_schema_version\tgtf_size_bytes\tgtf_sha256\tannotation_release\tpromoter_definition_id\tpromoter_upstream_bp\tpromoter_downstream_bp\tdownstream_definition_id\tdownstream_upstream_bp\tdownstream_downstream_bp\n'
    printf '0\t1\tMA0001.1\tband\tabc123\t9\t0\tnone\tensembl_113\ttss_upstream_2000_downstream_500_v1\t2000\t500\ttes_upstream_500_downstream_2000_v1\t500\t2000\n'
    printf '1\tX\tMA0002.1\tband\tabc123\t9\t0\tnone\tensembl_113\ttss_upstream_2000_downstream_500_v1\t2000\t500\ttes_upstream_500_downstream_2000_v1\t500\t2000\n'
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

# An anchor-only plan is a separate, normalized annotation layer: one summary
# package per chromosome, no fake cofactor, with the complete TSS/promoter
# relationship payload required by the finalization gate.
anchor_root="$temporary/anchor-run"
anchor_input="$anchor_root/inputs/task-0-chrom-1"
anchor_package="$anchor_root/packages/chrom-1/task-0"
mkdir -p "$anchor_root/plan" "$anchor_input"
cat > "$temporary/annotation.gtf" <<'EOF'
1	test	transcript	51	300	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENE1"; transcript_biotype "protein_coding";
1	test	exon	51	120	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENE1"; exon_number "1";
1	test	exon	181	300	.	+	.	gene_id "G1"; transcript_id "T1"; gene_name "GENE1"; exon_number "2";
EOF
gtf_size=$(wc -c < "$temporary/annotation.gtf" | tr -d '[:space:]')
if command -v sha256sum >/dev/null 2>&1; then
    gtf_sha256=$(sha256sum "$temporary/annotation.gtf" | awk '{print $1}')
else
    gtf_sha256=$(shasum -a 256 "$temporary/annotation.gtf" | awk '{print $1}')
fi
duckdb :memory: -bail -c "
COPY (
    SELECT * FROM (VALUES
        ('1', 100::BIGINT, 116::BIGINT, 'MA0861.2', 'TP73', '+',
         2.0::DOUBLE, 'log2_relative_risk', 1.0::DOUBLE, 0.8::DOUBLE,
         -5.0::DOUBLE,
         'homo_sapiens_grch38_ensembl113_primary',
         'jaspar2026_core_nonredundant', 'uniform_acgt_v1',
         'additive_per_base'),
        ('1', 130::BIGINT, 146::BIGINT, 'MA0861.2', 'TP73', '-',
         1.0::DOUBLE, 'log2_relative_risk', 1.0::DOUBLE, 0.7::DOUBLE,
         -5.0::DOUBLE,
         'homo_sapiens_grch38_ensembl113_primary',
         'jaspar2026_core_nonredundant', 'uniform_acgt_v1',
         'additive_per_base')
    ) AS t(chrom, start, \"end\", motif_id, motif_name, strand, score,
           score_mode, pseudocount, pwm_relative_score, minimum_score,
           genome_id, motif_set_id, background_model_id, pseudocount_scheme)
) TO '$temporary/anchor-hits.parquet' (FORMAT PARQUET);" >/dev/null

"$repository_root/scripts/build_motif_context.py" \
    --motif-hits "$temporary/anchor-hits.parquet" \
    --motif-hit-source-label "$anchor_input/**/*.parquet" \
    --input-uniqueness validated_scan_inventory \
    --gtf "$temporary/annotation.gtf" --gtf-size-bytes "$gtf_size" \
    --gtf-sha256 "$gtf_sha256" \
    --output "$anchor_package" --output-tier summary \
    --anchor-motif MA0861.2 --source-commit abc123 \
    --motif-set-id jaspar2026_core_nonredundant \
    --genome-id homo_sapiens_grch38_ensembl113_primary \
    --annotation-release ensembl_113 \
    --promoter-definition-id tss_upstream_2000_downstream_500_v1 \
    --promoter-upstream-bp 2000 --promoter-downstream-bp 500 \
    --anchor-minimum-score -1 --partner-minimum-score 0 \
    --anchor-selection-mode local_peak --anchor-local-peak-flank 150 \
    --score-mode log2_relative_risk --pseudocount 1 \
    --chrom 1 --capture-flank 150 --context-flank 150 \
    --tandem-flank 20 --cofactor-pair-flank 150 \
    --memory-limit 1GB --max-temp-size 1GB >/dev/null
cat > "$anchor_input/input_manifest.json" <<'EOF'
{
  "format_version": 1,
  "source_package": "/synthetic/scan",
  "source_manifest_sha256": "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
  "motifs": ["MA0861.2"],
  "chromosomes": ["1"],
  "files": []
}
EOF
cp "$anchor_input/input_manifest.json" "$anchor_package/input_manifest.json"
{
    printf 'task_index\tchrom\tcofactor_motif_ids\toutput_tier\tbuilder_source_commit\tcontext_schema_version\tgtf_size_bytes\tgtf_sha256\tannotation_release\tpromoter_definition_id\tpromoter_upstream_bp\tpromoter_downstream_bp\tdownstream_definition_id\tdownstream_upstream_bp\tdownstream_downstream_bp\ttask_kind\n'
    printf '0\t1\tnone\tsummary\tabc123\t9\t%s\t%s\tensembl_113\ttss_upstream_2000_downstream_500_v1\t2000\t500\ttes_upstream_500_downstream_2000_v1\t500\t2000\tanchor_annotation\n' \
        "$gtf_size" "$gtf_sha256"
} > "$anchor_root/plan/context_tasks.tsv"
"$repository_root/scripts/finalize_motif_context_run.py" \
    --run-root "$anchor_root" --task-file "$anchor_root/plan/context_tasks.tsv" \
    --duckdb duckdb >/dev/null
grep -Fq '"task_kind": "anchor_annotation"' \
    "$anchor_root/final/manifest.json"
duckdb "$anchor_root/final/context.duckdb" -readonly -bail -c "
SELECT CASE WHEN (SELECT task_kind FROM context_task_inventory)
                    <> 'anchor_annotation'
    THEN error('anchor task kind was not cataloged') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM context_file_inventory
    WHERE dataset = 'tp73_context_anchor'
) THEN error('anchor payload was not cataloged') END;" >/dev/null

echo "Motif-context finalizer tests passed."

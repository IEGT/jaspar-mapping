#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-genome-scan-manager.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

command -v duckdb >/dev/null 2>&1 || {
    echo "E: duckdb CLI is required." >&2
    exit 1
}
scanner="$repository_root/pssm_scan_parquet"
[[ -x $scanner ]] || {
    echo "E: Build pssm_scan_parquet before this test." >&2
    exit 1
}

cat > "$temporary/genome.fna" <<'EOF'
>1 synthetic chromosome
ACGTNACGT
EOF
python3 "$repository_root/scripts/build_fasta_index.py" \
    "$temporary/genome.fna" --output "$temporary/genome.fna.fai" >/dev/null

cat > "$temporary/motifs.jaspar" <<'EOF'
>MA0861.2 TP73_TEST
A [ 4 1 ]
C [ 1 4 ]
G [ 1 1 ]
T [ 1 1 ]
>MA0001.1 OTHER_TEST
A [ 1 ]
C [ 1 ]
G [ 4 ]
T [ 1 ]
EOF

source_commit=0123456789abcdef0123456789abcdef01234567
run_root="$temporary/run"
manager="$repository_root/scripts/manage_genome_scan.py"

task_count=$(
    "$manager" prepare \
        --run-root "$run_root" --run-id synthetic_jaspar_scan_v1 \
        --source-commit "$source_commit" \
        --genome "$temporary/genome.fna" \
        --fasta-index "$temporary/genome.fna.fai" \
        --jaspar "$temporary/motifs.jaspar" \
        --genome-id synthetic_genome_v1 --motif-set-id synthetic_jaspar_v1 \
        --taxon-id 9606 --species 'Homo sapiens synthetic' \
        --assembly-name SYN1 --assembly-accession GCA_SYNTHETIC \
        --ensembl-release 113 \
        --fasta-url https://example.invalid/genome.fna \
        --jaspar-url https://example.invalid/motifs.jaspar \
        --chrom 1 --motif-batch-size 1 --minimum-free-gib 0 \
        | awk -F= '$1 == "TASK_COUNT" { print $2 }'
)
[[ $task_count == 2 ]] || {
    echo "E: Expected two non-overlapping motif tasks, got $task_count." >&2
    exit 1
}

for task_index in 0 1; do
    "$manager" run-task --run-root "$run_root" --task-index "$task_index" \
        --scanner "$scanner" --duckdb duckdb --allow-local >/dev/null
done

# A retry must validate and reuse the promoted task rather than rescan it.
"$manager" run-task --run-root "$run_root" --task-index 0 \
    --scanner "$scanner" --duckdb duckdb --allow-local >/dev/null

"$manager" finalize --run-root "$run_root" --duckdb duckdb >/dev/null

(
    cd "$run_root/package"
    duckdb jaspar_genome_scan.duckdb -bail -c "
        SELECT CASE WHEN (SELECT COUNT(*) FROM scan_task) <> 2
            THEN error('scan_task inventory is incomplete') END;
        SELECT CASE WHEN (SELECT COUNT(*) FROM scan_file_inventory) <> 4
            THEN error('scan_file_inventory is incomplete') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM scan_task_completeness WHERE NOT complete
        ) THEN error('a promoted task is incomplete') END;
        SELECT CASE WHEN NOT EXISTS (
            SELECT 1 FROM scan_threshold_policy
            WHERE policy_id = 'tp73_calibrated' AND minimum_score = -5
        ) OR NOT EXISTS (
            SELECT 1 FROM scan_threshold_policy
            WHERE policy_id = 'default_uncalibrated' AND minimum_score = -1
        ) THEN error('threshold policies were not preserved') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM motif_hit
            WHERE run_id <> 'synthetic_jaspar_scan_v1'
               OR task_id IS NULL
               OR genome_id <> 'synthetic_genome_v1'
               OR motif_set_id <> 'synthetic_jaspar_v1'
               OR background_model_id <> 'uniform_acgt_v1'
               OR pseudocount_scheme <> 'additive_per_base'
               OR n_policy <> 'skip'
               OR coordinate_mode <> 'bed'
        ) THEN error('logical motif_hit lost scan configuration') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM scan_file_inventory
            WHERE state <> 'complete'
               OR expected_windows < emitted_hits
               OR expected_windows <> skipped_n_windows
                    + sentinel_score_windows + below_minimum_score_windows
                    + pwm_filtered_windows + emitted_hits
               OR bytes <= 0 OR length(sha256) <> 64
               OR source_commit <> '$source_commit'
               OR slurm_job_id IS NULL OR slurm_array_task_id IS NULL
               OR minimum_pwm_relative_score IS NOT NULL
               OR maximum_pwm_relative_score IS NOT NULL
               OR coordinate_mode <> 'bed' OR n_policy <> 'skip'
               OR matched_sequence_policy <> 'omitted'
        ) THEN error('file inventory lost validation or provenance fields') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM scan_run
            WHERE minimum_pwm_relative_score IS NOT NULL
               OR maximum_pwm_relative_score IS NOT NULL
               OR emitted_hit_count <= 0
        ) THEN error('scan_run lost explicit score configuration') END;
        SELECT CASE WHEN EXISTS (
            SELECT 1 FROM genome
            WHERE genome_id <> 'synthetic_genome_v1'
               OR taxon_id <> 9606
               OR assembly_accession <> 'GCA_SYNTHETIC'
               OR ensembl_release <> 113
               OR length(fasta_sha256) <> 64
               OR length(sequence_set_sha256) <> 64
        ) THEN error('genome identity or checksums are incomplete') END;
    " >/dev/null
)

status=$($manager status --run-root "$run_root")
grep -Fq $'tasks_complete\t2' <<< "$status"
grep -Fq $'finalized\ttrue' <<< "$status"

# A failed scanner attempt must remain in its unique staging directory and
# must never create a promoted task.
failed_run_root="$temporary/failed-run"
"$manager" prepare \
    --run-root "$failed_run_root" --run-id synthetic_failed_scan_v1 \
    --source-commit "$source_commit" \
    --genome "$temporary/genome.fna" \
    --fasta-index "$temporary/genome.fna.fai" \
    --jaspar "$temporary/motifs.jaspar" \
    --genome-id synthetic_genome_v1 --motif-set-id synthetic_jaspar_v1 \
    --taxon-id 9606 --species 'Homo sapiens synthetic' \
    --assembly-name SYN1 --assembly-accession GCA_SYNTHETIC \
    --ensembl-release 113 \
    --fasta-url https://example.invalid/genome.fna \
    --jaspar-url https://example.invalid/motifs.jaspar \
    --chrom 1 --motif-batch-size 1 --minimum-free-gib 0 >/dev/null
if "$manager" run-task --run-root "$failed_run_root" --task-index 0 \
    --scanner "$(type -P false)" --duckdb duckdb --allow-local \
    >/dev/null 2>&1
then
    echo "E: Deliberately failing scanner unexpectedly succeeded." >&2
    exit 1
fi
failure_count=$(find "$failed_run_root/staging" -name failure.json -type f | wc -l)
[[ $failure_count -eq 1 ]] || {
    echo "E: Failed attempt was not retained with exactly one failure marker." >&2
    exit 1
}
if find "$failed_run_root/package/task_data" -name task_result.json \
    -type f 2>/dev/null | grep -q .
then
    echo "E: Failed task was promoted." >&2
    exit 1
fi

echo "Genome-scan manager tests passed."

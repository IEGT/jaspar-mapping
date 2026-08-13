#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_tp73_genome_context_maxima_slurm.sh --run-root DIR
       --scan-package DIR --evidence-package DIR --threshold-registry FILE
       --runtime-prefix DIR [options]

Prepare and submit restart-safe autosomal TP73-context maxima for every
non-target JASPAR motif. Each array element stages the 22 small anchor files
once, processes a bounded motif batch, and publishes each completed motif
before advancing. Scan payloads are opened by exact inventory path.

Options:
  --run-root DIR             Dedicated output below /data/sm718
  --scan-package DIR         Finalized low-floor genome scan
  --evidence-package DIR     Finalized whole-genome TP73 evidence package
  --threshold-registry FILE  Chromosome-1 context-threshold Parquet
  --runtime-prefix DIR       Runtime containing DuckDB
  --source DIR               Repository root (default: script parent)
  --run-id ID                Immutable run identifier
  --threshold-set-id ID      Applied whole-genome threshold-set identifier
  --maximum-source-score-floor SCORE
                             Reject any non-target scan retained above this
                             floor (default: -1)
  --motifs-per-batch N       Motifs per restartable array element (default: 3)
  --scratch-root DIR         Node-local scratch parent (default: /scratch/$USER)
  --account NAME             Slurm account (default: cluster)
  --partition NAME           Slurm partition (default: requeue)
  --max-concurrent N         Maximum live tasks (default: 20)
  --array-size N             Maximum elements per array chunk (default: 1000)
  --cpus N                   CPUs per task (default: 4)
  --memory SIZE              Memory per task (default: 32G)
  --time LIMIT               Time per task (default: 02:00:00)
  --final-memory SIZE        Finalizer memory (default: 16G)
  --final-time LIMIT         Finalizer time (default: 01:00:00)
  --dry-run                  Print sbatch commands without submitting
  -h, --help                 Show this help

SIGUSR1 reports batch, motif, chromosome, phase, elapsed time, and child PID.
Array chunks are serialized so the declared concurrency remains global.
EOF
}

run_root=""
scan_package=""
evidence_package=""
threshold_registry=""
runtime_prefix=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
run_id=jaspar2026_grch38_tp73_context_maxima_autosomes_v1
threshold_set_id=jaspar2026_grch38_tp73_context_applied_scan_floor_v1
maximum_source_score_floor=-1
motifs_per_batch=3
scratch_root="/scratch/${USER:-sm718}"
account=cluster
partition=requeue
max_concurrent=20
array_size=1000
cpus=4
memory=32G
wall_time=02:00:00
final_memory=16G
final_time=01:00:00
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --scan-package) scan_package=${2:?}; shift 2 ;;
        --evidence-package) evidence_package=${2:?}; shift 2 ;;
        --threshold-registry) threshold_registry=${2:?}; shift 2 ;;
        --runtime-prefix) runtime_prefix=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --run-id) run_id=${2:?}; shift 2 ;;
        --threshold-set-id) threshold_set_id=${2:?}; shift 2 ;;
        --maximum-source-score-floor) maximum_source_score_floor=${2:?}; shift 2 ;;
        --motifs-per-batch) motifs_per_batch=${2:?}; shift 2 ;;
        --scratch-root) scratch_root=${2:?}; shift 2 ;;
        --account) account=${2:?}; shift 2 ;;
        --partition) partition=${2:?}; shift 2 ;;
        --max-concurrent) max_concurrent=${2:?}; shift 2 ;;
        --array-size) array_size=${2:?}; shift 2 ;;
        --cpus) cpus=${2:?}; shift 2 ;;
        --memory) memory=${2:?}; shift 2 ;;
        --time) wall_time=${2:?}; shift 2 ;;
        --final-memory) final_memory=${2:?}; shift 2 ;;
        --final-time) final_time=${2:?}; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scan_package && -n $evidence_package &&
   -n $threshold_registry && -n $runtime_prefix ]] || { usage >&2; exit 2; }
for value in "$motifs_per_batch" "$max_concurrent" "$array_size" "$cpus"; do
    [[ $value =~ ^[1-9][0-9]*$ ]] || {
        echo "E: Batch, concurrency, array-size, and CPU values must be positive." >&2
        exit 2
    }
done
for value in "$run_id" "$threshold_set_id"; do
    [[ $value =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
        echo "E: Run and threshold-set IDs contain unsafe characters." >&2
        exit 2
    }
done
for path in "$run_root" "$scan_package" "$evidence_package" \
            "$threshold_registry" "$runtime_prefix" "$source"; do
    case "$path" in
        /data/sm718/*) ;;
        *) echo "E: Production paths must be below /data/sm718: $path" >&2; exit 2 ;;
    esac
done

source=$(cd "$source" && pwd -P)
scan_package=$(cd "$scan_package" && pwd -P)
evidence_package=$(cd "$evidence_package" && pwd -P)
runtime_prefix=$(cd "$runtime_prefix" && pwd -P)
threshold_registry=$(cd "$(dirname "$threshold_registry")" && pwd -P)/$(basename "$threshold_registry")
mkdir -p "$run_root/plan" "$run_root/logs"
run_root=$(cd "$run_root" && pwd -P)
duckdb="$runtime_prefix/duckdb/bin/duckdb"
[[ -x $duckdb && -f $threshold_registry ]] || {
    echo "E: DuckDB runtime or threshold registry is missing." >&2
    exit 1
}
if ! git -C "$source" diff --quiet --ignore-submodules -- ||
   ! git -C "$source" diff --cached --quiet --ignore-submodules --; then
    echo "E: Production submission requires a tracked-clean source tree." >&2
    exit 1
fi

batch_count=$(
    python3 "$source/scripts/manage_tp73_genome_context_maxima.py" prepare \
        --run-root "$run_root" --scan-package "$scan_package" \
        --evidence-package "$evidence_package" \
        --threshold-registry "$threshold_registry" \
        --runtime-prefix "$runtime_prefix" --source "$source" \
        --duckdb "$duckdb" --run-id "$run_id" \
        --applied-threshold-set-id "$threshold_set_id" \
        --maximum-source-score-floor "$maximum_source_score_floor" \
        --motifs-per-batch "$motifs_per_batch" \
        --scratch-root "$scratch_root" --threads "$cpus" \
        --memory-limit 28GB --max-temp-size 100GB
)
[[ $batch_count =~ ^[1-9][0-9]*$ ]] || {
    echo "E: Invalid prepared batch count: $batch_count" >&2
    exit 1
}

submission_record="$run_root/plan/slurm_submission.tsv"
if [[ -e $submission_record && $dry_run -eq 0 ]]; then
    echo "E: Run already has a submission record: $submission_record" >&2
    exit 1
fi

array_jobs=()
dependency=()
array_start=0
array_chunk=0
while (( array_start < batch_count )); do
    chunk_tasks=$((batch_count - array_start))
    if (( chunk_tasks > array_size )); then
        chunk_tasks=$array_size
    fi
    array_end=$((chunk_tasks - 1))
    submission=(
        sbatch --parsable --account="$account" --partition="$partition" --requeue
        --job-name="tp73_ctx_${array_chunk}" --nodes=1 --ntasks=1
        --cpus-per-task="$cpus" --mem="$memory" --time="$wall_time"
        --array="0-${array_end}%${max_concurrent}" --signal=B:USR1@120
        --open-mode=append --chdir="$source"
        --output="$run_root/logs/context-%A_%a.out"
        --error="$run_root/logs/context-%A_%a.err"
    )
    if (( ${#dependency[@]} > 0 )); then
        submission+=("${dependency[@]}")
    fi
    submission+=(
        "$source/scripts/run_tp73_genome_context_maxima_batch.py"
        --run-root "$run_root" --batch-offset "$array_start"
    )
    if [[ $dry_run -eq 1 ]]; then
        printf '%q ' "${submission[@]}"; printf '\n'
        array_job="ARRAY_JOB_ID_${array_chunk}"
    else
        array_job=$("${submission[@]}")
    fi
    array_jobs+=("$array_job")
    dependency=(--dependency="afterok:$array_job")
    array_start=$((array_start + chunk_tasks))
    array_chunk=$((array_chunk + 1))
done
last_array_job=$array_job

final_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --dependency="afterok:$last_array_job" --job-name=tp73_ctx_final
    --nodes=1 --ntasks=1 --cpus-per-task=2 --mem="$final_memory"
    --time="$final_time" --open-mode=append --chdir="$source"
    --output="$run_root/logs/finalize-%j.out"
    --error="$run_root/logs/finalize-%j.err"
    "$source/scripts/run_tp73_genome_context_maxima_finalize.sh"
    --run-root "$run_root" --source "$source" --duckdb "$duckdb"
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${final_submission[@]}"; printf '\n'
    exit 0
fi
final_job=$("${final_submission[@]}")

array_jobs_text=$(IFS=,; printf '%s' "${array_jobs[*]}")
source_commit=$(git -C "$source" rev-parse HEAD)
printf 'array_jobs\tfinalizer_job\tbatch_count\tmotifs_per_batch\tmax_concurrent\tarray_chunks\tsource_commit\trun_id\n' \
    > "$submission_record"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$array_jobs_text" "$final_job" "$batch_count" "$motifs_per_batch" \
    "$max_concurrent" "$array_chunk" "$source_commit" "$run_id" \
    >> "$submission_record"
echo "I: Submitted context arrays $array_jobs_text and finalizer $final_job" >&2
printf '%s\t%s\n' "$array_jobs_text" "$final_job"

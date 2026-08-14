#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_tp73_cofactor_enrichment_slurm_batch.sh --run-root DIR
       --source-threshold-run DIR --task-file FILE --run-config FILE [options]

Run a bounded batch of enrichment motifs while staging the shared whole-genome
TP73 anchor evidence only once. Each motif is still validated and published
atomically by run_tp73_cofactor_enrichment_slurm_task.sh, so a requeue skips all
motifs already completed by an earlier attempt.

Options:
  --source-threshold-run DIR Completed context-maxima run (legacy threshold
                             calibration layouts remain supported)
  --source DIR               Repository root (default: script parent)
  --duckdb FILE              DuckDB CLI (default: duckdb)
  --rscript FILE             Rscript executable (default: Rscript)
  --motifs-per-job N         Consecutive motifs handled by one job (default: 2)
  --batch-offset N           Add N to SLURM_ARRAY_TASK_ID (default: 0)
  --block-size BP            Genomic uncertainty block (default: 5000000)
  --spline-df N              TP73 score spline degrees of freedom (default: 4)
  --minimum-class-fraction N Primary class-support floor (default: 0.01)
  -h, --help                 Show this help

SIGUSR1 prints exactly one batch-level progress line. It does not interrupt the
active motif evaluator.
EOF
}

run_root=""
source_threshold_run=""
task_file=""
run_config=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
duckdb=duckdb
rscript=Rscript
motifs_per_job=2
batch_offset=0
block_size=5000000
spline_df=4
minimum_class_fraction=0.01
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --source-threshold-run) source_threshold_run=${2:?}; shift 2 ;;
        --task-file) task_file=${2:?}; shift 2 ;;
        --run-config) run_config=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --rscript) rscript=${2:?}; shift 2 ;;
        --motifs-per-job) motifs_per_job=${2:?}; shift 2 ;;
        --batch-offset) batch_offset=${2:?}; shift 2 ;;
        --block-size) block_size=${2:?}; shift 2 ;;
        --spline-df) spline_df=${2:?}; shift 2 ;;
        --minimum-class-fraction) minimum_class_fraction=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $source_threshold_run && -n $task_file &&
   -n $run_config ]] || { usage >&2; exit 2; }
[[ ${SLURM_ARRAY_TASK_ID:-} =~ ^[0-9]+$ &&
   $motifs_per_job =~ ^[1-9][0-9]*$ && $batch_offset =~ ^[0-9]+$ ]] || {
    echo "E: Array index, motifs per job, or batch offset is invalid." >&2
    exit 2
}
for path in "$run_root" "$source_threshold_run" "$source"; do
    [[ -d $path ]] || { echo "E: Directory not found: $path" >&2; exit 1; }
done
for path in "$task_file" "$run_config"; do
    [[ -f $path ]] || { echo "E: File not found: $path" >&2; exit 1; }
done

read -r task_count anchor anchor_sha < <(python3 - "$run_config" <<'PY'
import json
import sys

config = json.load(open(sys.argv[1], encoding="utf-8"))
print(config["task_count"], config["anchor_evidence"],
      config["anchor_evidence_sha256"])
PY
)
[[ $task_count =~ ^[1-9][0-9]*$ && -f $anchor ]] || {
    echo "E: Immutable batch anchor contract is invalid." >&2
    exit 1
}

sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

batch_index=$((batch_offset + SLURM_ARRAY_TASK_ID))
first=$((batch_index * motifs_per_job))
(( first < task_count )) || {
    echo "E: Batch index $batch_index starts beyond $task_count motifs." >&2
    exit 2
}
last=$((first + motifs_per_job - 1))
(( last < task_count )) || last=$((task_count - 1))

start_epoch=$(date +%s)
phase=staging_anchor
current_task=none
child_pid=""
progress_report() {
    printf 'I: progress signal=SIGUSR1 phase=%s batch=%s motif_task=%s range=%s-%s elapsed_s=%s child_pid=%s\n' \
        "$phase" "$batch_index" "$current_task" "$first" "$last" \
        "$(($(date +%s) - start_epoch))" "${child_pid:-none}" >&2
}
trap progress_report USR1

run_child() {
    "$@" &
    child_pid=$!
    local status usr1_wait_status
    usr1_wait_status=$((128 + $(kill -l USR1)))
    while true; do
        set +e
        wait "$child_pid"
        status=$?
        set -e
        [[ $status -eq $usr1_wait_status ]] && continue
        kill -0 "$child_pid" 2>/dev/null && continue
        child_pid=""
        return "$status"
    done
}

scratch_base=${SLURM_TMPDIR:-/scratch/${USER:-sm718}}
scratch="$scratch_base/jaspar-enrichment-batch-${SLURM_JOB_ID:-manual}-${batch_index}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
mkdir -p "$scratch"
staged_anchor="$scratch/tp73_anchor_evidence.parquet"
run_child cp "$anchor" "$staged_anchor"
[[ $(sha256_file "$staged_anchor") == "$anchor_sha" ]] || {
    echo "E: Staged whole-genome anchor checksum differs." >&2
    exit 1
}

for ((task_index = first; task_index <= last; ++task_index)); do
    phase=evaluating_motif
    current_task=$task_index
    run_child env SLURM_ARRAY_TASK_ID="$task_index" \
        "$source/scripts/run_tp73_cofactor_enrichment_slurm_task.sh" \
        --run-root "$run_root" --source-threshold-run "$source_threshold_run" \
        --task-file "$task_file" --run-config "$run_config" \
        --source "$source" --duckdb "$duckdb" --rscript "$rscript" \
        --pre-staged-anchor "$staged_anchor" --task-offset 0 \
        --block-size "$block_size" --spline-df "$spline_df" \
        --minimum-class-fraction "$minimum_class_fraction"
done

phase=complete
current_task=none
echo "I: Completed enrichment batch $batch_index (tasks $first-$last)." >&2

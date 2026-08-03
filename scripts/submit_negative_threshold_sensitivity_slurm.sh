#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_negative_threshold_sensitivity_slurm.sh --run-root DIR
       --threshold-list FILE --jaspar FILE --genome FILE
       --anchor-evidence FILE --runtime-prefix DIR [options]

Prepare and submit a restart-safe chromosome-1 sensitivity analysis for every
motif whose tracked recommendation is exactly zero. Each short array task
rescans one motif only around TP73 anchors down to score -20 and evaluates the
inclusive integer thresholds -20..0. The final result is explicitly
non-normative and does not replace the tracked threshold list.

Options:
  --run-root DIR        Dedicated durable run below /data/sm718 (required)
  --threshold-list FILE Tracked two-column v1 threshold list (required)
  --jaspar FILE         Exact JASPAR 2026 matrix source (required)
  --genome FILE         Indexed GRCh38 FASTA (required)
  --fasta-index FILE    FASTA index (default: GENOME.fai)
  --anchor-evidence FILE
                        Existing TP73/CUT&RUN anchor Parquet (required)
  --runtime-prefix DIR  Existing calibration runtime with duckdb/ and r/
  --scanner FILE        Arrow scanner (default: SOURCE/pssm_scan_parquet)
  --source DIR          Repository root (default: script parent)
  --account NAME        Slurm account (default: cluster)
  --partition NAME      Slurm partition (default: requeue)
  --max-concurrent N    Live tasks (default: 20)
  --array-size N        Tasks per scheduler array (default: 1000)
  --memory SIZE         Memory per task (default: 16G)
  --time LIMIT          Time per task (default: 00:30:00)
  --minimum-free-gib N  Durable free-space reserve (default: 5)
  --dry-run             Print sbatch commands without submitting
  -h, --help            Show this help

The helper never installs software. Build the clean Arrow scanner first and
point --runtime-prefix at the already recorded threshold-calibration runtime.
EOF
}

run_root=""
threshold_list=""
jaspar=""
genome=""
fasta_index=""
anchor_evidence=""
runtime_prefix=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
scanner=""
account=cluster
partition=requeue
max_concurrent=20
array_size=1000
memory=16G
wall_time=00:30:00
minimum_free_gib=5
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --threshold-list) threshold_list=${2:?}; shift 2 ;;
        --jaspar) jaspar=${2:?}; shift 2 ;;
        --genome) genome=${2:?}; shift 2 ;;
        --fasta-index) fasta_index=${2:?}; shift 2 ;;
        --anchor-evidence) anchor_evidence=${2:?}; shift 2 ;;
        --runtime-prefix) runtime_prefix=${2:?}; shift 2 ;;
        --scanner) scanner=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --account) account=${2:?}; shift 2 ;;
        --partition) partition=${2:?}; shift 2 ;;
        --max-concurrent) max_concurrent=${2:?}; shift 2 ;;
        --array-size) array_size=${2:?}; shift 2 ;;
        --memory) memory=${2:?}; shift 2 ;;
        --time) wall_time=${2:?}; shift 2 ;;
        --minimum-free-gib) minimum_free_gib=${2:?}; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $threshold_list && -n $jaspar && -n $genome &&
   -n $anchor_evidence && -n $runtime_prefix ]] || {
    usage >&2
    exit 2
}
[[ $max_concurrent =~ ^[1-9][0-9]*$ &&
   $array_size =~ ^[1-9][0-9]*$ && $minimum_free_gib =~ ^[0-9]+$ ]] || {
    echo "E: Concurrency, array size, and free-space reserve are invalid." >&2
    exit 2
}
[[ $run_root == /data/sm718/* ]] || {
    echo "E: --run-root must be below /data/sm718." >&2
    exit 2
}

fasta_index=${fasta_index:-$genome.fai}
scanner=${scanner:-$source/pssm_scan_parquet}
duckdb="$runtime_prefix/duckdb/bin/duckdb"
rscript="$runtime_prefix/r/bin/Rscript"
for path in "$threshold_list" "$jaspar" "$genome" "$fasta_index" \
            "$anchor_evidence"; do
    [[ -f $path ]] || { echo "E: Required input not found: $path" >&2; exit 1; }
done
for executable in "$scanner" "$duckdb" "$rscript"; do
    [[ -x $executable ]] || { echo "E: Executable not found: $executable" >&2; exit 1; }
done

source=$(cd "$source" && pwd -P)
source_commit=$(git -C "$source" rev-parse HEAD)
if [[ -n $(git -C "$source" status --porcelain --untracked-files=normal) ]]; then
    echo "E: Sensitivity submission requires a clean source tree." >&2
    exit 1
fi
scanner_build=$($scanner --version-json)
python3 - "$scanner_build" "$source_commit" <<'PY'
import json
import sys

build = json.loads(sys.argv[1])
assert build["parquet_enabled"]
assert build["source_commit"] == sys.argv[2]
assert not build["source_dirty"]
PY

mkdir -p "$run_root/plan" "$run_root/logs" "$run_root/input" \
    "$run_root/tasks" "$run_root/final"
run_root=$(cd "$run_root" && pwd -P)
available_kib=$(df -Pk "$run_root" | awk 'NR == 2 { print $4 }')
required_kib=$((minimum_free_gib * 1024 * 1024))
if [[ ! $available_kib =~ ^[0-9]+$ ]] || (( available_kib < required_kib )); then
    echo "E: Insufficient or unknown free space at $run_root." >&2
    exit 1
fi
echo "I: Durable free space: $((available_kib / 1024 / 1024)) GiB." >&2

task_count=$(python3 "$source/scripts/manage_motif_threshold_calibration.py" \
    prepare-negative-sensitivity --run-root "$run_root" \
    --threshold-list "$threshold_list" --jaspar "$jaspar" \
    --anchor-evidence "$anchor_evidence" --genome "$genome" \
    --fasta-index "$fasta_index" --source "$source" --duckdb "$duckdb" \
    --selected-threshold 0 --source-minimum-score -20 --context-flank 150 \
    --minimum-class-fraction 0.01)
[[ $task_count =~ ^[1-9][0-9]*$ ]] || {
    echo "E: Invalid prepared task count: $task_count" >&2
    exit 1
}

submission_record="$run_root/plan/slurm_submission.tsv"
if [[ -e $submission_record && $dry_run -eq 0 ]]; then
    echo "E: Run already has a Slurm submission record: $submission_record" >&2
    exit 1
fi

array_jobs=()
array_start=0
array_chunk=0
dependency=()
while (( array_start < task_count )); do
    chunk_tasks=$((task_count - array_start))
    if (( chunk_tasks > array_size )); then
        chunk_tasks=$array_size
    fi
    array_end=$((chunk_tasks - 1))
    command=(
        sbatch --parsable --account="$account" --partition="$partition" --requeue
        --job-name="jaspar_neg_thr_${array_chunk}"
        --array="0-${array_end}%${max_concurrent}"
        --nodes=1 --ntasks=1 --cpus-per-task=1 --mem="$memory"
        --time="$wall_time" --signal=B:USR1@300 --chdir="$source"
        --output="$run_root/logs/sensitivity-%A_%a.out"
        --error="$run_root/logs/sensitivity-%A_%a.err"
        "${dependency[@]}"
        "$source/scripts/run_negative_threshold_sensitivity_slurm_task.sh"
        --run-root "$run_root" --scanner "$scanner" --duckdb "$duckdb"
        --rscript "$rscript" --source "$source" --task-offset "$array_start"
        --minimum-free-gib "$minimum_free_gib"
    )
    if [[ $dry_run -eq 1 ]]; then
        printf '%q ' "${command[@]}"; printf '\n'
        array_job="ARRAY_JOB_ID_${array_chunk}"
    else
        array_job=$("${command[@]}")
    fi
    array_jobs+=("$array_job")
    dependency=(--dependency="afterok:$array_job")
    array_start=$((array_start + chunk_tasks))
    array_chunk=$((array_chunk + 1))
done

final_command=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=jaspar_neg_thr_final --nodes=1 --ntasks=1 --cpus-per-task=1
    --mem=4G --time=00:30:00 --chdir="$source"
    --dependency="afterok:$array_job"
    --output="$run_root/logs/finalize-%j.out"
    --error="$run_root/logs/finalize-%j.err"
    "$source/scripts/run_negative_threshold_sensitivity_finalize.sh"
    --run-root "$run_root" --source "$source" --duckdb "$duckdb"
    --finalization-source-commit "$source_commit"
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${final_command[@]}"; printf '\n'
    exit 0
fi
final_job=$("${final_command[@]}")
array_jobs_text=$(IFS=,; printf '%s' "${array_jobs[*]}")
printf 'array_jobs\tfinalize_job\ttask_count\tmax_concurrent\tarray_size\tarray_chunks\tsource_commit\n' \
    > "$submission_record"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$array_jobs_text" "$final_job" "$task_count" "$max_concurrent" \
    "$array_size" "$array_chunk" "$source_commit" >> "$submission_record"
echo "I: Submitted arrays=$array_jobs_text finalizer=$final_job " \
     "tasks=$task_count live=$max_concurrent" >&2
printf '%s\n' "$array_jobs_text"

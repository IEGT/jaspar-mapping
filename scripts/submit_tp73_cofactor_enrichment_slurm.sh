#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_tp73_cofactor_enrichment_slurm.sh --run-root DIR
       --source-threshold-run DIR [options]

Prepare and submit a restart-safe all-JASPAR chromosome-1 TP73 cofactor
enrichment run. A setup job reconstructs the depth-bearing TP73 anchor table
from exact pinned bedGraphs. Each short array task reuses one cofactor-maxima
Parquet, stages it and the shared TP73/CUT&RUN anchors to node-local scratch,
and writes only compact statistics. Array chunks run serially with at most the
requested number of tasks live. A dependent finalizer validates all exact
outputs, performs the all-motif BH correction, and publishes Parquet + DuckDB.

Options:
  --run-root DIR             New durable run below /data/sm718
  --source-threshold-run DIR Completed threshold run below /data/sm718
  --source DIR               Repository root (default: script parent)
  --runtime-prefix DIR       Existing DuckDB/R runtime (default: source run)
  --account NAME             Slurm account (default: cluster)
  --partition NAME           Slurm partition (default: requeue)
  --max-concurrent N         Live tasks across the array (default: 20)
  --array-size N             Tasks per scheduler array (default: 1000)
  --cpus N                   CPUs per motif task (default: 1)
  --memory SIZE              Memory per motif task (default: 16G)
  --time LIMIT               Time per motif task (default: 00:30:00)
  --reuse-plan               Reuse an identical immutable prepared plan
  --dry-run                  Print sbatch commands without submitting
  -h, --help                 Show this help

The script never deletes source or durable task data. Requeued/repeated tasks
reuse a validated final task directory; interrupted attempts remain separately
named for inspection and never count as complete.
EOF
}

run_root=""
source_threshold_run=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
runtime_prefix=""
account=cluster
partition=requeue
max_concurrent=20
array_size=1000
cpus=1
memory=16G
wall_time=00:30:00
reuse_plan=0
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --source-threshold-run) source_threshold_run=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --runtime-prefix) runtime_prefix=${2:?}; shift 2 ;;
        --account) account=${2:?}; shift 2 ;;
        --partition) partition=${2:?}; shift 2 ;;
        --max-concurrent) max_concurrent=${2:?}; shift 2 ;;
        --array-size) array_size=${2:?}; shift 2 ;;
        --cpus) cpus=${2:?}; shift 2 ;;
        --memory) memory=${2:?}; shift 2 ;;
        --time) wall_time=${2:?}; shift 2 ;;
        --reuse-plan) reuse_plan=1; shift ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $source_threshold_run ]] || { usage >&2; exit 2; }
[[ $max_concurrent =~ ^[1-9][0-9]*$ &&
   $array_size =~ ^[1-9][0-9]*$ && $cpus =~ ^[1-9][0-9]*$ ]] || {
    echo "E: concurrency, array size, and CPUs must be positive integers." >&2
    exit 2
}
case "$run_root" in
    /data/sm718/*) ;;
    *) echo "E: --run-root must be below /data/sm718 on Haumea." >&2; exit 2 ;;
esac
case "$source_threshold_run" in
    /data/sm718/*) ;;
    *) echo "E: --source-threshold-run must be below /data/sm718." >&2; exit 2 ;;
esac

mkdir -p "$run_root/plan" "$run_root/logs" "$run_root/tasks" \
    "$run_root/final" "$run_root/input"
run_root=$(cd "$run_root" && pwd -P)
source_threshold_run=$(cd "$source_threshold_run" && pwd -P)
source=$(cd "$source" && pwd -P)
runtime_prefix=${runtime_prefix:-$source_threshold_run/runtime}
duckdb="$runtime_prefix/duckdb/bin/duckdb"
rscript="$runtime_prefix/r/bin/Rscript"
[[ -x $duckdb && -x $rscript ]] || {
    echo "E: Existing source-run DuckDB/R runtime is incomplete: $runtime_prefix" >&2
    exit 1
}
"$rscript" -e 'stopifnot(requireNamespace("data.table", quietly=TRUE))'

task_file="$run_root/plan/enrichment_tasks.tsv"
run_config="$run_root/plan/run_config.json"
setup_needed=0
if [[ $reuse_plan -eq 1 ]]; then
    [[ -s $task_file && -s $run_config ]] || {
        echo "E: --reuse-plan requires an existing immutable plan." >&2
        exit 1
    }
    plan_values=$(python3 - "$run_config" "$task_file" \
        "$source_threshold_run" "$source" <<'PY'
import csv
import json
import pathlib
import sys

config_path, tasks_path, source_run, source = map(pathlib.Path, sys.argv[1:])
config = json.loads(config_path.read_text())
expected = {
    "source_threshold_run": str(source_run.resolve()),
    "source": str(source.resolve()),
}
for key, value in expected.items():
    if config.get(key) != value:
        raise SystemExit(f"prepared plan {key} differs")
with tasks_path.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if len(rows) != int(config.get("task_count", -1)):
    raise SystemExit("prepared plan task count differs")
print(f'{len(rows)}\t{config["source_commit"]}')
PY
    )
    IFS=$'\t' read -r task_count plan_source_commit <<< "$plan_values"
    scientific_paths=(
        scripts/analyze_tp73_cofactor_enrichment.R
        scripts/build_tp73_anchor_evidence.py
        scripts/manage_tp73_cofactor_enrichment.py
        scripts/run_tp73_cofactor_enrichment_setup.sh
        scripts/run_tp73_cofactor_enrichment_slurm_task.sh
    )
    git -C "$source" cat-file -e "$plan_source_commit^{commit}" 2>/dev/null || {
        echo "E: Prepared plan source commit is unavailable." >&2
        exit 1
    }
    if ! git -C "$source" diff --quiet "$plan_source_commit..HEAD" -- \
        "${scientific_paths[@]}"; then
        echo "E: Scientific enrichment code changed since plan preparation." >&2
        exit 1
    fi
else
    [[ ! -e $task_file && ! -e $run_config ]] || {
        echo "E: Existing plan requires --reuse-plan." >&2
        exit 1
    }
    task_count=$(python3 - "$source_threshold_run/plan/calibration_tasks.tsv" <<'PY'
import csv
import sys
with open(sys.argv[1], newline="") as handle:
    print(sum(1 for _ in csv.DictReader(handle, delimiter="\t")))
PY
    )
    plan_source_commit=$(git -C "$source" rev-parse HEAD)
    setup_needed=1
fi
[[ $task_count =~ ^[1-9][0-9]*$ ]] || {
    echo "E: Invalid prepared task count: $task_count" >&2
    exit 1
}

execution_source_commit=$(git -C "$source" rev-parse HEAD)
execution_source_dirty=0
if ! git -C "$source" diff --quiet --ignore-submodules -- ||
   ! git -C "$source" diff --cached --quiet --ignore-submodules --; then
    execution_source_dirty=1
fi
if [[ $execution_source_dirty -eq 1 ]]; then
    echo "E: Scientific cluster execution requires a tracked-clean source tree." >&2
    exit 1
fi
finalization_provenance=(--finalization-source-commit "$execution_source_commit")

submission_record="$run_root/plan/slurm_submission.tsv"
if [[ -e $submission_record && $dry_run -eq 0 ]]; then
    echo "E: Run already has a Slurm submission record: $submission_record" >&2
    exit 1
fi

setup_job=reused
dependency=()
if [[ $setup_needed -eq 1 ]]; then
    setup_submission=(
        sbatch --parsable --account="$account" --partition="$partition" --requeue
        --job-name=jaspar_enrich_setup --nodes=1 --ntasks=1 --cpus-per-task=2
        --mem=16G --time=02:00:00 --signal=B:USR1@120 --chdir="$source"
        --output="$run_root/logs/setup-%j.out"
        --error="$run_root/logs/setup-%j.err"
        "$source/scripts/run_tp73_cofactor_enrichment_setup.sh"
        --run-root "$run_root" --source-threshold-run "$source_threshold_run"
        --source "$source" --source-commit "$execution_source_commit"
        --duckdb "$duckdb" --threads 2 --memory-limit 12GB
    )
    if [[ $dry_run -eq 1 ]]; then
        printf '%q ' "${setup_submission[@]}"; printf '\n'
        setup_job=SETUP_JOB_ID
    else
        setup_job=$("${setup_submission[@]}")
    fi
    dependency=(--dependency="afterok:$setup_job")
fi

array_jobs=()
array_start=0
array_chunk=0
while (( array_start < task_count )); do
    chunk_tasks=$((task_count - array_start))
    if (( chunk_tasks > array_size )); then
        chunk_tasks=$array_size
    fi
    array_end=$((chunk_tasks - 1))
    submission=(
        sbatch --parsable --account="$account" --partition="$partition" --requeue
        --job-name="jaspar_enrich_${array_chunk}"
        --array="0-${array_end}%${max_concurrent}"
        --nodes=1 --ntasks=1 --cpus-per-task="$cpus" --mem="$memory"
        --time="$wall_time" --signal=B:USR1@120 --chdir="$source"
        --output="$run_root/logs/enrichment-%A_%a.out"
        --error="$run_root/logs/enrichment-%A_%a.err"
    )
    if (( ${#dependency[@]} > 0 )); then
        submission+=("${dependency[@]}")
    fi
    submission+=(
        "$source/scripts/run_tp73_cofactor_enrichment_slurm_task.sh"
        --run-root "$run_root" --source-threshold-run "$source_threshold_run"
        --task-file "$task_file" --run-config "$run_config" --source "$source"
        --duckdb "$duckdb" --rscript "$rscript" --task-offset "$array_start"
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
array_jobs_text=$(IFS=,; printf '%s' "${array_jobs[*]}")

final_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=jaspar_enrich_final --nodes=1 --ntasks=1 --cpus-per-task=2
    --mem=16G --time=02:00:00 --chdir="$source"
    --dependency="afterok:$last_array_job"
    --output="$run_root/logs/finalize-%j.out"
    --error="$run_root/logs/finalize-%j.err"
    "$source/scripts/run_tp73_cofactor_enrichment_finalize.sh"
    --run-root "$run_root" --source "$source" --duckdb "$duckdb"
    "${finalization_provenance[@]}"
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${final_submission[@]}"; printf '\n'
    exit 0
fi
final_job=$("${final_submission[@]}")

printf 'setup_job\tarray_jobs\tfinalize_job\ttask_count\tmax_concurrent\tarray_size\tarray_chunks\tplan_source_commit\texecution_source_commit\treused_plan\n' \
    > "$submission_record"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$setup_job" "$array_jobs_text" "$final_job" "$task_count" "$max_concurrent" \
    "$array_size" "$array_chunk" "$plan_source_commit" \
    "$execution_source_commit" "$reuse_plan" >> "$submission_record"
echo "I: Submitted setup=$setup_job arrays=$array_jobs_text finalizer=$final_job tasks=$task_count live=$max_concurrent chunks=$array_chunk" >&2
printf '%s\n' "$array_jobs_text"

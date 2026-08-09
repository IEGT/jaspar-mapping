#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_motif_threshold_calibration_slurm_task.sh --run-root DIR
       --scan-package DIR --task-file FILE --anchor-evidence FILE [options]

Run one chromosome-1 motif-threshold task selected by SLURM_ARRAY_TASK_ID.
The task resolves exact immutable scan paths from the plan, builds TP73-local
maxima, and evaluates non-negative integer score thresholds. Scratch is local;
only a validated task directory is promoted to durable storage.

Options:
  --run-root DIR          Dedicated durable calibration run
  --scan-package DIR      Finalized JASPAR sparse-scan package
  --task-file FILE        Immutable calibration_tasks.tsv
  --anchor-evidence FILE  Fixed TP73/CUT&RUN anchor Parquet
  --source DIR            Repository root (default: script parent)
  --duckdb FILE           DuckDB CLI (default: duckdb)
  --rscript FILE          Rscript executable (default: Rscript)
  --task-offset N         Add N to the array task ID (default: 0)
  --threads N             DuckDB threads (default: 2)
  --memory-limit SIZE     DuckDB ceiling (default: 12GB)
  --max-temp-size SIZE    DuckDB scratch ceiling (default: 80GB)
  --minimum-class-fraction N  Fit floor for retained/absent anchors (default: 0.01)
  -h, --help              Show this help

Send SIGUSR1 to the batch process for one progress line. The signal does not
reach or interrupt the active DuckDB/R child.
EOF
}

run_root=""
scan_package=""
task_file=""
anchor_evidence=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
duckdb=duckdb
rscript=Rscript
task_offset=0
threads=2
memory_limit=12GB
max_temp_size=80GB
minimum_class_fraction=0.01
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --scan-package) scan_package=${2:?}; shift 2 ;;
        --task-file) task_file=${2:?}; shift 2 ;;
        --anchor-evidence) anchor_evidence=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --rscript) rscript=${2:?}; shift 2 ;;
        --task-offset) task_offset=${2:?}; shift 2 ;;
        --threads) threads=${2:?}; shift 2 ;;
        --memory-limit) memory_limit=${2:?}; shift 2 ;;
        --max-temp-size) max_temp_size=${2:?}; shift 2 ;;
        --minimum-class-fraction) minimum_class_fraction=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scan_package && -n $task_file && -n $anchor_evidence ]] || {
    usage >&2
    exit 2
}
[[ ${SLURM_ARRAY_TASK_ID:-} =~ ^[0-9]+$ ]] || {
    echo "E: SLURM_ARRAY_TASK_ID is required." >&2
    exit 2
}
[[ $task_offset =~ ^[0-9]+$ ]] || {
    echo "E: --task-offset must be a non-negative integer." >&2
    exit 2
}
[[ $threads =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --threads must be a positive integer." >&2
    exit 2
}
for path in "$run_root" "$scan_package"; do
    [[ -d $path ]] || { echo "E: Directory not found: $path" >&2; exit 1; }
done
for path in "$task_file" "$anchor_evidence"; do
    [[ -f $path ]] || { echo "E: File not found: $path" >&2; exit 1; }
done
[[ -x $duckdb ]] || duckdb=$(command -v "$duckdb" || true)
[[ -x $rscript ]] || rscript=$(command -v "$rscript" || true)
[[ -x $duckdb && -x $rscript ]] || {
    echo "E: DuckDB or Rscript executable is unavailable." >&2
    exit 1
}

global_task_index=$((task_offset + SLURM_ARRAY_TASK_ID))
task_row=$(awk -F '\t' -v task="$global_task_index" \
    'NR > 1 && $1 == task { print; found = 1; exit }
     END { if (!found) exit 1 }' "$task_file") || {
    echo "E: No task row for global task index $global_task_index." >&2
    exit 1
}
IFS=$'\t' read -r task_index motif_id motif_name plus_relative minus_relative \
    plus_bytes minus_bytes plus_sha minus_sha plus_hits minus_hits <<< "$task_row"
[[ $task_index == "$global_task_index" ]]
[[ $motif_id =~ ^[A-Za-z0-9._-]+$ ]]
for relative in "$plus_relative" "$minus_relative"; do
    [[ $relative != /* && $relative != *".."* ]] || {
        echo "E: Unsafe inventory-relative path: $relative" >&2
        exit 1
    }
done
plus="$scan_package/$plus_relative"
minus="$scan_package/$minus_relative"
[[ -f $plus && -f $minus ]] || {
    echo "E: Planned sparse scan input is missing for $motif_id." >&2
    exit 1
}
[[ $(stat -c %s "$plus") == "$plus_bytes" &&
   $(stat -c %s "$minus") == "$minus_bytes" ]] || {
    echo "E: Planned sparse scan size changed for $motif_id." >&2
    exit 1
}

final="$run_root/tasks/task-$(printf '%06d' "$task_index")-$motif_id"
if [[ -e $final ]]; then
    [[ -f $final/complete.json && -s $final/threshold_metrics.tsv ]] || {
        echo "E: Existing task output is incomplete: $final" >&2
        exit 1
    }
    recorded_motif=$(python3 -c \
        'import json,sys; print(json.load(open(sys.argv[1]))["motif_id"])' \
        "$final/complete.json")
    [[ $recorded_motif == "$motif_id" ]] || {
        echo "E: Existing task marker does not match $motif_id." >&2
        exit 1
    }
    echo "I: Reusing completed threshold task $task_index ($motif_id)." >&2
    exit 0
fi

start_epoch=$(date +%s)
phase=initializing
child_pid=""
progress_report() {
    local now elapsed
    now=$(date +%s)
    elapsed=$((now - start_epoch))
    printf 'I: progress signal=SIGUSR1 phase=%s task=%s array_task=%s task_offset=%s motif=%s elapsed_s=%s child_pid=%s maxima=%s metrics=%s\n' \
        "$phase" "$task_index" "$SLURM_ARRAY_TASK_ID" "$task_offset" \
        "$motif_id" "$elapsed" "${child_pid:-none}" \
        "$([[ -s ${maxima:-/nonexistent} ]] && echo ready || echo pending)" \
        "$([[ -s ${metrics:-/nonexistent} ]] && echo ready || echo pending)" >&2
}
trap progress_report USR1

run_child() {
    "$@" &
    child_pid=$!
    local status
    while true; do
        set +e
        wait "$child_pid"
        status=$?
        set -e
        if kill -0 "$child_pid" 2>/dev/null; then
            continue
        fi
        child_pid=""
        return "$status"
    done
}

scratch_base=${SLURM_TMPDIR:-/scratch/${USER:-sm718}}
scratch="$scratch_base/jaspar-threshold-${SLURM_JOB_ID:-manual}-${task_index}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
spill="$scratch/duckdb-spill"
mkdir -p "$spill"
maxima="$scratch/cofactor_maxima.parquet"
prefix="$scratch/evaluation"
metrics="${prefix}_threshold_metrics.tsv"

echo "I: Threshold task $task_index: $motif_id ($motif_name)" >&2
echo "I: Exact inputs: $plus_relative ; $minus_relative" >&2
echo "I: Node-local work: $scratch" >&2

phase=building_context_maxima
run_child "$source/scripts/build_sparse_context_maxima.py" \
    --anchor-parquet "$anchor_evidence" \
    --cofactor "$motif_id" "$plus" "$minus" \
    --output "$maxima" --flank 150 --source-score-floor -1 \
    --duckdb "$duckdb" --threads "$threads" --memory-limit "$memory_limit" \
    --max-temp-size "$max_temp_size" --temp-directory "$spill"

phase=fitting_threshold_models
run_child "$rscript" "$source/scripts/evaluate_tp73_cofactor_thresholds.R" \
    --anchor-evidence "$anchor_evidence" --cofactor-maxima "$maxima" \
    --output-prefix "$prefix" --thresholds auto --folds 5 \
    --comparison-mode fixed-negative-reference \
    --negative-reference-threshold -1 \
    --chrom-size 248956422 --spline-df 4 \
    --minimum-class-fraction "$minimum_class_fraction" --compact-output \
    --duckdb "$duckdb"

phase=validating
valid=$(
    "$duckdb" -csv -noheader :memory: -c "
SELECT count(*) > 0
   AND count(DISTINCT motif_id) = 1
   AND min(motif_id) = '$motif_id'
   AND min(threshold) = 0
   AND count(*) FILTER (WHERE evaluation_status IS NULL) = 0
FROM read_csv_auto('$metrics', delim='\\t', header=true, nullstr='NA');"
)
[[ $valid == true ]] || {
    echo "E: Threshold metrics failed validation for $motif_id." >&2
    exit 1
}

attempt="$run_root/tasks/.attempt-task-$(printf '%06d' "$task_index")-$motif_id-job-${SLURM_JOB_ID:-manual}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
[[ ! -e $attempt ]] || { echo "E: Durable attempt already exists: $attempt" >&2; exit 1; }
mkdir -p "$attempt"
cp "$maxima" "$attempt/cofactor_maxima.parquet"
cp "$maxima.run_config.tsv" "$attempt/cofactor_maxima.run_config.tsv"
cp "$metrics" "$attempt/threshold_metrics.tsv"
cp "${prefix}_baseline_metrics.tsv" "$attempt/baseline_metrics.tsv"
cp "${prefix}_fold_manifest.tsv" "$attempt/fold_manifest.tsv"
cp "${prefix}_run_config.tsv" "$attempt/evaluator_run_config.tsv"

phase=promoting
python3 -c '
import hashlib, json, pathlib, sys, time
output = pathlib.Path(sys.argv[1])
def digest(name):
    value = hashlib.sha256()
    with (output / name).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()
payload = {
    "schema_version": 1,
    "task_index": int(sys.argv[2]),
    "motif_id": sys.argv[3],
    "plus_inventory_sha256": sys.argv[4],
    "minus_inventory_sha256": sys.argv[5],
    "slurm_job_id": sys.argv[6],
    "completed_epoch": int(time.time()),
    "files": {name: {"bytes": (output / name).stat().st_size,
                      "sha256": digest(name)}
              for name in ("cofactor_maxima.parquet", "threshold_metrics.tsv",
                           "baseline_metrics.tsv", "fold_manifest.tsv",
                           "evaluator_run_config.tsv")},
}
(output / "complete.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
' "$attempt" "$task_index" "$motif_id" "$plus_sha" "$minus_sha" \
    "${SLURM_JOB_ID:-manual}"
[[ ! -e $final ]] || { echo "E: Final task appeared concurrently: $final" >&2; exit 1; }
mv "$attempt" "$final"
phase=complete
echo "I: Completed threshold task $task_index ($motif_id): $final" >&2

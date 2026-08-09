#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_negative_threshold_sensitivity_slurm_task.sh --run-root DIR
       --scanner FILE --duckdb FILE --rscript FILE [options]

Run one planned chromosome-1 negative-threshold sensitivity task selected by
SLURM_ARRAY_TASK_ID. The task stages chromosome 1 on node-local scratch, writes
one compact context-maximum Parquet directly from pssm_scan, evaluates every
integer threshold from the recorded source floor through zero, and promotes
only validated output to durable storage.

Options:
  --run-root DIR      Prepared sensitivity run below /data/sm718 (required)
  --scanner FILE      Arrow-enabled pssm_scan_parquet (required)
  --duckdb FILE       DuckDB CLI (required)
  --rscript FILE      Rscript executable (required)
  --source DIR        Repository root (default: script parent)
  --task-offset N     Add N to SLURM_ARRAY_TASK_ID (default: 0)
  --scratch-root DIR  Node-local parent (default: /scratch/USER)
  --minimum-free-gib N
                      Required durable free space before work (default: 5)
  -h, --help          Show this help

SIGUSR1 is forwarded to pssm_scan while scanning, producing its detailed
one-line report. During other phases the wrapper prints one task-level report.
EOF
}

run_root=""
scanner=""
duckdb=""
rscript=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
task_offset=0
minimum_free_gib=5
scratch_root="/scratch/${USER:-sm718}"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --scanner) scanner=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --rscript) rscript=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --task-offset) task_offset=${2:?}; shift 2 ;;
        --scratch-root) scratch_root=${2:?}; shift 2 ;;
        --minimum-free-gib) minimum_free_gib=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scanner && -n $duckdb && -n $rscript ]] || {
    usage >&2
    exit 2
}
[[ ${SLURM_ARRAY_TASK_ID:-} =~ ^[0-9]+$ ]] || {
    echo "E: SLURM_ARRAY_TASK_ID is required." >&2
    exit 2
}
[[ $task_offset =~ ^[0-9]+$ && $minimum_free_gib =~ ^[0-9]+$ ]] || {
    echo "E: task offset and free-space reserve must be non-negative integers." >&2
    exit 2
}
for executable in "$scanner" "$duckdb" "$rscript"; do
    [[ -x $executable ]] || { echo "E: Executable not found: $executable" >&2; exit 1; }
done
[[ -d $run_root && -d $source ]] || {
    echo "E: Run root or source repository is missing." >&2
    exit 1
}
mkdir -p "$scratch_root"
[[ -d $scratch_root && -w $scratch_root ]] || {
    echo "E: Writable scratch root is required: $scratch_root" >&2
    exit 1
}

task_file="$run_root/plan/calibration_tasks.tsv"
run_config="$run_root/plan/run_config.json"
[[ -f $task_file && -f $run_config ]] || {
    echo "E: Sensitivity plan is incomplete." >&2
    exit 1
}

global_task_index=$((10#$task_offset + 10#$SLURM_ARRAY_TASK_ID))
task_row=$(awk -F '\t' -v task="$global_task_index" \
    'NR > 1 && $1 == task { print; found = 1; exit }
     END { if (!found) exit 1 }' "$task_file") || {
    echo "E: No task row for global task index $global_task_index." >&2
    exit 1
}
IFS=$'\t' read -r task_index motif_id motif_name <<< "$task_row"
[[ $task_index == "$global_task_index" && $motif_id =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "E: Invalid task identity." >&2
    exit 1
}

config_values=$(python3 - "$run_config" <<'PY'
import json
import sys

config = json.load(open(sys.argv[1], encoding="utf-8"))
if config.get("analysis_mode") != "negative_threshold_sensitivity":
    raise SystemExit("run configuration is not a sensitivity plan")
keys = (
    "genome", "fasta_index", "jaspar", "anchor_evidence", "anchor_regions",
    "chrom", "chrom_size", "source_minimum_score", "context_flank_bp",
    "minimum_class_fraction", "anchor_count", "source_commit", "source_dirty",
)
print("\t".join(str(config[key]) for key in keys))
PY
)
IFS=$'\t' read -r genome fasta_index jaspar anchor_evidence anchor_regions \
    chrom chrom_size source_floor context_flank minimum_class_fraction \
    anchor_count source_commit source_dirty <<< "$config_values"
for path in "$genome" "$fasta_index" "$jaspar" "$anchor_evidence" "$anchor_regions"; do
    [[ -f $path ]] || { echo "E: Planned input is missing: $path" >&2; exit 1; }
done

scanner_build=$($scanner --version-json)
python3 - "$scanner_build" "$source_commit" "$source_dirty" <<'PY'
import json
import sys

build = json.loads(sys.argv[1])
if not build.get("parquet_enabled"):
    raise SystemExit("scanner lacks Parquet support")
expected_dirty = sys.argv[3].lower() == "true"
if (build.get("source_commit") != sys.argv[2]
        or bool(build.get("source_dirty")) != expected_dirty):
    raise SystemExit("scanner build does not match the planned source identity")
PY

available_kib=$(df -Pk "$run_root" | awk 'NR == 2 { print $4 }')
required_kib=$((minimum_free_gib * 1024 * 1024))
[[ $available_kib =~ ^[0-9]+$ ]] || {
    echo "E: Could not determine durable free space." >&2
    exit 1
}
if (( available_kib < required_kib )); then
    echo "E: Only $((available_kib / 1024 / 1024)) GiB free at $run_root; " \
         "$minimum_free_gib GiB is required before starting." >&2
    exit 1
fi

final="$run_root/tasks/task-$(printf '%06d' "$task_index")-$motif_id"
if [[ -e $final ]]; then
    [[ -f $final/complete.json && -s $final/threshold_metrics.tsv &&
       -s $final/cofactor_maxima.parquet ]] || {
        echo "E: Existing task output is incomplete: $final" >&2
        exit 1
    }
    recorded=$(python3 -c \
        'import json,sys; print(json.load(open(sys.argv[1]))["motif_id"])' \
        "$final/complete.json")
    [[ $recorded == "$motif_id" ]] || {
        echo "E: Existing task marker does not match $motif_id." >&2
        exit 1
    }
    echo "I: Reusing completed sensitivity task $task_index ($motif_id)." >&2
    exit 0
fi

start_epoch=$(date +%s)
phase=initializing
child_pid=""
child_accepts_usr1=0
progress_report() {
    local elapsed
    elapsed=$(($(date +%s) - start_epoch))
    if [[ $child_accepts_usr1 -eq 1 && -n $child_pid ]] &&
       kill -0 "$child_pid" 2>/dev/null; then
        kill -USR1 "$child_pid"
        return
    fi
    printf 'I: progress signal=SIGUSR1 phase=%s task=%s motif=%s elapsed_s=%s child_pid=%s\n' \
        "$phase" "$task_index" "$motif_id" "$elapsed" \
        "${child_pid:-none}" >&2
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

scratch="$scratch_root/jaspar-negative-threshold-${SLURM_JOB_ID:-manual}-${task_index}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
mkdir -p "$scratch"
export TMPDIR="$scratch"
staged_fasta="$scratch/chr${chrom}.fa"
maxima="$scratch/cofactor_maxima.parquet"
prefix="$scratch/evaluation"
metrics="${prefix}_threshold_metrics.tsv"

echo "I: Sensitivity task $task_index: $motif_id ($motif_name)" >&2
echo "I: Node-local work: $scratch" >&2

phase=staging_fasta
run_child python3 "$source/scripts/stage_fasta_region.py" \
    --genome "$genome" --fasta-index "$fasta_index" --sequence "$chrom" \
    --output "$staged_fasta" --expected-length "$chrom_size" \
    --metadata "$scratch/staged_fasta.json"

phase=scanning_context_maxima
child_accepts_usr1=1
run_child "$scanner" --context-maxima "$maxima" \
    --regions "$anchor_regions" --context-flank "$context_flank" \
    --genome "$staged_fasta" --fasta-index "$staged_fasta.fai" \
    --pssm "$jaspar" --motif "$motif_id" --strand both \
    --coordinate-mode bed --score-mode log2_relative_risk --pseudocount 1 \
    --threshold "$source_floor" --skip-N
child_accepts_usr1=0

phase=fitting_negative_thresholds
run_child "$rscript" "$source/scripts/evaluate_tp73_cofactor_thresholds.R" \
    --anchor-evidence "$anchor_evidence" --cofactor-maxima "$maxima" \
    --output-prefix "$prefix" --thresholds auto-source-floor --folds 5 \
    --comparison-mode threshold-complement \
    --chrom-size "$chrom_size" --spline-df 4 \
    --minimum-class-fraction "$minimum_class_fraction" --compact-output \
    --duckdb "$duckdb"

phase=validating
valid=$($duckdb -csv -noheader :memory: -c "
WITH maxima AS (SELECT * FROM read_parquet('$maxima')),
metrics AS (SELECT * FROM read_csv_auto(
    '$metrics', delim='\\t', header=true, nullstr='NA'))
SELECT (SELECT count(*) = $anchor_count
          AND count(DISTINCT motif_id) = 1
          AND min(motif_id) = '$motif_id'
          AND min(source_score_floor) = $source_floor
          AND max(source_score_floor) = $source_floor
          AND min(context_flank_bp) = $context_flank
          AND max(context_flank_bp) = $context_flank
          AND count(*) FILTER (
              WHERE context_score IS NOT NULL
                AND context_score < source_score_floor
          ) = 0
        FROM maxima)
   AND (SELECT count(*) = 1 - CAST($source_floor AS INTEGER)
          AND count(DISTINCT motif_id) = 1
          AND min(motif_id) = '$motif_id'
          AND min(threshold) = $source_floor
          AND max(threshold) = 0
          AND count(*) FILTER (WHERE evaluation_status IS NULL) = 0
        FROM metrics);"
)
[[ $valid == true ]] || {
    echo "E: Sensitivity task failed output validation for $motif_id." >&2
    exit 1
}

attempt="$run_root/tasks/.attempt-task-$(printf '%06d' "$task_index")-$motif_id-job-${SLURM_JOB_ID:-manual}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
[[ ! -e $attempt ]] || {
    echo "E: Durable attempt already exists: $attempt" >&2
    exit 1
}
mkdir -p "$attempt"
cp "$maxima" "$attempt/cofactor_maxima.parquet"
cp "$metrics" "$attempt/threshold_metrics.tsv"
cp "${prefix}_baseline_metrics.tsv" "$attempt/baseline_metrics.tsv"
cp "${prefix}_fold_manifest.tsv" "$attempt/fold_manifest.tsv"
cp "${prefix}_run_config.tsv" "$attempt/evaluator_run_config.tsv"
printf '%s\n' "$scanner_build" > "$attempt/scanner_version.json"
cp "$scratch/staged_fasta.json" "$attempt/staged_fasta.json"

phase=promoting
python3 - "$attempt" "$task_index" "$motif_id" "$source_commit" <<'PY'
import hashlib
import json
import pathlib
import sys
import time

output = pathlib.Path(sys.argv[1])

def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

files = {}
for path in sorted(output.iterdir()):
    files[path.name] = {"bytes": path.stat().st_size, "sha256": digest(path)}
payload = {
    "schema_version": 1,
    "analysis_mode": "negative_threshold_sensitivity",
    "task_index": int(sys.argv[2]),
    "motif_id": sys.argv[3],
    "source_commit": sys.argv[4],
    "slurm_job_id": str(__import__("os").environ.get("SLURM_JOB_ID", "manual")),
    "slurm_restart_count": int(__import__("os").environ.get("SLURM_RESTART_COUNT", 0)),
    "completed_epoch": int(time.time()),
    "files": files,
}
(output / "complete.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
PY
[[ ! -e $final ]] || {
    echo "E: Final task appeared concurrently: $final" >&2
    exit 1
}
mv "$attempt" "$final"
phase=complete
echo "I: Completed sensitivity task $task_index ($motif_id): $final" >&2

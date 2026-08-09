#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_tp73_cofactor_enrichment_slurm_task.sh --run-root DIR
       --source-threshold-run DIR --task-file FILE --run-config FILE [options]

Evaluate one planned cofactor motif against TP73 CUT&RUN evidence. The array
index plus --task-offset selects one exact immutable task row. The shared TP73
anchor evidence and that motif's existing cofactor-maxima Parquet are staged to
node-local scratch; only validated compact result tables are promoted.

Options:
  --run-root DIR             Dedicated durable enrichment run
  --source-threshold-run DIR Completed all-motif threshold-calibration run
  --task-file FILE           Immutable enrichment_tasks.tsv
  --run-config FILE          Immutable enrichment run_config.json
  --source DIR               Repository root (default: script parent)
  --duckdb FILE              DuckDB CLI (default: duckdb)
  --rscript FILE             Rscript executable (default: Rscript)
  --task-offset N            Add N to SLURM_ARRAY_TASK_ID (default: 0)
  --block-size BP            Genomic uncertainty block (default: 5000000)
  --spline-df N              TP73 score spline degrees of freedom (default: 4)
  --minimum-class-fraction N Primary class-support floor (default: 0.01)
  -h, --help                 Show this help

Send SIGUSR1 to the batch process for one progress line. It reports the phase,
global motif index, elapsed time, and active child without interrupting R.
EOF
}

run_root=""
source_threshold_run=""
task_file=""
run_config=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
duckdb=duckdb
rscript=Rscript
task_offset=0
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
        --task-offset) task_offset=${2:?}; shift 2 ;;
        --block-size) block_size=${2:?}; shift 2 ;;
        --spline-df) spline_df=${2:?}; shift 2 ;;
        --minimum-class-fraction) minimum_class_fraction=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $source_threshold_run && -n $task_file &&
   -n $run_config ]] || { usage >&2; exit 2; }
[[ ${SLURM_ARRAY_TASK_ID:-} =~ ^[0-9]+$ ]] || {
    echo "E: SLURM_ARRAY_TASK_ID is required." >&2
    exit 2
}
[[ $task_offset =~ ^[0-9]+$ && $block_size =~ ^[1-9][0-9]*$ &&
   $spline_df =~ ^[1-9][0-9]*$ ]] || {
    echo "E: task offset, block size, and spline df are invalid." >&2
    exit 2
}
for path in "$run_root" "$source_threshold_run" "$source"; do
    [[ -d $path ]] || { echo "E: Directory not found: $path" >&2; exit 1; }
done
for path in "$task_file" "$run_config"; do
    [[ -f $path ]] || { echo "E: File not found: $path" >&2; exit 1; }
done
[[ -x $duckdb ]] || duckdb=$(command -v "$duckdb" || true)
[[ -x $rscript ]] || rscript=$(command -v "$rscript" || true)
[[ -x $duckdb && -x $rscript ]] || {
    echo "E: DuckDB or Rscript executable is unavailable." >&2
    exit 1
}

planned_source=$(python3 - "$run_config" "$source" <<'PY'
import hashlib
import json
import pathlib
import sys

config = json.load(open(sys.argv[1], encoding="utf-8"))
source = pathlib.Path(sys.argv[2]).resolve()
expected = config.get("scientific_source_file_sha256")
if not isinstance(expected, dict) or not expected:
    raise SystemExit("immutable plan lacks scientific source checksums")
for relative_name, expected_digest in sorted(expected.items()):
    path = (source / relative_name).resolve()
    try:
        path.relative_to(source)
    except ValueError:
        raise SystemExit(f"unsafe scientific source path: {relative_name}")
    if not path.is_file():
        raise SystemExit(f"scientific source file is missing: {path}")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != expected_digest:
        raise SystemExit(f"scientific source checksum differs: {path}")
print("%s\t%s" % (
    config["source_commit"], str(bool(config["source_dirty"])).lower()
))
PY
)
IFS=$'\t' read -r planned_source_commit planned_source_dirty <<< "$planned_source"

sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

global_task_index=$((task_offset + SLURM_ARRAY_TASK_ID))
task_row=$(awk -F '\t' -v task="$global_task_index" \
    'NR > 1 && $1 == task { print; found = 1; exit }
     END { if (!found) exit 1 }' "$task_file") || {
    echo "E: No task row for global task index $global_task_index." >&2
    exit 1
}
IFS=$'\t' read -r task_index motif_id factor_name source_task_index \
    source_task_relative source_marker_sha source_maxima_bytes \
    source_maxima_sha _ positive_threshold _ positive_threshold_role _ _ \
    <<< "$task_row"
[[ $task_index == "$global_task_index" &&
   $motif_id =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "E: Invalid task identity in plan row." >&2
    exit 1
}
[[ $source_task_relative != /* && $source_task_relative != *".."* ]] || {
    echo "E: Unsafe source-task path: $source_task_relative" >&2
    exit 1
}
source_task="$source_threshold_run/$source_task_relative"
source_maxima="$source_task/cofactor_maxima.parquet"
source_marker="$source_task/complete.json"
anchor=$(python3 -c \
    'import json,sys; print(json.load(open(sys.argv[1]))["anchor_evidence"])' \
    "$run_config")
anchor_sha=$(python3 -c \
    'import json,sys; print(json.load(open(sys.argv[1]))["anchor_evidence_sha256"])' \
    "$run_config")
for path in "$source_maxima" "$source_marker" "$anchor"; do
    [[ -f $path ]] || { echo "E: Planned input not found: $path" >&2; exit 1; }
done
[[ $(wc -c < "$source_maxima" | tr -d ' ') == "$source_maxima_bytes" ]] || {
    echo "E: Source maxima size changed for $motif_id." >&2
    exit 1
}
[[ $(sha256_file "$source_marker") == "$source_marker_sha" &&
   $(sha256_file "$source_maxima") == "$source_maxima_sha" &&
   $(sha256_file "$anchor") == "$anchor_sha" ]] || {
    echo "E: Planned input checksum changed for $motif_id." >&2
    exit 1
}

final="$run_root/tasks/task-$(printf '%06d' "$task_index")-$motif_id"
if [[ -e $final ]]; then
    python3 - "$final/complete.json" "$task_index" "$motif_id" \
        "$source_maxima_sha" <<'PY'
import json
import pathlib
import sys

path = pathlib.Path(sys.argv[1])
value = json.loads(path.read_text())
if (value.get("task_index") != int(sys.argv[2])
        or value.get("motif_id") != sys.argv[3]
        or value.get("source_maxima_sha256") != sys.argv[4]):
    raise SystemExit("existing task marker has the wrong identity")
for name, record in value.get("files", {}).items():
    output = path.parent / name
    if not output.is_file() or output.stat().st_size != int(record["bytes"]):
        raise SystemExit(f"existing task output is incomplete: {output}")
PY
    echo "I: Reusing completed enrichment task $task_index ($motif_id)." >&2
    exit 0
fi

start_epoch=$(date +%s)
phase=initializing
child_pid=""
progress_report() {
    local elapsed
    elapsed=$(($(date +%s) - start_epoch))
    printf 'I: progress signal=SIGUSR1 phase=%s task=%s array_task=%s task_offset=%s motif=%s elapsed_s=%s child_pid=%s\n' \
        "$phase" "$task_index" "$SLURM_ARRAY_TASK_ID" "$task_offset" \
        "$motif_id" "$elapsed" "${child_pid:-none}" >&2
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
        # wait itself is interrupted when the batch shell handles SIGUSR1.
        # The child remains waitable even if it completed at that instant.
        if [[ $status -eq $usr1_wait_status ]]; then
            continue
        fi
        if kill -0 "$child_pid" 2>/dev/null; then
            continue
        fi
        child_pid=""
        return "$status"
    done
}

scratch_base=${SLURM_TMPDIR:-/scratch/${USER:-sm718}}
scratch="$scratch_base/jaspar-enrichment-${SLURM_JOB_ID:-manual}-${task_index}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
mkdir -p "$scratch"
staged_anchor="$scratch/tp73_anchor_evidence.parquet"
staged_maxima="$scratch/cofactor_maxima.parquet"
thresholds="$scratch/threshold.tsv"
prefix="$scratch/enrichment"

echo "I: Enrichment task $task_index: $motif_id ($factor_name)" >&2
echo "I: Source task $source_task_index: $source_task_relative" >&2
echo "I: Node-local work: $scratch" >&2

phase=staging_inputs
run_child cp "$anchor" "$staged_anchor"
run_child cp "$source_maxima" "$staged_maxima"
python3 - "$task_file" "$task_index" "$thresholds" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="") as handle:
    rows = csv.DictReader(handle, delimiter="\t")
    row = next(item for item in rows if item["task_index"] == sys.argv[2])
fields = [
    "motif_id", "positive_threshold", "factor_name",
    "positive_threshold_source", "selection_semantics",
]
with open(sys.argv[3], "w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t",
                            lineterminator="\n")
    writer.writeheader()
    writer.writerow({name: row[name] for name in fields})
PY

phase=evaluating
cd "$source"
run_child "$rscript" "$source/scripts/analyze_tp73_cofactor_enrichment.R" \
    --anchor-evidence "$staged_anchor" --cofactor-maxima "$staged_maxima" \
    --thresholds "$thresholds" --output-prefix "$prefix" \
    --negative-reference-thresholds "-1,0" \
    --primary-negative-reference -1 --block-size "$block_size" \
    --spline-df "$spline_df" \
    --minimum-class-fraction "$minimum_class_fraction" --duckdb "$duckdb" \
    --source-commit "$planned_source_commit" \
    --source-dirty "$planned_source_dirty"

phase=validating
valid=$(
    "$duckdb" -light-mode -csv -noheader :memory: -c "
WITH counts AS (
  SELECT
    (SELECT count(*) FROM read_csv_auto('${prefix}_class_counts.tsv', delim='\\t', header=true)) AS classes,
    (SELECT count(*) FROM read_csv_auto('${prefix}_depth_tier_manifest.tsv', delim='\\t', header=true)) AS depth,
    (SELECT count(*) FROM read_csv_auto('${prefix}_descriptive.tsv', delim='\\t', header=true)) AS descriptive,
    (SELECT count(*) FROM read_csv_auto('${prefix}_macro_summary.tsv', delim='\\t', header=true)) AS macro,
    (SELECT count(*) FROM read_csv_auto('${prefix}_primary_occupancy.tsv', delim='\\t', header=true)) AS primary_rows,
    (SELECT count(*) FROM read_csv_auto('${prefix}_run_config.tsv', delim='\\t', header=true)) AS config_rows,
    (SELECT count(DISTINCT motif_id) FROM read_csv_auto('${prefix}_primary_occupancy.tsv', delim='\\t', header=true)) AS motifs,
    (SELECT min(motif_id) FROM read_csv_auto('${prefix}_primary_occupancy.tsv', delim='\\t', header=true)) AS motif
)
SELECT classes = 18 AND depth > 0 AND depth % 12 = 0
   AND descriptive = depth * 9 AND macro = 108
   AND primary_rows = 1 AND config_rows = 1 AND motifs = 1
   AND motif = '$motif_id'
FROM counts;"
)
[[ $valid == true ]] || {
    echo "E: Enrichment outputs failed validation for $motif_id." >&2
    exit 1
}

attempt="$run_root/tasks/.attempt-task-$(printf '%06d' "$task_index")-$motif_id-job-${SLURM_JOB_ID:-manual}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
[[ ! -e $attempt ]] || {
    echo "E: Durable attempt already exists: $attempt" >&2
    exit 1
}
mkdir -p "$attempt"
cp "${prefix}_class_counts.tsv" "$attempt/class_counts.tsv"
cp "${prefix}_depth_tier_manifest.tsv" "$attempt/depth_tier_manifest.tsv"
cp "${prefix}_descriptive.tsv" "$attempt/descriptive.tsv"
cp "${prefix}_macro_summary.tsv" "$attempt/macro_summary.tsv"
cp "${prefix}_primary_occupancy.tsv" "$attempt/primary_occupancy.tsv"
cp "${prefix}_run_config.tsv" "$attempt/evaluator_run_config.tsv"

phase=promoting
python3 - "$attempt" "$task_index" "$motif_id" "$source_maxima_sha" \
    "$source_marker_sha" "$positive_threshold" "$positive_threshold_role" \
    "${SLURM_JOB_ID:-manual}" <<'PY'
import hashlib
import json
import pathlib
import sys
import time

output = pathlib.Path(sys.argv[1])
names = (
    "class_counts.tsv", "depth_tier_manifest.tsv", "descriptive.tsv",
    "macro_summary.tsv", "primary_occupancy.tsv", "evaluator_run_config.tsv",
)
def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()

payload = {
    "schema_version": 1,
    "task_index": int(sys.argv[2]),
    "motif_id": sys.argv[3],
    "source_maxima_sha256": sys.argv[4],
    "source_marker_sha256": sys.argv[5],
    "positive_threshold": float(sys.argv[6]),
    "positive_threshold_role": sys.argv[7],
    "slurm_job_id": sys.argv[8],
    "completed_epoch": int(time.time()),
    "files": {
        name: {"bytes": (output / name).stat().st_size,
               "sha256": digest(output / name)}
        for name in names
    },
}
(output / "complete.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n"
)
PY
[[ ! -e $final ]] || {
    echo "E: Final task appeared concurrently: $final" >&2
    exit 1
}
mv "$attempt" "$final"
phase=complete
echo "I: Completed enrichment task $task_index ($motif_id): $final" >&2

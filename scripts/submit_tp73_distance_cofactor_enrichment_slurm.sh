#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_tp73_distance_cofactor_enrichment_slurm.sh --run-root DIR
       --scan-package DIR --anchor-evidence FILE --thresholds FILE
       --threshold-set-id ID --jaspar-catalog DIR --duckdb FILE
       --run-id ID [options]

Prepare and submit the species-aware TP73 cofactor enrichment analysis by
isoform and exclusive interval-distance band. The primary task family defaults
to JASPAR vertebrate matrices. Each motif task processes the requested
chromosomes in order, staging only the current chromosome's exact Parquet files
on node-local scratch and publishing a durable checkpoint after each one.

Options:
  --run-root DIR          New durable run below /data/sm718
  --scan-package DIR      Finalized score-floor -1 whole-genome scan package
  --anchor-evidence FILE  Whole-genome TP73 CUT&RUN anchor-evidence Parquet
  --thresholds FILE       Motif operating-threshold registry Parquet
  --threshold-set-id ID   Registry threshold_set_id to apply
  --jaspar-catalog DIR    Official JASPAR metadata catalog directory
  --duckdb FILE           Absolute DuckDB CLI path available on compute nodes
  --run-id ID             Immutable analysis identifier
  --source DIR            Clean source checkout (default: script parent)
  --tax-group NAME        JASPAR tax_group or all (default: vertebrates)
  --chromosomes LIST      Comma-separated chromosome numbers (default: 1-22)
  --scratch-root DIR      Node-local parent (default: /scratch/$USER)
  --account NAME          Slurm account (default: cluster)
  --partition NAME        Slurm partition (default: requeue)
  --max-concurrent N      Maximum live motif tasks (default: 20)
  --array-size N          Maximum tasks per chained array (default: 1000)
  --cpus N                CPUs per task (default: 2)
  --memory SIZE           Memory per motif task (default: 32G)
  --memory-limit SIZE     DuckDB memory limit (default: 28GB)
  --max-temp-size SIZE    Per-task DuckDB spill ceiling (default: 64GB)
  --time LIMIT            Motif-task time limit (default: 04:00:00)
  --setup-memory SIZE     Anchor-split memory (default: 32G)
  --setup-time LIMIT      Anchor-split time limit (default: 02:00:00)
  --final-memory SIZE     Finalizer memory (default: 32G)
  --final-time LIMIT      Finalizer time limit (default: 02:00:00)
  --reuse-plan            Reuse an identical immutable prepared plan
  --dry-run               Prepare the plan and print sbatch commands
  -h, --help              Show this help

SIGUSR1 reports the current motif, chromosome, phase, and elapsed time. Slurm
requeue or a repeated task reuses validated chromosome checkpoints and final
motif outputs. Array chunks are serialized so concurrency remains globally
bounded. No source scan or checkpoint data are deleted.
EOF
}

run_root=""
scan_package=""
anchor_evidence=""
thresholds=""
threshold_set_id=""
jaspar_catalog=""
duckdb=""
run_id=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tax_group=vertebrates
chromosomes=$(seq -s, 1 22)
scratch_root="/scratch/${USER:-sm718}"
account=cluster
partition=requeue
max_concurrent=20
array_size=1000
cpus=2
memory=32G
memory_limit=28GB
max_temp_size=64GB
wall_time=04:00:00
setup_memory=32G
setup_time=02:00:00
final_memory=32G
final_time=02:00:00
reuse_plan=0
dry_run=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --scan-package) scan_package=${2:?}; shift 2 ;;
        --anchor-evidence) anchor_evidence=${2:?}; shift 2 ;;
        --thresholds) thresholds=${2:?}; shift 2 ;;
        --threshold-set-id) threshold_set_id=${2:?}; shift 2 ;;
        --jaspar-catalog) jaspar_catalog=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --run-id) run_id=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --tax-group) tax_group=${2:?}; shift 2 ;;
        --chromosomes) chromosomes=${2:?}; shift 2 ;;
        --scratch-root) scratch_root=${2:?}; shift 2 ;;
        --account) account=${2:?}; shift 2 ;;
        --partition) partition=${2:?}; shift 2 ;;
        --max-concurrent) max_concurrent=${2:?}; shift 2 ;;
        --array-size) array_size=${2:?}; shift 2 ;;
        --cpus) cpus=${2:?}; shift 2 ;;
        --memory) memory=${2:?}; shift 2 ;;
        --memory-limit) memory_limit=${2:?}; shift 2 ;;
        --max-temp-size) max_temp_size=${2:?}; shift 2 ;;
        --time) wall_time=${2:?}; shift 2 ;;
        --setup-memory) setup_memory=${2:?}; shift 2 ;;
        --setup-time) setup_time=${2:?}; shift 2 ;;
        --final-memory) final_memory=${2:?}; shift 2 ;;
        --final-time) final_time=${2:?}; shift 2 ;;
        --reuse-plan) reuse_plan=1; shift ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scan_package && -n $anchor_evidence &&
   -n $thresholds && -n $threshold_set_id && -n $jaspar_catalog &&
   -n $duckdb && -n $run_id ]] || { usage >&2; exit 2; }
for value in "$max_concurrent" "$array_size" "$cpus"; do
    [[ $value =~ ^[1-9][0-9]*$ ]] || {
        echo "E: concurrency, array size, and CPUs must be positive integers." >&2
        exit 2
    }
done
(( array_size <= 1000 )) || {
    echo "E: --array-size must not exceed 1000 on Haumea." >&2
    exit 2
}
for value in "$run_id" "$threshold_set_id" "$tax_group"; do
    [[ $value =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
        echo "E: run, threshold-set, and tax-group values contain unsafe characters." >&2
        exit 2
    }
done
[[ $chromosomes =~ ^[1-9][0-9]*(,[1-9][0-9]*)*$ ]] || {
    echo "E: --chromosomes must be comma-separated positive integers." >&2
    exit 2
}
case "$run_root" in
    /data/sm718/*) ;;
    *) echo "E: --run-root must be below /data/sm718." >&2; exit 2 ;;
esac
for path in "$scan_package" "$anchor_evidence" "$thresholds" \
            "$jaspar_catalog" "$source" "$duckdb"; do
    case "$path" in
        /data/sm718/*|/home/sm718/*) ;;
        *) echo "E: Production paths must be below /data/sm718 or /home/sm718: $path" >&2; exit 2 ;;
    esac
done

source=$(cd "$source" && pwd -P)
scan_package=$(cd "$scan_package" && pwd -P)
anchor_evidence=$(cd "$(dirname "$anchor_evidence")" && pwd -P)/$(basename "$anchor_evidence")
thresholds=$(cd "$(dirname "$thresholds")" && pwd -P)/$(basename "$thresholds")
jaspar_catalog=$(cd "$jaspar_catalog" && pwd -P)
duckdb=$(cd "$(dirname "$duckdb")" && pwd -P)/$(basename "$duckdb")
[[ -f $anchor_evidence && -f $thresholds && -x $duckdb ]] || {
    echo "E: Anchor evidence, thresholds, or DuckDB is missing." >&2
    exit 1
}
if ! git -C "$source" diff --quiet --ignore-submodules -- ||
   ! git -C "$source" diff --cached --quiet --ignore-submodules --; then
    echo "E: Production submission requires a tracked-clean source tree." >&2
    exit 1
fi
source_commit=$(git -C "$source" rev-parse HEAD)

if [[ $reuse_plan -eq 1 ]]; then
    config="$run_root/plan/run_config.json"
    tasks_file="$run_root/plan/tasks.tsv"
    [[ -s $config && -s $tasks_file ]] || {
        echo "E: --reuse-plan requires an existing plan." >&2
        exit 1
    }
    task_count=$(python3 - "$config" "$tasks_file" "$run_root" \
        "$scan_package" "$anchor_evidence" "$thresholds" "$threshold_set_id" \
        "$jaspar_catalog" "$source" "$source_commit" "$run_id" "$tax_group" \
        "$chromosomes" <<'PY'
import csv
import hashlib
import json
from pathlib import Path
import sys

(config_path, task_path, run_root, scan_package, anchors, thresholds,
 threshold_set_id, catalog, source, source_commit, run_id, tax_group,
 chromosomes) = sys.argv[1:]
config = json.loads(Path(config_path).read_text())
expected = {
    "scan_package": str(Path(scan_package).resolve()),
    "anchor_evidence": str(Path(anchors).resolve()),
    "thresholds": str(Path(thresholds).resolve()),
    "threshold_set_id": threshold_set_id,
    "jaspar_catalog": str(Path(catalog).resolve()),
    "source": str(Path(source).resolve()),
    "source_commit": source_commit,
    "run_id": run_id,
    "tax_group": tax_group,
    "chromosomes": chromosomes.split(","),
}
for key, value in expected.items():
    if config.get(key) != value:
        raise SystemExit(f"prepared plan {key} differs")
with Path(task_path).open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
if len(rows) != config.get("task_count"):
    raise SystemExit("prepared task count differs")
for filename, key in (("tasks.tsv", "tasks_sha256"),
                      ("scan_files.tsv", "scan_files_sha256")):
    path = Path(config_path).parent / filename
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != config.get(key):
        raise SystemExit(f"prepared {filename} checksum differs")
print(len(rows))
PY
    )
else
    task_count=$(
        "$source/scripts/manage_tp73_distance_cofactor_enrichment.py" prepare \
            --run-root "$run_root" --scan-package "$scan_package" \
            --anchor-evidence "$anchor_evidence" --thresholds "$thresholds" \
            --threshold-set-id "$threshold_set_id" \
            --jaspar-catalog "$jaspar_catalog" --source "$source" \
            --source-commit "$source_commit" --run-id "$run_id" \
            --tax-group "$tax_group" --chromosomes "$chromosomes" \
            --duckdb "$duckdb"
    )
fi
[[ $task_count =~ ^[1-9][0-9]*$ ]] || {
    echo "E: Invalid prepared task count: $task_count" >&2
    exit 1
}

submission_record="$run_root/plan/slurm_submission.tsv"
[[ ! -e $submission_record || $dry_run -eq 1 ]] || {
    echo "E: Run already has a submission record: $submission_record" >&2
    exit 1
}

setup_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=tp73_dist_setup --nodes=1 --ntasks=1 --cpus-per-task=2
    --mem="$setup_memory" --time="$setup_time" --chdir="$source"
    --output="$run_root/logs/setup-%j.out"
    --error="$run_root/logs/setup-%j.err"
    "$source/scripts/manage_tp73_distance_cofactor_enrichment.py"
    prepare-anchors --run-root "$run_root" --threads 2
    --memory-limit 28GB --duckdb "$duckdb"
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${setup_submission[@]}"; printf '\n'
    setup_job=SETUP_JOB_ID
else
    setup_job=$("${setup_submission[@]}")
fi

array_jobs=()
dependency=(--dependency="afterok:$setup_job")
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
        --job-name="tp73_dist_${array_chunk}"
        --array="0-${array_end}%${max_concurrent}"
        --nodes=1 --ntasks=1 --cpus-per-task="$cpus" --mem="$memory"
        --time="$wall_time" --signal=B:USR1@120 --open-mode=append
        --chdir="$source" "${dependency[@]}"
        --output="$run_root/logs/distance-%A_%a.out"
        --error="$run_root/logs/distance-%A_%a.err"
        "$source/scripts/manage_tp73_distance_cofactor_enrichment.py"
        run-task --run-root "$run_root" --task-offset "$array_start"
        --scratch "$scratch_root" --threads "$cpus"
        --memory-limit "$memory_limit" --max-temp-size "$max_temp_size"
        --duckdb "$duckdb"
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

final_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=tp73_dist_final --nodes=1 --ntasks=1 --cpus-per-task=2
    --mem="$final_memory" --time="$final_time" --chdir="$source"
    "${dependency[@]}" --output="$run_root/logs/finalize-%j.out"
    --error="$run_root/logs/finalize-%j.err"
    "$source/scripts/manage_tp73_distance_cofactor_enrichment.py"
    finalize --run-root "$run_root" --duckdb "$duckdb"
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${final_submission[@]}"; printf '\n'
    exit 0
fi
final_job=$("${final_submission[@]}")

array_jobs_text=$(IFS=,; printf '%s' "${array_jobs[*]}")
printf '%s\n' \
    $'setup_job\tarray_jobs\tfinalizer_job\ttask_count\tmax_concurrent\tarray_size\tarray_chunks\tsource_commit\trun_id\ttax_group\tchromosomes' \
    > "$submission_record"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$setup_job" "$array_jobs_text" "$final_job" "$task_count" \
    "$max_concurrent" "$array_size" "$array_chunk" "$source_commit" \
    "$run_id" "$tax_group" "$chromosomes" >> "$submission_record"
echo "I: Submitted setup=$setup_job arrays=$array_jobs_text " \
    "finalizer=$final_job tasks=$task_count tax_group=$tax_group" >&2
printf '%s\t%s\t%s\n' "$setup_job" "$array_jobs_text" "$final_job"

#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_tp73_genome_evidence_slurm.sh --run-root DIR
       --annotation-run DIR --scan-package DIR --track-root DIR
       --runtime-prefix DIR [options]

Submit one restart-safe task per sequence region to build TP73 and matched
negative-control CUT&RUN evidence. A dependent finalizer publishes autosomes as
the primary inference set, X/Y as a separate sensitivity set, and mitochondrial
evidence as a bystander-control set.

Options:
  --run-root DIR          Dedicated output below /data/sm718
  --annotation-run DIR    Finalized schema-7 TP73 annotation run
  --scan-package DIR      Finalized genome scan (chromosome identity source)
  --track-root DIR        Directory containing manifest-named BigWigs
  --runtime-prefix DIR    Runtime with DuckDB CLI and pyBigWig Python
  --source DIR            Repository root (default: script parent)
  --track-manifest FILE   Track manifest (default: source config)
  --scratch-root DIR      Node-local scratch parent (default: /scratch/$USER)
  --run-id ID             Immutable run identifier
  --account NAME          Slurm account (default: cluster)
  --partition NAME        Slurm partition (default: requeue)
  --max-concurrent N      Maximum live chromosome tasks (default: 20)
  --cpus N                CPUs per chromosome task (default: 4)
  --memory SIZE           Memory per chromosome task (default: 32G)
  --time LIMIT            Time per chromosome task (default: 02:00:00)
  --final-memory SIZE     Finalizer memory (default: 32G)
  --final-time LIMIT      Finalizer time (default: 02:00:00)
  --dry-run               Print sbatch commands without submitting
  -h, --help              Show this help

Completed chromosome packages are reused after contract validation. Interrupted
attempts remain immutable and never replace final output. SIGUSR1 reports the
current phase without terminating a chromosome worker.
EOF
}

run_root=""
annotation_run=""
scan_package=""
track_root=""
runtime_prefix=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
track_manifest=""
scratch_root="/scratch/${USER}"
run_id=jaspar2026_grch38_tp73_cutandrun_evidence_v1
account=cluster
partition=requeue
max_concurrent=20
cpus=4
memory=32G
wall_time=02:00:00
final_memory=32G
final_time=02:00:00
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --annotation-run) annotation_run=${2:?}; shift 2 ;;
        --scan-package) scan_package=${2:?}; shift 2 ;;
        --track-root) track_root=${2:?}; shift 2 ;;
        --runtime-prefix) runtime_prefix=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --track-manifest) track_manifest=${2:?}; shift 2 ;;
        --scratch-root) scratch_root=${2:?}; shift 2 ;;
        --run-id) run_id=${2:?}; shift 2 ;;
        --account) account=${2:?}; shift 2 ;;
        --partition) partition=${2:?}; shift 2 ;;
        --max-concurrent) max_concurrent=${2:?}; shift 2 ;;
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

[[ -n $run_root && -n $annotation_run && -n $scan_package &&
   -n $track_root && -n $runtime_prefix ]] || { usage >&2; exit 2; }
[[ $max_concurrent =~ ^[1-9][0-9]*$ && $cpus =~ ^[1-9][0-9]*$ ]] || {
    echo "E: concurrency and CPU values must be positive integers." >&2
    exit 2
}
[[ $run_id =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
    echo "E: --run-id contains unsafe characters." >&2
    exit 2
}
for path in "$run_root" "$annotation_run" "$scan_package" "$track_root" \
            "$runtime_prefix"; do
    case "$path" in
        /data/sm718/*) ;;
        *) echo "E: Production paths must be below /data/sm718: $path" >&2; exit 2 ;;
    esac
done

source=$(cd "$source" && pwd -P)
annotation_run=$(cd "$annotation_run" && pwd -P)
scan_package=$(cd "$scan_package" && pwd -P)
track_root=$(cd "$track_root" && pwd -P)
runtime_prefix=$(cd "$runtime_prefix" && pwd -P)
track_manifest=${track_manifest:-$source/config/h3k4me3_cutandrun_tracks_v1.tsv}
track_manifest=$(cd "$(dirname "$track_manifest")" && pwd -P)/$(basename "$track_manifest")
mkdir -p "$run_root/plan" "$run_root/logs"
run_root=$(cd "$run_root" && pwd -P)

duckdb="$runtime_prefix/duckdb/bin/duckdb"
bigwig_python="$runtime_prefix/duckdb/bin/python3"
[[ -x $duckdb && -x $bigwig_python ]] || {
    echo "E: Runtime lacks DuckDB or pyBigWig Python: $runtime_prefix" >&2
    exit 1
}
"$bigwig_python" -c 'import pyBigWig'
[[ -f $track_manifest && -f $annotation_run/final/manifest.json &&
   -f $scan_package/manifest.json ]] || {
    echo "E: A production manifest is missing." >&2
    exit 1
}
if ! git -C "$source" diff --quiet --ignore-submodules -- ||
   ! git -C "$source" diff --cached --quiet --ignore-submodules --; then
    echo "E: Production submission requires a tracked-clean source tree." >&2
    exit 1
fi

prepare=(
    python3 "$source/scripts/manage_tp73_genome_evidence.py" prepare
    --run-root "$run_root" --annotation-run "$annotation_run"
    --scan-package "$scan_package" --track-root "$track_root"
    --runtime-prefix "$runtime_prefix" --source "$source"
    --track-manifest "$track_manifest" --duckdb "$duckdb"
    --scratch-root "$scratch_root" --run-id "$run_id"
    --threads "$cpus" --memory-limit 28GB
)
task_count=$("${prepare[@]}")
[[ $task_count =~ ^[1-9][0-9]*$ ]] || {
    echo "E: Invalid prepared task count: $task_count" >&2
    exit 1
}
last_task=$((task_count - 1))
submission_record="$run_root/plan/slurm_submission.tsv"
if [[ -e $submission_record && $dry_run -eq 0 ]]; then
    echo "E: Run already has a submission record: $submission_record" >&2
    exit 1
fi

array_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=tp73_genome_evidence --nodes=1 --ntasks=1
    --cpus-per-task="$cpus" --mem="$memory" --time="$wall_time"
    --array="0-${last_task}%${max_concurrent}" --signal=B:USR1@120
    --open-mode=append --chdir="$source"
    --output="$run_root/logs/chromosome-%A_%a.out"
    --error="$run_root/logs/chromosome-%A_%a.err"
    "$source/scripts/manage_tp73_genome_evidence.py" run-task
    --run-root "$run_root"
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${array_submission[@]}"; printf '\n'
    array_job=ARRAY_JOB_ID
else
    array_job=$("${array_submission[@]}")
fi

final_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --dependency="afterok:${array_job}" --job-name=tp73_genome_finalize
    --nodes=1 --ntasks=1 --cpus-per-task=4 --mem="$final_memory"
    --time="$final_time" --open-mode=append --chdir="$source"
    --output="$run_root/logs/finalize-%j.out"
    --error="$run_root/logs/finalize-%j.err"
    "$source/scripts/run_tp73_genome_evidence_finalize.sh"
    --run-root "$run_root" --source "$source" --duckdb "$duckdb"
    --threads 4 --memory-limit 28GB
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${final_submission[@]}"; printf '\n'
    exit 0
fi
final_job=$("${final_submission[@]}")

source_commit=$(git -C "$source" rev-parse HEAD)
printf 'array_job_id\tfinalizer_job_id\ttask_count\tmax_concurrent\tsource_commit\trun_id\n' \
    > "$submission_record"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$array_job" "$final_job" \
    "$task_count" "$max_concurrent" "$source_commit" "$run_id" \
    >> "$submission_record"
echo "I: Submitted TP73 chromosome array $array_job and finalizer $final_job" >&2
printf '%s\t%s\n' "$array_job" "$final_job"

#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_motif_threshold_calibration_slurm.sh --run-root DIR
       --scan-package DIR --jaspar FILE --cutandrun-dir DIR [options]

Prepare and submit a requeueable chromosome-1 threshold calibration for every
non-target JASPAR motif. One array task reads exactly two inventory files for
one motif. A setup job reconstructs the shared TP73/CUT&RUN labels, and a final
dependency consolidates the finished metric tables into the threshold registry.
TP73 itself is deliberately excluded: its direct and tandem thresholds are
separate target policies, not a self-cofactor calibration.

Options:
  --run-root DIR          Dedicated durable run below /data/sm718
  --scan-package DIR      Finalized JASPAR 2026 sparse scan package
  --jaspar FILE           Exact JASPAR 2026 matrix source
  --cutandrun-dir DIR     Existing no-duplicates CUT&RUN resources
  --source DIR            Repository root (default: script parent)
  --runtime-prefix DIR    Micromamba runtime (default: RUN/runtime)
  --micromamba FILE       Micromamba executable (default: micromamba)
  --account NAME          Slurm account (default: cluster)
  --partition NAME        Slurm partition (default: requeue)
  --max-concurrent N      Live array tasks (default: 20)
  --cpus N                CPUs per motif task (default: 2)
  --memory SIZE           Memory per motif task (default: 16G)
  --time LIMIT            Time per motif task (default: 02:00:00)
  --dry-run               Print submissions without calling sbatch
  -h, --help              Show this help

The runtime is created once from conda-forge when absent and its explicit
package solution is recorded. Existing runtime, plan, anchor, and task outputs
are reused only after contract checks; source scan payloads are never changed.
EOF
}

run_root=""
scan_package=""
jaspar=""
cutandrun_dir=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
runtime_prefix=""
micromamba=${MAMBA_EXE:-micromamba}
account=cluster
partition=requeue
max_concurrent=20
cpus=2
memory=16G
wall_time=02:00:00
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --scan-package) scan_package=${2:?}; shift 2 ;;
        --jaspar) jaspar=${2:?}; shift 2 ;;
        --cutandrun-dir) cutandrun_dir=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --runtime-prefix) runtime_prefix=${2:?}; shift 2 ;;
        --micromamba) micromamba=${2:?}; shift 2 ;;
        --account) account=${2:?}; shift 2 ;;
        --partition) partition=${2:?}; shift 2 ;;
        --max-concurrent) max_concurrent=${2:?}; shift 2 ;;
        --cpus) cpus=${2:?}; shift 2 ;;
        --memory) memory=${2:?}; shift 2 ;;
        --time) wall_time=${2:?}; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scan_package && -n $jaspar && -n $cutandrun_dir ]] || {
    usage >&2
    exit 2
}
[[ $max_concurrent =~ ^[1-9][0-9]*$ && $cpus =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --max-concurrent and --cpus must be positive integers." >&2
    exit 2
}
mkdir -p "$run_root/plan" "$run_root/logs" "$run_root/input" "$run_root/tasks"
run_root=$(cd "$run_root" && pwd -P)
scan_package=$(cd "$scan_package" && pwd -P)
jaspar=$(cd "$(dirname "$jaspar")" && pwd -P)/$(basename "$jaspar")
cutandrun_dir=$(cd "$cutandrun_dir" && pwd -P)
source=$(cd "$source" && pwd -P)
runtime_prefix=${runtime_prefix:-$run_root/runtime}

if [[ ! -x $runtime_prefix/bin/duckdb || ! -x $runtime_prefix/bin/Rscript ||
      ! -x $runtime_prefix/bin/bigWigToBedGraph ]]; then
    [[ $dry_run -eq 0 ]] || {
        echo "E: Dry-run requires an existing runtime: $runtime_prefix" >&2
        exit 1
    }
    [[ $micromamba == */* ]] || micromamba=$(command -v "$micromamba" || true)
    [[ -x $micromamba ]] || { echo "E: Micromamba is unavailable." >&2; exit 1; }
    "$micromamba" create --yes --override-channels --channel conda-forge \
        --prefix "$runtime_prefix" \
        duckdb-cli=1.5.4 r-base=4.5.1 r-data.table=1.18.4 \
        ucsc-bigwigtobedgraph
    "$micromamba" list --prefix "$runtime_prefix" --explicit \
        > "$run_root/plan/runtime-explicit.txt"
fi
duckdb="$runtime_prefix/bin/duckdb"
rscript="$runtime_prefix/bin/Rscript"
bigwig_to_bedgraph="$runtime_prefix/bin/bigWigToBedGraph"
"$rscript" -e 'stopifnot(requireNamespace("data.table", quietly=TRUE))'

anchor="$run_root/input/tp73_chr1_anchor_evidence.parquet"
task_count=$(python3 "$source/scripts/manage_motif_threshold_calibration.py" prepare \
    --run-root "$run_root" --scan-package "$scan_package" --jaspar "$jaspar" \
    --anchor-evidence "$anchor" --source "$source" --duckdb "$duckdb")
[[ $task_count =~ ^[1-9][0-9]*$ ]] || {
    echo "E: Invalid prepared task count: $task_count" >&2
    exit 1
}
task_file="$run_root/plan/calibration_tasks.tsv"
submission_record="$run_root/plan/slurm_submission.tsv"
if [[ -e $submission_record && $dry_run -eq 0 ]]; then
    echo "E: Run already has a Slurm submission record: $submission_record" >&2
    exit 1
fi

setup_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=jaspar_thr_anchor --nodes=1 --ntasks=1 --cpus-per-task=4
    --mem=16G --time=02:00:00 --signal=B:USR1@300 --chdir="$source"
    --output="$run_root/logs/anchor-%j.out" --error="$run_root/logs/anchor-%j.err"
    "$source/scripts/run_motif_threshold_anchor_setup.sh"
    --run-root "$run_root" --scan-package "$scan_package"
    --cutandrun-dir "$cutandrun_dir" --output "$anchor" --source "$source"
    --duckdb "$duckdb" --bigwig-to-bedgraph "$bigwig_to_bedgraph"
    --threads 4 --memory-limit 12GB
)

dependency=()
setup_job="reused"
if [[ ! -s $anchor ]]; then
    if [[ $dry_run -eq 1 ]]; then
        printf '%q ' "${setup_submission[@]}"; printf '\n'
        setup_job=SETUP_JOB_ID
    else
        setup_job=$("${setup_submission[@]}")
    fi
    dependency=(--dependency="afterok:$setup_job")
fi

array_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=jaspar2026_thr --array="0-$((task_count - 1))%${max_concurrent}"
    --nodes=1 --ntasks=1 --cpus-per-task="$cpus" --mem="$memory"
    --time="$wall_time" --signal=B:USR1@300 --chdir="$source"
    --output="$run_root/logs/threshold-%A_%a.out"
    --error="$run_root/logs/threshold-%A_%a.err"
    "${dependency[@]}"
    "$source/scripts/run_motif_threshold_calibration_slurm_task.sh"
    --run-root "$run_root" --scan-package "$scan_package" --task-file "$task_file"
    --anchor-evidence "$anchor" --source "$source" --duckdb "$duckdb"
    --rscript "$rscript" --threads "$cpus" --memory-limit 12GB
    --max-temp-size 80GB --minimum-class-fraction 0.01
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${array_submission[@]}"; printf '\n'
    array_job=ARRAY_JOB_ID
else
    array_job=$("${array_submission[@]}")
fi

final_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=jaspar_thr_final --nodes=1 --ntasks=1 --cpus-per-task=2
    --mem=8G --time=01:00:00 --chdir="$source"
    --dependency="afterok:$array_job"
    --output="$run_root/logs/finalize-%j.out"
    --error="$run_root/logs/finalize-%j.err"
    "$source/scripts/run_motif_threshold_calibration_finalize.sh"
    --run-root "$run_root" --source "$source" --duckdb "$duckdb"
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${final_submission[@]}"; printf '\n'
    exit 0
fi
final_job=$("${final_submission[@]}")

printf 'setup_job\tarray_job\tfinalize_job\ttask_count\tmax_concurrent\n' \
    > "$submission_record"
printf '%s\t%s\t%s\t%s\t%s\n' \
    "$setup_job" "$array_job" "$final_job" "$task_count" "$max_concurrent" \
    >> "$submission_record"
echo "I: Submitted setup=$setup_job array=$array_job finalizer=$final_job tasks=$task_count live=$max_concurrent" >&2
printf '%s\n' "$array_job"

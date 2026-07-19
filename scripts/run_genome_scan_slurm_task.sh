#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_genome_scan_slurm_task.sh --run-root DIR --scanner FILE [OPTIONS]

Execute the task selected by SLURM_ARRAY_TASK_ID from an immutable genome-scan
plan. The shell is replaced by the Python coordinator so Slurm SIGUSR1 requests
are forwarded to the scanner and produce one progress line.

Options:
  --run-root DIR   Scan run created by manage_genome_scan.py prepare (required)
  --scanner FILE   Arrow-enabled pssm_scan_parquet executable (required)
  --manager FILE   Python scan coordinator (default: scripts/manage_genome_scan.py
                   beside the scanner's source tree)
  --duckdb FILE    DuckDB CLI executable or command name (default: duckdb)
  --task-offset N  Add N to SLURM_ARRAY_TASK_ID (default: 0)
  --duckdb-memory-limit SIZE
                    In-memory validation ceiling (default: 12GB)
  -h, --help       Show this help and exit
EOF
}

run_root=""
scanner=""
manager=""
duckdb=duckdb
task_offset=0
duckdb_memory_limit=12GB
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; run_root=$2; shift 2 ;;
        --scanner) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; scanner=$2; shift 2 ;;
        --manager) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; manager=$2; shift 2 ;;
        --duckdb) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; duckdb=$2; shift 2 ;;
        --task-offset) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; task_offset=$2; shift 2 ;;
        --duckdb-memory-limit) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; duckdb_memory_limit=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scanner ]] || { usage >&2; exit 2; }
[[ ${SLURM_ARRAY_TASK_ID:-} =~ ^[0-9]+$ ]] || {
    echo "E: SLURM_ARRAY_TASK_ID is not set; submit this script as a Slurm array." >&2
    exit 1
}
[[ $task_offset =~ ^[0-9]+$ ]] || {
    echo "E: --task-offset must be a non-negative integer." >&2
    exit 2
}
task_index=$((10#$task_offset + 10#$SLURM_ARRAY_TASK_ID))

if [[ -z $manager ]]; then
    scanner_directory=$(cd "$(dirname "$scanner")" && pwd)
    manager="$scanner_directory/scripts/manage_genome_scan.py"
fi
[[ -f $manager ]] || {
    echo "E: Scan coordinator does not exist: $manager" >&2
    exit 1
}
manager_arguments=(
    run-task
    --run-root "$run_root"
    --task-index "$task_index"
    --scanner "$scanner"
    --duckdb "$duckdb"
    --duckdb-memory-limit "$duckdb_memory_limit"
)

if [[ -d /scratch && -w /scratch ]]; then
    scratch_directory="/scratch/${USER:-unknown}/jaspar-${SLURM_JOB_ID:-local}-${SLURM_ARRAY_TASK_ID}-${SLURM_RESTART_COUNT:-0}-$$"
    mkdir -p "$scratch_directory/duckdb"
    export TMPDIR="$scratch_directory"
    manager_arguments+=(--duckdb-temp-directory "$scratch_directory/duckdb")
    echo "I: Using disposable job-local temporary storage: $scratch_directory" >&2
fi

exec python3 "$manager" \
    "${manager_arguments[@]}"

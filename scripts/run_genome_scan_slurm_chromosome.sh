#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_genome_scan_slurm_chromosome.sh --run-root DIR --scanner FILE [OPTIONS]

Stage the chromosome selected by SLURM_ARRAY_TASK_ID into this node's /scratch,
verify it against the immutable plan, and run every planned motif batch for it.
Completed batches remain independent durable tasks and are reused after requeue.

Options:
  --run-root DIR        Scan run created by manage_genome_scan.py prepare
  --scanner FILE        Arrow-enabled pssm_scan_parquet executable
  --manager FILE        Coordinator (default: beside scanner source tree)
  --duckdb FILE         DuckDB CLI executable or command name (default: duckdb)
  --chrom-offset N      Add N to SLURM_ARRAY_TASK_ID (default: 0)
  --duckdb-memory-limit SIZE
                        Validation ceiling (default: 12GB)
  --minimum-scratch-free-gib N
                        Scratch reserve checked before every batch (default: 5)
  --batch-workers N     Concurrent motif batches sharing staged FASTA (default: 1)
  --scratch-task-output Validate task output in scratch before durable copying;
                        copies remain until the cluster clears the job scratch
  -h, --help            Show this help and exit
EOF
}

run_root=""
scanner=""
manager=""
duckdb=duckdb
chrom_offset=0
duckdb_memory_limit=12GB
minimum_scratch_free_gib=5
batch_workers=1
scratch_task_output=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; run_root=$2; shift 2 ;;
        --scanner) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; scanner=$2; shift 2 ;;
        --manager) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; manager=$2; shift 2 ;;
        --duckdb) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; duckdb=$2; shift 2 ;;
        --chrom-offset) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; chrom_offset=$2; shift 2 ;;
        --duckdb-memory-limit) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; duckdb_memory_limit=$2; shift 2 ;;
        --minimum-scratch-free-gib) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; minimum_scratch_free_gib=$2; shift 2 ;;
        --batch-workers) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; batch_workers=$2; shift 2 ;;
        --scratch-task-output) scratch_task_output=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scanner ]] || { usage >&2; exit 2; }
[[ ${SLURM_ARRAY_TASK_ID:-} =~ ^[0-9]+$ ]] || {
    echo "E: SLURM_ARRAY_TASK_ID is not set." >&2
    exit 1
}
[[ $chrom_offset =~ ^[0-9]+$ ]] || {
    echo "E: --chrom-offset must be a non-negative integer." >&2
    exit 2
}
[[ -d /scratch && -w /scratch ]] || {
    echo "E: Writable node-local /scratch is required for chromosome staging." >&2
    exit 1
}

chrom_index=$((10#$chrom_offset + 10#$SLURM_ARRAY_TASK_ID))
if [[ -z $manager ]]; then
    scanner_directory=$(cd "$(dirname "$scanner")" && pwd)
    manager="$scanner_directory/scripts/manage_genome_scan.py"
fi
[[ -f $manager ]] || {
    echo "E: Scan coordinator does not exist: $manager" >&2
    exit 1
}

scratch_directory="/scratch/${USER:-unknown}/jaspar-${SLURM_JOB_ID:-local}-${SLURM_ARRAY_TASK_ID}-${SLURM_RESTART_COUNT:-0}-$$"
mkdir -p "$scratch_directory"
export TMPDIR="$scratch_directory"
arguments=(
    run-chromosome
    --run-root "$run_root"
    --chrom-index "$chrom_index"
    --scanner "$scanner"
    --scratch-directory "$scratch_directory/chromosome-worker"
    --duckdb "$duckdb"
    --duckdb-memory-limit "$duckdb_memory_limit"
    --minimum-scratch-free-gib "$minimum_scratch_free_gib"
    --batch-workers "$batch_workers"
)
if [[ $scratch_task_output -eq 1 ]]; then
    arguments+=(--scratch-task-output)
fi

echo "I: Using disposable node-local storage: $scratch_directory" >&2
exec python3 "$manager" "${arguments[@]}"

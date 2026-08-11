#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_h3k4me3_cofactor_analysis_finalize.sh --run-root DIR --source DIR
       --duckdb FILE [options]

Validate every motif checkpoint, combine compact result tables, recompute BH
correction across all planned motifs, and publish the DuckDB result catalog.

Options:
  --run-root DIR       Prepared all-motif H3K4me3 analysis run
  --source DIR         Repository root
  --duckdb FILE        DuckDB CLI
  --threads N          DuckDB threads (default: 4)
  --memory-limit SIZE  DuckDB memory ceiling (default: 28GB)
  -h, --help           Show this help
EOF
}

run_root=""
source=""
duckdb=""
threads=4
memory_limit=28GB
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --threads) threads=${2:?}; shift 2 ;;
        --memory-limit) memory_limit=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n $run_root && -n $source && -n $duckdb ]] || { usage >&2; exit 2; }
[[ $threads =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --threads must be a positive integer." >&2
    exit 2
}

scratch_base=${SLURM_TMPDIR:-/scratch/${USER:-sm718}}
scratch="$scratch_base/h3-cofactor-finalize-${SLURM_JOB_ID:-manual}-$$"
mkdir -p "$scratch"

exec python3 "$source/scripts/manage_h3k4me3_cofactor_analysis.py" finalize \
    --run-root "$run_root" --duckdb "$duckdb" --threads "$threads" \
    --memory-limit "$memory_limit" --temp-directory "$scratch"

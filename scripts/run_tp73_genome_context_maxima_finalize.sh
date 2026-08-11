#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_tp73_genome_context_maxima_finalize.sh --run-root DIR
       --source DIR --duckdb FILE [--verify-checksums]

Validate every autosomal motif-context task and publish the compact Parquet and
DuckDB file inventory. Large per-motif maxima remain in their atomic task
directories and are not copied by finalization.

Options:
  --run-root DIR       Prepared context-maxima run
  --source DIR         Repository root
  --duckdb FILE        DuckDB CLI
  --verify-checksums   Re-read every large maxima payload during finalization
  -h, --help           Show this help
EOF
}

run_root=""
source=""
duckdb=""
verify=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --verify-checksums) verify=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n $run_root && -n $source && -n $duckdb ]] || { usage >&2; exit 2; }

arguments=(
    python3 "$source/scripts/manage_tp73_genome_context_maxima.py" finalize
    --run-root "$run_root" --duckdb "$duckdb"
)
if [[ $verify -eq 1 ]]; then
    arguments+=(--verify-checksums)
fi
exec "${arguments[@]}"

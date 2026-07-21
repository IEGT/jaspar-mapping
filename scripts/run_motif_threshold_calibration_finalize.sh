#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_motif_threshold_calibration_finalize.sh --run-root DIR [options]

Consolidate every exact per-motif metric table and build the immutable,
provenance-rich motif_score_threshold Parquet registry. Finalization refuses
to proceed while any planned task is incomplete.

Options:
  --run-root DIR  Prepared calibration run (required)
  --source DIR    Repository root (default: script parent)
  --duckdb FILE   DuckDB CLI (default: duckdb)
  -h, --help      Show this help
EOF
}

run_root=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
duckdb=duckdb
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n $run_root ]] || { usage >&2; exit 2; }

exec python3 "$source/scripts/manage_motif_threshold_calibration.py" finalize \
    --run-root "$run_root" --duckdb "$duckdb"

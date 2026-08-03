#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_motif_threshold_calibration_finalize.sh --run-root DIR [options]

Consolidate every exact per-motif metric table and build the immutable,
provenance-rich motif_score_threshold Parquet registry. Finalization refuses
to proceed while any planned task is incomplete.

Options:
  --run-root DIR                  Prepared calibration run (required)
  --source DIR                    Repository root (default: script parent)
  --duckdb FILE                   DuckDB CLI (default: duckdb)
  --finalization-source-commit ID Source commit resolved before submission
                                    (required; 40 lowercase hexadecimal digits)
  --finalization-source-dirty     Record tracked source changes at submission
  -h, --help                      Show this help
EOF
}

run_root=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
duckdb=duckdb
finalization_source_commit=""
finalization_source_dirty=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --finalization-source-commit)
            finalization_source_commit=${2:?}; shift 2 ;;
        --finalization-source-dirty) finalization_source_dirty=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n $run_root && $finalization_source_commit =~ ^[0-9a-f]{40}$ ]] || {
    usage >&2
    exit 2
}

dirty_option=()
if [[ $finalization_source_dirty -eq 1 ]]; then
    dirty_option=(--finalization-source-dirty)
fi

exec python3 "$source/scripts/manage_motif_threshold_calibration.py" finalize \
    --run-root "$run_root" --duckdb "$duckdb" \
    --finalization-source-commit "$finalization_source_commit" \
    "${dirty_option[@]}"

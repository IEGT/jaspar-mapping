#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_negative_threshold_sensitivity_finalize.sh --run-root DIR [options]

Consolidate every completed per-motif negative-threshold sweep into one metric
table and a non-normative comparison of the best negative candidate with the
original zero threshold.

Options:
  --run-root DIR                  Prepared sensitivity run (required)
  --source DIR                    Repository root (default: script parent)
  --duckdb FILE                   DuckDB CLI (default: duckdb)
  --finalization-source-commit ID Full source commit resolved at submission
                                    (required)
  --finalization-source-dirty     Record tracked changes at submission
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
exec python3 "$source/scripts/manage_motif_threshold_calibration.py" \
    finalize-negative-sensitivity --run-root "$run_root" --duckdb "$duckdb" \
    --finalization-source-commit "$finalization_source_commit" \
    "${dirty_option[@]}"

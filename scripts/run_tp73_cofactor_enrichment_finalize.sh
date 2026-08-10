#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_tp73_cofactor_enrichment_finalize.sh --run-root DIR [options]

Validate every exact cofactor task, apply one Benjamini-Hochberg correction
across the complete planned motif family, and atomically publish the combined
Parquet tables plus their DuckDB query index.

Options:
  --run-root DIR                  Prepared enrichment run (required)
  --source DIR                    Repository root (default: script parent)
  --duckdb FILE                   DuckDB CLI (default: duckdb)
  --finalization-source-commit ID Full source commit fixed at submission
  --finalization-source-dirty     Record tracked changes at submission
  --publication-run-id ID         Corrected published identity (optional)
  --final-name NAME               Distinct directory below RUN/final
  -h, --help                      Show this help
EOF
}

run_root=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
duckdb=duckdb
finalization_source_commit=""
finalization_source_dirty=0
publication_run_id=""
final_name="cofactor_enrichment"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --finalization-source-commit)
            finalization_source_commit=${2:?}; shift 2 ;;
        --finalization-source-dirty) finalization_source_dirty=1; shift ;;
        --publication-run-id) publication_run_id=${2:?}; shift 2 ;;
        --final-name) final_name=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n $run_root && $finalization_source_commit =~ ^[0-9a-f]{40}$ ]] || {
    usage >&2
    exit 2
}

command=(
    python3 "$source/scripts/manage_tp73_cofactor_enrichment.py" finalize
    --run-root "$run_root" --duckdb "$duckdb"
    --finalization-source-commit "$finalization_source_commit"
    --final-name "$final_name"
)
if [[ -n $publication_run_id ]]; then
    command+=(--publication-run-id "$publication_run_id")
fi
if [[ $finalization_source_dirty -eq 1 ]]; then
    command+=(--finalization-source-dirty)
fi
exec "${command[@]}"

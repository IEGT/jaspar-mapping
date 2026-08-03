#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_motif_density_calibration_finalize.sh --run-root DIR [OPTIONS]

Finalize a complete motif-density calibration into an exact distribution
manifest and the per-motif TSV/JSON threshold registry.

Options:
  --run-root DIR  Prepared durable calibration run
  --source DIR    Repository root (default: script parent)
  -h, --help      Show this help and exit
EOF
}

run_root=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; run_root=$2; shift 2 ;;
        --source) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; source=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n $run_root ]] || { usage >&2; exit 2; }

exec python3 "$source/scripts/manage_motif_density_calibration.py" finalize \
    --run-root "$run_root"

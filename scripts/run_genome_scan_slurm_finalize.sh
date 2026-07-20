#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_genome_scan_slurm_finalize.sh --run-root DIR [OPTIONS]

Publish metadata and a metadata-only DuckDB catalog after all chromosome jobs
complete. This fast finalizer does not recompute payload checksums.

Options:
  --run-root DIR  Scan run created by manage_genome_scan.py prepare
  --manager FILE  Coordinator (otherwise resolve from script or working tree)
  --duckdb FILE   DuckDB CLI executable or command name (default: duckdb)
  -h, --help      Show this help and exit
EOF
}

run_root=""
script_directory=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
manager=""
duckdb=duckdb
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; run_root=$2; shift 2 ;;
        --manager) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; manager=$2; shift 2 ;;
        --duckdb) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; duckdb=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root ]] || { usage >&2; exit 2; }
if [[ -z $manager ]]; then
    if [[ -f $script_directory/manage_genome_scan.py ]]; then
        manager="$script_directory/manage_genome_scan.py"
    elif [[ -f $PWD/scripts/manage_genome_scan.py ]]; then
        # Slurm executes a spooled copy of this wrapper, while --chdir points
        # at the submitted source tree.
        manager="$PWD/scripts/manage_genome_scan.py"
    else
        echo "E: Cannot locate manage_genome_scan.py; pass --manager explicitly." >&2
        exit 1
    fi
fi
[[ -f $manager ]] || { echo "E: Scan coordinator does not exist: $manager" >&2; exit 1; }
exec python3 "$manager" finalize --run-root "$run_root" --duckdb "$duckdb"

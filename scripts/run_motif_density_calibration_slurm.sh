#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_motif_density_calibration_slurm.sh --run-root DIR --scanner FILE [OPTIONS]

Run a prepared chromosome-1 motif-density calibration. The chromosome is
copied once to this allocation's node-local /scratch directory; all remaining
motif batches reuse that copy. Completed batches remain durable below RUN_ROOT
and are reused after Slurm requeue.

Options:
  --run-root DIR        Prepared durable calibration run
  --scanner FILE        pssm_scan executable from the planned source commit
  --source DIR          Repository root (default: script parent)
  --batch-workers N     Concurrent scanner batches (default: 4)
  --minimum-scratch-free-gib N
                        Scratch reserve after chromosome staging (default: 5)
  -h, --help            Show this help and exit
EOF
}

run_root=""
scanner=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
batch_workers=4
minimum_scratch_free_gib=5
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; run_root=$2; shift 2 ;;
        --scanner) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; scanner=$2; shift 2 ;;
        --source) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; source=$2; shift 2 ;;
        --batch-workers) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; batch_workers=$2; shift 2 ;;
        --minimum-scratch-free-gib) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; minimum_scratch_free_gib=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scanner ]] || { usage >&2; exit 2; }
[[ ${SLURM_JOB_ID:-} =~ ^[0-9]+$ ]] || {
    echo "E: This wrapper requires a Slurm allocation." >&2
    exit 1
}
[[ -d /scratch && -w /scratch ]] || {
    echo "E: Writable node-local /scratch is required." >&2
    exit 1
}
[[ $batch_workers =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --batch-workers must be a positive integer." >&2
    exit 2
}

scratch_directory="/scratch/${USER:-unknown}/jaspar-density-${SLURM_JOB_ID}-${SLURM_RESTART_COUNT:-0}-$$"
mkdir -p "$scratch_directory"
export TMPDIR="$scratch_directory"
echo "I: Using disposable node-local storage: $scratch_directory" >&2

exec python3 "$source/scripts/manage_motif_density_calibration.py" run \
    --run-root "$run_root" --scanner "$scanner" \
    --scratch-directory "$scratch_directory/worker" \
    --batch-workers "$batch_workers" \
    --minimum-scratch-free-gib "$minimum_scratch_free_gib"

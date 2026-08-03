#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_motif_density_calibration_slurm.sh --run-root DIR --scanner FILE [OPTIONS]

Submit one requeue-safe chromosome-1 calibration worker and a dependent fast
finalizer. The worker stages chromosome 1 once and checkpoints every completed
motif batch below RUN_ROOT.

Options:
  --run-root DIR       Prepared run below /data/sm718
  --scanner FILE       pssm_scan executable from the planned source commit
  --source DIR         Repository root (default: script parent)
  --account NAME       Slurm account (default: cluster)
  --partition NAME     Slurm partition (default: requeue)
  --batch-workers N    Concurrent scanner batches (default: 4)
  --memory SIZE        Worker memory (default: 32G)
  --time LIMIT         Worker time limit (default: 0-08:00:00)
  --resume             Append a replacement worker after a prior submission;
                       completed batches will be reused
  --dry-run            Print sbatch commands without submitting
  -h, --help           Show this help and exit
EOF
}

run_root=""
scanner=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
account=cluster
partition=requeue
batch_workers=4
memory=32G
wall_time=0-08:00:00
resume=0
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; run_root=$2; shift 2 ;;
        --scanner) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; scanner=$2; shift 2 ;;
        --source) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; source=$2; shift 2 ;;
        --account) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; account=$2; shift 2 ;;
        --partition) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; partition=$2; shift 2 ;;
        --batch-workers) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; batch_workers=$2; shift 2 ;;
        --memory) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; memory=$2; shift 2 ;;
        --time) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; wall_time=$2; shift 2 ;;
        --resume) resume=1; shift ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scanner ]] || { usage >&2; exit 2; }
[[ $run_root == /data/sm718/* ]] || {
    echo "E: --run-root must be below /data/sm718." >&2
    exit 2
}
[[ -f $run_root/plan/density_plan.json ]] || {
    echo "E: Prepared density plan not found below $run_root." >&2
    exit 1
}
[[ -x $scanner ]] || { echo "E: Scanner is not executable: $scanner" >&2; exit 1; }
[[ -x $source/scripts/run_motif_density_calibration_slurm.sh &&
   -x $source/scripts/run_motif_density_calibration_finalize.sh ]] || {
    echo "E: Density-calibration wrappers are missing below $source." >&2
    exit 1
}
[[ $batch_workers =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --batch-workers must be a positive integer." >&2
    exit 2
}
mkdir -p "$run_root/logs"
submission_record="$run_root/plan/slurm_submission.tsv"
if [[ -e $submission_record && $dry_run -eq 0 && $resume -eq 0 ]]; then
    echo "E: Run already has a submission record; inspect it and use --resume: $submission_record" >&2
    exit 1
fi

plan_commit=$(python3 -c '
import json, sys
print(json.load(open(sys.argv[1], encoding="utf-8"))["source_commit"])
' "$run_root/plan/density_plan.json")
scanner_build=$("$scanner" --version-json)
python3 - "$scanner_build" "$plan_commit" <<'PY'
import json
import sys

build = json.loads(sys.argv[1])
assert build["program"] == "pssm_scan"
assert build["source_commit"] == sys.argv[2]
assert not build["source_dirty"]
PY

worker=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=jaspar_density_chr1 --nodes=1 --ntasks=1
    --cpus-per-task="$batch_workers" --mem="$memory" --time="$wall_time"
    --signal=B:USR1@300 --chdir="$source"
    --output="$run_root/logs/density-%j.out"
    --error="$run_root/logs/density-%j.err"
    "$source/scripts/run_motif_density_calibration_slurm.sh"
    --run-root "$run_root" --scanner "$scanner" --source "$source"
    --batch-workers "$batch_workers"
)

if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${worker[@]}"; printf '\n'
    worker_id=WORKER_JOB_ID
else
    worker_id=$("${worker[@]}")
fi

finalizer=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=jaspar_density_finalize --nodes=1 --ntasks=1
    --cpus-per-task=1 --mem=4G --time=00:30:00
    --dependency="afterok:$worker_id" --chdir="$source"
    --output="$run_root/logs/density-finalize-%j.out"
    --error="$run_root/logs/density-finalize-%j.err"
    "$source/scripts/run_motif_density_calibration_finalize.sh"
    --run-root "$run_root" --source "$source"
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${finalizer[@]}"; printf '\n'
    exit 0
fi
finalizer_id=$("${finalizer[@]}")
if [[ ! -e $submission_record ]]; then
    printf 'worker_job\tfinalizer_job\n' > "$submission_record"
fi
printf '%s\t%s\n' "$worker_id" "$finalizer_id" >> "$submission_record"
echo "I: Submitted density worker $worker_id and finalizer $finalizer_id." >&2
printf '%s\n' "$worker_id"

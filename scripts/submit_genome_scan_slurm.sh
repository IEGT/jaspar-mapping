#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_genome_scan_slurm.sh --run-root DIR --scanner FILE [OPTIONS]

Submit one requeue-safe Slurm array element per planned chromosome. Each element
copies its chromosome to node-local scratch once and runs all remaining motif
batches. The script submits jobs; use --dry-run to print the command only.

Options:
  --run-root DIR       Immutable prepared run
  --scanner FILE       Arrow-enabled pssm_scan_parquet executable
  --account NAME       Slurm account (default: cluster)
  --partition NAME     Slurm partition (default: requeue)
  --max-concurrent N   Maximum live chromosome workers (default: 20)
  --batch-workers N    Scanner processes per chromosome worker (default: 1)
  --memory SIZE        Memory per worker (default: 16G)
  --time D-HH:MM:SS    Time per worker (default: 0-08:00:00)
  --source DIR         Repository root (default: parent of this script)
  --duckdb FILE        DuckDB command (default: duckdb)
  --finalize           Submit a dependent fast finalizer
  --dry-run            Print the submission command without invoking sbatch
  -h, --help           Show this help and exit
EOF
}

run_root=""
scanner=""
account=cluster
partition=requeue
max_concurrent=20
batch_workers=1
memory=16G
wall_time=0-08:00:00
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
duckdb=duckdb
submit_finalize=0
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; run_root=$2; shift 2 ;;
        --scanner) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; scanner=$2; shift 2 ;;
        --account) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; account=$2; shift 2 ;;
        --partition) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; partition=$2; shift 2 ;;
        --max-concurrent) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; max_concurrent=$2; shift 2 ;;
        --batch-workers) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; batch_workers=$2; shift 2 ;;
        --memory) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; memory=$2; shift 2 ;;
        --time) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; wall_time=$2; shift 2 ;;
        --source) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; source=$2; shift 2 ;;
        --duckdb) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; duckdb=$2; shift 2 ;;
        --finalize) submit_finalize=1; shift ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scanner ]] || { usage >&2; exit 2; }
[[ -d $run_root ]] || { echo "E: Run root does not exist: $run_root" >&2; exit 1; }
[[ -d $source ]] || { echo "E: Source tree does not exist: $source" >&2; exit 1; }
[[ -x $scanner ]] || { echo "E: Scanner is not executable: $scanner" >&2; exit 1; }
run_root=$(cd "$run_root" && pwd)
source=$(cd "$source" && pwd)
scanner_directory=$(cd "$(dirname "$scanner")" && pwd)
scanner="$scanner_directory/$(basename "$scanner")"
[[ -x $source/scripts/run_genome_scan_slurm_chromosome.sh ]] || {
    echo "E: Chromosome Slurm wrapper is absent below source tree: $source" >&2
    exit 1
}
[[ -f $source/scripts/manage_genome_scan.py ]] || {
    echo "E: Scan coordinator is absent below source tree: $source" >&2
    exit 1
}
[[ $max_concurrent =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --max-concurrent must be a positive integer." >&2
    exit 2
}
[[ $batch_workers =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --batch-workers must be a positive integer." >&2
    exit 2
}
plan="$run_root/plan/scan_plan.json"
[[ -f $plan ]] || { echo "E: Scan plan not found: $plan" >&2; exit 1; }
mkdir -p "$run_root/logs"

chromosome_count=$(python3 -c '
import json, sys
plan = json.load(open(sys.argv[1], encoding="utf-8"))
print(sum(bool(row.get("included_in_scan", True)) for row in plan["sequence_regions"]))
' "$plan")
[[ $chromosome_count =~ ^[1-9][0-9]*$ ]] || {
    echo "E: Invalid chromosome count from scan plan: $chromosome_count" >&2
    exit 1
}

submission=(
    sbatch --parsable
    --account="$account" --partition="$partition" --requeue
    --array="0-$((chromosome_count - 1))%${max_concurrent}"
    --nodes=1 --ntasks=1 --cpus-per-task="$batch_workers" --mem="$memory" --time="$wall_time"
    --chdir="$source"
    --output="$run_root/logs/chromosome-%A_%a.out"
    --error="$run_root/logs/chromosome-%A_%a.err"
    "$source/scripts/run_genome_scan_slurm_chromosome.sh"
    --run-root "$run_root" --scanner "$scanner" --duckdb "$duckdb"
    --manager "$source/scripts/manage_genome_scan.py"
    --batch-workers "$batch_workers"
)

if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${submission[@]}"
    printf '\n'
    exit 0
fi

job_id=$("${submission[@]}")
echo "I: Submitted chromosome array $job_id with at most $max_concurrent workers." >&2
if [[ $submit_finalize -eq 1 ]]; then
    finalizer_submission=(
        sbatch --parsable
        --dependency="afterok:$job_id" --account="$account" --partition="$partition"
        --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G --time=0-02:00:00
        --chdir="$source"
        --output="$run_root/logs/finalize-%j.out"
        --error="$run_root/logs/finalize-%j.err"
        "$source/scripts/run_genome_scan_slurm_finalize.sh"
        --run-root "$run_root" --duckdb "$duckdb"
        --manager "$source/scripts/manage_genome_scan.py"
    )
    finalizer_id=$("${finalizer_submission[@]}")
    echo "I: Submitted dependent fast finalizer $finalizer_id." >&2
fi
printf '%s\n' "$job_id"

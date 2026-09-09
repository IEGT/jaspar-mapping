#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_regulatory_annotation_slurm.sh --run-root DIR
       --evidence-package DIR --annotation-catalog DIR --duckdb FILE [options]

Import the audited Ensembl source once, annotate all 22 autosomes with one
restart-safe task per chromosome, and finalize a zero-complete membership
catalog. This does not scan DNA or submit cofactor model fits.

Options:
  --source DIR           Clean repository checkout (default: script parent)
  --audit FILE           Coordinate audit JSON (default: tracked GRCh38 audit)
  --scratch-root DIR     Node-local staging (default: /scratch/sm718)
  --account NAME         Slurm account (default: cluster)
  --partition NAME       Slurm partition (default: requeue)
  --max-concurrent N     Concurrent chromosome tasks (default: 10)
  --memory SIZE          Memory per job (default: 8G)
  --time LIMIT           Time per job (default: 00:20:00)
  --dry-run              Pin inputs and print submission commands, but do not submit
  -h, --help             Show this help

Durable paths must be below /data/sm718. Requeue restages scratch and reuses
validated completed outputs. Source data and completed packages are never deleted.
SIGUSR1 requests a phase and elapsed-time report from the Python worker.
EOF
}

source=$(cd "$(dirname "$0")/.." && pwd)
run_root=""
evidence=""
catalog=""
duckdb=""
audit=""
scratch=/scratch/sm718
account=cluster
partition=requeue
concurrent=10
memory=8G
walltime=00:20:00
dry=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --evidence-package) evidence=${2:?}; shift 2 ;;
        --annotation-catalog) catalog=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --audit) audit=${2:?}; shift 2 ;;
        --scratch-root) scratch=${2:?}; shift 2 ;;
        --account) account=${2:?}; shift 2 ;;
        --partition) partition=${2:?}; shift 2 ;;
        --max-concurrent) concurrent=${2:?}; shift 2 ;;
        --memory) memory=${2:?}; shift 2 ;;
        --time) walltime=${2:?}; shift 2 ;;
        --dry-run) dry=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done
[[ -n $run_root && -n $evidence && -n $catalog && -n $duckdb ]] || { usage >&2; exit 2; }
[[ $concurrent =~ ^[1-9][0-9]*$ ]] || { echo "E: concurrency must be positive" >&2; exit 2; }
source=$(cd "$source" && pwd -P)
audit=${audit:-$source/resources/ensembl/regulatory_GRCh38_2025-05_coordinate_audit.json}
for path in "$source" "$run_root" "$evidence" "$catalog" "$duckdb" "$audit"; do
    case "$path" in /data/sm718/*) ;;
        *) echo "E: durable path must be below /data/sm718: $path" >&2; exit 2 ;;
    esac
done
case "$scratch" in /scratch/*) ;;
    *) echo "E: staging must use node-local /scratch" >&2; exit 2 ;;
esac
manager="$source/scripts/manage_regulatory_annotation.py"
python3 "$manager" prepare --run-root "$run_root" --source "$source" \
    --evidence-package "$evidence" --annotation-catalog "$catalog" \
    --duckdb "$duckdb" --audit "$audit"
printf -v setup 'exec python3 %q setup --run-root %q' "$manager" "$run_root"
printf -v worker 'exec python3 %q run-task --run-root %q --scratch-root %q' "$manager" "$run_root" "$scratch"
printf -v finalize 'exec python3 %q finalize --run-root %q' "$manager" "$run_root"
common=(--parsable --account="$account" --partition="$partition" --requeue
        --cpus-per-task=2 --mem="$memory" --time="$walltime" --chdir="$source")
if [[ $dry == 1 ]]; then
    printf 'sbatch'
    printf ' %q' "${common[@]}" --job-name=regulatory_setup "--wrap=$setup"
    printf '\n'
    printf 'sbatch'
    printf ' %q' "${common[@]}" '--dependency=afterok:SETUP_JOB' "--array=0-21%$concurrent" "--wrap=$worker"
    printf '\n'
    printf 'sbatch'
    printf ' %q' "${common[@]}" '--dependency=afterok:ARRAY_JOB' "--wrap=$finalize"
    printf '\n'
    exit 0
fi
setup_job=$(sbatch "${common[@]}" --job-name=regulatory_setup \
    --output="$run_root/logs/setup-%j.out" --wrap="$setup")
setup_job=${setup_job%%;*}
printf 'setup\t%s\n' "$setup_job" >> "$run_root/submissions.tsv"
array_job=$(sbatch "${common[@]}" --job-name=regulatory_chrom \
    --dependency="afterok:$setup_job" --array="0-21%$concurrent" \
    --output="$run_root/logs/chrom-%A_%a.out" --wrap="$worker")
array_job=${array_job%%;*}
printf 'chromosomes\t%s\n' "$array_job" >> "$run_root/submissions.tsv"
final_job=$(sbatch "${common[@]}" --job-name=regulatory_finalize \
    --dependency="afterok:$array_job" --output="$run_root/logs/final-%j.out" --wrap="$finalize")
final_job=${final_job%%;*}
printf 'finalize\t%s\n' "$final_job" >> "$run_root/submissions.tsv"
printf 'I: Submitted setup %s, chromosome array %s, finalizer %s\n' "$setup_job" "$array_job" "$final_job"

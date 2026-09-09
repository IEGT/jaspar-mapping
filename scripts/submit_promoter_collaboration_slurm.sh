#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_promoter_collaboration_slurm.sh --run-root DIR --analysis-run DIR
       --regulatory-features DIR --duckdb FILE --afterok JOB_ID [options]

Queue a compact collaborator export after a successful promoter finalizer.
The first job chooses the panel from completed common-score-zero results.
One array task per selected motif stages each chromosome's exact low-floor
Parquet inputs on /scratch; completed chromosome checkpoints survive requeue.
The finalizer publishes only a validated package below the byte limit.

Options:
  --run-root DIR           New dedicated directory below /data/sm718
  --analysis-run DIR       Completed/pending TP73 distance enrichment run
  --regulatory-features DIR Coordinate-audited Ensembl feature package
  --duckdb FILE            Absolute DuckDB CLI available on compute nodes
  --afterok JOB_ID         Existing promoter-result finalizer dependency
  --panel-size N           Detailed motif ceiling (default 100, maximum 200)
  --max-bytes N            Complete package ceiling (default 2000000000)
  --max-concurrent N       Concurrent motif tasks (default 4)
  --scratch-root DIR       Local scratch parent (default /scratch/sm718)
  --dry-run                Print submission commands without submitting
  -h, --help               Show this help

Fixed resources: cluster account, requeue partition, 2 CPUs; 32G/2h per motif,
8G/1h preparation, 8G/1h finalization. No input data are removed. Submission
IDs are journalled immediately; an existing journal blocks duplicate chains.
For failed jobs, inspect the journal and resubmit the affected worker with
the same arguments; workers reuse only checksum-validated checkpoints.
EOF
}

run_root= analysis_run= regulatory_features= duckdb= dependency= worker=
panel_size=100 max_bytes=2000000000 max_concurrent=4 scratch_root=/scratch/sm718 dry_run=0
while (($#)); do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2;;
        --analysis-run) analysis_run=${2:?}; shift 2;;
        --regulatory-features) regulatory_features=${2:?}; shift 2;;
        --duckdb) duckdb=${2:?}; shift 2;;
        --afterok) dependency=${2:?}; shift 2;;
        --panel-size) panel_size=${2:?}; shift 2;;
        --max-bytes) max_bytes=${2:?}; shift 2;;
        --max-concurrent) max_concurrent=${2:?}; shift 2;;
        --scratch-root) scratch_root=${2:?}; shift 2;;
        --worker) worker=${2:?}; shift 2;;
        --dry-run) dry_run=1; shift;;
        -h|--help) usage; exit 0;;
        *) echo "E: Unknown option $1" >&2; usage >&2; exit 2;;
    esac
done
[[ -n $run_root && -n $analysis_run && -n $regulatory_features && -x $duckdb ]] || { usage >&2; exit 2; }
for value in "$panel_size" "$max_bytes" "$max_concurrent"; do
    [[ $value =~ ^[1-9][0-9]*$ ]] || { echo 'E: Numeric options must be positive integers.' >&2; exit 2; }
done
((panel_size<=200)) || { echo 'E: Panel limit is 200.' >&2; exit 2; }
# Resolve before checking so ../ cannot move durable output outside the user's allocation.
run_root=$(realpath -m "$run_root")
case "$run_root" in /data/sm718/*) ;; *) echo 'E: Durable export must be below /data/sm718.' >&2; exit 2;; esac
case "$scratch_root" in /scratch/*) ;; *) echo 'E: Scratch must be node-local /scratch.' >&2; exit 2;; esac
source=$(cd "$(dirname "$0")/.." && pwd -P)
script="$source/scripts/export_promoter_collaboration.py"
self="$source/scripts/submit_promoter_collaboration_slurm.sh"

if [[ -n $worker ]]; then
    [[ -n ${SLURM_JOB_ID:-} ]] || { echo 'E: Workers require a Slurm allocation.' >&2; exit 2; }
    mkdir -p "$scratch_root"
    scratch=$(mktemp -d "$scratch_root/promoter-collaboration.${SLURM_JOB_ID}.XXXXXX")
    trap 'rm -rf -- "$scratch"' EXIT
    export TMPDIR="$scratch" TEMP="$scratch" TMP="$scratch" PYTHONDONTWRITEBYTECODE=1
    memory=6GB
    [[ $worker != motif ]] || memory=24GB
    common=(--duckdb "$duckdb" --memory-limit "$memory" --scratch-root "$scratch")
    case "$worker" in
        prepare) python3 "$script" "${common[@]}" prepare --run-root "$run_root" \
            --analysis-run "$analysis_run" --regulatory-features "$regulatory_features" \
            --panel-size "$panel_size" --max-bytes "$max_bytes";;
        motif) python3 "$script" "${common[@]}" motif --run-root "$run_root" --task-index "${SLURM_ARRAY_TASK_ID:?}";;
        finalize) python3 "$script" "${common[@]}" finalize --run-root "$run_root";;
        *) echo 'E: Unknown worker phase.' >&2; exit 2;;
    esac
    exit
fi

[[ $dependency =~ ^[1-9][0-9]*$ ]] || { echo 'E: --afterok must name an existing finalizer job.' >&2; exit 2; }
git -C "$source" diff --quiet
git -C "$source" diff --cached --quiet
commit=$(git -C "$source" rev-parse HEAD)
if [[ -e $run_root/submissions.tsv ]]; then
    echo 'E: Submission journal exists; inspect it rather than duplicating the chain.' >&2
    exit 1
fi
base=(sbatch --parsable --account=cluster --partition=requeue --requeue --cpus-per-task=2
      --output="$run_root/logs/%x-%A_%a.out" --error="$run_root/logs/%x-%A_%a.err")
worker_options=(--run-root "$run_root" --analysis-run "$analysis_run" --regulatory-features "$regulatory_features"
                --duckdb "$duckdb" --panel-size "$panel_size" --max-bytes "$max_bytes" --scratch-root "$scratch_root")
if (( ! dry_run )); then
    mkdir -p "$run_root/logs"
    # Noclobber also protects simultaneous submissions.
    (set -o noclobber; printf 'phase\tjob_id\tsource_commit\n' > "$run_root/submissions.tsv")
fi
job=$dependency
for phase in prepare motif finalize; do
    extra=(--mem=8G --time=01:00:00)
    if [[ $phase == motif ]]; then
        extra=(--mem=32G --time=02:00:00 --array="0-$((panel_size-1))%${max_concurrent}")
    fi
    command=("${base[@]}" "${extra[@]}" --job-name="glen_${phase}" --dependency="afterok:$job"
             "$self" --worker "$phase" "${worker_options[@]}")
    if ((dry_run)); then
        printf '%q ' "${command[@]}"; printf '\n'
        job=999999
    else
        job=$("${command[@]}")
        job=${job%%;*}
        [[ $job =~ ^[1-9][0-9]*$ ]] || { echo 'E: Ambiguous Slurm response; inspect queue before retry.' >&2; exit 1; }
        printf '%s\t%s\t%s\n' "$phase" "$job" "$commit" >> "$run_root/submissions.tsv"
        printf '%s: %s\n' "$phase" "$job"
    fi
done

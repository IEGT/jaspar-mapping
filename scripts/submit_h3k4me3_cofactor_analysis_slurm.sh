#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_h3k4me3_cofactor_analysis_slurm.sh --run-root DIR
       --h3-package DIR --evidence-package DIR --context-maxima-package DIR
       --annotation-catalog DIR --runtime-prefix DIR [options]

Submit restart-safe all-JASPAR H3K4me3 cofactor inference on autosomes 1-22.
A preflight job proves that H3K4me3 change, TP73 evidence, and schema-9
annotation use the same physical anchors. Batches stage those fixed inputs once
and publish each completed motif atomically. Finalization applies BH correction
over all planned motifs, rather than within one array task.

Options:
  --source DIR                 Repository root (default: script parent)
  --scratch-root DIR           Node-local scratch parent (default: /scratch/$USER)
  --run-id ID                  Immutable run identifier
  --motifs-per-batch N         Sequential motif checkpoints per job (default: 8)
  --account NAME               Slurm account (default: cluster)
  --partition NAME             Slurm partition (default: requeue)
  --max-concurrent N           Concurrent batches (default: 20)
  --array-size N               Maximum tasks per array submission (default: 500)
  --cpus N                     CPUs per batch (default: 4)
  --memory SIZE                Memory per batch (default: 32G)
  --time LIMIT                 Batch wall time (default: 04:00:00)
  --preflight-memory SIZE      Preflight memory (default: 32G)
  --preflight-time LIMIT       Preflight wall time (default: 02:00:00)
  --final-memory SIZE          Finalizer memory (default: 32G)
  --final-time LIMIT           Finalizer wall time (default: 02:00:00)
  --rscript FILE               Rscript executable (default: Rscript)
  --distance-bands CSV         Cofactor bands (default: all_150)
  --fixed-positive-threshold N Override motif thresholds (use 0 for score zero)
  --ranking-metadata FILE      JASPAR metadata for an afterok top-TF export
  --ranking-top N              TFs per criterion and band (default: 25)
  --adjust-gfp-baseline        Fit the GFP-baseline-adjusted sensitivity model
  --dry-run                    Prepare and print sbatch commands only
  -h, --help                   Show this help

Every durable path must be below /data/sm718. SIGUSR1 prints one progress line
without terminating a batch. Requeue repeats only unfinished motifs.
EOF
}

run_root=""
h3_package=""
evidence_package=""
context_maxima_package=""
annotation_catalog=""
runtime_prefix=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
scratch_root="/scratch/${USER:-sm718}"
run_id=jaspar2026_grch38_h3k4me3_cofactor_analysis_v1
motifs_per_batch=8
account=cluster
partition=requeue
max_concurrent=20
array_size=500
cpus=4
memory=32G
wall_time=04:00:00
preflight_memory=32G
preflight_time=02:00:00
final_memory=32G
final_time=02:00:00
rscript=Rscript
distance_bands=all_150
fixed_positive_threshold=""
ranking_metadata=""
ranking_top=25
adjust_gfp_baseline=0
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --h3-package) h3_package=${2:?}; shift 2 ;;
        --evidence-package) evidence_package=${2:?}; shift 2 ;;
        --context-maxima-package) context_maxima_package=${2:?}; shift 2 ;;
        --annotation-catalog) annotation_catalog=${2:?}; shift 2 ;;
        --runtime-prefix) runtime_prefix=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --scratch-root) scratch_root=${2:?}; shift 2 ;;
        --run-id) run_id=${2:?}; shift 2 ;;
        --motifs-per-batch) motifs_per_batch=${2:?}; shift 2 ;;
        --account) account=${2:?}; shift 2 ;;
        --partition) partition=${2:?}; shift 2 ;;
        --max-concurrent) max_concurrent=${2:?}; shift 2 ;;
        --array-size) array_size=${2:?}; shift 2 ;;
        --cpus) cpus=${2:?}; shift 2 ;;
        --memory) memory=${2:?}; shift 2 ;;
        --time) wall_time=${2:?}; shift 2 ;;
        --preflight-memory) preflight_memory=${2:?}; shift 2 ;;
        --preflight-time) preflight_time=${2:?}; shift 2 ;;
        --final-memory) final_memory=${2:?}; shift 2 ;;
        --final-time) final_time=${2:?}; shift 2 ;;
        --rscript) rscript=${2:?}; shift 2 ;;
        --distance-bands) distance_bands=${2:?}; shift 2 ;;
        --fixed-positive-threshold) fixed_positive_threshold=${2:?}; shift 2 ;;
        --ranking-metadata) ranking_metadata=${2:?}; shift 2 ;;
        --ranking-top) ranking_top=${2:?}; shift 2 ;;
        --adjust-gfp-baseline) adjust_gfp_baseline=1; shift ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $h3_package && -n $evidence_package &&
   -n $context_maxima_package && -n $annotation_catalog &&
   -n $runtime_prefix ]] || { usage >&2; exit 2; }
for value in "$motifs_per_batch" "$max_concurrent" "$array_size" "$cpus" \
             "$ranking_top"; do
    [[ $value =~ ^[1-9][0-9]*$ ]] || {
        echo "E: Batch, concurrency, array, and CPU values must be positive." >&2
        exit 2
    }
done
[[ $run_id =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] || {
    echo "E: --run-id contains unsafe characters." >&2
    exit 2
}
for path in "$run_root" "$h3_package" "$evidence_package" \
            "$context_maxima_package" "$annotation_catalog" "$runtime_prefix"; do
    case "$path" in
        /data/sm718/*) ;;
        *) echo "E: Production paths must be below /data/sm718: $path" >&2; exit 2 ;;
    esac
done

source=$(cd "$source" && pwd -P)
h3_package=$(cd "$h3_package" && pwd -P)
evidence_package=$(cd "$evidence_package" && pwd -P)
context_maxima_package=$(cd "$context_maxima_package" && pwd -P)
annotation_catalog=$(cd "$annotation_catalog" && pwd -P)
runtime_prefix=$(cd "$runtime_prefix" && pwd -P)
if [[ -n $ranking_metadata ]]; then
    ranking_metadata=$(cd "$(dirname "$ranking_metadata")" && pwd -P)/$(
        basename "$ranking_metadata"
    )
    [[ -f $ranking_metadata ]] || {
        echo "E: Ranking metadata is unavailable: $ranking_metadata" >&2
        exit 1
    }
    case "$ranking_metadata" in
        /data/sm718/*) ;;
        *) echo "E: Ranking metadata must be below /data/sm718." >&2; exit 2 ;;
    esac
fi
mkdir -p "$run_root/plan" "$run_root/logs"
run_root=$(cd "$run_root" && pwd -P)
duckdb="$runtime_prefix/duckdb/bin/duckdb"
[[ -x $duckdb ]] || { echo "E: DuckDB is unavailable: $duckdb" >&2; exit 1; }
[[ -x $rscript ]] || rscript=$(command -v "$rscript" || true)
[[ -x $rscript ]] || { echo "E: Rscript is unavailable." >&2; exit 1; }
if ! git -C "$source" diff --quiet --ignore-submodules -- ||
   ! git -C "$source" diff --cached --quiet --ignore-submodules --; then
    echo "E: Production submission requires a tracked-clean source tree." >&2
    exit 1
fi

prepare=(
    python3 "$source/scripts/manage_h3k4me3_cofactor_analysis.py" prepare
    --run-root "$run_root" --h3-package "$h3_package"
    --evidence-package "$evidence_package"
    --context-maxima-package "$context_maxima_package"
    --annotation-catalog "$annotation_catalog"
    --runtime-prefix "$runtime_prefix" --source "$source"
    --scratch-root "$scratch_root" --run-id "$run_id"
    --motifs-per-batch "$motifs_per_batch"
    --distance-bands "$distance_bands"
)
[[ -n $fixed_positive_threshold ]] &&
    prepare+=(--fixed-positive-threshold "$fixed_positive_threshold")
(( adjust_gfp_baseline == 1 )) && prepare+=(--adjust-gfp-baseline)
batch_count=$("${prepare[@]}")
[[ $batch_count =~ ^[1-9][0-9]*$ ]] || {
    echo "E: Invalid prepared batch count: $batch_count" >&2
    exit 1
}
submission_record="$run_root/plan/slurm_submission.tsv"
if [[ -e $submission_record && $dry_run -eq 0 ]]; then
    echo "E: Run already has a submission record: $submission_record" >&2
    exit 1
fi

preflight_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=h3_cofactor_preflight --nodes=1 --ntasks=1 --cpus-per-task=4
    --mem="$preflight_memory" --time="$preflight_time" --chdir="$source"
    --output="$run_root/logs/preflight-%j.out"
    --error="$run_root/logs/preflight-%j.err"
    "$source/scripts/manage_h3k4me3_cofactor_analysis.py" preflight
    --run-root "$run_root"
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${preflight_submission[@]}"; printf '\n'
    preflight_job=PREFLIGHT_JOB_ID
else
    preflight_job=$("${preflight_submission[@]}")
fi

array_jobs=()
batch_offset=0
dependency="afterok:${preflight_job}"
array_chunk=0
while (( batch_offset < batch_count )); do
    count=$((batch_count - batch_offset))
    (( count > array_size )) && count=$array_size
    submission=(
        sbatch --parsable --account="$account" --partition="$partition" --requeue
        --dependency="$dependency" --job-name="h3_cofactor_${array_chunk}"
        --nodes=1 --ntasks=1 --cpus-per-task="$cpus" --mem="$memory"
        --time="$wall_time" --array="0-$((count - 1))%${max_concurrent}"
        --signal=B:USR1@120 --open-mode=append --chdir="$source"
        --output="$run_root/logs/cofactor-%A_%a.out"
        --error="$run_root/logs/cofactor-%A_%a.err"
        "$source/scripts/manage_h3k4me3_cofactor_analysis.py" run-batch
        --run-root "$run_root" --batch-offset "$batch_offset"
        --rscript "$rscript"
    )
    if [[ $dry_run -eq 1 ]]; then
        printf '%q ' "${submission[@]}"; printf '\n'
        array_job="ARRAY_JOB_ID_${array_chunk}"
    else
        array_job=$("${submission[@]}")
    fi
    array_jobs+=("$array_job")
    dependency="afterok:${array_job}"
    batch_offset=$((batch_offset + count))
    array_chunk=$((array_chunk + 1))
done

final_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --dependency="$dependency" --job-name=h3_cofactor_finalize
    --nodes=1 --ntasks=1 --cpus-per-task=4 --mem="$final_memory"
    --time="$final_time" --chdir="$source"
    --output="$run_root/logs/finalize-%j.out"
    --error="$run_root/logs/finalize-%j.err"
    "$source/scripts/run_h3k4me3_cofactor_analysis_finalize.sh"
    --run-root "$run_root" --source "$source" --duckdb "$duckdb"
    --threads 4 --memory-limit 28GB
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${final_submission[@]}"; printf '\n'
    if [[ -n $ranking_metadata ]]; then
        printf '%q ' sbatch --parsable --dependency=afterok:FINALIZER_JOB_ID \
            --account="$account" --partition="$partition" --requeue \
            --job-name=h3_cofactor_rank --nodes=1 --ntasks=1 \
            --cpus-per-task=1 --mem=8G --time=01:00:00 --chdir="$source" \
            "$source/scripts/summarize_h3k4me3_isoform_rankings.py" \
            --analysis-database \
            "$run_root/final/h3k4me3_cofactor_analysis/h3k4me3_cofactor_analysis.duckdb" \
            --jaspar-metadata "$ranking_metadata" \
            --output-dir "$run_root/final/human_score0_isoform_rankings" \
            --duckdb "$duckdb" --top "$ranking_top"
        printf '\n'
    fi
    exit 0
fi
final_job=$("${final_submission[@]}")

ranking_job=""
if [[ -n $ranking_metadata ]]; then
    ranking_job=$(sbatch --parsable --account="$account" \
        --partition="$partition" --requeue --dependency="afterok:$final_job" \
        --job-name=h3_cofactor_rank --nodes=1 --ntasks=1 --cpus-per-task=1 \
        --mem=8G --time=01:00:00 --chdir="$source" \
        --output="$run_root/logs/ranking-%j.out" \
        --error="$run_root/logs/ranking-%j.err" \
        "$source/scripts/summarize_h3k4me3_isoform_rankings.py" \
        --analysis-database \
        "$run_root/final/h3k4me3_cofactor_analysis/h3k4me3_cofactor_analysis.duckdb" \
        --jaspar-metadata "$ranking_metadata" \
        --output-dir "$run_root/final/human_score0_isoform_rankings" \
        --duckdb "$duckdb" --top "$ranking_top")
fi

source_commit=$(git -C "$source" rev-parse HEAD)
array_jobs_text=$(IFS=,; printf '%s' "${array_jobs[*]}")
printf 'preflight_job_id\tarray_job_ids\tfinalizer_job_id\tranking_job_id\tbatch_count\tmotifs_per_batch\tsource_commit\trun_id\n' \
    > "$submission_record"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$preflight_job" \
    "$array_jobs_text" "$final_job" "$ranking_job" "$batch_count" "$motifs_per_batch" \
    "$source_commit" "$run_id" >> "$submission_record"
echo "I: Submitted H3K4me3 preflight $preflight_job, arrays $array_jobs_text, finalizer $final_job, ranking ${ranking_job:-not-requested}" >&2

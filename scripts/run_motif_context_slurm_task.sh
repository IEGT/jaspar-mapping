#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_motif_context_slurm_task.sh --run-root DIR --scan-package DIR
       --task-file FILE [OPTIONS]

Run one chromosome/motif-batch selected by SLURM_ARRAY_TASK_ID plus the optional
JASPAR_CONTEXT_TASK_OFFSET. Exact inventory payloads are staged on /data, copied
once to node-local scratch, and promoted atomically after validation.

Options:
  --run-root DIR          Dedicated context-run tree on durable storage
  --scan-package DIR      Finalized sparse genome-scan package
  --gtf FILE              GTF/GTF.gz (required for selected/summary tiers)
  --task-file FILE        TSV written by submit_motif_context_slurm.sh
  --source DIR            Repository root (default: parent of this script)
  --duckdb FILE           DuckDB executable (default: duckdb)
  --threads N             DuckDB threads (default: 4)
  --memory-limit SIZE     DuckDB memory ceiling (default: 24GB)
  --max-temp-size SIZE    DuckDB scratch ceiling (default: 100GB)
  -h, --help              Show this help and exit
EOF
}

run_root=""
scan_package=""
gtf=""
task_file=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
duckdb=duckdb
threads=4
memory_limit=24GB
max_temp_size=100GB
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; run_root=$2; shift 2 ;;
        --scan-package) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; scan_package=$2; shift 2 ;;
        --gtf) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; gtf=$2; shift 2 ;;
        --task-file) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; task_file=$2; shift 2 ;;
        --source) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; source=$2; shift 2 ;;
        --duckdb) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; duckdb=$2; shift 2 ;;
        --threads) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; threads=$2; shift 2 ;;
        --memory-limit) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; memory_limit=$2; shift 2 ;;
        --max-temp-size) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; max_temp_size=$2; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scan_package && -n $task_file ]] || { usage >&2; exit 2; }
[[ ${SLURM_ARRAY_TASK_ID:-} =~ ^[0-9]+$ ]] || {
    echo "E: SLURM_ARRAY_TASK_ID is required." >&2
    exit 2
}
task_offset=${JASPAR_CONTEXT_TASK_OFFSET:-0}
[[ $task_offset =~ ^[0-9]+$ ]] || {
    echo "E: JASPAR_CONTEXT_TASK_OFFSET must be a nonnegative integer." >&2
    exit 2
}
[[ $threads =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --threads must be a positive integer." >&2
    exit 2
}
[[ -d $run_root && -d $scan_package && -f $task_file ]] || {
    echo "E: Run root, scan package, or task file is missing." >&2
    exit 1
}
[[ -x $source/scripts/build_motif_context.py ]] || {
    echo "E: Context builder is not executable below $source." >&2
    exit 1
}
[[ -x $source/scripts/stage_motif_context_inputs.py ]] || {
    echo "E: Context input stager is not executable below $source." >&2
    exit 1
}

global_task_index=$((task_offset + SLURM_ARRAY_TASK_ID))
task_row=$(awk -F '\t' -v task="$global_task_index" \
    'NR > 1 && $1 == task { print; found = 1; exit } END { if (!found) exit 1 }' \
    "$task_file") || {
    echo "E: No task row for global task index $global_task_index." >&2
    exit 1
}
IFS=$'\t' read -r task_index chromosome cofactor_motif_ids output_tier \
    builder_source_commit <<< "$task_row"
[[ $task_index == "$global_task_index" ]]
[[ $chromosome =~ ^[A-Za-z0-9._-]+$ ]]
[[ $output_tier == selected || $output_tier == summary || $output_tier == band ]]
[[ $builder_source_commit =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "E: Task $task_index has an unsafe or absent builder source commit." >&2
    exit 2
}
if [[ $builder_source_commit != unknown ]]; then
    current_source_commit=$(git -C "$source" rev-parse --verify HEAD 2>/dev/null) || {
        echo "E: Cannot resolve the worker source commit below $source." >&2
        exit 1
    }
    [[ $current_source_commit == "$builder_source_commit" ]] || {
        echo "E: Worker source is $current_source_commit; task plan pins $builder_source_commit." >&2
        exit 1
    }
    if ! git -C "$source" diff --quiet -- \
       || ! git -C "$source" diff --cached --quiet --; then
        echo "E: Worker source has tracked changes: $source" >&2
        exit 1
    fi
fi
IFS=',' read -ra cofactor_motifs <<< "$cofactor_motif_ids"
[[ ${#cofactor_motifs[@]} -gt 0 ]]
for motif in "${cofactor_motifs[@]}"; do
    [[ $motif =~ ^[A-Za-z0-9._-]+$ && $motif != MA0861.2 ]] || {
        echo "E: Unsafe or anchor motif in task $task_index: $motif" >&2
        exit 2
    }
done
if [[ $output_tier != band ]]; then
    [[ -n $gtf && -f $gtf ]] || {
        echo "E: Task tier $output_tier requires --gtf." >&2
        exit 1
    }
fi

input="$run_root/inputs/task-$task_index-chrom-$chromosome"
output="$run_root/packages/chrom-$chromosome/task-$task_index"
mkdir -p "$(dirname "$input")" "$(dirname "$output")" "$run_root/staging/task-$task_index"

stage_arguments=(
    "$source/scripts/stage_motif_context_inputs.py"
    --package "$scan_package" --output "$input"
    --motif MA0861.2 --chrom "$chromosome" --duckdb "$duckdb"
)
for motif in "${cofactor_motifs[@]}"; do stage_arguments+=(--motif "$motif"); done
"${stage_arguments[@]}"

sql_motif_list=""
for motif in "${cofactor_motifs[@]}"; do
    if [[ -n $sql_motif_list ]]; then sql_motif_list+=","; fi
    sql_motif_list+="'$motif'"
done

validate_output() {
    local package=$1 valid unexpected
    [[ -f $package/context.duckdb && -f $package/input_manifest.json ]] || return 1
    cmp -s "$input/input_manifest.json" "$package/input_manifest.json" || return 1
    valid=$(
        cd "$package"
        "$duckdb" -light-mode -readonly -csv -noheader context.duckdb -c "
SELECT count(*)
FROM motif_context_run_config
WHERE schema_version = 5
  AND builder_source_commit = '$builder_source_commit'
  AND genome_id = 'homo_sapiens_grch38_ensembl113_primary'
  AND motif_set_id = 'jaspar2026_core_nonredundant'
  AND anchor_motif_id = 'MA0861.2'
  AND anchor_minimum_score = -1
  AND partner_minimum_score = 0
  AND anchor_selection_mode = 'local_peak'
  AND output_tier = '$output_tier'
  AND cofactor_pair_scope = 'at_least_one_member_is_a_tp73_context_locus'
  AND cofactor_motif_locus_scope = 'tp73_context_loci_plus_their_pair_partners'
  AND cofactor_locus_pair_feature_scope = 'tp73_context_loci_only';"
    ) || return 1
    [[ $valid == 1 ]] || return 1
    unexpected=$(
        cd "$package"
        "$duckdb" -light-mode -readonly -csv -noheader context.duckdb -c "
SELECT count(*) FROM anchor_motif_band_feature
WHERE neighbor_motif_id NOT IN ($sql_motif_list);"
    ) || return 1
    [[ $unexpected == 0 ]]
}

if [[ -e $output ]]; then
    validate_output "$output" || {
        echo "E: Existing context output is incomplete or incompatible: $output" >&2
        exit 1
    }
    echo "I: Reusing completed context task $task_index: $output" >&2
    exit 0
fi

scratch_base=${SLURM_TMPDIR:-/scratch/${USER:-sm718}}
scratch="$scratch_base/jaspar-context-${SLURM_JOB_ID:-manual}-${SLURM_ARRAY_TASK_ID}-global-$task_index-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
scratch_input="$scratch/input"
spill="$scratch/duckdb-spill"
mkdir -p "$scratch_input" "$spill"

input_kib=$(du -sk "$input" | awk '{print $1}')
available_kib=$(df -Pk "$scratch_base" | awk 'NR == 2 {print $4}')
case $max_temp_size in
    *KB) temp_number=${max_temp_size%KB}; temp_multiplier=1 ;;
    *MB) temp_number=${max_temp_size%MB}; temp_multiplier=1024 ;;
    *GB) temp_number=${max_temp_size%GB}; temp_multiplier=$((1024 * 1024)) ;;
    *TB) temp_number=${max_temp_size%TB}; temp_multiplier=$((1024 * 1024 * 1024)) ;;
    *)
        echo "E: --max-temp-size must use an integer KB, MB, GB, or TB suffix." >&2
        exit 2
        ;;
esac
[[ $temp_number =~ ^[0-9]+$ ]] || {
    echo "E: --max-temp-size must use an integer KB, MB, GB, or TB suffix." >&2
    exit 2
}
temp_kib=$((temp_number * temp_multiplier))
minimum_kib=$((input_kib + temp_kib + 10 * 1024 * 1024))
echo "I: Task $task_index: chromosome $chromosome, ${#cofactor_motifs[@]} cofactors, tier $output_tier" >&2
echo "I: Input size: $input_kib KiB; scratch available: $available_kib KiB" >&2
if (( available_kib < minimum_kib )); then
    echo "E: Scratch preflight failed; need at least $minimum_kib KiB before copying input." >&2
    exit 1
fi

echo "I: Copying exact chromosome/motif payloads from /data to $scratch_input" >&2
cp -R "$input/." "$scratch_input/"

attempt="$run_root/staging/task-$task_index/job-${SLURM_JOB_ID:-manual}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
build_arguments=(
    "$source/scripts/build_motif_context.py"
    --motif-hits "$scratch_input/**/*.parquet"
    --motif-hit-source-label "$input/**/*.parquet"
    --output "$attempt"
    --anchor-motif MA0861.2
    --source-commit "$builder_source_commit"
    --motif-set-id jaspar2026_core_nonredundant
    --genome-id homo_sapiens_grch38_ensembl113_primary
    --anchor-minimum-score -1 --tandem-minimum-score 0
    --anchor-selection-mode local_peak --anchor-local-peak-flank 150
    --score-mode log2_relative_risk --pseudocount 1
    --background-model-id uniform_acgt_v1
    --pseudocount-scheme additive_per_base
    --chrom "$chromosome"
    --capture-flank 150 --context-flank 150 --tandem-flank 20
    --cofactor-pair-flank 150 --output-tier "$output_tier"
    --threads "$threads" --memory-limit "$memory_limit"
    --max-temp-size "$max_temp_size" --temp-directory "$spill"
    --duckdb "$duckdb"
)
if [[ $output_tier != band ]]; then build_arguments+=(--gtf "$gtf"); fi

echo "I: Durable output after validation: $output" >&2
echo "I: Node-local DuckDB spill: $spill" >&2
"${build_arguments[@]}"
cp "$input/input_manifest.json" "$attempt/input_manifest.json"

validate_output "$attempt" || {
    echo "E: Newly built context package failed validation: $attempt" >&2
    exit 1
}
if [[ -e $output ]]; then
    validate_output "$output" || {
        echo "E: Output appeared concurrently but is incompatible: $output" >&2
        exit 1
    }
    echo "I: Another retry completed task $task_index first; retaining validated output." >&2
    exit 0
fi
mv "$attempt" "$output"
echo "I: Completed context task $task_index: $output" >&2

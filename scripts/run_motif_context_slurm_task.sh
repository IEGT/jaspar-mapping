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
  --annotation-release ID Annotation release (default: ensembl_113)
  --promoter-definition ID Versioned promoter definition
  --promoter-upstream-bp N Upstream promoter extent (default: 2000)
  --promoter-downstream-bp N Downstream promoter extent (default: 500)
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
annotation_release=ensembl_113
promoter_definition_id=tss_upstream_2000_downstream_500_v1
promoter_upstream_bp=2000
promoter_downstream_bp=500
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
        --annotation-release) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; annotation_release=$2; shift 2 ;;
        --promoter-definition) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; promoter_definition_id=$2; shift 2 ;;
        --promoter-upstream-bp) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; promoter_upstream_bp=$2; shift 2 ;;
        --promoter-downstream-bp) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; promoter_downstream_bp=$2; shift 2 ;;
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
[[ $annotation_release =~ ^[A-Za-z0-9._-]+$ \
   && $promoter_definition_id =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "E: Annotation release and promoter definition must be safe identifiers." >&2
    exit 2
}
[[ $promoter_upstream_bp =~ ^[0-9]+$ && $promoter_downstream_bp =~ ^[0-9]+$ ]] || {
    echo "E: Promoter extents must be nonnegative integers." >&2
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
    builder_source_commit context_schema_version gtf_size_bytes gtf_sha256 \
    task_annotation_release task_promoter_definition_id \
    task_promoter_upstream_bp task_promoter_downstream_bp task_kind \
    <<< "$task_row"
task_kind=${task_kind:-cofactor_context}
[[ $task_index == "$global_task_index" ]]
[[ $chromosome =~ ^[A-Za-z0-9._-]+$ ]]
[[ $output_tier == selected || $output_tier == summary || $output_tier == band ]]
[[ $context_schema_version == 7 ]] || {
    echo "E: Task $task_index requests unsupported context schema $context_schema_version." >&2
    exit 2
}
[[ $task_annotation_release == "$annotation_release" \
   && $task_promoter_definition_id == "$promoter_definition_id" \
   && $task_promoter_upstream_bp == "$promoter_upstream_bp" \
   && $task_promoter_downstream_bp == "$promoter_downstream_bp" ]] || {
    echo "E: Worker annotation/promoter arguments differ from task $task_index." >&2
    exit 2
}
[[ $builder_source_commit =~ ^[A-Za-z0-9._-]+$ ]] || {
    echo "E: Task $task_index has an unsafe or absent builder source commit." >&2
    exit 2
}
if [[ $builder_source_commit != unknown ]]; then
    [[ -f $source/source_commit.txt ]] || {
        echo "E: Worker source snapshot lacks source_commit.txt: $source" >&2
        exit 1
    }
    snapshot_source_commit=$(tr -d '\r\n' < "$source/source_commit.txt")
    [[ $snapshot_source_commit == "$builder_source_commit" ]] || {
        echo "E: Worker snapshot is $snapshot_source_commit; task plan pins $builder_source_commit." >&2
        exit 1
    }
fi
cofactor_motifs=(__jaspar_context_array_sentinel__)
case $task_kind in
    anchor_annotation)
        [[ $cofactor_motif_ids == none && $output_tier == summary ]] || {
            echo "E: Anchor-annotation task $task_index must use no cofactors and summary tier." >&2
            exit 2
        }
        ;;
    cofactor_context)
        IFS=',' read -ra parsed_cofactor_motifs <<< "$cofactor_motif_ids"
        [[ ${#parsed_cofactor_motifs[@]} -gt 0 ]]
        for motif in "${parsed_cofactor_motifs[@]}"; do
            [[ $motif =~ ^[A-Za-z0-9._-]+$ && $motif != MA0861.2 ]] || {
                echo "E: Unsafe or anchor motif in task $task_index: $motif" >&2
                exit 2
            }
            cofactor_motifs+=("$motif")
        done
        ;;
    *)
        echo "E: Unsupported task kind at task $task_index: $task_kind" >&2
        exit 2
        ;;
esac
cofactor_count=$((${#cofactor_motifs[@]} - 1))
if [[ $output_tier != band ]]; then
    [[ -n $gtf && -f $gtf ]] || {
        echo "E: Task tier $output_tier requires --gtf." >&2
        exit 1
    }
    [[ $gtf_size_bytes =~ ^[0-9]+$ && $gtf_sha256 =~ ^[0-9a-f]{64}$ ]] || {
        echo "E: Task $task_index lacks valid GTF size/SHA-256 provenance." >&2
        exit 2
    }
else
    [[ $gtf_size_bytes == 0 && $gtf_sha256 == none ]] || {
        echo "E: Band task $task_index unexpectedly carries GTF provenance." >&2
        exit 2
    }
fi

sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        echo "E: sha256sum or shasum is required." >&2
        return 1
    fi
}
if [[ $output_tier != band ]]; then
    actual_gtf_size=$(wc -c < "$gtf" | tr -d '[:space:]')
    actual_gtf_sha256=$(sha256_file "$gtf")
    [[ $actual_gtf_size == "$gtf_size_bytes" \
       && $actual_gtf_sha256 == "$gtf_sha256" ]] || {
        echo "E: GTF content no longer matches task $task_index provenance." >&2
        exit 1
    }
fi

input="$run_root/inputs/task-$task_index-chrom-$chromosome"
output="$run_root/packages/chrom-$chromosome/task-$task_index"
promotion_lock="$run_root/locks/context-task-$task_index.lock"
mkdir -p "$(dirname "$input")" "$(dirname "$output")" \
    "$run_root/staging/task-$task_index" "$(dirname "$promotion_lock")"

stage_arguments=(
    "$source/scripts/stage_motif_context_inputs.py"
    --package "$scan_package" --output "$input"
    --motif MA0861.2 --chrom "$chromosome" --duckdb "$duckdb"
)
for motif in "${cofactor_motifs[@]:1}"; do stage_arguments+=(--motif "$motif"); done
"${stage_arguments[@]}"

sql_motif_list=""
for motif in "${cofactor_motifs[@]:1}"; do
    if [[ -n $sql_motif_list ]]; then sql_motif_list+=","; fi
    sql_motif_list+="'$motif'"
done
if [[ $task_kind == anchor_annotation ]]; then
    unexpected_feature_sql="SELECT count(*) FROM anchor_motif_band_feature WHERE neighbor_motif_id <> 'MA0861.2';"
else
    unexpected_feature_sql="SELECT count(*) FROM anchor_motif_band_feature WHERE neighbor_motif_id NOT IN ($sql_motif_list);"
fi
if [[ $output_tier == band ]]; then
    gtf_validation_sql="AND gtf_source IS NULL AND gtf_sha256 IS NULL AND gtf_size_bytes IS NULL"
else
    gtf_validation_sql="AND gtf_source IS NOT NULL AND gtf_sha256 = '$gtf_sha256' AND gtf_size_bytes = $gtf_size_bytes"
fi

validate_output() {
    local package=$1 valid unexpected
    [[ -f $package/context.duckdb && -f $package/input_manifest.json ]] || return 1
    cmp -s "$input/input_manifest.json" "$package/input_manifest.json" || return 1
    valid=$(
        cd "$package"
        "$duckdb" -light-mode -readonly -csv -noheader context.duckdb -c "
SELECT count(*)
FROM motif_context_run_config
WHERE schema_version = $context_schema_version
  AND builder_source_commit = '$builder_source_commit'
  AND input_uniqueness = 'validated_scan_inventory'
  AND genome_id = 'homo_sapiens_grch38_ensembl113_primary'
  AND motif_set_id = 'jaspar2026_core_nonredundant'
  AND anchor_motif_id = 'MA0861.2'
  AND anchor_minimum_score = -1
  AND partner_minimum_score = 0
  AND anchor_selection_mode = 'local_peak'
  AND output_tier = '$output_tier'
  AND cofactor_pair_scope = 'at_least_one_member_is_a_tp73_context_locus'
  AND cofactor_motif_locus_scope = 'tp73_context_loci_plus_their_pair_partners'
  AND cofactor_locus_pair_feature_scope = 'tp73_context_loci_only'
  AND annotation_release = '$annotation_release'
  AND promoter_definition_id = '$promoter_definition_id'
  AND promoter_upstream_bp = $promoter_upstream_bp
  AND promoter_downstream_bp = $promoter_downstream_bp
  $gtf_validation_sql;"
    ) || return 1
    [[ $valid == 1 ]] || return 1
    unexpected=$(
        cd "$package"
        "$duckdb" -light-mode -readonly -csv -noheader context.duckdb -c "
$unexpected_feature_sql"
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
attempt_root="$run_root/staging/task-$task_index/job-${SLURM_JOB_ID:-manual}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
attempt="$attempt_root/package"

cleanup_owned_paths() {
    case ${attempt_root:-} in
        "$run_root/staging/task-$task_index/"job-*) rm -rf -- "$attempt_root" ;;
    esac
    case ${scratch:-} in
        "$scratch_base/"jaspar-context-*) rm -rf -- "$scratch" ;;
    esac
}
trap cleanup_owned_paths EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

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
echo "I: Task $task_index: chromosome $chromosome, kind $task_kind, $cofactor_count cofactors, tier $output_tier" >&2
echo "I: Input size: $input_kib KiB; scratch available: $available_kib KiB" >&2
if (( available_kib < minimum_kib )); then
    echo "E: Scratch preflight failed; need at least $minimum_kib KiB before copying input." >&2
    exit 1
fi

echo "I: Copying exact chromosome/motif payloads from /data to $scratch_input" >&2
cp -R "$input/." "$scratch_input/"

mkdir -p "$attempt_root"
build_arguments=(
    "$source/scripts/build_motif_context.py"
    --motif-hits "$scratch_input/**/*.parquet"
    --motif-hit-source-label "$input/**/*.parquet"
    --input-uniqueness validated_scan_inventory
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
    --annotation-release "$annotation_release"
    --promoter-definition-id "$promoter_definition_id"
    --promoter-upstream-bp "$promoter_upstream_bp"
    --promoter-downstream-bp "$promoter_downstream_bp"
    --chrom "$chromosome"
    --capture-flank 150 --context-flank 150 --tandem-flank 20
    --cofactor-pair-flank 150 --output-tier "$output_tier"
    --threads "$threads" --memory-limit "$memory_limit"
    --max-temp-size "$max_temp_size" --temp-directory "$spill"
    --duckdb "$duckdb"
)
if [[ $output_tier != band ]]; then
    build_arguments+=(
        --gtf "$gtf" --gtf-size-bytes "$gtf_size_bytes"
        --gtf-sha256 "$gtf_sha256"
    )
fi

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
promote_no_replace() {
    python3 - "$1" "$2" "$3" <<'PY'
import fcntl
import os
import sys

source, destination, lock_path = sys.argv[1:]
with open(lock_path, "a", encoding="utf-8") as lock:
    fcntl.flock(lock, fcntl.LOCK_EX)
    if os.path.lexists(destination):
        raise SystemExit(3)
    try:
        os.rename(source, destination)
    except OSError as error:
        print(f"E: atomic context promotion failed: {error}", file=sys.stderr)
        raise SystemExit(1)
PY
}
if promote_no_replace "$attempt" "$output" "$promotion_lock"; then
    :
else
    promotion_status=$?
    if [[ $promotion_status -eq 3 && -e $output ]] \
       && validate_output "$output"; then
        echo "I: Another retry promoted task $task_index first; retaining validated output." >&2
        exit 0
    fi
    echo "E: Could not promote context task $task_index without replacing output." >&2
    exit 1
fi
echo "I: Completed context task $task_index: $output" >&2

#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_motif_context_slurm.sh --run-root DIR --scan-package DIR
       (--motif ID [--motif ID ...] | --motif-file FILE) [OPTIONS]

Create an immutable chromosome/motif-batch plan and submit requeue-enabled
Slurm arrays. Each task builds one independently reusable context package.

Options:
  --run-root DIR          New or matching context-run tree under /data
  --scan-package DIR      Finalized sparse genome-scan package
  --gtf FILE              GTF/GTF.gz (required for selected/summary tiers)
  --motif ID              Cofactor motif accession; repeat or comma-separate
  --motif-file FILE       Cofactor accessions, one per line; comments allowed
  --chrom NAME            Chromosome; repeat/comma-separate (default: 1)
  --chrom-file FILE       Chromosomes, one per line; comments allowed
  --motifs-per-task N     Cofactors built together per chromosome (default: 20)
  --array-chunk-size N    Tasks per chained array, at most 1000 (default: 1000)
  --output-tier TIER      selected, summary, or band (default: selected)
  --source DIR            Repository root (default: parent of this script)
  --duckdb FILE           DuckDB executable (default: duckdb)
  --account NAME          Slurm account (default: cluster)
  --partition NAME        Slurm partition (default: requeue)
  --max-concurrent N      Maximum live tasks across chained arrays (default: 8)
  --cpus N                CPUs and DuckDB threads per task (default: 4)
  --memory SIZE           Slurm memory per task (default: 32G)
  --memory-limit SIZE     DuckDB memory ceiling (default: 24GB)
  --max-temp-size SIZE    DuckDB scratch ceiling (default: 100GB)
  --time D-HH:MM:SS       Wall time per task (default: 0-08:00:00)
  --dry-run               Render all sbatch commands without submitting
  -h, --help              Show this help and exit
EOF
}

run_root=""
scan_package=""
gtf=""
# Bash 3 treats an empty array expansion as unbound under `set -u`. Sentinels
# keep the command usable on the macOS system Bash used by the local tests.
motif_values=(__jaspar_context_array_sentinel__)
motif_files=(__jaspar_context_array_sentinel__)
chrom_values=(__jaspar_context_array_sentinel__)
chrom_files=(__jaspar_context_array_sentinel__)
motifs_per_task=20
array_chunk_size=1000
output_tier=selected
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
duckdb=duckdb
account=cluster
partition=requeue
max_concurrent=8
cpus=4
memory=32G
memory_limit=24GB
max_temp_size=100GB
wall_time=0-08:00:00
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; run_root=$2; shift 2 ;;
        --scan-package) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; scan_package=$2; shift 2 ;;
        --gtf) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; gtf=$2; shift 2 ;;
        --motif) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; motif_values+=("$2"); shift 2 ;;
        --motif-file) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; motif_files+=("$2"); shift 2 ;;
        --chrom) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; chrom_values+=("$2"); shift 2 ;;
        --chrom-file) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; chrom_files+=("$2"); shift 2 ;;
        --motifs-per-task) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; motifs_per_task=$2; shift 2 ;;
        --array-chunk-size) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; array_chunk_size=$2; shift 2 ;;
        --output-tier) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; output_tier=$2; shift 2 ;;
        --source) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; source=$2; shift 2 ;;
        --duckdb) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; duckdb=$2; shift 2 ;;
        --account) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; account=$2; shift 2 ;;
        --partition) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; partition=$2; shift 2 ;;
        --max-concurrent) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; max_concurrent=$2; shift 2 ;;
        --cpus) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; cpus=$2; shift 2 ;;
        --memory) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; memory=$2; shift 2 ;;
        --memory-limit) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; memory_limit=$2; shift 2 ;;
        --max-temp-size) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; max_temp_size=$2; shift 2 ;;
        --time) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; wall_time=$2; shift 2 ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scan_package ]] || { usage >&2; exit 2; }
[[ $output_tier == selected || $output_tier == summary || $output_tier == band ]] || {
    echo "E: --output-tier must be selected, summary, or band." >&2
    exit 2
}
for value in "$max_concurrent" "$cpus" "$motifs_per_task" "$array_chunk_size"; do
    [[ $value =~ ^[1-9][0-9]*$ ]] || {
        echo "E: Concurrency, CPU, batching, and chunk values must be positive integers." >&2
        exit 2
    }
done
(( array_chunk_size <= 1000 )) || {
    echo "E: --array-chunk-size must not exceed Haumea's safe limit of 1000." >&2
    exit 2
}
[[ -d $scan_package && -f $scan_package/manifest.json && -d $source ]] || {
    echo "E: Scan package, finalized manifest, or source tree is missing." >&2
    exit 1
}
if [[ $output_tier != band ]]; then
    [[ -n $gtf && -f $gtf ]] || {
        echo "E: --gtf is required for selected and summary output tiers." >&2
        exit 1
    }
elif [[ -n $gtf && ! -f $gtf ]]; then
    echo "E: GTF does not exist: $gtf" >&2
    exit 1
fi
[[ -x $source/scripts/run_motif_context_slurm_task.sh ]] || {
    echo "E: Context Slurm worker is not executable below $source." >&2
    exit 1
}

append_file_values() {
    local destination=$1 file=$2 line item
    [[ -f $file ]] || { echo "E: List file does not exist: $file" >&2; return 1; }
    while IFS= read -r line || [[ -n $line ]]; do
        line=${line%%#*}
        for item in $line; do
            if [[ $destination == motif ]]; then
                motif_values+=("$item")
            else
                chrom_values+=("$item")
            fi
        done
    done < "$file"
}
for file in "${motif_files[@]:1}"; do append_file_values motif "$file"; done
for file in "${chrom_files[@]:1}"; do append_file_values chrom "$file"; done
motif_values=("${motif_values[@]:1}")
chrom_values=("${chrom_values[@]:1}")
[[ ${#motif_values[@]} -gt 0 ]] || {
    echo "E: At least one --motif or --motif-file entry is required." >&2
    exit 2
}

split_unique() {
    local value item existing duplicate
    local -a pieces
    collected=(__jaspar_context_array_sentinel__)
    for value in "$@"; do
        IFS=',' read -ra pieces <<< "$value"
        for item in "${pieces[@]}"; do
            [[ -n $item && $item =~ ^[A-Za-z0-9._-]+$ ]] || {
                echo "E: Invalid empty or unsafe motif/chromosome value: $item" >&2
                return 1
            }
            duplicate=0
            for existing in "${collected[@]}"; do
                if [[ $existing == "$item" ]]; then duplicate=1; break; fi
            done
            if [[ $duplicate -eq 0 ]]; then collected+=("$item"); fi
        done
    done
    collected=("${collected[@]:1}")
}

split_unique "${motif_values[@]}"
motifs=("${collected[@]}")
if [[ ${#chrom_values[@]} -eq 0 ]]; then chrom_values=(1); fi
split_unique "${chrom_values[@]}"
chromosomes=("${collected[@]}")
for motif in "${motifs[@]}"; do
    [[ $motif != MA0861.2 ]] || {
        echo "E: Cofactor lists must omit TP73 MA0861.2; it is included automatically." >&2
        exit 2
    }
done

mkdir -p "$run_root/plan" "$run_root/logs" "$run_root/inputs" \
    "$run_root/packages" "$run_root/staging"
run_root=$(cd "$run_root" && pwd)
scan_package=$(cd "$scan_package" && pwd)
source=$(cd "$source" && pwd)
source_commit=$(git -C "$source" rev-parse --verify HEAD 2>/dev/null) || {
    echo "E: Source tree has no resolvable Git commit: $source" >&2
    exit 1
}
if [[ $dry_run -eq 0 ]] && {
       ! git -C "$source" diff --quiet -- \
       || ! git -C "$source" diff --cached --quiet --;
}; then
    echo "E: Refusing production submission from a source tree with tracked changes." >&2
    exit 1
fi
source_snapshot="$run_root/source"
mkdir -p "$source_snapshot/scripts"
for relative in scripts/build_motif_context.py scripts/stage_motif_context_inputs.py; do
    snapshot_target="$source_snapshot/$relative"
    if [[ -e $snapshot_target ]]; then
        cmp -s "$source/$relative" "$snapshot_target" || {
            echo "E: Existing immutable source snapshot differs: $snapshot_target" >&2
            exit 1
        }
    else
        cp -p "$source/$relative" "$snapshot_target"
    fi
done
if [[ -e $source_snapshot/source_commit.txt ]]; then
    [[ $(tr -d '\r\n' < "$source_snapshot/source_commit.txt") == "$source_commit" ]] || {
        echo "E: Existing source snapshot pins a different commit." >&2
        exit 1
    }
else
    printf '%s\n' "$source_commit" > "$source_snapshot/source_commit.txt"
fi
if [[ -n $gtf ]]; then
    gtf_directory=$(cd "$(dirname "$gtf")" && pwd)
    gtf="$gtf_directory/$(basename "$gtf")"
fi

task_file="$run_root/plan/context_tasks.tsv"
candidate=$(mktemp "$run_root/plan/.context_tasks.XXXXXX")
trap 'rm -f "$candidate"' EXIT HUP INT TERM
printf 'task_index\tchrom\tcofactor_motif_ids\toutput_tier\tbuilder_source_commit\n' \
    > "$candidate"
task_index=0
motif_offset=0
while (( motif_offset < ${#motifs[@]} )); do
    batch=("${motifs[@]:motif_offset:motifs_per_task}")
    joined=""
    for motif in "${batch[@]}"; do
        if [[ -n $joined ]]; then joined+=","; fi
        joined+=$motif
    done
    for chromosome in "${chromosomes[@]}"; do
        printf '%d\t%s\t%s\t%s\t%s\n' \
            "$task_index" "$chromosome" "$joined" "$output_tier" \
            "$source_commit" >> "$candidate"
        task_index=$((task_index + 1))
    done
    motif_offset=$((motif_offset + motifs_per_task))
done
if [[ -f $task_file ]]; then
    cmp -s "$candidate" "$task_file" || {
        echo "E: Existing immutable task plan differs: $task_file" >&2
        exit 1
    }
else
    mv "$candidate" "$task_file"
fi

worker_arguments=(
    "$source/scripts/run_motif_context_slurm_task.sh"
    --run-root "$run_root" --scan-package "$scan_package"
    --task-file "$task_file" --source "$source_snapshot" --duckdb "$duckdb"
    --threads "$cpus" --memory-limit "$memory_limit"
    --max-temp-size "$max_temp_size"
)
if [[ -n $gtf ]]; then worker_arguments+=(--gtf "$gtf"); fi

job_ids=()
task_offset=0
chunk_number=0
previous_job_id=""
while (( task_offset < task_index )); do
    remaining=$((task_index - task_offset))
    chunk_tasks=$array_chunk_size
    if (( remaining < chunk_tasks )); then chunk_tasks=$remaining; fi
    submission=(
        sbatch --parsable --account="$account" --partition="$partition" --requeue
        --job-name=tp73_context_v5
        --array="0-$((chunk_tasks - 1))%${max_concurrent}"
        --export="ALL,JASPAR_CONTEXT_TASK_OFFSET=$task_offset"
        --nodes=1 --ntasks=1 --cpus-per-task="$cpus" --mem="$memory" --time="$wall_time"
        --chdir="$source"
        --output="$run_root/logs/context-%A_%a.out"
        --error="$run_root/logs/context-%A_%a.err"
    )
    if [[ -n $previous_job_id ]]; then
        submission+=(--dependency="afterany:$previous_job_id")
    fi
    submission+=("${worker_arguments[@]}")

    if [[ $dry_run -eq 1 ]]; then
        printf '%q ' "${submission[@]}"
        printf '\n'
        previous_job_id="DRY_RUN_CHUNK_$chunk_number"
    else
        job_id=$("${submission[@]}")
        job_ids+=("$job_id")
        previous_job_id=$job_id
    fi
    task_offset=$((task_offset + chunk_tasks))
    chunk_number=$((chunk_number + 1))
done

if [[ $dry_run -eq 1 ]]; then exit 0; fi
echo "I: Submitted $task_index context tasks in $chunk_number chained arrays; at most $max_concurrent run at once." >&2
printf '%s\n' "${job_ids[@]}"

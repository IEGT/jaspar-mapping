#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_motif_context_slurm.sh --run-root DIR --scan-package DIR
       --gtf FILE --motif ID [--motif ID ...] [OPTIONS]

Create an immutable motif/chromosome task plan and submit a requeue-enabled
Slurm array. Each task builds one independently reusable context package.

Options:
  --run-root DIR          New or matching context-run tree under /data
  --scan-package DIR      Finalized sparse genome-scan package
  --gtf FILE              Ensembl-compatible GTF or GTF.gz
  --motif ID              Cofactor motif accession; repeat or comma-separate
  --chrom NAME            Chromosome; repeat/comma-separate (default: 1)
  --output-tier TIER      selected or summary (default: selected)
  --source DIR            Repository root (default: parent of this script)
  --duckdb FILE           DuckDB executable (default: duckdb)
  --account NAME          Slurm account (default: cluster)
  --partition NAME        Slurm partition (default: requeue)
  --max-concurrent N      Maximum live tasks (default: 8)
  --cpus N                CPUs and DuckDB threads per task (default: 4)
  --memory SIZE           Slurm memory per task (default: 32G)
  --memory-limit SIZE     DuckDB memory ceiling (default: 24GB)
  --max-temp-size SIZE    DuckDB scratch ceiling (default: 100GB)
  --time D-HH:MM:SS       Wall time per task (default: 0-12:00:00)
  --dry-run               Print the sbatch command without submitting
  -h, --help              Show this help and exit
EOF
}

run_root=""
scan_package=""
gtf=""
motif_values=()
chrom_values=()
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
wall_time=0-12:00:00
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; run_root=$2; shift 2 ;;
        --scan-package) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; scan_package=$2; shift 2 ;;
        --gtf) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; gtf=$2; shift 2 ;;
        --motif) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; motif_values+=("$2"); shift 2 ;;
        --chrom) [[ $# -ge 2 ]] || { usage >&2; exit 2; }; chrom_values+=("$2"); shift 2 ;;
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

[[ -n $run_root && -n $scan_package && -n $gtf && ${#motif_values[@]} -gt 0 ]] || {
    usage >&2
    exit 2
}
[[ $output_tier == selected || $output_tier == summary ]]
for value in "$max_concurrent" "$cpus"; do
    [[ $value =~ ^[1-9][0-9]*$ ]] || {
        echo "E: Concurrency and CPU values must be positive integers." >&2
        exit 2
    }
done
[[ -d $scan_package && -f $scan_package/manifest.json && -f $gtf && -d $source ]] || {
    echo "E: Scan package, finalized manifest, GTF, or source tree is missing." >&2
    exit 1
}
[[ -x $source/scripts/run_motif_context_slurm_task.sh ]] || {
    echo "E: Context Slurm worker is not executable below $source." >&2
    exit 1
}

split_unique() {
    local value item existing duplicate
    local -a pieces
    # Bash 3 with `set -u` treats an empty array expansion as unbound. Keep a
    # private sentinel while collecting and strip it before returning.
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
                if [[ $existing == "$item" ]]; then
                    duplicate=1
                    break
                fi
            done
            if [[ $duplicate -eq 0 ]]; then
                collected+=("$item")
            fi
        done
    done
    collected=("${collected[@]:1}")
}

split_unique "${motif_values[@]}"
motifs=("${collected[@]}")
if [[ ${#chrom_values[@]} -eq 0 ]]; then
    chrom_values=(1)
fi
split_unique "${chrom_values[@]}"
chromosomes=("${collected[@]}")
for motif in "${motifs[@]}"; do
    [[ $motif != MA0861.2 ]] || {
        echo "E: --motif lists cofactors; TP73 MA0861.2 is included automatically." >&2
        exit 2
    }
done

mkdir -p "$run_root/plan" "$run_root/logs" "$run_root/inputs" "$run_root/packages"
run_root=$(cd "$run_root" && pwd)
scan_package=$(cd "$scan_package" && pwd)
source=$(cd "$source" && pwd)
gtf_directory=$(cd "$(dirname "$gtf")" && pwd)
gtf="$gtf_directory/$(basename "$gtf")"
task_file="$run_root/plan/context_tasks.tsv"
candidate=$(mktemp "$run_root/plan/.context_tasks.XXXXXX")
trap 'rm -f "$candidate"' EXIT HUP INT TERM
printf 'task_index\tchrom\tcofactor_motif_id\toutput_tier\n' > "$candidate"
task_index=0
for motif in "${motifs[@]}"; do
    for chromosome in "${chromosomes[@]}"; do
        printf '%d\t%s\t%s\t%s\n' \
            "$task_index" "$chromosome" "$motif" "$output_tier" >> "$candidate"
        task_index=$((task_index + 1))
    done
done
if [[ -f $task_file ]]; then
    cmp -s "$candidate" "$task_file" || {
        echo "E: Existing immutable task plan differs: $task_file" >&2
        exit 1
    }
else
    mv "$candidate" "$task_file"
fi

submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=tp73_context_v4
    --array="0-$((task_index - 1))%${max_concurrent}"
    --nodes=1 --ntasks=1 --cpus-per-task="$cpus" --mem="$memory" --time="$wall_time"
    --chdir="$source"
    --output="$run_root/logs/context-%A_%a.out"
    --error="$run_root/logs/context-%A_%a.err"
    "$source/scripts/run_motif_context_slurm_task.sh"
    --run-root "$run_root" --scan-package "$scan_package" --gtf "$gtf"
    --task-file "$task_file" --source "$source" --duckdb "$duckdb"
    --threads "$cpus" --memory-limit "$memory_limit"
    --max-temp-size "$max_temp_size"
)

if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${submission[@]}"
    printf '\n'
    exit 0
fi
job_id=$("${submission[@]}")
echo "I: Submitted context array $job_id with $task_index tasks and at most $max_concurrent live." >&2
printf '%s\n' "$job_id"

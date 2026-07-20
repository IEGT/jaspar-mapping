#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_motif_context_slurm_task.sh --run-root DIR --scan-package DIR
       --gtf FILE --task-file FILE [OPTIONS]

Run one motif/chromosome context task selected by SLURM_ARRAY_TASK_ID. Exact
scan payloads are staged as hard links on durable storage; DuckDB spill uses a
unique node-local scratch directory. Completed matching outputs are reused.

Options:
  --run-root DIR          Dedicated context-run tree on durable storage
  --scan-package DIR      Finalized sparse genome-scan package
  --gtf FILE              Ensembl-compatible GTF or GTF.gz
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

[[ -n $run_root && -n $scan_package && -n $gtf && -n $task_file ]] || {
    usage >&2
    exit 2
}
[[ ${SLURM_ARRAY_TASK_ID:-} =~ ^[0-9]+$ ]] || {
    echo "E: SLURM_ARRAY_TASK_ID is required." >&2
    exit 2
}
[[ $threads =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --threads must be a positive integer." >&2
    exit 2
}
[[ -d $run_root && -d $scan_package && -f $gtf && -f $task_file ]] || {
    echo "E: Run root, scan package, GTF, or task file is missing." >&2
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

task_row=$(awk -F '\t' -v task="$SLURM_ARRAY_TASK_ID" \
    'NR > 1 && $1 == task { print; found = 1; exit } END { if (!found) exit 1 }' \
    "$task_file") || {
    echo "E: No task row for array index $SLURM_ARRAY_TASK_ID." >&2
    exit 1
}
IFS=$'\t' read -r task_index chromosome cofactor_motif output_tier <<< "$task_row"
[[ $task_index == "$SLURM_ARRAY_TASK_ID" ]]
[[ $chromosome =~ ^[A-Za-z0-9._-]+$ && $cofactor_motif =~ ^[A-Za-z0-9._-]+$ ]]
[[ $output_tier == selected || $output_tier == summary ]]

# This wrapper must not look like a Hive partition. DuckDB scans partition
# labels from every path component; an outer motif_id=... would silently
# replace the genuine TP73/cofactor labels carried by the linked source paths.
input="$run_root/inputs/task-$task_index-cofactor-$cofactor_motif-chrom-$chromosome"
output="$run_root/packages/motif_id=$cofactor_motif/chrom=$chromosome"
mkdir -p "$(dirname "$input")" "$(dirname "$output")"

"$source/scripts/stage_motif_context_inputs.py" \
    --package "$scan_package" --output "$input" \
    --motif MA0861.2 --motif "$cofactor_motif" --chrom "$chromosome" \
    --duckdb "$duckdb"

if [[ -e $output ]]; then
    [[ -f $output/context.duckdb ]] || {
        echo "E: Existing context output is incomplete: $output" >&2
        exit 1
    }
    valid=$(
        cd "$output"
        "$duckdb" -readonly -csv -noheader context.duckdb -c "
SELECT count(*)
FROM motif_context_run_config
WHERE schema_version = 4
  AND genome_id = 'homo_sapiens_grch38_ensembl113_primary'
  AND motif_set_id = 'jaspar2026_core_nonredundant'
  AND anchor_motif_id = 'MA0861.2'
  AND anchor_minimum_score = -1
  AND partner_minimum_score = 0
  AND anchor_selection_mode = 'local_peak'
  AND output_tier = '$output_tier';"
    )
    [[ $valid == 1 ]] || {
        echo "E: Existing context output has incompatible provenance: $output" >&2
        exit 1
    }
    echo "I: Reusing completed context task $task_index: $output" >&2
    exit 0
fi

scratch_base=${SLURM_TMPDIR:-/scratch/${USER:-sm718}}
scratch="$scratch_base/jaspar-context-${SLURM_JOB_ID:-manual}-${SLURM_ARRAY_TASK_ID}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
spill="$scratch/duckdb-spill"
mkdir -p "$spill"

echo "I: Task $task_index: chromosome $chromosome, cofactor $cofactor_motif, tier $output_tier" >&2
echo "I: Durable output: $output" >&2
echo "I: Node-local spill: $spill" >&2
"$source/scripts/build_motif_context.py" \
    --motif-hits "$input/**/*.parquet" \
    --gtf "$gtf" --output "$output" \
    --anchor-motif MA0861.2 \
    --motif-set-id jaspar2026_core_nonredundant \
    --genome-id homo_sapiens_grch38_ensembl113_primary \
    --anchor-minimum-score -1 --tandem-minimum-score 0 \
    --anchor-selection-mode local_peak --anchor-local-peak-flank 150 \
    --score-mode log2_relative_risk --pseudocount 1 \
    --background-model-id uniform_acgt_v1 \
    --pseudocount-scheme additive_per_base \
    --chrom "$chromosome" \
    --capture-flank 150 --context-flank 150 --tandem-flank 20 \
    --cofactor-pair-flank 150 --output-tier "$output_tier" \
    --threads "$threads" --memory-limit "$memory_limit" \
    --max-temp-size "$max_temp_size" --temp-directory "$spill" \
    --duckdb "$duckdb"

echo "I: Completed context task $task_index: $output" >&2

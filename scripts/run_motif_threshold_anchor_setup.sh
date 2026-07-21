#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_motif_threshold_anchor_setup.sh --run-root DIR --scan-package DIR
       --cutandrun-dir DIR --output FILE [options]

Build the fixed chromosome-1 TP73 anchor/CUT&RUN evidence table used by every
motif-threshold task. Exact TP73 sparse files come from the immutable task plan.
Anti-p73 bedGraphs are used directly; six matched control BigWigs are exported
to chromosome-local bedGraph before strict-immersion labels are calculated.

Options:
  --run-root DIR          Prepared threshold-calibration run
  --scan-package DIR      Finalized sparse-scan package
  --cutandrun-dir DIR     CUT&RUN resource directory
  --output FILE           New TP73 anchor-evidence Parquet
  --source DIR            Repository root (default: script parent)
  --duckdb FILE           DuckDB CLI (default: duckdb)
  --bigwig-to-bedgraph FILE  UCSC converter (default: bigWigToBedGraph)
  --threads N             DuckDB threads (default: 4)
  --memory-limit SIZE     DuckDB ceiling (default: 12GB)
  -h, --help              Show this help
EOF
}

run_root=""
scan_package=""
cutandrun_dir=""
output=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
duckdb=duckdb
bigwig_to_bedgraph=bigWigToBedGraph
threads=4
memory_limit=12GB
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --scan-package) scan_package=${2:?}; shift 2 ;;
        --cutandrun-dir) cutandrun_dir=${2:?}; shift 2 ;;
        --output) output=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --bigwig-to-bedgraph) bigwig_to_bedgraph=${2:?}; shift 2 ;;
        --threads) threads=${2:?}; shift 2 ;;
        --memory-limit) memory_limit=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scan_package && -n $cutandrun_dir && -n $output ]] || {
    usage >&2
    exit 2
}
for path in "$run_root" "$scan_package" "$cutandrun_dir"; do
    [[ -d $path ]] || { echo "E: Directory not found: $path" >&2; exit 1; }
done
target_plan="$run_root/plan/target_anchor_files.tsv"
[[ -f $target_plan ]] || { echo "E: Target anchor plan is missing: $target_plan" >&2; exit 1; }
[[ -x $duckdb ]] || duckdb=$(command -v "$duckdb" || true)
[[ -x $bigwig_to_bedgraph ]] || bigwig_to_bedgraph=$(command -v "$bigwig_to_bedgraph" || true)
[[ -x $duckdb && -x $bigwig_to_bedgraph ]] || {
    echo "E: DuckDB or bigWigToBedGraph is unavailable." >&2
    exit 1
}

if [[ -e $output ]]; then
    valid=$("$duckdb" -csv -noheader :memory: -c "
SELECT count(*) = 310782
FROM read_parquet('$output');")
    [[ $valid == true ]] || { echo "E: Existing anchor evidence is incompatible." >&2; exit 1; }
    echo "I: Reusing completed TP73 anchor evidence: $output" >&2
    exit 0
fi

plus_relative=$(awk -F '\t' '$1 == "+" {print $2; found=1} END {if (!found) exit 1}' "$target_plan")
minus_relative=$(awk -F '\t' '$1 == "-" {print $2; found=1} END {if (!found) exit 1}' "$target_plan")
plus="$scan_package/$plus_relative"
minus="$scan_package/$minus_relative"
[[ -f $plus && -f $minus ]] || { echo "E: Planned TP73 sparse files are missing." >&2; exit 1; }

start_epoch=$(date +%s)
phase=exporting_controls
sample=none
progress_report() {
    printf 'I: progress signal=SIGUSR1 phase=%s sample=%s elapsed_s=%s output=%s\n' \
        "$phase" "$sample" "$(( $(date +%s) - start_epoch ))" \
        "$([[ -s $output ]] && echo ready || echo pending)" >&2
}
trap progress_report USR1

control_dir="$run_root/input/control_bedgraph"
mkdir -p "$control_dir" "$(dirname "$output")"
samples=(saos2_TA saos2_DN skmel29_1_TA skmel29_1_DN skmel29_2_TA skmel29_2_DN)
coverage_arguments=()
for sample in "${samples[@]}"; do
    anti="$cutandrun_dir/tp73_${sample}_R1.clipped.bedGraph"
    control_bigwig="$cutandrun_dir/neg_${sample}_R1.bigWig"
    control="$control_dir/neg_${sample}_chr1.bedGraph"
    [[ -s $anti && -s $control_bigwig ]] || {
        echo "E: CUT&RUN anti/control input is missing for $sample." >&2
        exit 1
    }
    if [[ ! -e $control ]]; then
        staging="$control.export-job-${SLURM_JOB_ID:-manual}-pid-$$"
        [[ ! -e $staging ]] || { echo "E: Control staging exists: $staging" >&2; exit 1; }
        "$bigwig_to_bedgraph" -chrom=1 "$control_bigwig" "$staging"
        [[ -s $staging ]] || { echo "E: Empty BigWig export for $sample." >&2; exit 1; }
        [[ ! -e $control ]] || { echo "E: Control appeared concurrently: $control" >&2; exit 1; }
        mv "$staging" "$control"
    fi
    coverage_arguments+=(--coverage "supported_anti_${sample}=$anti")
    coverage_arguments+=(--coverage "supported_control_${sample}=$control")
done

phase=building_anchor_evidence
sample=all
"$source/scripts/build_tp73_anchor_evidence.py" \
    --anchor-plus "$plus" --anchor-minus "$minus" \
    "${coverage_arguments[@]}" --output "$output" --chrom 1 \
    --minimum-anchor-score 0 --duckdb "$duckdb" \
    --threads "$threads" --memory-limit "$memory_limit"

phase=validating
valid=$("$duckdb" -csv -noheader :memory: -c "
DESCRIBE SELECT * FROM read_parquet('$output');
SELECT count(*) = 310782
   AND min(anchor_score) >= 0
   AND count(*) FILTER (
       WHERE supported_anti_saos2_TA IS NULL
          OR supported_control_saos2_TA IS NULL
          OR supported_anti_skmel29_2_DN IS NULL
          OR supported_control_skmel29_2_DN IS NULL
   ) = 0
FROM read_parquet('$output');" | tail -1)
[[ $valid == true ]] || { echo "E: TP73 anchor evidence failed validation." >&2; exit 1; }
phase=complete
echo "I: Completed TP73 anchor evidence: $output" >&2

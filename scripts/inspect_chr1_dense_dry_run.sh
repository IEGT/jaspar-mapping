#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
package_dir="$repository_root/dry_runs/chr1_patz1_tp73_from-3600000-to-3700000"
database_name=jaspar2026_chr1_patz1_tp73.duckdb

usage() {
    cat <<'EOF'
Usage: inspect_chr1_dense_dry_run.sh [--package DIR] COMMAND [ARGS]

Commands:
  overview
  files
  shell
  examples
  region MOTIF MODE STRAND START END [PSEUDOCOUNT]
  summary MOTIF MODE STRAND START END [PSEUDOCOUNT]
  histogram MOTIF MODE STRAND BIN_WIDTH START END [PSEUDOCOUNT]
  bins MOTIF MODE STRAND START END [PSEUDOCOUNT]

MODE is log2_relative_risk or log_odds; STRAND is + or -.
All coordinates are BED-style, 0-based, half-open alignment-start ranges.
EOF
}

if [[ ${1:-} == "--package" ]]; then
    [[ $# -ge 2 ]] || { usage >&2; exit 2; }
    package_dir=$2
    shift 2
fi

command=${1:-overview}
[[ $# -eq 0 ]] || shift

command -v duckdb >/dev/null 2>&1 || {
    echo "E: duckdb CLI is required." >&2
    exit 1
}

[[ -d $package_dir ]] || {
    echo "E: Package directory not found: $package_dir" >&2
    exit 1
}
package_dir=$(cd "$package_dir" && pwd -P)
[[ -f "$package_dir/$database_name" ]] || {
    echo "E: DuckDB file not found: $package_dir/$database_name" >&2
    exit 1
}

validate_motif() {
    [[ $1 =~ ^MA[0-9]{4}\.[0-9]+$ ]] || {
        echo "E: Invalid motif ID: $1" >&2
        exit 2
    }
}

validate_mode() {
    [[ $1 == "log2_relative_risk" || $1 == "log_odds" ]] || {
        echo "E: Invalid score mode: $1" >&2
        exit 2
    }
}

validate_strand() {
    [[ $1 == "+" || $1 == "-" ]] || {
        echo "E: Strand must be + or -." >&2
        exit 2
    }
}

validate_integer() {
    [[ $1 =~ ^[0-9]+$ ]] || {
        echo "E: Expected a non-negative integer, got: $1" >&2
        exit 2
    }
}

validate_number() {
    [[ $1 =~ ^[0-9]+([.][0-9]+)?$ ]] || {
        echo "E: Expected a non-negative number, got: $1" >&2
        exit 2
    }
}

run_query() {
    local sql=$1
    (
        cd "$package_dir"
        duckdb -readonly -box "$database_name" -c "$sql"
    )
}

read_region_arguments() {
    [[ $# -ge 5 && $# -le 6 ]] || { usage >&2; exit 2; }
    motif=$1
    mode=$2
    strand=$3
    start=$4
    end=$5
    pseudocount=${6:-1}
    validate_motif "$motif"
    validate_mode "$mode"
    validate_strand "$strand"
    validate_integer "$start"
    validate_integer "$end"
    validate_number "$pseudocount"
    (( end > start )) || { echo "E: END must be greater than START." >&2; exit 2; }
}

case "$command" in
    overview)
        run_query "SELECT * FROM run_manifest;
                   SELECT * FROM motif_metadata ORDER BY motif_id;
                   SELECT * FROM dense_run_inventory ORDER BY motif_id, score_mode, strand;"
        ;;
    files)
        run_query "SELECT motif_id, score_mode, pseudocount, chrom, strand,
                          block_start, len(scores) AS n_windows, source_file
                   FROM motif_score_dense_block
                   ORDER BY motif_id, score_mode, strand, block_start;"
        ;;
    shell)
        cd "$package_dir"
        exec duckdb -readonly "$database_name"
        ;;
    examples)
        cd "$package_dir"
        exec duckdb -readonly -box "$database_name" \
            -f "$repository_root/sql/chr1_dense_dry_run_examples.sql"
        ;;
    region|summary|bins)
        read_region_arguments "$@"
        macro=dense_scores_region
        [[ $command == summary ]] && macro=dense_score_summary
        [[ $command == bins ]] && macro=dense_score_calibration_bins
        run_query "SELECT * FROM $macro(
                       '$motif', '$mode', $pseudocount, '1', '$strand', $start, $end
                   ) ORDER BY ALL;"
        ;;
    histogram)
        [[ $# -ge 6 && $# -le 7 ]] || { usage >&2; exit 2; }
        motif=$1
        mode=$2
        strand=$3
        bin_width=$4
        start=$5
        end=$6
        pseudocount=${7:-1}
        validate_motif "$motif"
        validate_mode "$mode"
        validate_strand "$strand"
        validate_number "$bin_width"
        validate_integer "$start"
        validate_integer "$end"
        validate_number "$pseudocount"
        (( end > start )) || { echo "E: END must be greater than START." >&2; exit 2; }
        run_query "SELECT * FROM dense_score_histogram(
                       '$motif', '$mode', $pseudocount, '1', '$strand',
                       $start, $end, $bin_width
                   );"
        ;;
    -h|--help|help)
        usage
        ;;
    *)
        echo "E: Unknown command: $command" >&2
        usage >&2
        exit 2
        ;;
esac

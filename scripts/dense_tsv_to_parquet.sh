#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage: dense_tsv_to_parquet.sh INPUT [OUTPUT.parquet]

Convert dense score blocks from tab-separated text to ZSTD-compressed Parquet.
INPUT may be one .tsv file or a directory, which is searched recursively for
.tsv files. For one file, OUTPUT defaults to the input basename with .parquet;
directory mode writes each .parquet beside its corresponding .tsv file.

Options:
  -h, --help  Show this help and exit

Requirements:
  The DuckDB command-line client must be available in PATH. Input TSV files
  must have `block_start` and `scores` columns, where scores is a FLOAT list.
EOF
}

sql_quote() {
    printf "%s" "$1" | sed "s/'/''/g"
}

convert_one() {
    local input_tsv="$1"
    local output_parquet="${2:-${input_tsv%.tsv}.parquet}"
    local input_sql
    local output_sql

    input_sql=$(sql_quote "$input_tsv")
    output_sql=$(sql_quote "$output_parquet")
    mkdir -p "$(dirname "$output_parquet")"

    duckdb -c "COPY (
        SELECT CAST(block_start AS BIGINT) AS block_start, scores
        FROM read_csv('${input_sql}',
                      delim = '\t',
                      header = true,
                      columns = {'block_start': 'BIGINT', 'scores': 'FLOAT[]'})
    ) TO '${output_sql}' (FORMAT PARQUET, COMPRESSION ZSTD);"
}

case ${1:-} in
    -h|--help)
        usage
        exit 0
        ;;
esac

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    usage >&2
    exit 2
fi

if ! command -v duckdb >/dev/null 2>&1; then
    echo "E: duckdb CLI not found in PATH." >&2
    exit 1
fi

input="$1"
if [ -d "$input" ]; then
    if [ "$#" -ne 1 ]; then
        echo "E: explicit output path is only supported for a single TSV input." >&2
        exit 2
    fi
    find "$input" -name '*.tsv' -print0 | while IFS= read -r -d '' tsv; do
        convert_one "$tsv"
    done
elif [ -f "$input" ]; then
    convert_one "$input" "${2:-}"
else
    echo "E: input not found: $input" >&2
    exit 1
fi

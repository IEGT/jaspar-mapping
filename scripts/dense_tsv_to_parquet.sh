#!/usr/bin/env bash
set -euo pipefail

usage() {
    echo "Usage: $0 <dense-score.tsv | dense-score-directory> [output.parquet]" >&2
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

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    usage
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

#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: bedtools_map_serial.sh REGIONS.bed VALUES1.bed [VALUES2.bed ...]

Map column 5 from each VALUES file onto REGIONS with `bedtools map -o max`.
The first result is written to mapped_max_step1.bed in the current directory;
each additional VALUES file produces mapped_max_stepN.bed using the preceding
result as its input.

Options:
  -h, --help  Show this help and exit

Requirements:
  bedtools must be available in PATH. Existing mapped_max_stepN.bed files are
  replaced by the corresponding bedtools invocation.
EOF
}

case ${1:-} in
    -h|--help)
        usage
        exit 0
        ;;
esac

if [ "$#" -lt 2 ]; then
    echo "E: Provide a regions BED and at least one values BED." >&2
    usage >&2
    exit 2
fi

command -v bedtools >/dev/null 2>&1 || {
    echo "E: bedtools is required in PATH." >&2
    exit 1
}

regions="$1"
shift

output="mapped_max_step1.bed"
bedtools map -a "$regions" -b "$1" -null 0 -o max -c 5 > "$output"
echo "Mapped $1 to $output"
shift

step=2
while [ "$#" -gt 0 ]; do
    prev_output="$output"
    output="mapped_max_step${step}.bed"
    # Append new column to previous output
    bedtools map -a "$prev_output" -b "$1" -null 0 -o max -c 5 > "$output"
    #paste "$prev_output" - | cut -f1-$(awk '{print NF}' "$prev_output" | head -1),$(awk '{print NF+1}' "$prev_output" | head -1) > "$output"
    echo "Mapped $1 to $output"
    shift
    step=$((step+1))
done

#!/bin/bash
set -e

usage() {
    cat <<'EOF'
Usage: CHR=CHROMOSOME run_what_is_missing.sh

Ask the Makefile for the expected bidirectional BED archives on CHROMOSOME and
print the paths that are absent from output_ChrCHROMOSOME. No scans are run and
no files are removed.

Options:
  -h, --help  Show this help and exit

Environment:
  CHR  Required chromosome label, for example 1, X, or Y

Run this command from the repository root so that it can invoke `make`.
EOF
}

case ${1:-} in
    -h|--help)
        usage
        exit 0
        ;;
esac
if [ "$#" -ne 0 ]; then
    echo "E: This command accepts no positional arguments." >&2
    usage >&2
    exit 2
fi

if [ -z "${CHR:-}" ]; then
    export CHR=unset
    echo "E: Set 'CHR' environment variable."
    exit 1
fi

outDir="output_Chr$CHR"
if [ -d "${outDir}_done" ]; then
    echo "E: Directory '${outDir}_done' already exists."
    exit 1
fi

mkdir -p "$outDir"

for i in $(make CHR="$CHR" echo_bidirect)
do
    f="$outDir/$i"
    if [ -f "$f" ]; then
        # echo "I: Found $i - skipped"
        echo -n
    else
        # echo "I: Need to make $f"
        # echo make $f
        echo "$f"
    fi
done

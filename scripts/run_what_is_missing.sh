#!/bin/bash
set -e

if [ -z "$CHR" ]; then
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

for i in $(make CHR=$CHR echo_bidirect)
do
    f="$outDir/$i"
    if [ -f "$f" ]; then
        # echo "I: Found $i - skipped"
        echo -n
    else
        # echo "I: Need to make $f"
        # echo make $f
        echo $f
    fi
done

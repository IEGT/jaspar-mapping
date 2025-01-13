#!/bin/bash
set -e

if [ -z "$CHR" ]; then
    export CHR=unset
    echo "E: Set 'CHR' environment variable."
    exit 1
fi

mkdir -p output_Chr$CHR

for i in $(make CHR=$CHR echo_bidirect)
do
    f="output_Chr$CHR/$i"
    if [ -f "$f" ]; then
        # echo "I: Found $i - skipped"
        echo -n
    else
        # echo "I: Need to make $f"
        # echo make $f
        echo $f
    fi
done

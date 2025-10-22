#!/bin/bash

# Usage: ./run_bedtools_map.sh regions.bed values1.bed values2.bed ...

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 regions.bed values1.bed [values2.bed ...]"
    exit 1
fi

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

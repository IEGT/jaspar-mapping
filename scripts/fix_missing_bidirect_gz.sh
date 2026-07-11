#!/bin/bash

set -e

usage() {
    cat <<'EOF'
Usage: CHR=CHROMOSOME fix_missing_bidirect_gz.sh

Inspect output_ChrCHROMOSOME for matching positive- and negative-strand
.bed.gz files whose bidirectional archive is missing. Print shell commands for
performing those merges, inserting `wait` after each group of ten commands.
This script reports commands only; it does not execute merges or remove files.

Options:
  -h, --help  Show this help and exit

Environment:
  CHR  Required chromosome label, for example 1, X, or Y
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
    CHR=_unset
fi

cd "output_Chr$CHR"
prefixes=$(find . -name "*_negative_$CHR.bed.gz" | sed "s/_negative_$CHR.bed.gz//g" |sort | uniq)
#echo $prefixes

#echo "Now in: $(pwd)"

number=0
for i in $prefixes
do
    negbed="${i}_negative_$CHR.bed.gz"
    if [ ! -f "$negbed" ]; then
        echo "Missing file: '$negbed' - some fiddling in parallel with this directory?"
        exit 1
    fi

    posbed="${i}_positive_$CHR.bed.gz"
    if [ ! -f "$posbed" ]; then
        echo "Missing file: '$posbed' - please have a look."
        exit 1
    fi

    bibed="${i}_bidirect_$CHR.bed.gz"

    if [ -f "$bibed" ]; then
        echo "File '$bibed' already exists. Not merging."
    else    
        number=$((number+1))
        echo "Merging files for $i: $negbed and $posbed > $bibed"
        #zcat $negbed $posbed | grep -v ^Chromo | sort  -S 5G -k 1,1 -k2,2n | gzip -c > "${i}_bidirect_$CHR.bed.gz"
        echo "cd output_Chr$CHR && nice zcat $negbed $posbed | grep -v ^Chromo | sort -S 5G -k 1,1 -k2,2n | gzip -c > \"${i}_bidirect_$CHR.bed.gz\" &"
    fi

    if [ 10 -eq "$number" ]; then
	echo "wait"
        number=0
    fi

done

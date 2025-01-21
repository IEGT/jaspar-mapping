#!/bin/bash

set -e

if [ -z "$CHR" ]; then
	CHR=unset
fi


cd output_Chr$CHR
prefixes=$(find . -name "*_negative_$CHR.bed" | sed 's/_negative_$CHR.bed//g' |sort | uniq)

for i in $prefixes
do
    negbed=$i"_negative_$CHR.bed"
    if [ ! -f "$negbed" ]; then
        echo "Missing file: '$negbed' - some fiddling in parallel with this directory?"
        exit 1
    fi

    posbed=$i"_positive_$CHR.bed"
    if [ ! -f "$posbed" ]; then
        echo "Missing file: '$posbed' - please have a look."
        exit 1
    fi

    bibed="${i}_bidirect_$CHR.bed.gz"

    if [ -f "$bibed" ]; then
        echo "File '$bibed' already exists. Not merging."
    else
        echo "Merging files for $i: $negbed and $posbed > $bibed"
        cat $negbed $posbed | grep -v ^Chromo | sort  -k 1,1 -k2,2n | gzip -c > "${i}_bidirect_$CHR.bed.gz"
    fi

    echo "Merge successful, now removing individual files"
    #if [ -f "$negbed.gz" ]; then
    #    echo "File '$negbed.gz' already exists. Removing unpacked version."
        rm "$negbed"
    #else
    #    gzip "$negbed"
    #fi
    #if [ -f "$posbed.gz" ]; then
    #    echo "File '$posbed.gz' already exists. Removing unpacked version."
        rm "$posbed"
    #else
    #    gzip "$posbed"
    #fi

done

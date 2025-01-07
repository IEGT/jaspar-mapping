#!/bin/bash

set -e

cd output_Chr1
prefixes=$(find . -name "*_negative_1.bed" | sed 's/_negative_1.bed//g' |sort | uniq)

for i in $prefixes 
do
    negbed=$i"_negative_1.bed"
    if [ ! -f "$negbed" ]; then
        echo "Missing file: '$negbed' - some fiddling in parallel with this directory?"
        exit 1
    fi

    posbed=$i"_positive_1.bed"
    if [ ! -f "$posbed" ]; then
        echo "Missing file: '$posbed' - please have a look."
        exit 1
    fi

    bibed="${i}_bidirect_1.bed.gz"

    if [ -f "$bibed" ]; then
        echo "File '$bibed' already exists. Not merging."
    else    
        echo "Merging files for $i: $negbed and $posbed > $bibed"
        cat $negbed $posbed | grep -v ^Chromo | sort  -k 1,1 -k2,2n | gzip -c > "${i}_bidirect_1.bed.gz"    
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

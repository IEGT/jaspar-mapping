#!/bin/bash

set -e

cd output_Chr1
prefixes=$(find . -name "*_negative_1.bed.gz" | sed 's/_negative_1.bed.gz//g' |sort | uniq)
#echo $prefixes

#echo "Now in: $(pwd)"

for i in $prefixes 
do
    negbed="${i}_negative_1.bed.gz"
    if [ ! -f "$negbed" ]; then
        echo "Missing file: '$negbed' - some fiddling in parallel with this directory?"
        exit 1
    fi

    posbed="${i}_positive_1.bed.gz"
    if [ ! -f "$posbed" ]; then
        echo "Missing file: '$posbed' - please have a look."
        exit 1
    fi

    bibed="${i}_bidirect_1.bed.gz"

    if [ -f "$bibed" ]; then
        echo "File '$bibed' already exists. Not merging."
    else    
        echo "Merging files for $i: $negbed and $posbed > $bibed"
        #zcat $negbed $posbed | grep -v ^Chromo | sort  -k 1,1 -k2,2n | gzip -c > "${i}_bidirect_1.bed.gz"    
        echo "zcat $negbed $posbed | grep -v ^Chromo | sort  -k 1,1 -k2,2n | gzip -c > \"${i}_bidirect_1.bed.gz\" &"
    fi

done

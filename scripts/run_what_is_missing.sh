#!/bin/bash
set -e

for i in $(make echo_bidirect)
do
    f="output_Chr1/$i"
    if [ -f "$f" ]; then
        # echo "I: Found $i - skipped"
        echo 
    else
        # echo "I: Need to make $f"
        # echo make $f
        echo $f
    fi
done

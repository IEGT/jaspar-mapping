#!/usr/bin/awk -f
# shift_bed.awk: Shift BED regions 500 bp downstream based on strand

{
    chrom = $1
    start = $2
    end   = $3
    strand = ($6 == "+" || $6 == "-") ? $6 : "+"

    if (strand == "+") {
        start += 500
        end   += 500
    } else if (strand == "-") {
        start -= 500
        end   -= 500
        if (start < 0) start = 0
        if (end < 0) end = 0
    } else {
        print "Error: Invalid strand value '" strand "' in line: " $0 > "/dev/stderr"
        next
    }

    $2 = start
    $3 = end
    OFS = "\t"
    print $0
}
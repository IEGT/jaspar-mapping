#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: shift_bed.awk [BED ...]

Shift every BED interval 500 bp in its strand direction. Plus-strand intervals
move toward larger coordinates; minus-strand intervals move toward zero and
are clipped at zero. A missing or unrecognized column 6 is treated as plus.
Results are written as tab-separated BED rows to standard output.

With no BED arguments, input is read from standard input. This executable is a
small shell launcher around its embedded AWK transformation.

Options:
  -h, --help  Show this help and exit
EOF
}

case ${1:-} in
    -h|--help)
        usage
        exit 0
        ;;
esac

awk '
BEGIN {
    OFS = "\t"
}

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
        print "Error: Invalid strand value \047" strand "\047 in line: " $0 > "/dev/stderr"
        next
    }

    $2 = start
    $3 = end
    print $0
}
' "$@"

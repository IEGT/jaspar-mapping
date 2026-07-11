#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: CHR=CHROMOSOME fix_missing_bidirect.sh

Create missing bidirectional .bed.gz archives in output_ChrCHROMOSOME by
combining matching uncompressed positive- and negative-strand BED files.
Each new archive is sorted, gzip-validated, and row-count-validated before its
two source BED files are removed. Existing archives are never replaced and
their source files are retained.

Options:
  -h, --help  Show this help and exit

Environment:
  CHR  Required chromosome label, for example 1, X, or Y

Run this command from the repository root. This command deletes source BED
files only after a newly created replacement archive passes all validations.
EOF
}

case ${1:-} in
    -h|--help)
        usage
        exit 0
        ;;
esac
if [[ $# -ne 0 ]]; then
    echo "E: This command accepts no positional arguments." >&2
    usage >&2
    exit 2
fi

CHR=${CHR:-unset}
output_directory="output_Chr${CHR}"

if [ ! -d "$output_directory" ]; then
    echo "E: Output directory '$output_directory' does not exist." >&2
    exit 1
fi
cd "$output_directory"

# This script removes the strand-specific inputs only after it has created and
# validated a new bidirectional archive. Existing archives are never treated
# as proof that the current source files were merged into them.
prefixes=$(find . -type f -name "*_negative_${CHR}.bed" |
    sed "s/_negative_${CHR}\\.bed$//" |
    LC_ALL=C sort -u)

temporary_output=""
cleanup_temporary_output() {
    if [ -n "${temporary_output:-}" ] && [ -e "$temporary_output" ]; then
        rm -f "$temporary_output"
    fi
}
trap cleanup_temporary_output EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

for prefix in $prefixes; do
    negbed="${prefix}_negative_${CHR}.bed"
    posbed="${prefix}_positive_${CHR}.bed"
    bibed="${prefix}_bidirect_${CHR}.bed.gz"

    if [ ! -f "$negbed" ]; then
        echo "E: Missing negative-strand file '$negbed'." >&2
        exit 1
    fi
    if [ ! -f "$posbed" ]; then
        echo "E: Missing positive-strand file '$posbed'." >&2
        exit 1
    fi

    if [ -e "$bibed" ]; then
        if gzip -t "$bibed"; then
            echo "I: '$bibed' already exists and is a valid gzip file; leaving both source BEDs untouched."
        else
            echo "E: '$bibed' already exists but fails gzip validation; leaving both source BEDs untouched." >&2
            exit 1
        fi
        continue
    fi

    expected_rows=$(awk '!/^Chromo/ { rows++ } END { print rows + 0 }' \
        "$negbed" "$posbed")
    if [ "$expected_rows" -eq 0 ]; then
        echo "E: '$negbed' and '$posbed' contain no data rows; source BEDs retained." >&2
        exit 1
    fi

    temporary_output=$(mktemp "${bibed}.tmp.XXXXXX")
    echo "I: Merging '$negbed' and '$posbed' into temporary archive '$temporary_output'."
    if ! awk '!/^Chromo/' "$negbed" "$posbed" |
        LC_ALL=C sort -k 1,1 -k2,2n |
        gzip -n -c > "$temporary_output"; then
        echo "E: Merge pipeline failed; source BEDs retained." >&2
        exit 1
    fi

    if [ ! -s "$temporary_output" ]; then
        echo "E: Merge produced an empty archive; source BEDs retained." >&2
        exit 1
    fi
    if ! gzip -t "$temporary_output"; then
        echo "E: Merge produced an invalid gzip archive; source BEDs retained." >&2
        exit 1
    fi

    actual_rows=$(gzip -cd "$temporary_output" | awk 'END { print NR + 0 }')
    if [ "$actual_rows" -ne "$expected_rows" ]; then
        echo "E: Merge row-count mismatch: expected $expected_rows, got $actual_rows; source BEDs retained." >&2
        exit 1
    fi

    # Do not replace a file that appeared while the temporary output was being
    # built. This can happen if another cleanup process is operating in the
    # same directory.
    if [ -e "$bibed" ]; then
        echo "E: '$bibed' appeared during the merge; source BEDs retained." >&2
        exit 1
    fi

    mv "$temporary_output" "$bibed"
    temporary_output=""
    echo "I: Validated '$bibed' with $actual_rows rows. Removing its two source BEDs."
    rm "$negbed" "$posbed"
done

#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-fasta-index.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

printf '>chr1 description\nACGT\nAC\n>chr2\nNNN\n' > "$temporary/tiny.fasta"
python3 "$repository_root/scripts/build_fasta_index.py" \
    "$temporary/tiny.fasta" >/dev/null 2>&1

printf 'chr1\t6\t18\t4\t5\nchr2\t3\t32\t3\t4\n' > "$temporary/expected.fai"
cmp "$temporary/expected.fai" "$temporary/tiny.fasta.fai"

printf '>crlf\r\nAC\r\nG\r\n' > "$temporary/crlf.fasta"
python3 "$repository_root/scripts/build_fasta_index.py" \
    "$temporary/crlf.fasta" >/dev/null 2>&1
printf 'crlf\t3\t7\t2\t4\n' > "$temporary/expected-crlf.fai"
cmp "$temporary/expected-crlf.fai" "$temporary/crlf.fasta.fai"

printf '>full\nAC\nGT' > "$temporary/full-final-line.fasta"
python3 "$repository_root/scripts/build_fasta_index.py" \
    "$temporary/full-final-line.fasta" >/dev/null 2>&1
printf 'full\t4\t6\t2\t3\n' > "$temporary/expected-full-final-line.fai"
cmp "$temporary/expected-full-final-line.fai" \
    "$temporary/full-final-line.fasta.fai"

printf '>bad\nACGT\nACGTA\n' > "$temporary/invalid.fasta"
if python3 "$repository_root/scripts/build_fasta_index.py" \
    "$temporary/invalid.fasta" >/dev/null 2>&1; then
    echo "E: FASTA indexer accepted inconsistent line widths." >&2
    exit 1
fi
[[ ! -e "$temporary/invalid.fasta.fai" ]] || {
    echo "E: FASTA indexer retained output after a validation error." >&2
    exit 1
}

echo "FASTA index builder tests passed."

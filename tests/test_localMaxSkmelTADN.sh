#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

input="$tmp_dir/input.tsv"
actual="$tmp_dir/actual.tsv"
expected="$tmp_dir/expected.tsv"

make_row() {
    local marker=$1
    local motif_strand=$2
    local dn_activity=$3
    local ta_activity=$4
    local gene=$5
    local promoter_strand=$6

    printf '1\t%s\t%s\t%s\t0\t%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t%s\t0\t%s\t1\t0\t1\t%s\t0\t%s\t1\n' \
        "$marker" "$marker" "$marker" "$motif_strand" \
        "$dn_activity" "$ta_activity" "$gene" "$promoter_strand"
}

{
    # Higher activity must win even if the lower-activity row has strand agreement.
    make_row lower-activity + 10 0 GENE1 +
    make_row higher-activity + 11 0 GENE1 -

    # Strand agreement breaks an exact activity tie.
    make_row tied-mismatch - 4 3 GENE2 +
    make_row tied-match + 5 2 GENE2 +

    # The final group exercises the EOF flush.
    make_row final-gene - 6 3 CACNA1B -
} > "$input"

{
    make_row higher-activity + 11 0 GENE1 -
    make_row tied-match + 5 2 GENE2 +
    make_row final-gene - 6 3 CACNA1B -
} > "$expected"

perl "$repo_dir/OverlapTfPromoters/localMaxSkmelTADN.pl" "$input" > "$actual"
diff -u "$expected" "$actual"

echo "localMaxSkmelTADN regression test passed"

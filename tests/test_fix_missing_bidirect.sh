#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
script="$repository_root/scripts/fix_missing_bidirect.sh"
test_root=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-bidirect-test.XXXXXX")
trap 'rm -rf "$test_root"' EXIT HUP INT TERM

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

existing_case="$test_root/existing"
mkdir -p "$existing_case/output_Chr1"
printf 'Chromosome\tFrom\tTo\n1\t20\t30\n' > \
    "$existing_case/output_Chr1/EXISTING_negative_1.bed"
printf 'Chromosome\tFrom\tTo\n1\t10\t15\n' > \
    "$existing_case/output_Chr1/EXISTING_positive_1.bed"
printf 'previous output\n' | gzip -n -c > \
    "$existing_case/output_Chr1/EXISTING_bidirect_1.bed.gz"

(
    cd "$existing_case"
    CHR=1 bash "$script"
)

[ -f "$existing_case/output_Chr1/EXISTING_negative_1.bed" ] ||
    fail "negative source was removed when combined output already existed"
[ -f "$existing_case/output_Chr1/EXISTING_positive_1.bed" ] ||
    fail "positive source was removed when combined output already existed"
existing_contents=$(gzip -cd \
    "$existing_case/output_Chr1/EXISTING_bidirect_1.bed.gz")
[ "$existing_contents" = "previous output" ] ||
    fail "existing combined output was replaced"

success_case="$test_root/success"
mkdir -p "$success_case/output_Chr1"
printf 'Chromosome\tFrom\tTo\n1\t20\t30\n' > \
    "$success_case/output_Chr1/MERGE_negative_1.bed"
printf 'Chromosome\tFrom\tTo\n1\t10\t15\n' > \
    "$success_case/output_Chr1/MERGE_positive_1.bed"

(
    cd "$success_case"
    CHR=1 bash "$script"
)

[ ! -e "$success_case/output_Chr1/MERGE_negative_1.bed" ] ||
    fail "negative source remained after a validated merge"
[ ! -e "$success_case/output_Chr1/MERGE_positive_1.bed" ] ||
    fail "positive source remained after a validated merge"
gzip -t "$success_case/output_Chr1/MERGE_bidirect_1.bed.gz" ||
    fail "combined output is not a valid gzip file"
merged_contents=$(gzip -cd "$success_case/output_Chr1/MERGE_bidirect_1.bed.gz")
expected_contents=$(printf '1\t10\t15\n1\t20\t30')
[ "$merged_contents" = "$expected_contents" ] ||
    fail "combined output is incomplete or incorrectly sorted"

empty_case="$test_root/empty"
mkdir -p "$empty_case/output_Chr1"
printf 'Chromosome\tFrom\tTo\n' > \
    "$empty_case/output_Chr1/EMPTY_negative_1.bed"
printf 'Chromosome\tFrom\tTo\n' > \
    "$empty_case/output_Chr1/EMPTY_positive_1.bed"

if (
    cd "$empty_case"
    CHR=1 bash "$script"
); then
    fail "header-only inputs unexpectedly produced a combined output"
fi

[ -f "$empty_case/output_Chr1/EMPTY_negative_1.bed" ] ||
    fail "negative source was removed after failed validation"
[ -f "$empty_case/output_Chr1/EMPTY_positive_1.bed" ] ||
    fail "positive source was removed after failed validation"
[ ! -e "$empty_case/output_Chr1/EMPTY_bidirect_1.bed.gz" ] ||
    fail "invalid combined output was installed"

echo "Bidirectional cleanup tests passed."

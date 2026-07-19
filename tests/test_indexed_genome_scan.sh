#!/usr/bin/env bash
set -euo pipefail

repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
pssm_scan_bin=${PSSM_SCAN_BIN:-"$repo_dir/pssm_scan"}
tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

fasta="$tmp_dir/genome.fa"
fasta_index="$fasta.fai"
pssm="$tmp_dir/motif.pfm"
output_dir="$tmp_dir/output"

printf '>chr1\nAAAA\n>chr2\nCCCC\n>chr10\nGGGG\n' > "$fasta"
printf 'chr1\t4\t6\t4\t5\nchr2\t4\t17\t4\t5\nchr10\t4\t29\t4\t5\n' > "$fasta_index"
printf '>TEST.1 Any\nA [ 1 ]\nC [ 1 ]\nG [ 1 ]\nT [ 1 ]\n' > "$pssm"

"$pssm_scan_bin" \
    --genome "$fasta" \
    --pssm "$pssm" \
    --motif TEST.1 \
    --motif-set-id synthetic_test_motifs \
    --genome-id synthetic_multichrom_v1 \
    --outdir "$output_dir" \
    --threshold 0 \
    --strand both \
    --coordinate-mode bed \
    > "$tmp_dir/stdout.log" 2> "$tmp_dir/stderr.log"

positive_file=$(find "$output_dir" -type f -name 'Any_TEST.1_positive.bed' -print)
negative_file=$(find "$output_dir" -type f -name 'Any_TEST.1_negative.bed' -print)
[[ -f "$positive_file" ]] || fail "positive-strand BED was not created"
[[ -f "$negative_file" ]] || fail "negative-strand BED was not created"

expected_order=$(printf 'chr1\nchr2\nchr10')
for bed_file in "$positive_file" "$negative_file"; do
    header_count=$(awk '$1 == "Chromosome" { count++ } END { print count + 0 }' "$bed_file")
    [[ "$header_count" -eq 1 ]] || fail "$bed_file has $header_count headers"

    row_count=$(awk 'NR > 1 { count++ } END { print count + 0 }' "$bed_file")
    [[ "$row_count" -eq 12 ]] || fail "$bed_file has $row_count data rows, expected 12"

    actual_order=$(awk 'NR > 1 && !seen[$1]++ { print $1 }' "$bed_file")
    [[ "$actual_order" == "$expected_order" ]] ||
        fail "$bed_file chromosome order is '$actual_order', expected '$expected_order'"
done

# A new run must truncate on the first indexed chromosome, then append later ones.
"$pssm_scan_bin" \
    --genome "$fasta" \
    --pssm "$pssm" \
    --motif TEST.1 \
    --outdir "$output_dir" \
    --threshold 0 \
    --strand both \
    --coordinate-mode bed \
    > "$tmp_dir/rerun.stdout" 2> "$tmp_dir/rerun.stderr"

for bed_file in "$positive_file" "$negative_file"; do
    row_count=$(awk 'NR > 1 { count++ } END { print count + 0 }' "$bed_file")
    [[ "$row_count" -eq 12 ]] ||
        fail "$bed_file retained stale rows after rerun ($row_count rows)"
done

distribution_dir="$tmp_dir/distribution"
mkdir -p "$distribution_dir"
"$pssm_scan_bin" \
    --genome "$fasta" \
    --fasta-index "$fasta_index" \
    --pssm "$pssm" \
    --motif TEST.1 \
    --outdir "$distribution_dir" \
    --score-distribution \
    --distribution-bin-width 1 \
    --strand both \
    > "$tmp_dir/distribution.stdout" 2> "$tmp_dir/distribution.stderr"

distribution_file=$(find "$distribution_dir" -type f -name 'Any_TEST.1_score_distribution_*.tsv' -print)
[[ -f "$distribution_file" ]] || fail "score-distribution output was not created"
distribution_count=$(awk -F '\t' 'NR > 1 { count += $16 } END { print count + 0 }' "$distribution_file")
[[ "$distribution_count" -eq 24 ]] ||
    fail "score distribution counted $distribution_count windows, expected 24"

dense_dir="$tmp_dir/dense"
"$pssm_scan_bin" \
    --genome "$fasta" \
    --fasta-index "$fasta_index" \
    --pssm "$pssm" \
    --motif TEST.1 \
    --chr chr2 \
    --outdir "$dense_dir" \
    --dense-scores \
    --motif-set-id synthetic_test_motifs \
    --genome-id synthetic_multichrom_v1 \
    --strand both \
    > "$tmp_dir/dense.stdout" 2> "$tmp_dir/dense.stderr"

dense_file_count=$(find "$dense_dir" -type f -name '*.tsv' | wc -l | tr -d ' ')
[[ "$dense_file_count" -eq 2 ]] ||
    fail "dense scan created $dense_file_count files, expected 2"

regions="$tmp_dir/regions.tsv"
printf 'chr10\t0\t4\tlate\nchr1\t0\t4\tearly\n' > "$regions"
regions_dir="$tmp_dir/regions-output"
"$pssm_scan_bin" \
    --genome "$fasta" \
    --fasta-index "$fasta_index" \
    --pssm "$pssm" \
    --motif TEST.1 \
    --regions "$regions" \
    --outdir "$regions_dir" \
    --threshold 0 \
    --strand + \
    --coordinate-mode bed \
    > "$tmp_dir/regions.stdout" 2> "$tmp_dir/regions.stderr"

regions_file=$(find "$regions_dir" -type f -name 'Any_TEST.1_positive.bed' -print)
[[ -f "$regions_file" ]] || fail "regions BED was not created"
regions_order=$(awk 'NR > 1 && !seen[$1]++ { print $1 }' "$regions_file")
expected_regions_order=$(printf 'chr1\nchr10')
[[ "$regions_order" == "$expected_regions_order" ]] ||
    fail "regions BED chromosome order is '$regions_order', expected '$expected_regions_order'"

if "$pssm_scan_bin" \
    --genome "$fasta" \
    --fasta-index "$tmp_dir/missing.fai" \
    --pssm "$pssm" \
    --motif TEST.1 \
    --outdir "$tmp_dir/missing-index-output" \
    --strand + \
    > "$tmp_dir/missing-index.stdout" 2> "$tmp_dir/missing-index.stderr"; then
    fail "scan without a readable FASTA index unexpectedly succeeded"
fi

grep -q 'Indexed scanning requires a readable FASTA index' "$tmp_dir/missing-index.stderr" ||
    fail "missing-index failure did not explain the indexed-scan requirement"

echo "Indexed chromosome-at-a-time scan tests passed"

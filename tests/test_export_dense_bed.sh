#!/usr/bin/env bash

set -euo pipefail

command -v duckdb >/dev/null 2>&1 || {
    echo "E: duckdb CLI is required for this test." >&2
    exit 1
}

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-export-bed.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

metadata_dir="$temporary/tables/jaspar2026/motif_metadata"
dense_root="$temporary/tables/jaspar2026/motif_score_dense"
motif_set_id=synthetic_dense_v1
genome_id=synthetic_genome_v1
partition_prefix="motif_set_id=$motif_set_id/genome_id=$genome_id"
scoring_suffix="background_model_id=uniform_acgt_v1/pseudocount_scheme=additive_per_base"
mkdir -p "$metadata_dir"

for partition in \
    "motif_id=MA0001.1/score_mode=log_odds/pseudocount=1/$scoring_suffix/chrom=1/strand=plus" \
    "motif_id=MA0001.1/score_mode=log_odds/pseudocount=1/$scoring_suffix/chrom=1/strand=minus" \
    "motif_id=MA0001.1/score_mode=log_odds/pseudocount=1/$scoring_suffix/chrom=2/strand=plus" \
    "motif_id=MA0001.1/score_mode=log_odds/pseudocount=1/$scoring_suffix/chrom=2/strand=minus" \
    "motif_id=MA0001.1/score_mode=log2_relative_risk/pseudocount=1/$scoring_suffix/chrom=1/strand=plus" \
    "motif_id=MA0002.1/score_mode=log_odds/pseudocount=1/$scoring_suffix/chrom=1/strand=plus"
do
    mkdir -p "$dense_root/$partition_prefix/$partition"
done

(
    cd "$temporary"
    duckdb :memory: -bail -c "
        COPY (
            SELECT * FROM (
                VALUES
                    ('$motif_set_id'::VARCHAR, 'MA0001.1'::VARCHAR,
                     'TEST'::VARCHAR, 4::INTEGER,
                     2026::INTEGER, 'fixture'::VARCHAR),
                    ('$motif_set_id'::VARCHAR, 'MA0002.1'::VARCHAR,
                     'OTHER'::VARCHAR, 2::INTEGER,
                     2026::INTEGER, 'fixture'::VARCHAR)
            ) AS t(motif_set_id, motif_id, motif_name, motif_length,
                   jaspar_version, source_sha256)
        ) TO 'tables/jaspar2026/motif_metadata/part-000000.parquet'
          (FORMAT PARQUET, COMPRESSION ZSTD);

        COPY (SELECT 100::BIGINT AS block_start,
                     [1.25::FLOAT, NULL::FLOAT, -1.5::FLOAT, 2.75::FLOAT] AS scores)
          TO 'tables/jaspar2026/motif_score_dense/$partition_prefix/motif_id=MA0001.1/score_mode=log_odds/pseudocount=1/$scoring_suffix/chrom=1/strand=plus/part-000000.parquet'
          (FORMAT PARQUET, COMPRESSION ZSTD);
        COPY (SELECT 100::BIGINT AS block_start,
                     [0.5::FLOAT, 4.0::FLOAT, -2.0::FLOAT, NULL::FLOAT] AS scores)
          TO 'tables/jaspar2026/motif_score_dense/$partition_prefix/motif_id=MA0001.1/score_mode=log_odds/pseudocount=1/$scoring_suffix/chrom=1/strand=minus/part-000000.parquet'
          (FORMAT PARQUET, COMPRESSION ZSTD);
        COPY (SELECT 200::BIGINT AS block_start,
                     [3.5::FLOAT, -3.0::FLOAT] AS scores)
          TO 'tables/jaspar2026/motif_score_dense/$partition_prefix/motif_id=MA0001.1/score_mode=log_odds/pseudocount=1/$scoring_suffix/chrom=2/strand=plus/part-000000.parquet'
          (FORMAT PARQUET, COMPRESSION ZSTD);
        COPY (SELECT 200::BIGINT AS block_start,
                     [5.0::FLOAT, 0.0::FLOAT] AS scores)
          TO 'tables/jaspar2026/motif_score_dense/$partition_prefix/motif_id=MA0001.1/score_mode=log_odds/pseudocount=1/$scoring_suffix/chrom=2/strand=minus/part-000000.parquet'
          (FORMAT PARQUET, COMPRESSION ZSTD);
        COPY (SELECT 100::BIGINT AS block_start, [9.0::FLOAT] AS scores)
          TO 'tables/jaspar2026/motif_score_dense/$partition_prefix/motif_id=MA0001.1/score_mode=log2_relative_risk/pseudocount=1/$scoring_suffix/chrom=1/strand=plus/part-000000.parquet'
          (FORMAT PARQUET, COMPRESSION ZSTD);
        COPY (SELECT 100::BIGINT AS block_start, [8.0::FLOAT] AS scores)
          TO 'tables/jaspar2026/motif_score_dense/$partition_prefix/motif_id=MA0002.1/score_mode=log_odds/pseudocount=1/$scoring_suffix/chrom=1/strand=plus/part-000000.parquet'
          (FORMAT PARQUET, COMPRESSION ZSTD);" >/dev/null

    duckdb fixture.duckdb -bail \
        -f "$repository_root/sql/chr1_dense_dry_run_schema.sql" >/dev/null
)

exporter=(
    python3 "$repository_root/scripts/export_dense_bed.py"
    --package "$temporary"
    --genome-id "$genome_id"
    --motif-set-id "$motif_set_id"
)

assert_equal() {
    local label=$1
    local expected=$2
    local actual=$3
    if [[ $actual != "$expected" ]]; then
        echo "E: $label" >&2
        diff -u <(printf '%s\n' "$expected") <(printf '%s\n' "$actual") >&2 || true
        exit 1
    fi
}

configs=$("${exporter[@]}" --list-configs)
grep -F "$genome_id"$'\t'"$motif_set_id"$'\tMA0001.1\tTEST\tlog_odds\t1.0\t2\t-' <<<"$configs" >/dev/null || {
    echo "E: --list-configs omitted a stored partition." >&2
    exit 1
}

default_export=$("${exporter[@]}" \
    --motif MA0001.1 --score-mode log_odds \
    --chrom 1 --from 100 --to 104)
assert_equal "default score threshold or BED ordering is wrong" \
    $'Chromosome\tFrom\tTo\tName\tScore\tStrand\n1\t100\t104\tTEST\t1.25\t+\n1\t100\t104\tTEST\t0.5\t-\n1\t101\t105\tTEST\t4.0\t-\n1\t103\t107\tTEST\t2.75\t+' \
    "$default_export"

negative_export=$("${exporter[@]}" \
    --motif MA0001.1 --score-mode log_odds \
    --chrom 1 --strand + --from 102 --to 103 \
    --all-scores --max-score -1 --no-header)
assert_equal "all-score or upper-bound filtering is wrong" \
    $'1\t102\t106\tTEST\t-1.5\t+' "$negative_export"

provenance_export=$("${exporter[@]}" \
    --motif MA0001.1 --score-mode log_odds \
    --chrom 2 --min-score 2 --columns provenance --score-decimals 3)
assert_equal "provenance columns or score formatting are wrong" \
    $'Chromosome\tFrom\tTo\tName\tScore\tStrand\tGenomeID\tMotifSetID\tMotifID\tScoreMode\tPseudocount\tBackgroundModelID\tPseudocountScheme\tCoordinateMode\n2\t200\t204\tTEST\t3.500\t+\t'$genome_id$'\t'$motif_set_id$'\tMA0001.1\tlog_odds\t1.0\tuniform_acgt_v1\tadditive_per_base\tbed\n2\t200\t204\tTEST\t5.000\t-\t'$genome_id$'\t'$motif_set_id$'\tMA0001.1\tlog_odds\t1.0\tuniform_acgt_v1\tadditive_per_base\tbed' \
    "$provenance_export"

count=$("${exporter[@]}" \
    --motif MA0001.1 --score-mode log_odds \
    --chrom 1,2 --all-scores --count-only)
assert_equal "multi-chromosome count-only result is wrong" $'Matches\n10' "$count"

gzip_output="$temporary/export.bed.gz"
"${exporter[@]}" \
    --motif MA0001.1 --score-mode log_odds \
    --chrom 2 --min-score 3 --no-header --output "$gzip_output" >/dev/null
gzip_rows=$(gzip -cd "$gzip_output")
assert_equal "direct gzip export is wrong" \
    $'2\t200\t204\tTEST\t3.5\t+\n2\t200\t204\tTEST\t5.0\t-' "$gzip_rows"

if "${exporter[@]}" \
    --motif MA0001.1 --score-mode log_odds \
    --chrom 2 --output "$gzip_output" >/dev/null 2>&1
then
    echo "E: exporter overwrote an existing file without --force." >&2
    exit 1
fi

"${exporter[@]}" \
    --motif MA0001.1 --score-mode log_odds \
    --chrom 2 --min-score 4 --no-header --force \
    --output "$gzip_output" >/dev/null
gzip_rows=$(gzip -cd "$gzip_output")
assert_equal "--force did not replace the compressed export" \
    $'2\t200\t204\tTEST\t5.0\t-' "$gzip_rows"

if "${exporter[@]}" \
    --motif MA0001.1 --score-mode log2_relative_risk \
    --chrom 1 --strand both >/dev/null 2>&1
then
    echo "E: exporter accepted a missing chromosome/strand partition." >&2
    exit 1
fi

echo "Dense score BED export tests passed."

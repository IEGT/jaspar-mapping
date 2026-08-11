#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-genome-evidence.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping TP73 genome-evidence manager test." >&2
    exit 0
}
duckdb=$(command -v "$duckdb")

scan="$temporary/scan-package"
annotation="$temporary/annotation-run"
run_root="$temporary/run"
mkdir -p "$scan" "$annotation/final" "$temporary/tracks" \
    "$temporary/runtime" "$temporary/finalizer-tmp"

chromosomes=({1..22} X Y 25)
sequence_values=""
annotation_values=""
for index in "${!chromosomes[@]}"; do
    chrom=${chromosomes[$index]}
    separator=,
    [[ $index -eq $((${#chromosomes[@]} - 1)) ]] && separator=""
    sequence_values+="($index, '$chrom', $((100000 + index)), repeat('a', 64))$separator"
    annotation_chrom=$chrom
    [[ $chrom == 25 ]] && annotation_chrom=MT
    annotation_values+="('tp73_context_anchor', '$annotation_chrom', true, '/synthetic/chr${annotation_chrom}.parquet')$separator"
done

"$duckdb" -batch "$scan/scan.duckdb" >/dev/null <<SQL
CREATE TABLE sequence_region AS
SELECT * FROM (VALUES $sequence_values)
    AS v(sequence_order, chrom, length, sequence_sha256);
SQL
printf '{"database":"scan.duckdb"}\n' > "$scan/manifest.json"

"$duckdb" -batch "$annotation/final/context.duckdb" >/dev/null <<SQL
CREATE TABLE context_file_inventory AS
SELECT * FROM (VALUES $annotation_values)
    AS v(dataset, chrom, is_parquet, absolute_path);
SQL
printf '{}\n' > "$annotation/final/manifest.json"

task_count=$(
    "$repository_root/scripts/manage_tp73_genome_evidence.py" prepare \
        --run-root "$run_root" --annotation-run "$annotation" \
        --scan-package "$scan" --track-root "$temporary/tracks" \
        --runtime-prefix "$temporary/runtime" --source "$repository_root" \
        --track-manifest \
            "$repository_root/config/h3k4me3_cutandrun_tracks_v1.tsv" \
        --duckdb "$duckdb" --scratch-root "$temporary/scratch" \
        --run-id synthetic_tp73_genome_evidence --threads 1 \
        --memory-limit 512MB --minimum-free-run-gb 0 \
        --minimum-free-scratch-gb 0
)
[[ $task_count -eq 25 ]] || {
    echo "E: Expected 25 sequence-region tasks, found $task_count." >&2
    exit 1
}

read -r -d '' evidence_columns <<'SQL' || true
false::BOOLEAN AS supported_tp73_saos2_GFP,
0::DOUBLE AS depth_tp73_saos2_GFP,
false::BOOLEAN AS supported_negative_control_saos2_GFP,
0::DOUBLE AS depth_negative_control_saos2_GFP,
false::BOOLEAN AS supported_tp73_saos2_TA,
0::DOUBLE AS depth_tp73_saos2_TA,
false::BOOLEAN AS supported_negative_control_saos2_TA,
0::DOUBLE AS depth_negative_control_saos2_TA,
false::BOOLEAN AS supported_tp73_saos2_DN,
0::DOUBLE AS depth_tp73_saos2_DN,
false::BOOLEAN AS supported_negative_control_saos2_DN,
0::DOUBLE AS depth_negative_control_saos2_DN,
false::BOOLEAN AS supported_tp73_skmel29_2_GFP,
0::DOUBLE AS depth_tp73_skmel29_2_GFP,
false::BOOLEAN AS supported_negative_control_skmel29_2_GFP,
0::DOUBLE AS depth_negative_control_skmel29_2_GFP,
false::BOOLEAN AS supported_tp73_skmel29_2_TA,
0::DOUBLE AS depth_tp73_skmel29_2_TA,
false::BOOLEAN AS supported_negative_control_skmel29_2_TA,
0::DOUBLE AS depth_negative_control_skmel29_2_TA,
false::BOOLEAN AS supported_tp73_skmel29_2_DN,
0::DOUBLE AS depth_tp73_skmel29_2_DN,
false::BOOLEAN AS supported_negative_control_skmel29_2_DN,
0::DOUBLE AS depth_negative_control_skmel29_2_DN
SQL

sha256_file() {
    shasum -a 256 "$1" | awk '{print $1}'
}

tail -n +2 "$run_root/plan/chromosome_tasks.tsv" |
while IFS=$'\t' read -r _task_index _sequence_order chrom _annotation_chrom \
        chrom_length _sequence_sha256 analysis_partition _primary_inference; do
    child="$run_root/chromosomes/chrom-$chrom/final"
    mkdir -p "$child"
    evidence="$child/tp73_anchor_evidence.parquet"
    "$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT '$chrom'::VARCHAR AS chrom, 10::BIGINT AS anchor_start,
           26::BIGINT AS anchor_end, 0.5::FLOAT AS anchor_score,
           $evidence_columns
) TO '$evidence' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL
    bytes=$(wc -c < "$evidence" | tr -d ' ')
    digest=$(sha256_file "$evidence")
    {
        printf '{"state":"complete","chrom":"%s","chrom_length":%s,' \
            "$chrom" "$chrom_length"
        printf '"analysis_partition":"%s","outputs":[{' \
            "$analysis_partition"
        printf '"relative_path":"tp73_anchor_evidence.parquet",'
        printf '"bytes":%s,"sha256":"%s"}]}\n' "$bytes" "$digest"
    } > "$child/manifest.json"
done

status=$(
    "$repository_root/scripts/manage_tp73_genome_evidence.py" status \
        --run-root "$run_root"
)
grep -Fq $'planned\t25' <<< "$status"
grep -Fq $'complete\t25' <<< "$status"
grep -Fq $'complete_autosome\t22' <<< "$status"
grep -Fq $'complete_sex_chromosome\t2' <<< "$status"
grep -Fq $'complete_mitochondrial_bystander_control\t1' <<< "$status"

finalize=(
    "$repository_root/scripts/manage_tp73_genome_evidence.py" finalize
    --run-root "$run_root" --duckdb "$duckdb" --threads 1
    --memory-limit 512MB --temp-directory "$temporary/finalizer-tmp"
)
"${finalize[@]}"
"${finalize[@]}"

final="$run_root/final/genome_evidence"
(
    cd "$final"
    "$duckdb" -batch tp73_genome_evidence.duckdb >/dev/null <<'SQL'
SELECT CASE WHEN (SELECT count(*) FROM tp73_anchor_evidence_autosome) <> 22
    THEN error('autosome evidence is incomplete') END;
SELECT CASE WHEN (SELECT count(*) FROM tp73_anchor_evidence_sex_chromosome) <> 2
    THEN error('sex-chromosome evidence is incomplete') END;
SELECT CASE WHEN (SELECT count(*) FROM
    tp73_anchor_evidence_mitochondrial_bystander_control) <> 1
    THEN error('mitochondrial bystander evidence is incomplete') END;
SELECT CASE WHEN (SELECT count(*) FROM chromosome_file_inventory) <> 25
    THEN error('chromosome inventory is incomplete') END;
SELECT CASE WHEN (SELECT count(*) FROM tp73_coverage_summary_by_chromosome) <> 300
    THEN error('coverage summary is not chromosome x support-track complete') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM tp73_coverage_summary_by_chromosome
    WHERE supported_anchors <> 0 OR supported_fraction <> 0
) THEN error('zero-support fixture changed during aggregation') END;
SQL
)

python3 - "$final/manifest.json" <<'PY'
import json
import pathlib
import sys

manifest = json.loads(pathlib.Path(sys.argv[1]).read_text())
assert manifest["partition_counts"] == {
    "autosome": 22,
    "mitochondrial_bystander_control": 1,
    "sex_chromosome": 2,
}
assert manifest["partition_chromosomes"]["mitochondrial_bystander_control"] == ["25"]
assert manifest["validation"]["chromosomes"] == 25
assert manifest["validation"]["support_depth_mismatches"] == 0
PY

echo "TP73 genome-evidence manager tests passed."

#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-h3-genome.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping H3K4me3 genome manager test." >&2
    exit 0
}
duckdb=$(command -v "$duckdb")
python=$(command -v python3)

evidence="$temporary/evidence"
run_root="$temporary/run"
runtime="$temporary/runtime"
tracks="$temporary/tracks"
mkdir -p "$evidence/chromosomes" "$runtime/duckdb/bin" "$tracks" \
    "$temporary/scratch" "$temporary/finalizer-tmp"
ln -s "$duckdb" "$runtime/duckdb/bin/duckdb"
ln -s "$python" "$runtime/duckdb/bin/python3"

track_manifest="$repository_root/config/h3k4me3_cutandrun_tracks_v1.tsv"
awk -F '\t' 'NR > 1 && $7 == "true" && ($4 == "h3k4me3" || $4 == "input") {
    print $6
}' "$track_manifest" |
while IFS= read -r filename; do
    touch "$tracks/$filename"
done

sha256_file() {
    shasum -a 256 "$1" | awk '{print $1}'
}

printf 'sequence_order\tchrom\tannotation_chrom\tchrom_length\tanalysis_partition\tprimary_inference\trelative_path\tbytes\tsha256\n' \
    > "$evidence/chromosome_file_inventory.tsv"
chromosomes=({1..22} X Y 25)
for index in "${!chromosomes[@]}"; do
    chrom=${chromosomes[$index]}
    partition=autosome
    primary=true
    if [[ $chrom == X || $chrom == Y ]]; then
        partition=sex_chromosome
        primary=false
    elif [[ $chrom == 25 ]]; then
        partition=mitochondrial_bystander_control
        primary=false
    fi
    relative="chromosomes/chrom-$chrom/tp73_anchor_evidence.parquet"
    output="$evidence/$relative"
    mkdir -p "$(dirname "$output")"
    "$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
  SELECT '$chrom'::VARCHAR AS chrom, 100::BIGINT AS anchor_start,
         116::BIGINT AS anchor_end, 0.5::FLOAT AS anchor_score
) TO '$output' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL
    bytes=$(wc -c < "$output" | tr -d ' ')
    digest=$(sha256_file "$output")
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$index" "$chrom" "$chrom" "$((100000 + index))" "$partition" \
        "$primary" "$relative" "$bytes" "$digest" \
        >> "$evidence/chromosome_file_inventory.tsv"
done
printf '{"state":"complete","database":"synthetic.duckdb"}\n' \
    > "$evidence/manifest.json"

task_count=$(
    "$repository_root/scripts/manage_h3k4me3_genome_signal.py" prepare \
        --run-root "$run_root" --evidence-package "$evidence" \
        --track-root "$tracks" --runtime-prefix "$runtime" \
        --source "$repository_root" --track-manifest "$track_manifest" \
        --scratch-root "$temporary/scratch" --run-id synthetic_h3_genome \
        --threads 1 --memory-limit 512MB --minimum-free-run-gb 0 \
        --minimum-free-scratch-gb 0
)
[[ $task_count -eq 24 ]] || {
    echo "E: Expected 24 nuclear chromosome tasks, found $task_count." >&2
    exit 1
}
[[ $(wc -l < "$run_root/plan/track_file_inventory.tsv" | tr -d ' ') -eq 13 ]] || {
    echo "E: Expected one header and 12 pinned H3K4me3/input tracks." >&2
    exit 1
}
reprepared=$(
    "$repository_root/scripts/manage_h3k4me3_genome_signal.py" prepare \
        --run-root "$run_root" --evidence-package "$evidence" \
        --track-root "$tracks" --runtime-prefix "$runtime" \
        --source "$repository_root" --track-manifest "$track_manifest" \
        --scratch-root "$temporary/scratch" --run-id synthetic_h3_genome \
        --threads 1 --memory-limit 512MB --minimum-free-run-gb 0 \
        --minimum-free-scratch-gb 0
)
[[ $reprepared -eq 24 ]] || {
    echo "E: Immutable H3K4me3 preparation was not idempotent." >&2
    exit 1
}
if grep -Fq $'\t25\t' "$run_root/plan/chromosome_tasks.tsv"; then
    echo "E: Mitochondrial evidence entered the H3K4me3 task plan." >&2
    exit 1
fi

tail -n +2 "$run_root/plan/chromosome_tasks.tsv" |
while IFS=$'\t' read -r _task_index _sequence_order chrom chrom_length \
        partition primary evidence_path _evidence_bytes _evidence_sha; do
    child="$run_root/chromosomes/chrom-$chrom/final"
    mkdir -p "$child"
    signal_path="$child/h3k4me3_anchor_signal.parquet"
    change_path="$child/h3k4me3_anchor_change.parquet"
    "$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
  WITH series AS (
    SELECT * FROM (VALUES
      ('saos2', 'SaOS-2'), ('skmel29_2', 'SK-Mel-29')
    ) AS t(series_id, cell_line)
  ), conditions AS (
    SELECT * FROM (VALUES ('GFP'), ('TA'), ('DN')) AS t(condition)
  ), windows AS (
    SELECT * FROM (VALUES
      ('motif_span', 0, 16, 1, 16),
      ('central_150', 0, 150, 1, 300),
      ('central_500', 0, 500, 1, 1000),
      ('central_1000', 0, 1000, 1, 2000),
      ('flank_150_1000', 150, 1000, 2, 1700)
    ) AS t(window_name, inner_bp, outer_bp, segment_count, effective_window_bp)
  )
  SELECT '$chrom'::VARCHAR AS chrom, 100::BIGINT AS anchor_start,
         116::BIGINT AS anchor_end, 0.5::FLOAT AS anchor_score,
         s.series_id, s.cell_line, c.condition, 'R1'::VARCHAR AS replicate,
         w.window_name, w.inner_bp::BIGINT AS inner_bp,
         w.outer_bp::BIGINT AS outer_bp,
         w.segment_count::INTEGER AS segment_count,
         w.effective_window_bp::BIGINT AS effective_window_bp,
         CASE c.condition WHEN 'TA' THEN 8.0 WHEN 'DN' THEN 1.0 ELSE 2.0 END
           AS h3k4me3_area,
         0.01::DOUBLE AS h3k4me3_mean,
         CASE c.condition WHEN 'TA' THEN 8.0 WHEN 'DN' THEN 1.0 ELSE 2.0 END
           AS h3k4me3_max,
         10::BIGINT AS h3k4me3_positive_bp,
         0.1::DOUBLE AS h3k4me3_positive_fraction,
         1.0::DOUBLE AS input_area, 0.005::DOUBLE AS input_mean,
         1.0::DOUBLE AS input_max, 10::BIGINT AS input_positive_bp,
         0.1::DOUBLE AS input_positive_fraction
  FROM series s CROSS JOIN conditions c CROSS JOIN windows w
  ORDER BY s.series_id, c.condition, w.window_name
) TO '$signal_path' (FORMAT PARQUET, COMPRESSION ZSTD);

COPY (
  SELECT '$chrom'::VARCHAR AS chrom, 100::BIGINT AS anchor_start,
         116::BIGINT AS anchor_end, 0.5::FLOAT AS anchor_score,
         series_id, cell_line, 'R1'::VARCHAR AS replicate,
         'flank_150_1000'::VARCHAR AS window_name,
         150::BIGINT AS inner_bp, 1000::BIGINT AS outer_bp,
         2::INTEGER AS segment_count, 1700::BIGINT AS effective_window_bp,
         2.0::DOUBLE AS gfp_h3k4me3_area, 1.0::DOUBLE AS gfp_input_area,
         2.0::DOUBLE AS gfp_h3k4me3_max, 1.0::DOUBLE AS gfp_input_max,
         8.0::DOUBLE AS ta_h3k4me3_area, 1.0::DOUBLE AS ta_input_area,
         8.0::DOUBLE AS ta_h3k4me3_max, 1.0::DOUBLE AS ta_input_max,
         1.0::DOUBLE AS dn_h3k4me3_area, 1.0::DOUBLE AS dn_input_area,
         1.0::DOUBLE AS dn_h3k4me3_max, 1.0::DOUBLE AS dn_input_max,
         log2(3.0 / 2.0)::DOUBLE AS gfp_log2_h3k4me3_input_ratio,
         log2(9.0 / 2.0)::DOUBLE AS ta_log2_h3k4me3_input_ratio,
         0.0::DOUBLE AS dn_log2_h3k4me3_input_ratio,
         log2(3.0)::DOUBLE AS delta_ta_vs_gfp,
         -log2(3.0 / 2.0)::DOUBLE AS delta_dn_vs_gfp,
         true AS has_any_h3k4me3_signal, true AS has_any_input_signal,
         false AS uninformative_all_h3k4me3_zero
  FROM (VALUES ('saos2', 'SaOS-2'), ('skmel29_2', 'SK-Mel-29'))
       AS t(series_id, cell_line)
) TO '$change_path' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL
    signal_bytes=$(wc -c < "$signal_path" | tr -d ' ')
    change_bytes=$(wc -c < "$change_path" | tr -d ' ')
    signal_sha=$(sha256_file "$signal_path")
    change_sha=$(sha256_file "$change_path")
    cat > "$child/manifest.json" <<JSON
{
  "state": "complete",
  "chrom": "$chrom",
  "chrom_length": $chrom_length,
  "analysis_partition": "$partition",
  "primary_inference": $primary,
  "validation": {"anchors": 1, "signal_rows": 30, "change_rows": 2},
  "outputs": [
    {"relative_path": "h3k4me3_anchor_signal.parquet", "bytes": $signal_bytes, "sha256": "$signal_sha"},
    {"relative_path": "h3k4me3_anchor_change.parquet", "bytes": $change_bytes, "sha256": "$change_sha"}
  ]
}
JSON
done

# Exercise the production change derivation and its full signal/evidence
# cardinality checks on one transparent chromosome fixture.
python3 - "$repository_root" "$duckdb" "$temporary" "$evidence" \
    "$run_root" <<'PY'
import importlib.util
import pathlib
import sys

root, duckdb, temporary, evidence, run_root = map(pathlib.Path, sys.argv[1:])
spec = importlib.util.spec_from_file_location(
    "h3_manager", root / "scripts" / "manage_h3k4me3_genome_signal.py"
)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)
signal = run_root / "chromosomes" / "chrom-1" / "final" / \
    "h3k4me3_anchor_signal.parquet"
source = evidence / "chromosomes" / "chrom-1" / \
    "tp73_anchor_evidence.parquet"
derived = temporary / "derived-change.parquet"
tmp = temporary / "derived-tmp"
tmp.mkdir()
module.build_change_table(duckdb, signal, derived, 1, "512MB", tmp)
row = {
    "chrom": "1", "chrom_length": "100000",
    "analysis_partition": "autosome",
}
validation = module.validate_outputs(
    duckdb, row, source, signal, derived, ["saos2", "skmel29_2"]
)
assert validation["anchors"] == 1
assert validation["signal_rows"] == 30
assert validation["change_rows"] == 2
PY
"$duckdb" -batch :memory: >/dev/null <<SQL
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM read_parquet('$temporary/derived-change.parquet')
  WHERE abs(delta_ta_vs_gfp - log2(3.0)) > 1e-12
     OR abs(delta_dn_vs_gfp + log2(1.5)) > 1e-12
) THEN error('GFP-referenced change arithmetic is incorrect') END;
SQL

status=$(
    "$repository_root/scripts/manage_h3k4me3_genome_signal.py" status \
        --run-root "$run_root"
)
grep -Fq $'planned\t24' <<< "$status"
grep -Fq $'complete\t24' <<< "$status"
grep -Fq $'complete_autosome\t22' <<< "$status"
grep -Fq $'complete_sex_chromosome\t2' <<< "$status"

finalize=(
    "$repository_root/scripts/manage_h3k4me3_genome_signal.py" finalize
    --run-root "$run_root" --duckdb "$duckdb" --threads 1
    --memory-limit 512MB --temp-directory "$temporary/finalizer-tmp"
)
"${finalize[@]}"
"${finalize[@]}"

final="$run_root/final/genome_h3k4me3_signal"
for provenance_file in run_config.json chromosome_tasks.tsv \
        track_file_inventory.tsv track_manifest.tsv; do
    [[ -f "$final/provenance/$provenance_file" ]] || {
        echo "E: Missing finalized provenance file: $provenance_file" >&2
        exit 1
    }
done
(
    cd "$final"
    "$duckdb" -batch genome_h3k4me3_signal.duckdb >/dev/null <<'SQL'
SELECT CASE WHEN (SELECT count(*) FROM h3k4me3_anchor_signal) <> 720
  THEN error('nuclear signal view is incomplete') END;
SELECT CASE WHEN (SELECT count(*) FROM h3k4me3_anchor_change) <> 48
  THEN error('nuclear change view is incomplete') END;
SELECT CASE WHEN (SELECT count(*) FROM h3k4me3_anchor_signal_autosome) <> 660
  THEN error('autosome signal view is incomplete') END;
SELECT CASE WHEN (SELECT count(*) FROM h3k4me3_anchor_signal_sex_chromosome) <> 60
  THEN error('sex-chromosome signal view is incomplete') END;
SELECT CASE WHEN (SELECT count(*) FROM chromosome_file_inventory) <> 48
  THEN error('exact chromosome/dataset inventory is incomplete') END;
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM chromosome_file_inventory WHERE chrom = '25'
) THEN error('mitochondrial evidence entered the H3K4me3 package') END;
SELECT CASE WHEN (SELECT count(*) FROM h3k4me3_change_summary_by_chromosome) <> 48
  THEN error('chromosome-by-series summary is incomplete') END;
SQL
)

python3 - "$final/manifest.json" <<'PY'
import json
import pathlib
import sys

manifest = json.loads(pathlib.Path(sys.argv[1]).read_text())
assert manifest["partition_counts"] == {
    "autosome": 22,
    "sex_chromosome": 2,
}
assert manifest["mitochondrial_policy"] == (
    "excluded_from_histone_mark_inference"
)
assert manifest["validation"] == {
    "anchors": 24,
    "change_rows": 48,
    "chromosomes": 24,
    "signal_rows": 720,
}
PY

echo "H3K4me3 genome-signal manager tests passed."

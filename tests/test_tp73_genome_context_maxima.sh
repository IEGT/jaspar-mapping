#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-genome-context.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping genome-context test." >&2
    exit 0
}
duckdb=$(command -v "$duckdb")

scan="$temporary/scan-package"
evidence="$temporary/evidence-package"
runtime="$temporary/runtime"
run_root="$temporary/context-run"
legacy_run_root="$temporary/legacy-context-run"
unsafe_run_root="$temporary/unsafe-context-run"
enrichment_root="$temporary/enrichment-run"
mkdir -p "$scan/task_data/task_id=1" "$scan/task_data/task_id=2" \
    "$evidence/tables/by-chromosome" "$runtime/duckdb/bin"
ln -s "$duckdb" "$runtime/duckdb/bin/duckdb"

sql_file="$temporary/build-fixture.sql"
inventory_values=""
sequence_values=""
anchor_paths=""
for chrom in {1..22}; do
    m1_plus="$scan/task_data/task_id=1/chr${chrom}-plus.parquet"
    m1_minus="$scan/task_data/task_id=1/chr${chrom}-minus.parquet"
    m2_plus="$scan/task_data/task_id=2/chr${chrom}-plus.parquet"
    m2_minus="$scan/task_data/task_id=2/chr${chrom}-minus.parquet"
    anchor="$evidence/tables/by-chromosome/chr${chrom}.parquet"
    cat >> "$sql_file" <<SQL
COPY (SELECT 120::BIGINT AS start, 130::BIGINT AS "end", 2.0::FLOAT AS score)
TO '$m1_plus' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT NULL::BIGINT AS start, NULL::BIGINT AS "end",
             NULL::FLOAT AS score WHERE false)
TO '$m1_minus' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT NULL::BIGINT AS start, NULL::BIGINT AS "end",
             NULL::FLOAT AS score WHERE false)
TO '$m2_plus' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT NULL::BIGINT AS start, NULL::BIGINT AS "end",
             NULL::FLOAT AS score WHERE false)
TO '$m2_minus' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (
  SELECT '$chrom'::VARCHAR AS chrom, 100::BIGINT AS anchor_start,
         116::BIGINT AS anchor_end, 0.5::FLOAT AS anchor_score,
         true::BOOLEAN AS supported_tp73_s1,
         3.0::DOUBLE AS depth_tp73_s1,
         false::BOOLEAN AS supported_negative_control_s1,
         0.0::DOUBLE AS depth_negative_control_s1
) TO '$anchor' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL
done
"$duckdb" -batch :memory: < "$sql_file" >/dev/null

sha256_file() {
    shasum -a 256 "$1" | awk '{print $1}'
}

for chrom in {1..22}; do
    separator=,
    [[ $chrom -eq 22 ]] && separator=""
    sequence_values+="($((chrom - 1)), '$chrom', 1000)$separator"
    anchor="$evidence/tables/by-chromosome/chr${chrom}.parquet"
    anchor_bytes=$(wc -c < "$anchor" | tr -d ' ')
    anchor_sha=$(sha256_file "$anchor")
    anchor_paths+="$((chrom - 1))\t$chrom\t1000\t"
    anchor_paths+="autosome\ttrue\t"
    anchor_paths+="tables/by-chromosome/chr${chrom}.parquet\t"
    anchor_paths+="$anchor_bytes\t$anchor_sha\n"
    for specification in \
        "M1 1 + chr${chrom}-plus.parquet -1" \
        "M1 1 - chr${chrom}-minus.parquet -1" \
        "M2 2 + chr${chrom}-plus.parquet -1" \
        "M2 2 - chr${chrom}-minus.parquet -1"; do
        read -r motif task strand relative floor <<< "$specification"
        path="$scan/task_data/task_id=$task/$relative"
        bytes=$(wc -c < "$path" | tr -d ' ')
        digest=$(sha256_file "$path")
        emitted=0
        [[ $motif == M1 && $strand == + ]] && emitted=1
        [[ -n $inventory_values ]] && inventory_values+=,
        inventory_values+="('$motif','$chrom','$strand',$task,'$relative',"
        inventory_values+="$bytes,'$digest',$emitted,$floor)"
    done
done

"$duckdb" -batch "$scan/scan.duckdb" >/dev/null <<SQL
CREATE TABLE sequence_region AS
SELECT * FROM (VALUES $sequence_values)
    AS v(sequence_order, chrom, length);
CREATE TABLE motif_metadata AS
SELECT * FROM (VALUES ('M1', 10), ('M2', 12))
    AS v(motif_id, motif_length);
CREATE TABLE scan_motif_threshold AS
SELECT * FROM (VALUES
    ('M1', -1.0, 'scan-density200-v1', false),
    ('M2', -1.0, 'scan-low-floor-v1', false),
    ('MA0861.2', -1.0, 'scan-density200-v1', false)
) AS v(motif_id, final_minimum_score, threshold_set_id, density_limited);
CREATE TABLE scan_file_inventory AS
SELECT * FROM (VALUES $inventory_values) AS v(
    motif_id, chrom, strand, task_id, output_relative_path,
    bytes, sha256, emitted_hits, minimum_score
);
SQL
printf '{"database":"scan.duckdb"}\n' > "$scan/manifest.json"

anchor_list=""
for chrom in {1..22}; do
    [[ -n $anchor_list ]] && anchor_list+=,
    anchor_list+="'$evidence/tables/by-chromosome/chr${chrom}.parquet'"
done
"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
  SELECT * FROM read_parquet([$anchor_list], hive_partitioning=false)
  ORDER BY CAST(chrom AS INTEGER), anchor_start
) TO '$evidence/tables/tp73_anchor_evidence_autosome.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
SQL
printf 'sequence_order\tchrom\tchrom_length\tanalysis_partition\tprimary_inference\trelative_path\tbytes\tsha256\n%b' \
    "$anchor_paths" > "$evidence/chromosome_file_inventory.tsv"
printf '{"state":"complete","primary_inference_partition":"autosome"}\n' \
    > "$evidence/manifest.json"

thresholds="$temporary/source-thresholds.parquet"
"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
  SELECT * FROM (VALUES
    ('M1', 'Motif one',  'chr1-context-v1', 'GRCh38', 'JASPAR2026',
     'convenient_context_threshold', 'MA0861.2', 'log2_relative_risk',
     0.0, 'uniform', 'none', 'chr1_context_150bp', 1.0,
     'informative', NULL::BIGINT, 150::BIGINT, 'any', true,
     'signed_interval_edge_distance', -20.0),
    ('M2', 'Motif two',  'chr1-context-v1', 'GRCh38', 'JASPAR2026',
     'convenient_context_threshold', 'MA0861.2', 'log2_relative_risk',
     0.0, 'uniform', 'none', 'chr1_context_150bp', NULL::DOUBLE,
     'not_estimable', NULL::BIGINT, 150::BIGINT, 'any', true,
     'signed_interval_edge_distance', -20.0)
  ) AS v(
    motif_id, motif_name, threshold_set_id, genome_id, motif_set_id,
    threshold_role, target_motif_id, score_mode, pseudocount,
    background_model_id, pseudocount_scheme, calibration_stratum_id,
    recommended_threshold, calibration_status,
    context_min_interval_distance_bp, context_max_interval_distance_bp,
    context_relation_filter, threshold_inclusive, context_distance_metric,
    source_minimum_score
  )
) TO '$thresholds' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL

legacy_scan="$temporary/legacy-scan-package"
cp -R "$scan" "$legacy_scan"
"$duckdb" -batch "$legacy_scan/scan.duckdb" >/dev/null <<'SQL'
DROP TABLE scan_motif_threshold;
SQL
legacy_batch_count=$(
    "$repository_root/scripts/manage_tp73_genome_context_maxima.py" prepare \
        --run-root "$legacy_run_root" --scan-package "$legacy_scan" \
        --evidence-package "$evidence" --threshold-registry "$thresholds" \
        --runtime-prefix "$runtime" --source "$repository_root" \
        --duckdb "$duckdb" --run-id synthetic_legacy_genome_context \
        --applied-threshold-set-id legacy-context-applied-v1 \
        --motifs-per-batch 2 --threads 1 --memory-limit 512MB \
        --max-temp-size 1GB --scratch-root "$temporary/legacy-scratch" \
        --minimum-free-run-gb 0 --minimum-free-scratch-gb 0
)
[[ $legacy_batch_count -eq 1 ]] || {
    echo "E: Legacy inventory-floor catalog was not accepted." >&2
    exit 1
}
python3 - "$legacy_run_root/plan/run_config.json" <<'PY'
import json
import pathlib
import sys

config = json.loads(pathlib.Path(sys.argv[1]).read_text())
assert config["scan_threshold_source"] == "scan_file_inventory_unique_floor"
assert config["maximum_source_score_floor"] == -1
PY

unsafe_scan="$temporary/unsafe-scan-package"
cp -R "$scan" "$unsafe_scan"
"$duckdb" -batch "$unsafe_scan/scan.duckdb" >/dev/null <<'SQL'
UPDATE scan_motif_threshold
SET final_minimum_score = 2.0
WHERE motif_id = 'M1';
UPDATE scan_file_inventory
SET minimum_score = 2.0
WHERE motif_id = 'M1';
SQL
if "$repository_root/scripts/manage_tp73_genome_context_maxima.py" prepare \
    --run-root "$unsafe_run_root" --scan-package "$unsafe_scan" \
    --evidence-package "$evidence" --threshold-registry "$thresholds" \
    --runtime-prefix "$runtime" --source "$repository_root" \
    --duckdb "$duckdb" --run-id unsafe_genome_context \
    --applied-threshold-set-id unsafe-context-applied-v1 \
    --motifs-per-batch 2 --threads 1 --memory-limit 512MB \
    --max-temp-size 1GB --scratch-root "$temporary/unsafe-scratch" \
    --minimum-free-run-gb 0 --minimum-free-scratch-gb 0 \
    >"$temporary/unsafe-prepare.out" 2>"$temporary/unsafe-prepare.err"; then
    echo "E: A scan source floor above -1 passed context preflight." >&2
    exit 1
fi
grep -Fq 'applied threshold registry is invalid' \
    "$temporary/unsafe-prepare.err"

batch_count=$(
    "$repository_root/scripts/manage_tp73_genome_context_maxima.py" prepare \
        --run-root "$run_root" --scan-package "$scan" \
        --evidence-package "$evidence" --threshold-registry "$thresholds" \
        --runtime-prefix "$runtime" --source "$repository_root" \
        --duckdb "$duckdb" --run-id synthetic_genome_context \
        --applied-threshold-set-id genome-context-applied-v1 \
        --motifs-per-batch 2 --threads 1 --memory-limit 512MB \
        --max-temp-size 1GB --scratch-root "$temporary/scratch" \
        --minimum-free-run-gb 0 --minimum-free-scratch-gb 0
)
[[ $batch_count -eq 1 ]] || {
    echo "E: Expected one context batch, found $batch_count." >&2
    exit 1
}

"$repository_root/scripts/run_tp73_genome_context_maxima_batch.py" \
    --run-root "$run_root" --batch-index 0
"$repository_root/scripts/run_tp73_genome_context_maxima_batch.py" \
    --run-root "$run_root" --batch-index 0

status=$(
    "$repository_root/scripts/manage_tp73_genome_context_maxima.py" status \
        --run-root "$run_root"
)
grep -Fq $'planned_motifs\t2' <<< "$status"
grep -Fq $'complete_motifs\t2' <<< "$status"

"$repository_root/scripts/manage_tp73_genome_context_maxima.py" finalize \
    --run-root "$run_root" --duckdb "$duckdb" --verify-checksums
"$repository_root/scripts/manage_tp73_genome_context_maxima.py" finalize \
    --run-root "$run_root" --duckdb "$duckdb"

m1="$run_root/tasks/task-000000-M1/cofactor_maxima.parquet"
m2="$run_root/tasks/task-000001-M2/cofactor_maxima.parquet"
final="$run_root/final/context_maxima"
"$duckdb" -batch :memory: >/dev/null <<SQL
SELECT CASE WHEN (SELECT count(*) FROM read_parquet('$m1')) <> 22
    THEN error('M1 autosomal context is incomplete') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM read_parquet('$m1')
    WHERE schema_version <> 2 OR context_score <> 2
       OR n_neighbor_loci_at_source_floor <> 1
       OR n_neighbor_loci_at_or_above_zero <> 1
       OR n_neighbor_loci_above_threshold <> 1
) THEN error('M1 nearby hit was not retained and counted') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM read_parquet('$m2')
    WHERE schema_version <> 2 OR context_score IS NOT NULL
       OR n_neighbor_loci_at_source_floor <> 0
       OR n_neighbor_loci_at_or_above_zero <> 0
       OR n_neighbor_loci_above_threshold <> 0
) THEN error('empty M2 partitions did not produce explicit absence') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM read_parquet(
      '$final/tables/jaspar2026/motif_score_threshold/part-000000.parquet'
    )
    WHERE motif_id = 'M2' AND recommended_threshold = 0
      AND source_recommended_threshold IS NULL
      AND calibration_status = 'fallback_zero_no_source_recommendation'
      AND used_fallback AND NOT raised_to_scan_floor
) THEN error('fallback threshold did not preserve the low source floor') END;
SQL

"$duckdb" -batch "$final/tp73_genome_context_maxima.duckdb" >/dev/null <<SQL
SELECT CASE WHEN (
    SELECT count(*) FROM tp73_motif_threshold_count_files(['$m1'])
) <> 22 THEN error('threshold/count file macro did not expose all autosomes') END;
SQL

"$repository_root/scripts/manage_tp73_cofactor_enrichment.py" prepare \
    --run-root "$enrichment_root" --source-threshold-run "$run_root" \
    --source "$repository_root" --duckdb "$duckdb" \
    --run-id synthetic_genome_enrichment >/dev/null
python3 - "$enrichment_root/plan/run_config.json" <<'PY'
import json
import pathlib
import sys

config = json.loads(pathlib.Path(sys.argv[1]).read_text())
assert config["analysis_scope"] == "all_non_target_jaspar2026_motifs_vs_tp73_autosomes"
assert "autosomes 1-22" in config["multiple_testing_scope_label"]
PY

echo "TP73 genome-context manager tests passed."

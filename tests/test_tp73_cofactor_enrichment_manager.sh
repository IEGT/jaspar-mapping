#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-enrichment-manager.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping enrichment-manager test." >&2
    exit 0
}
command -v Rscript >/dev/null 2>&1 || {
    echo "I: Rscript unavailable; skipping enrichment-manager test." >&2
    exit 0
}
Rscript -e 'library(data.table)' >/dev/null 2>&1 || {
    echo "I: data.table unavailable; skipping enrichment-manager test." >&2
    exit 0
}

source_run="$temporary/source-threshold-run"
run_root="$temporary/enrichment-run"
mkdir -p "$source_run/plan" "$source_run/input" \
    "$source_run/final/threshold_calibration/tables/jaspar2026/motif_score_threshold"
printf 'task_index\tmotif_id\tmotif_name\n0\tM_ENRICHED\tSynthetic enriched\n1\tM_DEPLETED\tSynthetic depleted\n' \
    > "$source_run/plan/calibration_tasks.tsv"

"$duckdb" -light-mode -batch :memory: >/dev/null <<SQL
COPY (
    SELECT
        '1'::VARCHAR AS chrom,
        (i * 1000)::BIGINT AS anchor_start,
        (i * 1000 + 16)::BIGINT AS anchor_end,
        (((i * 7) % 17) / 2.0)::FLOAT AS anchor_score,
        ((i % 8 IN (0, 1) AND floor(i / 8) % 5 <> 0)
          OR (i % 8 IN (4, 5) AND floor(i / 8) % 5 = 0)
          OR i % 8 = 6) AS supported_anti_s1,
        ((i % 8 IN (4, 5) AND floor(i / 8) % 5 <> 0)
          OR (i % 8 IN (0, 1) AND floor(i / 8) % 5 = 0)
          OR i % 8 = 2) AS supported_control_s1,
        CASE WHEN supported_anti_s1 THEN 1 + (i % 10) ELSE 0 END::FLOAT
            AS depth_anti_s1,
        CASE WHEN supported_control_s1 THEN 1 + (i % 7) ELSE 0 END::FLOAT
            AS depth_control_s1,
        ((i % 8 IN (0, 1) AND floor(i / 8) % 5 <> 1)
          OR (i % 8 IN (4, 5) AND floor(i / 8) % 5 = 1)
          OR i % 8 = 7) AS supported_anti_s2,
        ((i % 8 IN (4, 5) AND floor(i / 8) % 5 <> 1)
          OR (i % 8 IN (0, 1) AND floor(i / 8) % 5 = 1)
          OR i % 8 = 3) AS supported_control_s2,
        CASE WHEN supported_anti_s2 THEN 1 ELSE 0 END::FLOAT
            AS depth_anti_s2,
        CASE WHEN supported_control_s2 THEN 1 ELSE 0 END::FLOAT
            AS depth_control_s2
    FROM range(400) AS r(i)
) TO '$source_run/input/tp73_chr1_anchor_evidence.parquet'
  (FORMAT PARQUET);

COPY (
    SELECT '1'::VARCHAR AS chrom, (i * 1000)::BIGINT AS anchor_start,
           (i * 1000 + 16)::BIGINT AS anchor_end,
           'M_ENRICHED'::VARCHAR AS motif_id,
           CASE WHEN i % 8 IN (0, 1) THEN 5.0
                WHEN i % 8 = 2 THEN -0.5
                WHEN i % 8 = 3 THEN -1.0 ELSE NULL END::FLOAT AS context_score,
           -1.0::DOUBLE AS source_score_floor, 150::BIGINT AS context_flank_bp,
           'signed_interval_edge_distance'::VARCHAR AS context_distance_metric
    FROM range(400) AS r(i)
) TO '$temporary/enriched.parquet' (FORMAT PARQUET);
COPY (
    SELECT '1'::VARCHAR AS chrom, (i * 1000)::BIGINT AS anchor_start,
           (i * 1000 + 16)::BIGINT AS anchor_end,
           'M_DEPLETED'::VARCHAR AS motif_id,
           CASE WHEN i % 8 IN (4, 5) THEN 5.0
                WHEN i % 8 = 6 THEN -0.5
                WHEN i % 8 = 7 THEN -1.0 ELSE NULL END::FLOAT AS context_score,
           -1.0::DOUBLE AS source_score_floor, 150::BIGINT AS context_flank_bp,
           'signed_interval_edge_distance'::VARCHAR AS context_distance_metric
    FROM range(400) AS r(i)
) TO '$temporary/depleted.parquet' (FORMAT PARQUET);
COPY (
    SELECT * FROM (VALUES
      ('M_ENRICHED', 'Synthetic enriched', 4.0::DOUBLE, 'recommended'),
      ('M_DEPLETED', 'Synthetic depleted', NULL::DOUBLE, 'no_eligible_threshold')
    ) AS t(motif_id, motif_name, recommended_threshold, calibration_status)
) TO '$source_run/final/threshold_calibration/tables/jaspar2026/motif_score_threshold/part-000000.parquet'
  (FORMAT PARQUET);
SQL

python3 - "$source_run" "$temporary/enriched.parquet" \
    "$temporary/depleted.parquet" <<'PY'
import hashlib
import json
import pathlib
import shutil
import sys

root = pathlib.Path(sys.argv[1])
anchor = root / "input/tp73_chr1_anchor_evidence.parquet"
(root / "plan/run_config.json").write_text(json.dumps({
    "anchor_evidence": str(anchor.resolve()),
    "source_minimum_score": -1,
    "context_flank_bp": 150,
}, indent=2) + "\n")
for index, (motif, source) in enumerate((
    ("M_ENRICHED", pathlib.Path(sys.argv[2])),
    ("M_DEPLETED", pathlib.Path(sys.argv[3])),
)):
    directory = root / "tasks" / f"task-{index:06d}-{motif}"
    directory.mkdir(parents=True)
    target = directory / "cofactor_maxima.parquet"
    shutil.copyfile(source, target)
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    (directory / "complete.json").write_text(json.dumps({
        "schema_version": 1,
        "task_index": index,
        "motif_id": motif,
        "files": {"cofactor_maxima.parquet": {
            "bytes": target.stat().st_size, "sha256": digest,
        }},
    }, indent=2, sort_keys=True) + "\n")
PY

modern_source_run="$temporary/modern-context-run"
modern_package="$modern_source_run/final/context_maxima"
modern_run="$temporary/modern-enrichment-run"
mkdir -p "$modern_source_run/plan" \
    "$modern_package/tables/jaspar2026/motif_score_threshold"
python3 - "$source_run" "$modern_source_run" "$modern_package" <<'PY'
import csv
import hashlib
import json
import pathlib
import shutil
import sys

legacy, modern, package = map(pathlib.Path, sys.argv[1:])
anchor = legacy / "input/tp73_chr1_anchor_evidence.parquet"
(modern / "plan/run_config.json").write_text(json.dumps({
    "analysis_scope": "all_non_target_jaspar2026_motifs_vs_tp73_autosomes",
    "anchor_evidence": str(anchor.resolve()),
    "context_flank_bp": 150,
}, indent=2) + "\n")
rows = []
for index, (motif, name, threshold) in enumerate((
    ("M_ENRICHED", "Synthetic enriched", 6.0),
    ("M_DEPLETED", "Synthetic depleted", 3.0),
)):
    source = legacy / "tasks" / f"task-{index:06d}-{motif}"
    target = modern / "tasks" / source.name
    shutil.copytree(source, target)
    marker_path = target / "complete.json"
    marker = json.loads(marker_path.read_text())
    marker["scan_minimum_score"] = -1.0
    marker_path.write_text(json.dumps(marker, indent=2, sort_keys=True) + "\n")
    maxima = target / "cofactor_maxima.parquet"
    rows.append({
        "task_index": index,
        "motif_id": motif,
        "motif_name": name,
        "motif_length": 10,
        "analysis_partition": "autosome",
        "chromosome_count": 1,
        "anchor_count": 400,
        "scan_minimum_score": -1.0,
        "applied_context_threshold": threshold,
        "relative_path": str(maxima.relative_to(modern)),
        "absolute_path": str(maxima.resolve()),
        "bytes": maxima.stat().st_size,
        "sha256": hashlib.sha256(maxima.read_bytes()).hexdigest(),
    })
inventory = package / "context_maxima_file_inventory.tsv"
with inventory.open("w", newline="") as stream:
    writer = csv.DictWriter(
        stream, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
PY

"$duckdb" -light-mode -batch :memory: >/dev/null <<SQL
COPY (
  SELECT motif_id, motif_name, recommended_threshold,
         recommended_threshold AS source_recommended_threshold,
         'recommended'::VARCHAR AS calibration_status,
         'source_recommendation_reachable'::VARCHAR
           AS threshold_application_status,
         'synthetic_lowfloor_v2'::VARCHAR AS threshold_set_id
  FROM (VALUES
    ('M_ENRICHED', 'Synthetic enriched', 6.0::DOUBLE),
    ('M_DEPLETED', 'Synthetic depleted', 3.0::DOUBLE)
  ) AS t(motif_id, motif_name, recommended_threshold)
) TO '$modern_package/tables/jaspar2026/motif_score_threshold/part-000000.parquet'
  (FORMAT PARQUET);
SQL
"$duckdb" -light-mode -batch "$modern_package/tp73_genome_context_maxima.duckdb" \
    >/dev/null <<SQL
CREATE TABLE context_maxima_file_inventory AS SELECT * FROM read_csv_auto(
  '$modern_package/context_maxima_file_inventory.tsv', delim='\t', header=true);
CREATE TABLE motif_score_threshold AS SELECT * FROM read_parquet(
  '$modern_package/tables/jaspar2026/motif_score_threshold/part-000000.parquet');
SQL
printf '{"state":"complete","database":"tp73_genome_context_maxima.duckdb"}\n' \
    > "$modern_package/manifest.json"

modern_count=$(python3 \
    "$repository_root/scripts/manage_tp73_cofactor_enrichment.py" prepare \
    --run-root "$modern_run" --source-threshold-run "$modern_source_run" \
    --anchor-evidence "$source_run/input/tp73_chr1_anchor_evidence.parquet" \
    --source "$repository_root" --duckdb "$duckdb" \
    --run-id synthetic_modern_context_v2)
[[ $modern_count == 2 ]]
python3 - "$modern_run" "$modern_package" <<'PY'
import csv
import json
import pathlib
import sys

root, package = map(pathlib.Path, sys.argv[1:])
config = json.loads((root / "plan/run_config.json").read_text())
assert config["source_input_contract"] == "finalized_context_maxima_inventory"
assert config["source_context_maxima_package"] == str(package.resolve())
assert config["source_threshold_task_plan"].endswith(
    "context_maxima_file_inventory.tsv"
)
assert config["source_threshold_registry"].endswith(
    "motif_score_threshold/part-000000.parquet"
)
with (root / "plan/enrichment_tasks.tsv").open(newline="") as stream:
    rows = {row["motif_id"]: row for row in csv.DictReader(stream, delimiter="\t")}
assert rows["M_ENRICHED"]["positive_threshold"] == "6"
assert rows["M_DEPLETED"]["positive_threshold"] == "3"
assert rows["M_ENRICHED"]["source_task_relative_path"].startswith("tasks/")
PY

mkdir -p "$temporary/scratch-modern"
SLURM_ARRAY_TASK_ID=0 SLURM_TMPDIR="$temporary/scratch-modern" \
    bash "$repository_root/scripts/run_tp73_cofactor_enrichment_slurm_task.sh" \
    --run-root "$modern_run" --source-threshold-run "$modern_source_run" \
    --task-file "$modern_run/plan/enrichment_tasks.tsv" \
    --run-config "$modern_run/plan/run_config.json" \
    --source "$repository_root" --duckdb "$duckdb" --rscript Rscript \
    --block-size 50000 --spline-df 3

"$duckdb" -light-mode -batch :memory: >/dev/null <<SQL
COPY (
    SELECT * EXCLUDE (
        depth_anti_s1, depth_control_s1, depth_anti_s2, depth_control_s2
    ) FROM read_parquet('$source_run/input/tp73_chr1_anchor_evidence.parquet')
) TO '$temporary/anchors-without-depth.parquet' (FORMAT PARQUET);
SQL
if python3 "$repository_root/scripts/manage_tp73_cofactor_enrichment.py" \
    prepare --run-root "$temporary/invalid-run" \
    --source-threshold-run "$source_run" \
    --anchor-evidence "$temporary/anchors-without-depth.parquet" \
    --source "$repository_root" --duckdb "$duckdb" --run-id invalid_test_v1 \
    > "$temporary/invalid.out" 2> "$temporary/invalid.err"; then
    echo "E: manager accepted anchor evidence without depth" >&2
    exit 1
fi
grep -Fq 'lacks matched support/depth columns' "$temporary/invalid.err"

"$duckdb" -light-mode -batch :memory: >/dev/null <<SQL
COPY (
    SELECT chrom, anchor_start, anchor_end, anchor_score,
           supported_anti_s1 AS supported_tp73_s1,
           supported_control_s1 AS supported_negative_control_s1,
           depth_anti_s1 AS depth_tp73_s1,
           depth_control_s1 AS depth_negative_control_s1,
           supported_anti_s2 AS supported_tp73_s2,
           supported_control_s2 AS supported_negative_control_s2,
           depth_anti_s2 AS depth_tp73_s2,
           depth_control_s2 AS depth_negative_control_s2
    FROM read_parquet('$source_run/input/tp73_chr1_anchor_evidence.parquet')
) TO '$temporary/manifest-anchors.parquet' (FORMAT PARQUET);
COPY (
    SELECT * REPLACE (
        CASE WHEN motif_id = 'M_DEPLETED' THEN 3.0::DOUBLE
             ELSE recommended_threshold END AS recommended_threshold
    )
    FROM read_parquet(
      '$source_run/final/threshold_calibration/tables/jaspar2026/motif_score_threshold/part-000000.parquet'
    )
) TO '$temporary/fixed-registry.parquet' (FORMAT PARQUET);
SQL
external_count=$(python3 \
    "$repository_root/scripts/manage_tp73_cofactor_enrichment.py" prepare \
    --run-root "$temporary/external-run" --source-threshold-run "$source_run" \
    --anchor-evidence "$temporary/manifest-anchors.parquet" \
    --threshold-registry "$temporary/fixed-registry.parquet" \
    --source "$repository_root" --duckdb "$duckdb" \
    --run-id synthetic_external_v1)
[[ $external_count == 2 ]]
python3 - "$temporary/external-run" <<'PY'
import csv
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
config = json.loads((root / "plan/run_config.json").read_text())
assert config["evidence_column_scheme"] == "supported_tp73_and_negative_control"
assert config["sample_ids"] == ["s1", "s2"]
assert config["sample_count"] == 2
assert config["source_threshold_registry"].endswith("fixed-registry.parquet")
with (root / "plan/enrichment_tasks.tsv").open(newline="") as handle:
    rows = {row["motif_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
assert rows["M_DEPLETED"]["positive_threshold"] == "3"
PY

task_count=$(python3 "$repository_root/scripts/manage_tp73_cofactor_enrichment.py" \
    prepare --run-root "$run_root" --source-threshold-run "$source_run" \
    --source "$repository_root" --duckdb "$duckdb" --block-size 50000 \
    --spline-df 3 --run-id synthetic_enrichment_v1)
[[ $task_count == 2 ]]

real_rscript=$(command -v Rscript)
cat > "$temporary/slow-rscript" <<EOF
#!/usr/bin/env bash
sleep 2
exec "$real_rscript" "\$@"
EOF
chmod +x "$temporary/slow-rscript"
mkdir -p "$temporary/scratch-0"
progress_log="$temporary/task-0.err"
SLURM_ARRAY_TASK_ID=0 SLURM_TMPDIR="$temporary/scratch-0" \
    bash "$repository_root/scripts/run_tp73_cofactor_enrichment_slurm_task.sh" \
    --run-root "$run_root" --source-threshold-run "$source_run" \
    --task-file "$run_root/plan/enrichment_tasks.tsv" \
    --run-config "$run_root/plan/run_config.json" \
    --source "$repository_root" --duckdb "$duckdb" \
    --rscript "$temporary/slow-rscript" --block-size 50000 --spline-df 3 \
    2> "$progress_log" &
worker_pid=$!
for _ in $(seq 1 100); do
    grep -q 'Node-local work:' "$progress_log" 2>/dev/null && break
    kill -0 "$worker_pid" 2>/dev/null || {
        cat "$progress_log" >&2
        echo "E: enrichment worker exited before SIGUSR1 test" >&2
        exit 1
    }
    sleep 0.05
done
kill -USR1 "$worker_pid"
kill -0 "$worker_pid"
wait "$worker_pid"
[[ $(grep -c 'progress signal=SIGUSR1' "$progress_log") == 1 ]]

mkdir -p "$temporary/scratch-batch"
SLURM_ARRAY_TASK_ID=0 SLURM_TMPDIR="$temporary/scratch-batch" \
    bash "$repository_root/scripts/run_tp73_cofactor_enrichment_slurm_batch.sh" \
    --run-root "$run_root" --source-threshold-run "$source_run" \
    --task-file "$run_root/plan/enrichment_tasks.tsv" \
    --run-config "$run_root/plan/run_config.json" \
    --source "$repository_root" --duckdb "$duckdb" --rscript Rscript \
    --motifs-per-job 2 --block-size 50000 --spline-df 3

reuse=$(SLURM_ARRAY_TASK_ID=0 SLURM_TMPDIR="$temporary/scratch-0" \
    bash "$repository_root/scripts/run_tp73_cofactor_enrichment_slurm_task.sh" \
    --run-root "$run_root" --source-threshold-run "$source_run" \
    --task-file "$run_root/plan/enrichment_tasks.tsv" \
    --run-config "$run_root/plan/run_config.json" \
    --source "$repository_root" --duckdb "$duckdb" --rscript Rscript \
    --block-size 50000 --spline-df 3 2>&1)
[[ $reuse == *"Reusing completed enrichment task"* ]]

status=$(python3 "$repository_root/scripts/manage_tp73_cofactor_enrichment.py" \
    status --run-root "$run_root")
[[ $status == *$'complete\t2'* && $status == *$'pending\t0'* ]]

commit=$(git -C "$repository_root" rev-parse HEAD)
bash "$repository_root/scripts/run_tp73_cofactor_enrichment_finalize.sh" \
    --run-root "$run_root" --source "$repository_root" --duckdb "$duckdb" \
    --finalization-source-commit "$commit"

database="$run_root/final/cofactor_enrichment/tp73_cofactor_enrichment.duckdb"
"$duckdb" -light-mode -batch "$database" >/dev/null <<'SQL'
SELECT CASE WHEN (SELECT count(*) FROM cofactor_overview) <> 2
    THEN error('overview cardinality is wrong') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM cofactor_task
    WHERE motif_id = 'M_DEPLETED' AND recommended_threshold IS NULL
      AND positive_threshold = 0
      AND positive_threshold_role = 'descriptive_fallback'
) THEN error('null recommendation fallback was not explicit') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM cofactor_primary_occupancy
    WHERE requested_motifs_in_multiple_testing_scope <> 2
      OR multiple_testing_scope NOT LIKE 'all planned%'
) THEN error('all-motif multiple-testing scope is wrong') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM pragma_table_info('cofactor_primary_occupancy')
    WHERE name = 'q_value_bh_panel'
) THEN error('per-task q value leaked into final output') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM cofactor_primary_occupancy
    WHERE q_value_bh_all_jaspar IS NOT NULL
) THEN error('all-motif BH values were not calculated') END;
SQL

# A corrected publication identity can reuse the exact validated task outputs,
# but must be published beside rather than over the original final package.
bash "$repository_root/scripts/run_tp73_cofactor_enrichment_finalize.sh" \
    --run-root "$run_root" --source "$repository_root" --duckdb "$duckdb" \
    --finalization-source-commit "$commit" \
    --publication-run-id synthetic_enrichment_corrected_v2 \
    --final-name cofactor_enrichment_corrected
python3 - "$run_root" <<'PY'
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
original = json.loads(
    (root / "final/cofactor_enrichment/manifest.json").read_text()
)
corrected = json.loads(
    (root / "final/cofactor_enrichment_corrected/manifest.json").read_text()
)
assert original["run_id"] == original["plan_run_id"] == "synthetic_enrichment_v1"
assert not original["publication_identity_overrides_plan"]
assert corrected["run_id"] == "synthetic_enrichment_corrected_v2"
assert corrected["plan_run_id"] == "synthetic_enrichment_v1"
assert corrected["publication_identity_overrides_plan"]
assert corrected["row_counts"] == original["row_counts"]
PY

[[ ! -e "$run_root/tasks/task-000000-M_ENRICHED/cofactor_maxima.parquet" ]]
echo "I: TP73 cofactor enrichment manager test passed." >&2

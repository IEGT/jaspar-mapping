#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-enrichment-setup.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping enrichment-setup test." >&2
    exit 0
}

source_run="$temporary/source-run"
run_root="$temporary/enrichment-run"
base_anchor="$source_run/input/tp73_chr1_anchor_evidence.parquet"
mkdir -p "$source_run/input" "$source_run/plan" \
    "$source_run/final/threshold_calibration/tables/jaspar2026/motif_score_threshold" \
    "$source_run/tasks/task-000000-M_TEST" "$run_root" "$temporary/scratch"

"$duckdb" -light-mode -batch :memory: >/dev/null <<SQL
COPY (
  SELECT * FROM (VALUES
    (10::BIGINT, 16::BIGINT, 2.0::FLOAT),
    (20::BIGINT, 26::BIGINT, 3.0::FLOAT)
  ) AS v(start, "end", score)
) TO '$temporary/plus.parquet' (FORMAT PARQUET);
COPY (
  SELECT * FROM (VALUES
    (10::BIGINT, 16::BIGINT, 4.0::FLOAT)
  ) AS v(start, "end", score)
) TO '$temporary/minus.parquet' (FORMAT PARQUET);
SQL
printf '%s\n' $'chr1\t5\t17\t2' $'chr1\t19\t27\t5' \
    > "$temporary/coverage.bedGraph"

coverage_arguments=()
for sample in s1 s2 s3 s4 s5 s6; do
    coverage_arguments+=(
        --coverage "supported_anti_${sample}=$temporary/coverage.bedGraph"
        --coverage "supported_control_${sample}=$temporary/coverage.bedGraph"
    )
done
"$repository_root/scripts/build_tp73_anchor_evidence.py" \
    --anchor-plus "$temporary/plus.parquet" \
    --anchor-minus "$temporary/minus.parquet" \
    "${coverage_arguments[@]}" --output "$base_anchor" --chrom 1 \
    --minimum-anchor-score 0 --duckdb "$duckdb"

printf 'task_index\tmotif_id\tmotif_name\n0\tM_TEST\tSynthetic factor\n' \
    > "$source_run/plan/calibration_tasks.tsv"
python3 - "$source_run" <<'PY'
import hashlib
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
anchor = root / "input/tp73_chr1_anchor_evidence.parquet"
(root / "plan/run_config.json").write_text(json.dumps({
    "anchor_evidence": str(anchor.resolve()),
}, indent=2) + "\n")
PY

"$duckdb" -light-mode -batch :memory: >/dev/null <<SQL
COPY (
  SELECT '1'::VARCHAR AS chrom, anchor_start, anchor_end,
         'M_TEST'::VARCHAR AS motif_id, 2.0::FLOAT AS context_score,
         -1.0::DOUBLE AS source_score_floor, 150::BIGINT AS context_flank_bp,
         'signed_interval_edge_distance'::VARCHAR AS context_distance_metric
  FROM read_parquet('$base_anchor')
) TO '$source_run/tasks/task-000000-M_TEST/cofactor_maxima.parquet'
  (FORMAT PARQUET);
COPY (
  SELECT 'M_TEST'::VARCHAR AS motif_id, 'Synthetic factor'::VARCHAR AS motif_name,
         1.0::DOUBLE AS recommended_threshold, 'recommended'::VARCHAR AS calibration_status
) TO '$source_run/final/threshold_calibration/tables/jaspar2026/motif_score_threshold/part-000000.parquet'
  (FORMAT PARQUET);
SQL

python3 - "$source_run/tasks/task-000000-M_TEST" <<'PY'
import hashlib
import json
import pathlib
import sys

directory = pathlib.Path(sys.argv[1])
path = directory / "cofactor_maxima.parquet"
(directory / "complete.json").write_text(json.dumps({
    "schema_version": 1,
    "task_index": 0,
    "motif_id": "M_TEST",
    "files": {"cofactor_maxima.parquet": {
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }},
}, indent=2, sort_keys=True) + "\n")
PY

SLURM_TMPDIR="$temporary/scratch" \
    bash "$repository_root/scripts/run_tp73_cofactor_enrichment_setup.sh" \
    --run-root "$run_root" --source-threshold-run "$source_run" \
    --source "$repository_root" \
    --source-commit 0123456789abcdef0123456789abcdef01234567 \
    --duckdb "$duckdb" --threads 1 \
    --memory-limit 1GB --expected-anchor-count 2

depth_anchor="$run_root/input/depth_anchor/tp73_chr1_anchor_evidence.parquet"
"$duckdb" -light-mode -batch :memory: >/dev/null <<SQL
SELECT CASE WHEN (SELECT count(*) FROM read_parquet('$depth_anchor')) <> 2
    THEN error('setup anchor cardinality differs') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM read_parquet('$depth_anchor')
    WHERE NOT supported_anti_s1 OR depth_anti_s1 NOT IN (2, 5)
      OR supported_anti_s1 <> (depth_anti_s1 > 0)
) THEN error('setup effective depths are wrong') END;
SQL
[[ -s $run_root/input/depth_anchor/complete.json &&
   -s $run_root/plan/enrichment_tasks.tsv &&
   -s $run_root/plan/run_config.json ]]
[[ $(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["task_count"])' \
    "$run_root/plan/run_config.json") == 1 ]]

SLURM_TMPDIR="$temporary/scratch" \
    bash "$repository_root/scripts/run_tp73_cofactor_enrichment_setup.sh" \
    --run-root "$run_root" --source-threshold-run "$source_run" \
    --source "$repository_root" \
    --source-commit 0123456789abcdef0123456789abcdef01234567 \
    --duckdb "$duckdb" --threads 1 \
    --memory-limit 1GB --expected-anchor-count 2 >/dev/null
echo "I: TP73 cofactor enrichment setup test passed." >&2

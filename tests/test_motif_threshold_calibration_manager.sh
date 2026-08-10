#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-threshold-manager.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM
duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping threshold-manager test." >&2
    exit 0
}

git_source="$temporary/git-source"
mkdir "$git_source"
git -C "$git_source" init -q
git -C "$git_source" config user.name Test
git -C "$git_source" config user.email test@example.invalid
printf 'tracked\n' > "$git_source/tracked.txt"
git -C "$git_source" add tracked.txt
git -C "$git_source" -c commit.gpgsign=false commit -q -m initial
printf 'untracked\n' > "$git_source/untracked.txt"
python3 - "$repository_root/scripts/manage_motif_threshold_calibration.py" \
    "$git_source" false <<'PY'
import importlib.util
import pathlib
import sys

spec = importlib.util.spec_from_file_location("threshold_manager", sys.argv[1])
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)
_, dirty = module.git_identity(pathlib.Path(sys.argv[2]))
assert dirty is (sys.argv[3] == "true")
PY
printf 'changed\n' > "$git_source/tracked.txt"
python3 - "$repository_root/scripts/manage_motif_threshold_calibration.py" \
    "$git_source" true <<'PY'
import importlib.util
import pathlib
import sys

spec = importlib.util.spec_from_file_location("threshold_manager", sys.argv[1])
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)
_, dirty = module.git_identity(pathlib.Path(sys.argv[2]))
assert dirty is (sys.argv[3] == "true")
PY

package="$temporary/package"
mkdir -p "$package/task_data/task_id=0" "$package/task_data/task_id=1" \
    "$package/task_data/task_id=2"
for task in 0 1 2; do
    for strand in plus minus; do
        "$duckdb" -batch :memory: -c "
COPY (SELECT 1::BIGINT AS start, 4::BIGINT AS \"end\", 1.0::FLOAT AS score)
TO '$package/task_data/task_id=$task/$strand.parquet' (FORMAT PARQUET);" >/dev/null
    done
done

cat > "$temporary/jaspar.txt" <<'EOF'
>MA0861.2 TP73
A [ 1 0 ]
>MA0001.1 TEST_ONE
A [ 1 0 ]
>MA0002.1 TEST_TWO
A [ 1 0 ]
EOF

python3 - "$package" > "$temporary/inventory.tsv" <<'PY'
import hashlib
import pathlib
import sys

package = pathlib.Path(sys.argv[1])
motifs = [(0, "MA0861.2"), (1, "MA0001.1"), (2, "MA0002.1")]
print("motif_id\tstrand\ttask_id\toutput_relative_path\tbytes\tsha256\temitted_hits")
for task, motif in motifs:
    for strand, filename in (("+", "plus.parquet"), ("-", "minus.parquet")):
        path = package / "task_data" / f"task_id={task}" / filename
        print(motif, strand, task, filename, path.stat().st_size,
              hashlib.sha256(path.read_bytes()).hexdigest(), 1, sep="\t")
PY

mkdir -p "$package/tables/jaspar2026/scan_file_inventory"
pushd "$package" >/dev/null
"$duckdb" jaspar_genome_scan.duckdb -batch <<SQL >/dev/null
COPY (
    SELECT *, '1'::VARCHAR AS chrom
    FROM read_csv_auto('$temporary/inventory.tsv', delim='\t', header=true)
) TO '$package/tables/jaspar2026/scan_file_inventory/part-000000.parquet'
  (FORMAT PARQUET);
CREATE VIEW scan_file_inventory AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/scan_file_inventory/**/*.parquet'
);
SQL
popd >/dev/null
printf '%s\n' '{"database":"jaspar_genome_scan.duckdb"}' > "$package/manifest.json"

run="$temporary/run"
mkdir -p "$run/input"
"$duckdb" -batch :memory: -c "
COPY (
    SELECT * FROM (VALUES
        ('1', 10::BIGINT, 26::BIGINT, -0.75::DOUBLE),
        ('1', 40::BIGINT, 56::BIGINT, 2.5::DOUBLE)
    ) AS anchors(chrom, anchor_start, anchor_end, anchor_score)
) TO '$run/input/anchors.parquet' (FORMAT PARQUET);" >/dev/null
cat > "$run/input/anchors.parquet.run_config.tsv" <<'EOF'
schema_version	chrom	anchor_motif_id	anchor_source_mode	anchor_score_rule	minimum_anchor_score	strict_immersion_rule	effective_depth_rule	anchor_count	observed_minimum_score	observed_maximum_score
3	1	MA0861.2	schema7_local_peak_context_anchor	maximum_score_across_preselected_orientation_records_at_same_span	-1	strict	maximum	2	-0.75	2.5
EOF
task_count=$(python3 "$repository_root/scripts/manage_motif_threshold_calibration.py" prepare \
    --run-root "$run" --scan-package "$package" --jaspar "$temporary/jaspar.txt" \
    --anchor-evidence "$run/input/anchors.parquet" --source "$repository_root" \
    --run-id synthetic_schema7_local_peak_v1 \
    --threshold-set-id synthetic_schema7_local_peak_thresholds_v1 \
    --duckdb "$duckdb")
[[ $task_count == 2 && $(wc -l < "$run/plan/calibration_tasks.tsv") -eq 3 ]] || {
    echo "E: Calibration manager prepared the wrong task count." >&2
    exit 1
}
grep -q '^MA0001.1$' "$run/plan/motifs.txt"
grep -q '^MA0002.1$' "$run/plan/motifs.txt"
! grep -q '^MA0861.2$' "$run/plan/motifs.txt"
python3 - "$run/plan/run_config.json" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    config = json.load(handle)
assert config["schema_version"] == 2
assert config["run_id"] == "synthetic_schema7_local_peak_v1"
assert config["threshold_set_id"] == "synthetic_schema7_local_peak_thresholds_v1"
assert config["anchor_provenance"]["anchor_source_mode"] == \
    "schema7_local_peak_context_anchor"
assert config["anchor_provenance"]["anchor_minimum_score"] == -1
assert config["anchor_provenance"]["anchor_count"] == 2
PY

python3 - "$run" <<'PY'
import csv
import json
import pathlib
import sys

run = pathlib.Path(sys.argv[1])
with (run / "plan" / "calibration_tasks.tsv").open(newline="") as handle:
    tasks = list(csv.DictReader(handle, delimiter="\t"))

fields = [
    "motif_id", "threshold", "anchors_total", "anchors_retained",
    "retained_fraction", "discordant_observations",
    "baseline_macro_roc_auc", "augmented_macro_roc_auc",
    "delta_macro_roc_auc", "baseline_macro_average_precision",
    "augmented_macro_average_precision", "delta_macro_average_precision",
    "baseline_macro_log_loss", "augmented_macro_log_loss",
    "delta_macro_log_loss", "baseline_macro_brier_score",
    "augmented_macro_brier_score", "delta_macro_brier_score",
    "median_adjusted_odds_ratio", "minimum_adjusted_odds_ratio",
    "maximum_adjusted_odds_ratio", "median_raw_sample_odds_ratio",
    "samples_with_raw_odds_ratio_below_one", "samples_total",
    "sample_fold_cells", "sample_fold_cells_with_roc_auc_gain",
    "sample_fold_cells_with_log_loss_gain", "evaluation_status",
    "evaluation_note",
]
for task in tasks:
    index = int(task["task_index"])
    motif = task["motif_id"]
    directory = run / "tasks" / f"task-{index:06d}-{motif}"
    directory.mkdir(parents=True)
    row = dict.fromkeys(fields, "")
    row.update({
        "motif_id": motif, "threshold": 0, "anchors_total": 100,
        "anchors_retained": 50, "retained_fraction": 0.5,
        "discordant_observations": 120, "baseline_macro_roc_auc": 0.6,
        "baseline_macro_average_precision": 0.6,
        "baseline_macro_log_loss": 0.7, "baseline_macro_brier_score": 0.2,
        "samples_with_raw_odds_ratio_below_one": 0, "samples_total": 6,
        "sample_fold_cells": 30, "sample_fold_cells_with_roc_auc_gain": 0,
        "sample_fold_cells_with_log_loss_gain": 0,
    })
    if index == 0:
        row.update({
            "augmented_macro_roc_auc": 0.7, "delta_macro_roc_auc": 0.1,
            "augmented_macro_average_precision": 0.7,
            "delta_macro_average_precision": 0.1,
            "augmented_macro_log_loss": 0.6, "delta_macro_log_loss": -0.1,
            "augmented_macro_brier_score": 0.18,
            "delta_macro_brier_score": -0.02,
            "median_adjusted_odds_ratio": 1.2,
            "minimum_adjusted_odds_ratio": 1.1,
            "maximum_adjusted_odds_ratio": 1.3,
            "median_raw_sample_odds_ratio": 1.2,
            "sample_fold_cells_with_roc_auc_gain": 25,
            "sample_fold_cells_with_log_loss_gain": 25,
            "evaluation_status": "evaluated",
        })
    else:
        row.update({
            "evaluation_status": "fit_failed",
            "evaluation_note": "no finite held-out metric",
        })
    with (directory / "threshold_metrics.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(row)
    (directory / "complete.json").write_text(
        json.dumps({"motif_id": motif}) + "\n", encoding="utf-8"
    )
PY

python3 "$repository_root/scripts/manage_motif_threshold_calibration.py" status \
    --run-root "$run" | grep -q $'^complete\t2$'
finalization_commit=$(git -C "$repository_root" rev-parse HEAD)
python_executable=$(command -v python3)
duckdb_executable=$(command -v "$duckdb")
mkdir "$temporary/no-git-path"
PATH="$temporary/no-git-path" \
    "$python_executable" \
    "$repository_root/scripts/manage_motif_threshold_calibration.py" finalize \
    --run-root "$run" --duckdb "$duckdb_executable" \
    --finalization-source-commit "$finalization_commit"

registry="$run/final/threshold_calibration/tables/jaspar2026/motif_score_threshold/part-000000.parquet"
"$duckdb" -batch :memory: >/dev/null <<SQL
SELECT CASE WHEN (SELECT count(*) FROM read_parquet('$registry')) <> 2
    THEN error('final threshold registry has the wrong motif count') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM read_parquet('$registry')
    WHERE motif_id = 'MA0001.1'
      AND calibration_status = 'exploratory_positive_gain'
) THEN error('synthetic positive threshold was not selected') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM read_parquet('$registry')
    WHERE motif_id = 'MA0002.1'
      AND recommended_threshold IS NULL
      AND selected_metric_gain IS NULL
      AND calibration_status = 'no_finite_metric'
) THEN error('blank non-evaluable metrics did not retain null semantics') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM read_parquet('$registry')
    WHERE calibration_stratum_id <> 'schema7_local_peak_tp73_anchors'
       OR CAST(json_extract(calibration_stratum,
              '$.anchor_minimum_score') AS DOUBLE) <> -1
       OR json_extract_string(calibration_stratum,
              '$.anchor_source_mode') <> 'schema7_local_peak_context_anchor'
) THEN error('schema-7 local-peak calibration stratum was not retained') END;
SQL

python3 - "$run/final/threshold_calibration/manifest.json" \
    "$finalization_commit" <<'PY'
import json
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    manifest = json.load(handle)
assert manifest["schema_version"] == 3
assert manifest["source_commit"] == manifest["metric_source_commit"]
assert manifest["finalization_source_commit"] == sys.argv[2]
assert isinstance(manifest["metric_source_dirty"], bool)
assert isinstance(manifest["finalization_source_dirty"], bool)
assert manifest["calibration_stratum_id"] == \
    "schema7_local_peak_tp73_anchors"
assert manifest["calibration_stratum"]["anchor_minimum_score"] == -1
assert manifest["calibration_stratum"]["anchor_count"] == 2
assert len(manifest["anchor_evidence_sha256"]) == 64
assert len(manifest["anchor_run_config_sha256"]) == 64
PY

echo "Motif threshold-calibration manager tests passed."

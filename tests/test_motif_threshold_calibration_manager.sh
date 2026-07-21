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
task_count=$(python3 "$repository_root/scripts/manage_motif_threshold_calibration.py" prepare \
    --run-root "$run" --scan-package "$package" --jaspar "$temporary/jaspar.txt" \
    --anchor-evidence "$run/input/anchors.parquet" --source "$repository_root" \
    --duckdb "$duckdb")
[[ $task_count == 2 && $(wc -l < "$run/plan/calibration_tasks.tsv") -eq 3 ]] || {
    echo "E: Calibration manager prepared the wrong task count." >&2
    exit 1
}
grep -q '^MA0001.1$' "$run/plan/motifs.txt"
grep -q '^MA0002.1$' "$run/plan/motifs.txt"
! grep -q '^MA0861.2$' "$run/plan/motifs.txt"

header=$'motif_id\tthreshold\tanchors_total\tanchors_retained\tretained_fraction\tdiscordant_observations\tbaseline_macro_roc_auc\taugmented_macro_roc_auc\tdelta_macro_roc_auc\tbaseline_macro_average_precision\taugmented_macro_average_precision\tdelta_macro_average_precision\tbaseline_macro_log_loss\taugmented_macro_log_loss\tdelta_macro_log_loss\tbaseline_macro_brier_score\taugmented_macro_brier_score\tdelta_macro_brier_score\tmedian_adjusted_odds_ratio\tminimum_adjusted_odds_ratio\tmaximum_adjusted_odds_ratio\tmedian_raw_sample_odds_ratio\tsamples_with_raw_odds_ratio_below_one\tsamples_total\tsample_fold_cells\tsample_fold_cells_with_roc_auc_gain\tsample_fold_cells_with_log_loss_gain\tevaluation_status\tevaluation_note'
for task in 0 1; do
    motif=$(awk -F '\t' -v task="$task" 'NR > 1 && $1 == task {print $2}' \
        "$run/plan/calibration_tasks.tsv")
    directory="$run/tasks/task-$(printf '%06d' "$task")-$motif"
    mkdir -p "$directory"
    printf '%s\n' "$header" > "$directory/threshold_metrics.tsv"
    printf '%s\n' "$motif"$'\t0\t100\t50\t0.5\t120\t0.6\t0.7\t0.1\t0.6\t0.7\t0.1\t0.7\t0.6\t-0.1\t0.2\t0.18\t-0.02\t1.2\t1.1\t1.3\t1.2\t0\t6\t30\t25\t25\tevaluated\tNA' \
        >> "$directory/threshold_metrics.tsv"
    printf '{"motif_id":"%s"}\n' "$motif" > "$directory/complete.json"
done

python3 "$repository_root/scripts/manage_motif_threshold_calibration.py" status \
    --run-root "$run" | grep -q $'^complete\t2$'
python3 "$repository_root/scripts/manage_motif_threshold_calibration.py" finalize \
    --run-root "$run" --duckdb "$duckdb"

registry="$run/final/threshold_calibration/tables/jaspar2026/motif_score_threshold/part-000000.parquet"
"$duckdb" -batch :memory: >/dev/null <<SQL
SELECT CASE WHEN (SELECT count(*) FROM read_parquet('$registry')) <> 2
    THEN error('final threshold registry has the wrong motif count') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM read_parquet('$registry')
    WHERE calibration_status <> 'exploratory_positive_gain'
) THEN error('synthetic positive threshold was not selected') END;
SQL

echo "Motif threshold-calibration manager tests passed."

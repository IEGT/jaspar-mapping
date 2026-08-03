#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-negative-sensitivity.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

scanner="$repository_root/pssm_scan_parquet"
duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "E: duckdb CLI is required." >&2
    exit 1
}
command -v Rscript >/dev/null 2>&1 || {
    echo "I: Rscript unavailable; skipping negative-sensitivity test." >&2
    exit 0
}
Rscript -e 'library(data.table)' >/dev/null 2>&1 || {
    echo "I: R data.table unavailable; skipping negative-sensitivity test." >&2
    exit 0
}
[[ -x $scanner ]] || {
    echo "E: Build the Arrow scanner with 'make pssm_scan_parquet'." >&2
    exit 1
}

genome="$temporary/genome.fa"
awk 'BEGIN {
    print ">1 synthetic"
    line = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    for (i = 0; i < 1100; ++i) print line
}' > "$genome"
python3 "$repository_root/scripts/build_fasta_index.py" "$genome" \
    --output "$genome.fai" >/dev/null

cat > "$temporary/motif.jaspar" <<'EOF'
>TEST.1	SENSITIVITY_TEST
A [ 1 ]
C [ 1 ]
G [ 1 ]
T [ 1 ]
EOF
printf 'motif_id\trecommended_threshold\nTEST.1\t0\n' \
    > "$temporary/thresholds.tsv"

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT
        '1'::VARCHAR AS chrom,
        (i * 1000)::BIGINT AS anchor_start,
        (i * 1000 + 16)::BIGINT AS anchor_end,
        ((i % 7) / 2.0)::FLOAT AS anchor_score,
        (i % 2 = 0) AS supported_anti_saos2_TA,
        (i % 2 = 0) AS supported_anti_saos2_DN,
        (i % 2 = 0) AS supported_anti_skmel29_1_TA,
        (i % 2 = 0) AS supported_anti_skmel29_1_DN,
        (i % 2 = 0) AS supported_anti_skmel29_2_TA,
        (i % 2 = 0) AS supported_anti_skmel29_2_DN,
        (i % 2 = 1) AS supported_control_saos2_TA,
        (i % 2 = 1) AS supported_control_saos2_DN,
        (i % 2 = 1) AS supported_control_skmel29_1_TA,
        (i % 2 = 1) AS supported_control_skmel29_1_DN,
        (i % 2 = 1) AS supported_control_skmel29_2_TA,
        (i % 2 = 1) AS supported_control_skmel29_2_DN
    FROM range(100) AS r(i)
) TO '$temporary/anchors.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL

run="$temporary/run"
task_count=$(
    "$repository_root/scripts/manage_motif_threshold_calibration.py" \
        prepare-negative-sensitivity --run-root "$run" \
        --threshold-list "$temporary/thresholds.tsv" \
        --jaspar "$temporary/motif.jaspar" \
        --anchor-evidence "$temporary/anchors.parquet" \
        --genome "$genome" --fasta-index "$genome.fai" \
        --source "$repository_root" --duckdb "$duckdb" \
        --chrom-size 110000 --source-minimum-score -20
)
[[ $task_count -eq 1 ]] || {
    echo "E: Synthetic sensitivity plan has the wrong task count." >&2
    exit 1
}

SLURM_ARRAY_TASK_ID=0 \
    "$repository_root/scripts/run_negative_threshold_sensitivity_slurm_task.sh" \
    --run-root "$run" --scanner "$scanner" \
    --duckdb "$(command -v "$duckdb")" --rscript "$(command -v Rscript)" \
    --source "$repository_root" --scratch-root "$temporary/scratch" \
    --minimum-free-gib 0 >/dev/null

task="$run/tasks/task-000000-TEST.1"
[[ -s $task/complete.json && -s $task/cofactor_maxima.parquet &&
   -s $task/threshold_metrics.tsv ]] || {
    echo "E: Synthetic sensitivity task was not promoted." >&2
    exit 1
}

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW maxima AS SELECT * FROM read_parquet('$task/cofactor_maxima.parquet');
CREATE VIEW metric AS SELECT * FROM read_csv_auto(
    '$task/threshold_metrics.tsv', delim='\t', header=true, nullstr='NA');
SELECT CASE WHEN (SELECT count(*) FROM maxima) <> 100
    THEN error('sensitivity maxima lost anchors') END;
SELECT CASE WHEN (SELECT count(*) FROM metric) <> 21
                  OR (SELECT min(threshold) FROM metric) <> -20
                  OR (SELECT max(threshold) FROM metric) <> 0
    THEN error('sensitivity metric grid is incomplete') END;
SQL

# A repeated/requeued task must validate and reuse the immutable completion.
SLURM_ARRAY_TASK_ID=0 \
    "$repository_root/scripts/run_negative_threshold_sensitivity_slurm_task.sh" \
    --run-root "$run" --scanner "$scanner" \
    --duckdb "$(command -v "$duckdb")" --rscript "$(command -v Rscript)" \
    --source "$repository_root" --scratch-root "$temporary/scratch" \
    --minimum-free-gib 0 2> "$temporary/reuse.log"
grep -Fq 'Reusing completed sensitivity task' "$temporary/reuse.log"

source_commit=$(git -C "$repository_root" rev-parse HEAD)
dirty=()
if ! git -C "$repository_root" diff --quiet --ignore-submodules -- ||
   ! git -C "$repository_root" diff --cached --quiet --ignore-submodules --; then
    dirty=(--finalization-source-dirty)
fi
"$repository_root/scripts/manage_motif_threshold_calibration.py" \
    finalize-negative-sensitivity --run-root "$run" --duckdb "$duckdb" \
    --finalization-source-commit "$source_commit" "${dirty[@]}" >/dev/null

summary="$run/final/negative_threshold_sensitivity/negative_threshold_sensitivity.parquet"
[[ -s $summary ]] || {
    echo "E: Sensitivity finalizer did not publish its summary." >&2
    exit 1
}
"$duckdb" -batch :memory: >/dev/null <<SQL
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM read_parquet('$summary')
    WHERE motif_id = 'TEST.1'
      AND candidate_count = 21
      AND tested_minimum_threshold = -20
      AND tested_maximum_threshold = 0
      AND sensitivity_conclusion = 'zero_not_evaluable'
) THEN error('sensitivity summary is inconsistent') END;
SQL

echo "Negative-threshold sensitivity tests passed."

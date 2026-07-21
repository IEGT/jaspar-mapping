#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-slurm-offset.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

mkdir -p "$temporary/bin"
cat > "$temporary/bin/python3" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$@" > "$CAPTURE_ARGUMENTS"
EOF
chmod +x "$temporary/bin/python3"

source_tree="$temporary/source"
mkdir -p "$source_tree/scripts" "$temporary/slurm-spool"
source_tree=$(cd "$source_tree" && pwd)
touch "$source_tree/pssm_scan_parquet" "$source_tree/scripts/manage_genome_scan.py"
cp "$repository_root/scripts/run_genome_scan_slurm_task.sh" \
    "$temporary/slurm-spool/submitted-script"
chmod +x "$temporary/slurm-spool/submitted-script"

capture="$temporary/arguments.txt"
CAPTURE_ARGUMENTS="$capture" \
PATH="$temporary/bin:/usr/bin:/bin" \
SLURM_ARRAY_TASK_ID=7 \
"$temporary/slurm-spool/submitted-script" \
    --run-root /data/example/run --scanner "$source_tree/pssm_scan_parquet" \
    --duckdb /data/example/duckdb --task-offset 100

manager=$(head -n 1 "$capture")
[[ $manager == "$source_tree/scripts/manage_genome_scan.py" ]] || {
    echo "E: Spool copy resolved the manager as '$manager'." >&2
    exit 1
}

task_index=$(awk 'previous == "--task-index" { print; exit } { previous = $0 }' "$capture")
[[ $task_index == 107 ]] || {
    echo "E: Expected global task index 107, got '$task_index'." >&2
    exit 1
}

if SLURM_ARRAY_TASK_ID=0 \
    "$temporary/slurm-spool/submitted-script" \
        --run-root /data/example/run --scanner "$source_tree/pssm_scan_parquet" \
        --task-offset -1 >/dev/null 2>&1; then
    echo "E: A negative Slurm task offset was accepted." >&2
    exit 1
fi

# Threshold calibration uses the same offset scheme when a motif set is larger
# than Slurm's maximum local array index. A completed synthetic task lets this
# exercise index resolution without invoking DuckDB or R.
mkdir -p "$temporary/threshold-bin"
cat > "$temporary/threshold-bin/stat" <<'EOF'
#!/usr/bin/env bash
[[ $1 == -c && $2 == %s && $# -eq 3 ]]
bytes=$(wc -c < "$3")
printf '%d\n' "$bytes"
EOF
chmod +x "$temporary/threshold-bin/stat"
threshold_run="$temporary/threshold-run"
threshold_scan="$temporary/threshold-scan"
mkdir -p "$threshold_run/plan" "$threshold_run/input" \
    "$threshold_scan/task_data/task_id=example"
touch "$threshold_run/input/anchors.parquet"
printf 'plus\n' > "$threshold_scan/task_data/task_id=example/plus.parquet"
printf 'minus\n' > "$threshold_scan/task_data/task_id=example/minus.parquet"
plus_bytes=$(( $(wc -c < "$threshold_scan/task_data/task_id=example/plus.parquet") ))
minus_bytes=$(( $(wc -c < "$threshold_scan/task_data/task_id=example/minus.parquet") ))
printf 'task_index\tmotif_id\tmotif_name\tplus_relative_path\tminus_relative_path\tplus_bytes\tminus_bytes\tplus_sha256\tminus_sha256\tplus_emitted_hits\tminus_emitted_hits\n' \
    > "$threshold_run/plan/calibration_tasks.tsv"
printf '1007\tMTEST.1\tTEST\ttask_data/task_id=example/plus.parquet\ttask_data/task_id=example/minus.parquet\t%s\t%s\tplus-sha\tminus-sha\t1\t1\n' \
    "$plus_bytes" "$minus_bytes" >> "$threshold_run/plan/calibration_tasks.tsv"
threshold_final="$threshold_run/tasks/task-001007-MTEST.1"
mkdir -p "$threshold_final"
printf '{"motif_id":"MTEST.1"}\n' > "$threshold_final/complete.json"
printf 'header\n' > "$threshold_final/threshold_metrics.tsv"
threshold_output=$(
    PATH="$temporary/threshold-bin:/usr/bin:/bin" SLURM_ARRAY_TASK_ID=7 \
    "$repository_root/scripts/run_motif_threshold_calibration_slurm_task.sh" \
        --run-root "$threshold_run" --scan-package "$threshold_scan" \
        --task-file "$threshold_run/plan/calibration_tasks.tsv" \
        --anchor-evidence "$threshold_run/input/anchors.parquet" \
        --source "$repository_root" --duckdb /usr/bin/true \
        --rscript /usr/bin/true --task-offset 1000 2>&1
)
grep -Fq 'Reusing completed threshold task 1007 (MTEST.1).' \
    <<< "$threshold_output" || {
    echo "E: Threshold task offset did not resolve global task 1007." >&2
    exit 1
}

mkdir -p "$temporary/prepared-run/plan"
cat > "$temporary/prepared-run/plan/scan_plan.json" <<'EOF'
{
  "sequence_regions": [
    {"chrom": "1", "included_in_scan": true},
    {"chrom": "2", "included_in_scan": true},
    {"chrom": "MT", "included_in_scan": false}
  ]
}
EOF
submission=$(
    "$repository_root/scripts/submit_genome_scan_slurm.sh" \
        --run-root "$temporary/prepared-run" \
        --scanner "$repository_root/pssm_scan" \
        --source "$repository_root" --max-concurrent 2 --dry-run
)
grep -Fq -- '--array=0-1%2' <<< "$submission" || {
    echo "E: Chromosome submission did not derive its array from the plan." >&2
    exit 1
}
grep -Fq 'run_genome_scan_slurm_chromosome.sh' <<< "$submission" || {
    echo "E: Chromosome submission omitted the scratch-aware wrapper." >&2
    exit 1
}

# A manually submitted finalizer also runs from a Slurm spool copy. With the
# source tree as its submitted working directory, it must still find the real
# coordinator rather than looking beside the spool file.
cp "$repository_root/scripts/run_genome_scan_slurm_finalize.sh" \
    "$temporary/slurm-spool/finalizer-script"
chmod +x "$temporary/slurm-spool/finalizer-script"
(
    cd "$source_tree"
    CAPTURE_ARGUMENTS="$capture" PATH="$temporary/bin:/usr/bin:/bin" \
        "$temporary/slurm-spool/finalizer-script" \
        --run-root /data/example/run --duckdb /data/example/duckdb
)
manager=$(head -n 1 "$capture")
[[ $manager == "$source_tree/scripts/manage_genome_scan.py" ]] || {
    echo "E: Spool-copy finalizer resolved the manager as '$manager'." >&2
    exit 1
}

echo "Slurm task-offset tests passed."

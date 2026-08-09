#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_tp73_cofactor_enrichment_setup.sh --run-root DIR
       --source-threshold-run DIR [options]

Build the shared depth-bearing TP73/CUT&RUN anchor table and then prepare the
immutable all-motif enrichment plan. Exact TP73 sparse inputs and all 12
bedGraphs are recovered from the completed threshold run's provenance. Inputs
are copied once to node-local scratch and verified there before use.

Options:
  --run-root DIR             New durable enrichment run
  --source-threshold-run DIR Completed threshold-calibration run
  --source DIR               Repository root (default: script parent)
  --duckdb FILE              DuckDB CLI (default: duckdb)
  --threads N                DuckDB threads (default: 2)
  --memory-limit SIZE        DuckDB memory ceiling (default: 12GB)
  --expected-anchor-count N  Required anchor rows (default: 310782)
  -h, --help                 Show this help

SIGUSR1 prints one phase/input/elapsed-time line without interrupting the
active copy, DuckDB build, or plan-preparation child.
EOF
}

run_root=""
source_threshold_run=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
duckdb=duckdb
threads=2
memory_limit=12GB
expected_anchor_count=310782
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --source-threshold-run) source_threshold_run=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --threads) threads=${2:?}; shift 2 ;;
        --memory-limit) memory_limit=${2:?}; shift 2 ;;
        --expected-anchor-count) expected_anchor_count=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $source_threshold_run ]] || { usage >&2; exit 2; }
[[ $threads =~ ^[1-9][0-9]*$ &&
   $expected_anchor_count =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --threads and --expected-anchor-count must be positive integers." >&2
    exit 2
}
for path in "$run_root" "$source_threshold_run" "$source"; do
    [[ -d $path ]] || { echo "E: Directory not found: $path" >&2; exit 1; }
done
[[ -x $duckdb ]] || duckdb=$(command -v "$duckdb" || true)
[[ -x $duckdb ]] || { echo "E: DuckDB executable is unavailable." >&2; exit 1; }

sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

anchor_package="$run_root/input/depth_anchor"
anchor="$anchor_package/tp73_chr1_anchor_evidence.parquet"
if [[ -e $anchor_package ]]; then
    [[ -s $anchor && -s $anchor_package/complete.json ]] || {
        echo "E: Existing depth-anchor package is incomplete: $anchor_package" >&2
        exit 1
    }
    recorded=$(python3 -c \
        'import json,sys; print(json.load(open(sys.argv[1]))["anchor_sha256"])' \
        "$anchor_package/complete.json")
    [[ $(sha256_file "$anchor") == "$recorded" ]] || {
        echo "E: Existing depth-anchor checksum differs." >&2
        exit 1
    }
    echo "I: Reusing completed depth-bearing TP73 anchors: $anchor" >&2
    python3 "$source/scripts/manage_tp73_cofactor_enrichment.py" prepare \
        --run-root "$run_root" --source-threshold-run "$source_threshold_run" \
        --anchor-evidence "$anchor" --source "$source" --duckdb "$duckdb"
    exit 0
fi

base_anchor="$source_threshold_run/input/tp73_chr1_anchor_evidence.parquet"
base_config="$base_anchor.run_config.tsv"
coverage_manifest="$base_anchor.coverage_manifest.tsv"
for path in "$base_anchor" "$base_config" "$coverage_manifest"; do
    [[ -f $path ]] || { echo "E: Source anchor provenance is missing: $path" >&2; exit 1; }
done

start_epoch=$(date +%s)
phase=initializing
current_input=none
child_pid=""
progress_report() {
    printf 'I: progress signal=SIGUSR1 phase=%s input=%s elapsed_s=%s child_pid=%s\n' \
        "$phase" "$current_input" "$(( $(date +%s) - start_epoch ))" \
        "${child_pid:-none}" >&2
}
trap progress_report USR1

run_child() {
    "$@" &
    child_pid=$!
    local status usr1_wait_status
    usr1_wait_status=$((128 + $(kill -l USR1)))
    while true; do
        set +e
        wait "$child_pid"
        status=$?
        set -e
        if [[ $status -eq $usr1_wait_status ]]; then
            continue
        fi
        if kill -0 "$child_pid" 2>/dev/null; then
            continue
        fi
        child_pid=""
        return "$status"
    done
}

scratch_base=${SLURM_TMPDIR:-/scratch/${USER:-sm718}}
scratch="$scratch_base/jaspar-enrichment-setup-${SLURM_JOB_ID:-manual}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
mkdir -p "$scratch/inputs"
staging_plan="$scratch/staging_inputs.tsv"
python3 - "$base_config" "$coverage_manifest" "$staging_plan" <<'PY'
import csv
import pathlib
import sys

with open(sys.argv[1], newline="") as handle:
    config = next(csv.DictReader(handle, delimiter="\t"))
rows = [
    ("anchor_plus", "anchor_plus", config["anchor_plus"], config["anchor_plus_sha256"]),
    ("anchor_minus", "anchor_minus", config["anchor_minus"], config["anchor_minus_sha256"]),
]
with open(sys.argv[2], newline="") as handle:
    coverage = list(csv.DictReader(handle, delimiter="\t"))
if len(coverage) != 12:
    raise SystemExit(f"expected 12 coverage tracks, found {len(coverage)}")
for row in coverage:
    rows.append(("coverage", row["support_column"], row["path"], row["sha256"]))
with open(sys.argv[3], "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(("role", "support_column", "source", "sha256", "staged"))
    for index, (role, support, source, digest) in enumerate(rows):
        suffix = pathlib.Path(source).suffix or ".input"
        writer.writerow((role, support, source, digest,
                         f"input-{index:02d}{suffix}"))
PY

input_bytes=0
while IFS=$'\t' read -r role _ input _ _; do
    [[ $role != role ]] || continue
    [[ -f $input ]] || { echo "E: Planned setup input is missing: $input" >&2; exit 1; }
    input_bytes=$((input_bytes + $(wc -c < "$input" | tr -d ' ')))
done < "$staging_plan"
available_kib=$(df -Pk "$scratch" | awk 'NR == 2 {print $4}')
required_kib=$((input_bytes / 1024 + 2 * 1024 * 1024))
[[ $available_kib =~ ^[0-9]+$ ]] || {
    echo "E: Could not determine node-local scratch capacity." >&2
    exit 1
}
if (( available_kib < required_kib )); then
    echo "E: Setup needs staged inputs plus a 2 GiB reserve; scratch has only $((available_kib / 1024 / 1024)) GiB free." >&2
    exit 1
fi
echo "I: Staging $((input_bytes / 1024 / 1024)) MiB of pinned inputs; scratch has $((available_kib / 1024 / 1024)) GiB free." >&2

phase=staging_inputs
while IFS=$'\t' read -r role support input expected_sha staged_name; do
    [[ $role != role ]] || continue
    current_input=${support:-$role}
    [[ -f $input ]] || { echo "E: Planned setup input is missing: $input" >&2; exit 1; }
    staged="$scratch/inputs/$staged_name"
    run_child cp "$input" "$staged"
    [[ $(sha256_file "$staged") == "$expected_sha" ]] || {
        echo "E: Staged input checksum differs: $input" >&2
        exit 1
    }
done < "$staging_plan"

plus=$(awk -F '\t' '$1 == "anchor_plus" {print $5}' "$staging_plan")
minus=$(awk -F '\t' '$1 == "anchor_minus" {print $5}' "$staging_plan")
coverage_arguments=()
while IFS=$'\t' read -r role support _ _ staged_name; do
    [[ $role == coverage ]] || continue
    coverage_arguments+=(--coverage "$support=$scratch/inputs/$staged_name")
done < "$staging_plan"

phase=building_depth_anchor
current_input=all_tracks
scratch_anchor="$scratch/tp73_chr1_anchor_evidence.parquet"
run_child "$source/scripts/build_tp73_anchor_evidence.py" \
    --anchor-plus "$scratch/inputs/$plus" --anchor-minus "$scratch/inputs/$minus" \
    "${coverage_arguments[@]}" --output "$scratch_anchor" --chrom 1 \
    --minimum-anchor-score 0 --duckdb "$duckdb" --threads "$threads" \
    --memory-limit "$memory_limit"

phase=validating_depth_anchor
current_input=output
support_columns=$(python3 - "$coverage_manifest" <<'PY'
import csv
import sys
with open(sys.argv[1], newline="") as handle:
    print(",".join(row["support_column"] for row in csv.DictReader(handle, delimiter="\t")))
PY
)
valid=$("$duckdb" -light-mode -csv -noheader :memory: -c "
WITH old AS (
  SELECT chrom, anchor_start, anchor_end, anchor_score, $support_columns
  FROM read_parquet('$base_anchor')
), new AS (
  SELECT chrom, anchor_start, anchor_end, anchor_score, $support_columns
  FROM read_parquet('$scratch_anchor')
)
SELECT (SELECT count(*) FROM old) = $expected_anchor_count
   AND (SELECT count(*) FROM new) = (SELECT count(*) FROM old)
   AND NOT EXISTS (SELECT * FROM old EXCEPT SELECT * FROM new)
   AND NOT EXISTS (SELECT * FROM new EXCEPT SELECT * FROM old);"
)
[[ $valid == true ]] || {
    "$duckdb" -light-mode -table :memory: -c "
WITH old AS (
  SELECT chrom, anchor_start, anchor_end, anchor_score, $support_columns
  FROM read_parquet('$base_anchor')
), new AS (
  SELECT chrom, anchor_start, anchor_end, anchor_score, $support_columns
  FROM read_parquet('$scratch_anchor')
)
SELECT (SELECT count(*) FROM old) AS old_rows,
       (SELECT count(*) FROM new) AS new_rows,
       (SELECT count(*) FROM (SELECT * FROM old EXCEPT SELECT * FROM new))
           AS old_only_rows,
       (SELECT count(*) FROM (SELECT * FROM new EXCEPT SELECT * FROM old))
           AS new_only_rows;" >&2
    echo "E: Depth-anchor build changed the established anchor/support table." >&2
    exit 1
}

attempt="$run_root/input/.attempt-depth-anchor-job-${SLURM_JOB_ID:-manual}-restart-${SLURM_RESTART_COUNT:-0}-pid-$$"
[[ ! -e $attempt ]] || { echo "E: Durable anchor attempt exists: $attempt" >&2; exit 1; }
mkdir -p "$attempt"
cp "$scratch_anchor" "$attempt/tp73_chr1_anchor_evidence.parquet"
cp "$scratch_anchor.run_config.tsv" \
    "$attempt/tp73_chr1_anchor_evidence.parquet.run_config.tsv"
cp "$scratch_anchor.coverage_manifest.tsv" \
    "$attempt/tp73_chr1_anchor_evidence.parquet.coverage_manifest.tsv"
cp "$base_config" "$attempt/source_anchor_run_config.tsv"
cp "$coverage_manifest" "$attempt/source_coverage_manifest.tsv"
python3 - "$attempt" "$base_anchor" "$coverage_manifest" \
    "${SLURM_JOB_ID:-manual}" <<'PY'
import hashlib
import json
import pathlib
import sys
import time

directory = pathlib.Path(sys.argv[1])
anchor = directory / "tp73_chr1_anchor_evidence.parquet"
digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
(directory / "complete.json").write_text(json.dumps({
    "schema_version": 1,
    "anchor_sha256": digest(anchor),
    "anchor_bytes": anchor.stat().st_size,
    "source_anchor": str(pathlib.Path(sys.argv[2]).resolve()),
    "source_anchor_sha256": digest(pathlib.Path(sys.argv[2])),
    "source_coverage_manifest_sha256": digest(pathlib.Path(sys.argv[3])),
    "slurm_job_id": sys.argv[4],
    "completed_epoch": int(time.time()),
}, indent=2, sort_keys=True) + "\n")
PY
[[ ! -e $anchor_package ]] || {
    echo "E: Final anchor package appeared concurrently: $anchor_package" >&2
    exit 1
}
mv "$attempt" "$anchor_package"

phase=preparing_enrichment_plan
run_child python3 "$source/scripts/manage_tp73_cofactor_enrichment.py" prepare \
    --run-root "$run_root" --source-threshold-run "$source_threshold_run" \
    --anchor-evidence "$anchor" --source "$source" --duckdb "$duckdb"
phase=complete
echo "I: Completed enrichment setup and immutable plan: $run_root" >&2

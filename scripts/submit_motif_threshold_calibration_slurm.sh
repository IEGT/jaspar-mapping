#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_motif_threshold_calibration_slurm.sh --run-root DIR
       --scan-package DIR --jaspar FILE
       (--cutandrun-dir DIR | --anchor-evidence FILE) [options]

Prepare and submit a requeueable chromosome-1 threshold calibration for every
non-target JASPAR motif. One array task reads exactly two inventory files for
one motif. A setup job reconstructs the shared TP73/CUT&RUN labels, and a final
dependency consolidates the finished metric tables into the threshold registry.
TP73 itself is deliberately excluded: its direct and tandem thresholds are
separate target policies, not a self-cofactor calibration.

Options:
  --run-root DIR          Dedicated durable run below /data/sm718
  --scan-package DIR      Finalized JASPAR 2026 sparse scan package
  --jaspar FILE           Exact JASPAR 2026 matrix source
  --cutandrun-dir DIR     Existing no-duplicates CUT&RUN resources
  --anchor-evidence FILE  Prebuilt TP73/CUT&RUN evidence; skips anchor setup
  --source DIR            Repository root (default: script parent)
  --runtime-prefix DIR    Micromamba runtime (default: RUN/runtime)
  --micromamba FILE       Micromamba executable (default: micromamba)
  --account NAME          Slurm account (default: cluster)
  --partition NAME        Slurm partition (default: requeue)
  --max-concurrent N      Live array tasks (default: 20)
  --array-size N          Tasks per scheduler array (default: 1000)
  --cpus N                CPUs per motif task (default: 2)
  --memory SIZE           Memory per motif task (default: 16G)
  --time LIMIT            Time per motif task (default: 02:00:00)
  --reuse-plan            Reuse a compatible immutable prepared plan
  --dry-run               Print submissions without calling sbatch
  -h, --help              Show this help

Separate DuckDB/Python and R runtimes are created once from conda-forge when
absent, avoiding incompatible ICU constraints; both explicit solutions are
recorded. Existing runtime, plan, anchor, and task outputs
are reused only after contract checks; source scan payloads are never changed.
Plan reuse also verifies that scientific scripts match the recorded plan commit.
EOF
}

run_root=""
scan_package=""
jaspar=""
cutandrun_dir=""
anchor_evidence=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
runtime_prefix=""
micromamba=${MAMBA_EXE:-micromamba}
account=cluster
partition=requeue
max_concurrent=20
array_size=1000
cpus=2
memory=16G
wall_time=02:00:00
reuse_plan=0
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --scan-package) scan_package=${2:?}; shift 2 ;;
        --jaspar) jaspar=${2:?}; shift 2 ;;
        --cutandrun-dir) cutandrun_dir=${2:?}; shift 2 ;;
        --anchor-evidence) anchor_evidence=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --runtime-prefix) runtime_prefix=${2:?}; shift 2 ;;
        --micromamba) micromamba=${2:?}; shift 2 ;;
        --account) account=${2:?}; shift 2 ;;
        --partition) partition=${2:?}; shift 2 ;;
        --max-concurrent) max_concurrent=${2:?}; shift 2 ;;
        --array-size) array_size=${2:?}; shift 2 ;;
        --cpus) cpus=${2:?}; shift 2 ;;
        --memory) memory=${2:?}; shift 2 ;;
        --time) wall_time=${2:?}; shift 2 ;;
        --reuse-plan) reuse_plan=1; shift ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $scan_package && -n $jaspar ]] || {
    usage >&2
    exit 2
}
if [[ -n $anchor_evidence && -n $cutandrun_dir ]] ||
   [[ -z $anchor_evidence && -z $cutandrun_dir ]]; then
    echo "E: Provide exactly one of --cutandrun-dir or --anchor-evidence." >&2
    exit 2
fi
[[ $max_concurrent =~ ^[1-9][0-9]*$ &&
   $array_size =~ ^[1-9][0-9]*$ && $cpus =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --max-concurrent, --array-size, and --cpus must be positive integers." >&2
    exit 2
}
mkdir -p "$run_root/plan" "$run_root/logs" "$run_root/input" "$run_root/tasks"
run_root=$(cd "$run_root" && pwd -P)
scan_package=$(cd "$scan_package" && pwd -P)
jaspar=$(cd "$(dirname "$jaspar")" && pwd -P)/$(basename "$jaspar")
if [[ -n $cutandrun_dir ]]; then
    cutandrun_dir=$(cd "$cutandrun_dir" && pwd -P)
else
    anchor_evidence=$(cd "$(dirname "$anchor_evidence")" && pwd -P)/$(basename "$anchor_evidence")
    [[ -f $anchor_evidence ]] || {
        echo "E: Prebuilt anchor evidence is missing: $anchor_evidence" >&2
        exit 1
    }
fi
source=$(cd "$source" && pwd -P)
runtime_prefix=${runtime_prefix:-$run_root/runtime}

duckdb_prefix="$runtime_prefix/duckdb"
r_prefix="$runtime_prefix/r"
if [[ ! -x $duckdb_prefix/bin/duckdb ||
      ! -x $duckdb_prefix/bin/python3 ]]; then
    [[ $dry_run -eq 0 ]] || {
        echo "E: Dry-run requires an existing runtime: $runtime_prefix" >&2
        exit 1
    }
    [[ $micromamba == */* ]] || micromamba=$(command -v "$micromamba" || true)
    [[ -x $micromamba ]] || { echo "E: Micromamba is unavailable." >&2; exit 1; }
    "$micromamba" create --yes --override-channels --channel conda-forge \
        --channel bioconda \
        --prefix "$duckdb_prefix" duckdb-cli=1.5.4 pybigwig python=3.12
fi
if [[ ! -x $r_prefix/bin/Rscript ]]; then
    [[ $dry_run -eq 0 ]] || {
        echo "E: Dry-run requires an existing R runtime: $r_prefix" >&2
        exit 1
    }
    [[ $micromamba == */* ]] || micromamba=$(command -v "$micromamba" || true)
    [[ -x $micromamba ]] || { echo "E: Micromamba is unavailable." >&2; exit 1; }
    "$micromamba" create --yes --override-channels --channel conda-forge \
        --prefix "$r_prefix" r-base=4.5.1 r-data.table=1.18.4
fi
if [[ ! -s $run_root/plan/runtime-duckdb-explicit.txt ||
      ! -s $run_root/plan/runtime-r-explicit.txt ]]; then
    [[ $micromamba == */* ]] || micromamba=$(command -v "$micromamba" || true)
    [[ -x $micromamba ]] || { echo "E: Micromamba is needed to record runtimes." >&2; exit 1; }
    "$micromamba" list --prefix "$duckdb_prefix" --explicit \
        > "$run_root/plan/runtime-duckdb-explicit.txt"
    "$micromamba" list --prefix "$r_prefix" --explicit \
        > "$run_root/plan/runtime-r-explicit.txt"
fi
duckdb="$duckdb_prefix/bin/duckdb"
rscript="$r_prefix/bin/Rscript"
bigwig_python="$duckdb_prefix/bin/python3"
"$rscript" -e 'stopifnot(requireNamespace("data.table", quietly=TRUE))'
"$bigwig_python" -c 'import pyBigWig'

anchor=${anchor_evidence:-$run_root/input/tp73_chr1_anchor_evidence.parquet}
task_file="$run_root/plan/calibration_tasks.tsv"
run_config="$run_root/plan/run_config.json"
if [[ $reuse_plan -eq 1 ]]; then
    [[ -s $task_file && -s $run_config ]] || {
        echo "E: --reuse-plan requires an existing task file and run configuration." >&2
        exit 1
    }
    plan_values=$(python3 - "$run_config" "$task_file" "$scan_package" \
        "$jaspar" "$source" "$anchor" <<'PY'
import csv
import json
import pathlib
import sys

config_path, task_path, scan, jaspar, source, anchor = map(pathlib.Path, sys.argv[1:])
config = json.loads(config_path.read_text())
expected = {
    "scan_package": str(scan.resolve()),
    "jaspar": str(jaspar.resolve()),
    "source": str(source.resolve()),
    "anchor_evidence": str(anchor.resolve()),
}
for key, value in expected.items():
    if config.get(key) != value:
        raise SystemExit(f"prepared plan {key} differs: {config.get(key)!r} != {value!r}")
with task_path.open(newline="") as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
count = int(config.get("cofactor_task_count", -1))
if count <= 0 or len(rows) != count:
    raise SystemExit("prepared plan task count is inconsistent")
if [int(row["task_index"]) for row in rows] != list(range(count)):
    raise SystemExit("prepared plan task indices are not contiguous")
print(f'{count}\t{config["source_commit"]}')
PY
    )
    IFS=$'\t' read -r task_count plan_source_commit <<< "$plan_values"
    scientific_paths=(
        scripts/build_tp73_anchor_evidence.py
        scripts/build_sparse_context_maxima.py
        scripts/evaluate_tp73_cofactor_thresholds.R
        scripts/build_motif_score_thresholds.py
        scripts/manage_motif_threshold_calibration.py
        scripts/export_bigwig_chrom_bedgraph.py
    )
    git -C "$source" cat-file -e "$plan_source_commit^{commit}" 2>/dev/null || {
        echo "E: Prepared plan source commit is unavailable: $plan_source_commit" >&2
        exit 1
    }
    if ! git -C "$source" diff --quiet "$plan_source_commit..HEAD" -- \
        "${scientific_paths[@]}"; then
        echo "E: Scientific calibration code changed since the prepared plan." >&2
        exit 1
    fi
    if ! git -C "$source" diff --quiet -- "${scientific_paths[@]}" ||
       ! git -C "$source" diff --cached --quiet -- "${scientific_paths[@]}"; then
        echo "E: Scientific calibration code has uncommitted changes." >&2
        exit 1
    fi
else
    task_count=$(python3 "$source/scripts/manage_motif_threshold_calibration.py" prepare \
        --run-root "$run_root" --scan-package "$scan_package" --jaspar "$jaspar" \
        --anchor-evidence "$anchor" --source "$source" --duckdb "$duckdb")
    plan_source_commit=$(python3 -c \
        'import json,sys; print(json.load(open(sys.argv[1]))["source_commit"])' \
        "$run_config")
fi
[[ $task_count =~ ^[1-9][0-9]*$ ]] || {
    echo "E: Invalid prepared task count: $task_count" >&2
    exit 1
}
execution_source_commit=$(git -C "$source" rev-parse HEAD)
execution_source_dirty=0
if ! git -C "$source" diff --quiet --ignore-submodules -- ||
   ! git -C "$source" diff --cached --quiet --ignore-submodules --; then
    execution_source_dirty=1
fi
finalization_provenance=(
    --finalization-source-commit "$execution_source_commit"
)
if [[ $execution_source_dirty -eq 1 ]]; then
    finalization_provenance+=(--finalization-source-dirty)
fi
submission_record="$run_root/plan/slurm_submission.tsv"
if [[ -e $submission_record && $dry_run -eq 0 ]]; then
    echo "E: Run already has a Slurm submission record: $submission_record" >&2
    exit 1
fi

dependency=()
setup_job="reused"
if [[ -z $anchor_evidence && ! -s $anchor ]]; then
    setup_submission=(
        sbatch --parsable --account="$account" --partition="$partition" --requeue
        --job-name=jaspar_thr_anchor --nodes=1 --ntasks=1 --cpus-per-task=4
        --mem=16G --time=02:00:00 --signal=B:USR1@300 --chdir="$source"
        --output="$run_root/logs/anchor-%j.out" --error="$run_root/logs/anchor-%j.err"
        "$source/scripts/run_motif_threshold_anchor_setup.sh"
        --run-root "$run_root" --scan-package "$scan_package"
        --cutandrun-dir "$cutandrun_dir" --output "$anchor" --source "$source"
        --duckdb "$duckdb" --bigwig-python "$bigwig_python"
        --bigwig-exporter "$source/scripts/export_bigwig_chrom_bedgraph.py"
        --threads 4 --memory-limit 12GB
    )
    if [[ $dry_run -eq 1 ]]; then
        printf '%q ' "${setup_submission[@]}"; printf '\n'
        setup_job=SETUP_JOB_ID
    else
        setup_job=$("${setup_submission[@]}")
    fi
    dependency=(--dependency="afterok:$setup_job")
fi

array_jobs=()
array_start=0
array_chunk=0
while (( array_start < task_count )); do
    chunk_tasks=$((task_count - array_start))
    if (( chunk_tasks > array_size )); then
        chunk_tasks=$array_size
    fi
    array_end=$((chunk_tasks - 1))
    array_submission=(
        sbatch --parsable --account="$account" --partition="$partition" --requeue
        --job-name="jaspar2026_thr_${array_chunk}"
        --array="0-${array_end}%${max_concurrent}"
        --nodes=1 --ntasks=1 --cpus-per-task="$cpus" --mem="$memory"
        --time="$wall_time" --signal=B:USR1@300 --chdir="$source"
        --output="$run_root/logs/threshold-%A_%a.out"
        --error="$run_root/logs/threshold-%A_%a.err"
        "${dependency[@]}"
        "$source/scripts/run_motif_threshold_calibration_slurm_task.sh"
        --run-root "$run_root" --scan-package "$scan_package" --task-file "$task_file"
        --anchor-evidence "$anchor" --source "$source" --duckdb "$duckdb"
        --rscript "$rscript" --task-offset "$array_start"
        --threads "$cpus" --memory-limit 12GB
        --max-temp-size 80GB --minimum-class-fraction 0.01
    )
    if [[ $dry_run -eq 1 ]]; then
        printf '%q ' "${array_submission[@]}"; printf '\n'
        array_job="ARRAY_JOB_ID_${array_chunk}"
    else
        array_job=$("${array_submission[@]}")
    fi
    array_jobs+=("$array_job")
    dependency=(--dependency="afterok:$array_job")
    array_start=$((array_start + chunk_tasks))
    array_chunk=$((array_chunk + 1))
done
last_array_job=$array_job
array_jobs_text=$(IFS=,; printf '%s' "${array_jobs[*]}")

final_submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name=jaspar_thr_final --nodes=1 --ntasks=1 --cpus-per-task=2
    --mem=8G --time=01:00:00 --chdir="$source"
    --dependency="afterok:$last_array_job"
    --output="$run_root/logs/finalize-%j.out"
    --error="$run_root/logs/finalize-%j.err"
    "$source/scripts/run_motif_threshold_calibration_finalize.sh"
    --run-root "$run_root" --source "$source" --duckdb "$duckdb"
    "${finalization_provenance[@]}"
)
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${final_submission[@]}"; printf '\n'
    exit 0
fi
final_job=$("${final_submission[@]}")

printf 'setup_job\tarray_jobs\tfinalize_job\ttask_count\tmax_concurrent\tarray_size\tarray_chunks\tplan_source_commit\texecution_source_commit\treused_plan\n' \
    > "$submission_record"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$setup_job" "$array_jobs_text" "$final_job" "$task_count" \
    "$max_concurrent" "$array_size" "$array_chunk" "$plan_source_commit" \
    "$execution_source_commit" "$reuse_plan" \
    >> "$submission_record"
echo "I: Submitted setup=$setup_job arrays=$array_jobs_text finalizer=$final_job tasks=$task_count live=$max_concurrent chunks=$array_chunk" >&2
printf '%s\n' "$array_jobs_text"

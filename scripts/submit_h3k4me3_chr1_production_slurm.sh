#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: submit_h3k4me3_chr1_production_slurm.sh --run-root DIR
       --annotation-run DIR --scan-package DIR --track-root DIR
       --runtime-prefix DIR [options]

Submit one restart-safe single-chromosome production job. The job resolves the
schema-7 TP73 local-peak anchor file through its finalized DuckDB inventory,
stages the nine prespecified cofactor motifs by exact scan-inventory path,
copies motif inputs to node-local scratch, builds TP73/control and H3K4me3/input
evidence, and fits the GFP-referenced H3K4me3 change analysis.

Options:
  --run-root DIR          Dedicated durable output below /data/sm718
  --annotation-run DIR    Finalized schema-7 TP73 annotation run
  --scan-package DIR      Finalized sparse-v3 scan package (score floor -1)
  --track-root DIR        Directory containing manifest-named BigWig files
  --runtime-prefix DIR    Runtime with duckdb/, r/, and pyBigWig Python
  --chrom NAME            Chromosome identifier (default: 1)
  --chrom-length BP       Chromosome length (required unless chromosome 1)
  --source DIR            Repository root (default: script parent)
  --scratch-root DIR      Node-local scratch parent (default: /scratch/$USER)
  --account NAME          Slurm account (default: cluster)
  --partition NAME        Slurm partition (default: requeue)
  --cpus N                CPUs (default: 4)
  --memory SIZE           Memory (default: 32G)
  --time LIMIT            Wall time (default: 08:00:00)
  --include-occupancy-analysis
                          Also fit matched anti-p73/control occupancy effects
  --dry-run               Print the exact sbatch command only
  -h, --help              Show this help

The source tree must be tracked-clean. The job listens to SIGUSR1 and reports
its phase, elapsed time, durable bytes, and scratch bytes without terminating.
An interrupted attempt is retained; only a completely validated attempt is
atomically promoted to final/.
EOF
}

run_root=""
annotation_run=""
scan_package=""
track_root=""
runtime_prefix=""
chrom=1
chrom_length=""
source=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
scratch_root="/scratch/${USER}"
account=cluster
partition=requeue
cpus=4
memory=32G
wall_time=08:00:00
include_occupancy_analysis=0
dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root=${2:?}; shift 2 ;;
        --annotation-run) annotation_run=${2:?}; shift 2 ;;
        --scan-package) scan_package=${2:?}; shift 2 ;;
        --track-root) track_root=${2:?}; shift 2 ;;
        --runtime-prefix) runtime_prefix=${2:?}; shift 2 ;;
        --chrom) chrom=${2:?}; shift 2 ;;
        --chrom-length) chrom_length=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --scratch-root) scratch_root=${2:?}; shift 2 ;;
        --account) account=${2:?}; shift 2 ;;
        --partition) partition=${2:?}; shift 2 ;;
        --cpus) cpus=${2:?}; shift 2 ;;
        --memory) memory=${2:?}; shift 2 ;;
        --time) wall_time=${2:?}; shift 2 ;;
        --include-occupancy-analysis) include_occupancy_analysis=1; shift ;;
        --dry-run) dry_run=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

[[ -n $run_root && -n $annotation_run && -n $scan_package &&
   -n $track_root && -n $runtime_prefix ]] || { usage >&2; exit 2; }
[[ $cpus =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --cpus must be a positive integer." >&2
    exit 2
}
[[ $chrom =~ ^[A-Za-z0-9_.-]+$ ]] || {
    echo "E: --chrom contains unsafe characters." >&2
    exit 2
}
if [[ -z $chrom_length ]]; then
    [[ $chrom == 1 ]] || {
        echo "E: --chrom-length is required when --chrom is not 1." >&2
        exit 2
    }
    chrom_length=248956422
fi
[[ $chrom_length =~ ^[1-9][0-9]*$ ]] || {
    echo "E: --chrom-length must be a positive integer." >&2
    exit 2
}
for path in "$run_root" "$annotation_run" "$scan_package" "$track_root" \
            "$runtime_prefix"; do
    case "$path" in
        /data/sm718/*) ;;
        *) echo "E: Production paths must be below /data/sm718: $path" >&2; exit 2 ;;
    esac
done

source=$(cd "$source" && pwd -P)
annotation_run=$(cd "$annotation_run" && pwd -P)
scan_package=$(cd "$scan_package" && pwd -P)
track_root=$(cd "$track_root" && pwd -P)
runtime_prefix=$(cd "$runtime_prefix" && pwd -P)
mkdir -p "$run_root/plan" "$run_root/logs"
run_root=$(cd "$run_root" && pwd -P)

duckdb="$runtime_prefix/duckdb/bin/duckdb"
bigwig_python="$runtime_prefix/duckdb/bin/python3"
rscript="$runtime_prefix/r/bin/Rscript"
[[ -x $duckdb && -x $bigwig_python && -x $rscript ]] || {
    echo "E: Runtime lacks DuckDB, pyBigWig Python, or Rscript: $runtime_prefix" >&2
    exit 1
}
"$bigwig_python" -c 'import pyBigWig'
"$rscript" -e 'stopifnot(requireNamespace("data.table", quietly=TRUE))'
[[ -f $annotation_run/final/manifest.json && -f $scan_package/manifest.json ]] || {
    echo "E: Annotation or scan package is not finalized." >&2
    exit 1
}
if ! git -C "$source" diff --quiet --ignore-submodules -- ||
   ! git -C "$source" diff --cached --quiet --ignore-submodules --; then
    echo "E: Production execution requires a tracked-clean source tree." >&2
    exit 1
fi

submission_record="$run_root/plan/slurm_submission.tsv"
if [[ -e $submission_record && $dry_run -eq 0 ]]; then
    echo "E: Run already has a Slurm submission record: $submission_record" >&2
    exit 1
fi
source_commit=$(git -C "$source" rev-parse HEAD)
submission=(
    sbatch --parsable --account="$account" --partition="$partition" --requeue
    --job-name="tp73_h3k4me3_chr${chrom}" --nodes=1 --ntasks=1
    --cpus-per-task="$cpus" --mem="$memory" --time="$wall_time"
    --signal=B:USR1@300 --open-mode=append --chdir="$source"
    --output="$run_root/logs/production-%j.out"
    --error="$run_root/logs/production-%j.err"
    "$source/scripts/run_h3k4me3_chr1_production.py"
    --run-root "$run_root" --annotation-run "$annotation_run"
    --scan-package "$scan_package" --track-root "$track_root"
    --source "$source" --source-commit "$source_commit"
    --chrom "$chrom" --chrom-length "$chrom_length"
    --duckdb "$duckdb" --rscript "$rscript"
    --bigwig-python "$bigwig_python" --scratch-root "$scratch_root"
    --threads "$cpus" --memory-limit 28GB --minimum-free-scratch-gb 30
)
if [[ $include_occupancy_analysis -eq 1 ]]; then
    submission+=(--include-occupancy-analysis)
fi
if [[ $dry_run -eq 1 ]]; then
    printf '%q ' "${submission[@]}"
    printf '\n'
    exit 0
fi

job_id=$("${submission[@]}")
printf 'job_id\tsource_commit\tannotation_run\tscan_package\ttrack_root\truntime_prefix\tchrom\tchrom_length\tinclude_occupancy_analysis\n' \
    > "$submission_record"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$job_id" "$source_commit" \
    "$annotation_run" "$scan_package" "$track_root" "$runtime_prefix" \
    "$chrom" "$chrom_length" "$include_occupancy_analysis" \
    >> "$submission_record"
echo "I: Submitted H3K4me3 production job $job_id" >&2
printf '%s\n' "$job_id"

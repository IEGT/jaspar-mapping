#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-context-worker.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

source_tree="$temporary/source"
run_root="$temporary/run"
scan_package="$temporary/scan"
mkdir -p "$source_tree/scripts" "$run_root/plan" "$scan_package" \
    "$run_root/packages/chrom-X/task-107"
printf '{}\n' > "$scan_package/manifest.json"

cat > "$source_tree/scripts/stage_motif_context_inputs.py" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
output=""
while [[ $# -gt 0 ]]; do
    if [[ $1 == --output ]]; then output=$2; shift 2; else shift; fi
done
mkdir -p "$output"
printf '{"motifs":["MA0861.2","MA0001.1","MA0002.1"],"chromosomes":["X"]}\n' \
    > "$output/input_manifest.json"
EOF
cat > "$source_tree/scripts/build_motif_context.py" <<'EOF'
#!/usr/bin/env bash
exit 99
EOF
chmod +x "$source_tree/scripts/stage_motif_context_inputs.py" \
    "$source_tree/scripts/build_motif_context.py"

cat > "$temporary/duckdb" <<'EOF'
#!/usr/bin/env bash
case "$*" in
    *motif_context_run_config*) printf '1\n' ;;
    *anchor_motif_band_feature*) printf '0\n' ;;
    *) exit 2 ;;
esac
EOF
chmod +x "$temporary/duckdb"

printf 'task_index\tchrom\tcofactor_motif_ids\toutput_tier\tbuilder_source_commit\n' \
    > "$run_root/plan/context_tasks.tsv"
printf '107\tX\tMA0001.1,MA0002.1\tband\tunknown\n' \
    >> "$run_root/plan/context_tasks.tsv"
touch "$run_root/packages/chrom-X/task-107/context.duckdb"
printf '{"motifs":["MA0861.2","MA0001.1","MA0002.1"],"chromosomes":["X"]}\n' \
    > "$run_root/packages/chrom-X/task-107/input_manifest.json"

output=$(
    SLURM_ARRAY_TASK_ID=7 JASPAR_CONTEXT_TASK_OFFSET=100 \
    "$repository_root/scripts/run_motif_context_slurm_task.sh" \
        --run-root "$run_root" --scan-package "$scan_package" \
        --task-file "$run_root/plan/context_tasks.tsv" \
        --source "$source_tree" --duckdb "$temporary/duckdb" 2>&1
)
grep -Fq 'Reusing completed context task 107' <<< "$output"
grep -Fq 'packages/chrom-X/task-107' <<< "$output"

if SLURM_ARRAY_TASK_ID=7 JASPAR_CONTEXT_TASK_OFFSET=-1 \
    "$repository_root/scripts/run_motif_context_slurm_task.sh" \
        --run-root "$run_root" --scan-package "$scan_package" \
        --task-file "$run_root/plan/context_tasks.tsv" \
        --source "$source_tree" --duckdb "$temporary/duckdb" \
        >/dev/null 2>&1; then
    echo "E: Negative context task offset was accepted." >&2
    exit 1
fi

echo "Motif-context Slurm worker tests passed."

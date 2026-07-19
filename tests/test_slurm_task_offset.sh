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

echo "Slurm task-offset tests passed."

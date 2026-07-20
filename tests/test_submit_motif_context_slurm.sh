#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-context-submit.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM
mkdir -p "$temporary/scan"
printf '{}\n' > "$temporary/scan/manifest.json"
printf 'gtf\n' > "$temporary/annotation.gtf"

rendered=$("$repository_root/scripts/submit_motif_context_slurm.sh" \
    --run-root "$temporary/run" --scan-package "$temporary/scan" \
    --gtf "$temporary/annotation.gtf" \
    --motif MA1961.2 --motif MA0024.3 --chrom 1 --chrom X \
    --max-concurrent 3 --cpus 2 --memory 20G \
    --memory-limit 14GB --max-temp-size 80GB --dry-run)

grep -Fq -- '--array=0-3%3' <<< "$rendered"
grep -Fq -- '--cpus-per-task=2' <<< "$rendered"
grep -Fq -- '--mem=20G' <<< "$rendered"
grep -Fq -- '--memory-limit 14GB' <<< "$rendered"
grep -Fq -- '--max-temp-size 80GB' <<< "$rendered"
[[ $(wc -l < "$temporary/run/plan/context_tasks.tsv") -eq 5 ]]

# An identical dry run reuses the plan; changing it is rejected.
"$repository_root/scripts/submit_motif_context_slurm.sh" \
    --run-root "$temporary/run" --scan-package "$temporary/scan" \
    --gtf "$temporary/annotation.gtf" \
    --motif MA1961.2 --motif MA0024.3 --chrom 1 --chrom X --dry-run >/dev/null
if "$repository_root/scripts/submit_motif_context_slurm.sh" \
    --run-root "$temporary/run" --scan-package "$temporary/scan" \
    --gtf "$temporary/annotation.gtf" \
    --motif MA1961.2 --chrom 1 --dry-run >/dev/null 2>&1; then
    echo "E: changed context task plan was accepted" >&2
    exit 1
fi

echo "Motif-context Slurm submission tests passed."

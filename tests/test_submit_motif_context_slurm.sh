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
    --motifs-per-task 1 \
    --max-concurrent 3 --cpus 2 --memory 20G \
    --memory-limit 14GB --max-temp-size 80GB --dry-run)

grep -Fq -- '--array=0-3%3' <<< "$rendered"
grep -Fq -- '--cpus-per-task=2' <<< "$rendered"
grep -Fq -- '--mem=20G' <<< "$rendered"
grep -Fq -- '--memory-limit 14GB' <<< "$rendered"
grep -Fq -- '--max-temp-size 80GB' <<< "$rendered"
[[ $(wc -l < "$temporary/run/plan/context_tasks.tsv") -eq 5 ]]
cmp -s "$repository_root/scripts/build_motif_context.py" \
    "$temporary/run/source/scripts/build_motif_context.py"
cmp -s "$repository_root/scripts/stage_motif_context_inputs.py" \
    "$temporary/run/source/scripts/stage_motif_context_inputs.py"
[[ $(tr -d '\r\n' < "$temporary/run/source/source_commit.txt") == \
   "$(git -C "$repository_root" rev-parse HEAD)" ]]
snapshot=$(cd "$temporary/run/source" && pwd)
grep -Fq -- "--source $snapshot" <<< "$rendered"

# An identical dry run reuses the plan; changing it is rejected.
"$repository_root/scripts/submit_motif_context_slurm.sh" \
    --run-root "$temporary/run" --scan-package "$temporary/scan" \
    --gtf "$temporary/annotation.gtf" \
    --motif MA1961.2 --motif MA0024.3 --chrom 1 --chrom X \
    --motifs-per-task 1 --dry-run >/dev/null
if "$repository_root/scripts/submit_motif_context_slurm.sh" \
    --run-root "$temporary/run" --scan-package "$temporary/scan" \
    --gtf "$temporary/annotation.gtf" \
    --motif MA1961.2 --chrom 1 --dry-run >/dev/null 2>&1; then
    echo "E: changed context task plan was accepted" >&2
    exit 1
fi

cat > "$temporary/motifs.txt" <<'EOF'
# Five synthetic cofactor accessions.
MA0001.1
MA0002.1
MA0003.1
MA0004.1
MA0005.1
EOF
printf '1\nX\n' > "$temporary/chromosomes.txt"
batched=$("$repository_root/scripts/submit_motif_context_slurm.sh" \
    --run-root "$temporary/batched" --scan-package "$temporary/scan" \
    --motif-file "$temporary/motifs.txt" \
    --chrom-file "$temporary/chromosomes.txt" \
    --motifs-per-task 2 --array-chunk-size 4 --output-tier band \
    --max-concurrent 3 --dry-run)

[[ $(wc -l < "$temporary/batched/plan/context_tasks.tsv") -eq 7 ]]
[[ $(wc -l <<< "$batched") -eq 2 ]]
grep -Fq -- '--array=0-3%3' <<< "$batched"
grep -Fq -- '--array=0-1%3' <<< "$batched"
grep -Fq -- 'JASPAR_CONTEXT_TASK_OFFSET=0' <<< "$batched"
grep -Fq -- 'JASPAR_CONTEXT_TASK_OFFSET=4' <<< "$batched"
grep -Fq -- 'afterany:DRY_RUN_CHUNK_0' <<< "$batched"
grep -Fq $'0\t1\tMA0001.1,MA0002.1\tband' \
    "$temporary/batched/plan/context_tasks.tsv"
grep -Fq $'1\tX\tMA0001.1,MA0002.1\tband' \
    "$temporary/batched/plan/context_tasks.tsv"
grep -Fq $'2\t1\tMA0003.1,MA0004.1\tband' \
    "$temporary/batched/plan/context_tasks.tsv"
grep -Fq $'5\tX\tMA0005.1\tband' \
    "$temporary/batched/plan/context_tasks.tsv"

echo "Motif-context Slurm submission tests passed."

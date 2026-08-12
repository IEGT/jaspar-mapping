#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-context-submit.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

expected_nuclear_chromosomes=$(printf '%s\n' {1..22} X Y)
actual_nuclear_chromosomes=$(
    sed '/^[[:space:]]*#/d; /^[[:space:]]*$/d' \
        "$repository_root/config/grch38_primary_nuclear_chromosomes.txt"
)
[[ $actual_nuclear_chromosomes == "$expected_nuclear_chromosomes" ]] || {
    echo "E: GRCh38 primary nuclear chromosome registry differs." >&2
    exit 1
}

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
grep -Fq -- '--job-name=tp73_context_v9' <<< "$rendered"
[[ $(wc -l < "$temporary/run/plan/context_tasks.tsv") -eq 5 ]]
cmp -s "$repository_root/scripts/build_motif_context.py" \
    "$temporary/run/source/scripts/build_motif_context.py"
cmp -s "$repository_root/scripts/stage_motif_context_inputs.py" \
    "$temporary/run/source/scripts/stage_motif_context_inputs.py"
cmp -s "$repository_root/scripts/finalize_motif_context_run.py" \
    "$temporary/run/source/scripts/finalize_motif_context_run.py"
[[ $(tr -d '\r\n' < "$temporary/run/source/source_commit.txt") == \
   "$(git -C "$repository_root" rev-parse HEAD)" ]]
snapshot=$(cd "$temporary/run/source" && pwd)
grep -Fq -- "--source $snapshot" <<< "$rendered"
grep -Fq -- 'tp73_context_finalize' <<< "$rendered"
grep -Eq $'\t9\t[0-9]+\t[0-9a-f]{64}\tensembl_113\ttss_upstream_2000_downstream_500_v1\t2000\t500\ttes_upstream_500_downstream_2000_v1\t500\t2000\tcofactor_context$' \
    "$temporary/run/plan/context_tasks.tsv"

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
[[ $(wc -l <<< "$batched") -eq 3 ]]
grep -Fq -- '--array=0-3%3' <<< "$batched"
grep -Fq -- '--array=0-1%3' <<< "$batched"
grep -Fq -- 'JASPAR_CONTEXT_TASK_OFFSET=0' <<< "$batched"
grep -Fq -- 'JASPAR_CONTEXT_TASK_OFFSET=4' <<< "$batched"
grep -Fq -- 'afterany:DRY_RUN_CHUNK_0' <<< "$batched"
grep -Fq -- 'afterany:DRY_RUN_CHUNK_1' <<< "$batched"
grep -Fq -- 'finalize_motif_context_run.py' <<< "$batched"
grep -Fq $'0\t1\tMA0001.1,MA0002.1\tband' \
    "$temporary/batched/plan/context_tasks.tsv"
grep -Fq $'1\tX\tMA0001.1,MA0002.1\tband' \
    "$temporary/batched/plan/context_tasks.tsv"
grep -Fq $'2\t1\tMA0003.1,MA0004.1\tband' \
    "$temporary/batched/plan/context_tasks.tsv"
grep -Fq $'5\tX\tMA0005.1\tband' \
    "$temporary/batched/plan/context_tasks.tsv"
awk -F '\t' 'NR > 1 && ($6 != 9 || $7 != 0 || $8 != "none" \
    || $9 != "ensembl_113" \
    || $10 != "tss_upstream_2000_downstream_500_v1" \
    || $11 != 2000 || $12 != 500 \
    || $13 != "tes_upstream_500_downstream_2000_v1" \
    || $14 != 500 || $15 != 2000 || $16 != "cofactor_context") { exit 1 }' \
    "$temporary/batched/plan/context_tasks.tsv"

anchor_rendered=$("$repository_root/scripts/submit_motif_context_slurm.sh" \
    --run-root "$temporary/anchor" --scan-package "$temporary/scan" \
    --gtf "$temporary/annotation.gtf" --anchor-only \
    --chrom 1 --chrom X --output-tier summary --dry-run)
[[ $(wc -l < "$temporary/anchor/plan/context_tasks.tsv") -eq 3 ]]
grep -Fq $'0\t1\tnone\tsummary' "$temporary/anchor/plan/context_tasks.tsv"
grep -Fq $'1\tX\tnone\tsummary' "$temporary/anchor/plan/context_tasks.tsv"
awk -F '\t' 'NR > 1 && $16 != "anchor_annotation" { exit 1 }' \
    "$temporary/anchor/plan/context_tasks.tsv"
grep -Fq -- '--array=0-1%8' <<< "$anchor_rendered"

if "$repository_root/scripts/submit_motif_context_slurm.sh" \
    --run-root "$temporary/invalid-anchor" --scan-package "$temporary/scan" \
    --gtf "$temporary/annotation.gtf" --anchor-only --motif MA0001.1 \
    --output-tier summary --dry-run >/dev/null 2>&1; then
    echo "E: Anchor-only submission accepted a cofactor motif." >&2
    exit 1
fi
if "$repository_root/scripts/submit_motif_context_slurm.sh" \
    --run-root "$temporary/invalid-anchor-tier" --scan-package "$temporary/scan" \
    --gtf "$temporary/annotation.gtf" --anchor-only \
    --output-tier selected --dry-run >/dev/null 2>&1; then
    echo "E: Anchor-only submission accepted a non-summary tier." >&2
    exit 1
fi

echo "Motif-context Slurm submission tests passed."

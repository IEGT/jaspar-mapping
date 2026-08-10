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
chrom=""
motifs=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --output) output=$2; shift 2 ;;
        --chrom) chrom=$2; shift 2 ;;
        --motif) motifs+=("$2"); shift 2 ;;
        *) shift ;;
    esac
done
mkdir -p "$output"
joined=""
for motif in "${motifs[@]}"; do
    [[ -z $joined ]] || joined+=","
    joined+="\"$motif\""
done
printf '{"motifs":[%s],"chromosomes":["%s"]}\n' "$joined" "$chrom" \
    > "$output/input_manifest.json"
EOF
cat > "$source_tree/scripts/build_motif_context.py" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
output=""
while [[ $# -gt 0 ]]; do
    if [[ $1 == --output ]]; then output=$2; shift 2; else shift; fi
done
mkdir -p "$output"
touch "$output/context.duckdb"
EOF
chmod +x "$source_tree/scripts/stage_motif_context_inputs.py" \
    "$source_tree/scripts/build_motif_context.py"
printf 'abc123\n' > "$source_tree/source_commit.txt"
printf 'gtf\n' > "$temporary/annotation.gtf"
printf 'bad\n' > "$temporary/annotation-wrong.gtf"
gtf_size=$(wc -c < "$temporary/annotation.gtf" | tr -d '[:space:]')
if command -v sha256sum >/dev/null 2>&1; then
    gtf_sha256=$(sha256sum "$temporary/annotation.gtf" | awk '{print $1}')
else
    gtf_sha256=$(shasum -a 256 "$temporary/annotation.gtf" | awk '{print $1}')
fi

cat > "$temporary/duckdb" <<'EOF'
#!/usr/bin/env bash
case "$*" in
    *motif_context_run_config*) printf '%s\n' "${FAKE_CONFIG_COUNT:-1}" ;;
    *anchor_motif_band_feature*) printf '0\n' ;;
    *) exit 2 ;;
esac
EOF
chmod +x "$temporary/duckdb"

{
    printf 'task_index\tchrom\tcofactor_motif_ids\toutput_tier\tbuilder_source_commit\tcontext_schema_version\tgtf_size_bytes\tgtf_sha256\tannotation_release\tpromoter_definition_id\tpromoter_upstream_bp\tpromoter_downstream_bp\ttask_kind\n'
    printf '107\tX\tMA0001.1,MA0002.1\tband\tabc123\t7\t0\tnone\tensembl_113\ttss_upstream_2000_downstream_500_v1\t2000\t500\tcofactor_context\n'
    printf '108\tX\tMA0001.1,MA0002.1\tband\tabc123\t7\t0\tnone\tensembl_113\ttss_upstream_2000_downstream_500_v1\t2000\t500\tcofactor_context\n'
    printf '109\tX\tMA0001.1,MA0002.1\tband\tabc123\t7\t0\tnone\tensembl_113\ttss_upstream_2000_downstream_500_v1\t2000\t500\tcofactor_context\n'
    printf '110\tX\tMA0001.1,MA0002.1\tselected\tabc123\t7\t%s\t%s\tensembl_113\ttss_upstream_2000_downstream_500_v1\t2000\t500\tcofactor_context\n' \
        "$gtf_size" "$gtf_sha256"
    printf '111\tX\tnone\tsummary\tabc123\t7\t%s\t%s\tensembl_113\ttss_upstream_2000_downstream_500_v1\t2000\t500\tanchor_annotation\n' \
        "$gtf_size" "$gtf_sha256"
} > "$run_root/plan/context_tasks.tsv"
touch "$run_root/packages/chrom-X/task-107/context.duckdb"
printf '{"motifs":["MA0861.2","MA0001.1","MA0002.1"],"chromosomes":["X"]}\n' \
    > "$run_root/packages/chrom-X/task-107/input_manifest.json"
mkdir -p "$run_root/packages/chrom-X/task-110"
touch "$run_root/packages/chrom-X/task-110/context.duckdb"
cp "$run_root/packages/chrom-X/task-107/input_manifest.json" \
    "$run_root/packages/chrom-X/task-110/input_manifest.json"

output=$(
    SLURM_ARRAY_TASK_ID=7 JASPAR_CONTEXT_TASK_OFFSET=100 \
    "$repository_root/scripts/run_motif_context_slurm_task.sh" \
        --run-root "$run_root" --scan-package "$scan_package" \
        --task-file "$run_root/plan/context_tasks.tsv" \
        --source "$source_tree" --duckdb "$temporary/duckdb" 2>&1
)
grep -Fq 'Reusing completed context task 107' <<< "$output"
grep -Fq 'packages/chrom-X/task-107' <<< "$output"

# Exercise scratch copy, fresh build, validation, atomic promotion, and the
# task-owned cleanup trap rather than returning through the reuse branch.
mkdir -p "$temporary/scratch-fresh"
fresh_output=$(
    SLURM_ARRAY_TASK_ID=8 JASPAR_CONTEXT_TASK_OFFSET=100 \
    SLURM_TMPDIR="$temporary/scratch-fresh" \
    "$repository_root/scripts/run_motif_context_slurm_task.sh" \
        --run-root "$run_root" --scan-package "$scan_package" \
        --task-file "$run_root/plan/context_tasks.tsv" \
        --source "$source_tree" --duckdb "$temporary/duckdb" \
        --max-temp-size 1MB 2>&1
)
grep -Fq 'Completed context task 108' <<< "$fresh_output"
[[ -f $run_root/packages/chrom-X/task-108/context.duckdb ]]
[[ -f $run_root/packages/chrom-X/task-108/input_manifest.json ]]
[[ -f $run_root/locks/context-task-108.lock ]]
[[ -z $(find "$run_root/staging/task-108" -mindepth 1 -print -quit) ]]
[[ -z $(find "$temporary/scratch-fresh" -mindepth 1 -print -quit) ]]

# A newly built package that fails validation must not be promoted, and both
# durable and node-local task-owned staging must be removed.
mkdir -p "$temporary/scratch-failure"
if FAKE_CONFIG_COUNT=0 SLURM_ARRAY_TASK_ID=9 JASPAR_CONTEXT_TASK_OFFSET=100 \
    SLURM_TMPDIR="$temporary/scratch-failure" \
    "$repository_root/scripts/run_motif_context_slurm_task.sh" \
        --run-root "$run_root" --scan-package "$scan_package" \
        --task-file "$run_root/plan/context_tasks.tsv" \
        --source "$source_tree" --duckdb "$temporary/duckdb" \
        --max-temp-size 1MB >/dev/null 2>&1; then
    echo "E: Invalid fresh context package was accepted." >&2
    exit 1
fi
[[ ! -e $run_root/packages/chrom-X/task-109 ]]
[[ -z $(find "$run_root/staging/task-109" -mindepth 1 -print -quit) ]]
[[ -z $(find "$temporary/scratch-failure" -mindepth 1 -print -quit) ]]

selected_output=$(
    SLURM_ARRAY_TASK_ID=10 JASPAR_CONTEXT_TASK_OFFSET=100 \
    "$repository_root/scripts/run_motif_context_slurm_task.sh" \
        --run-root "$run_root" --scan-package "$scan_package" \
        --task-file "$run_root/plan/context_tasks.tsv" \
        --source "$source_tree" --duckdb "$temporary/duckdb" \
        --gtf "$temporary/annotation.gtf" 2>&1
)
grep -Fq 'Reusing completed context task 110' <<< "$selected_output"

mkdir -p "$temporary/scratch-anchor"
anchor_output=$(
    SLURM_ARRAY_TASK_ID=11 JASPAR_CONTEXT_TASK_OFFSET=100 \
    SLURM_TMPDIR="$temporary/scratch-anchor" \
    "$repository_root/scripts/run_motif_context_slurm_task.sh" \
        --run-root "$run_root" --scan-package "$scan_package" \
        --task-file "$run_root/plan/context_tasks.tsv" \
        --source "$source_tree" --duckdb "$temporary/duckdb" \
        --gtf "$temporary/annotation.gtf" --max-temp-size 1MB 2>&1
)
grep -Fq 'kind anchor_annotation, 0 cofactors' <<< "$anchor_output"
grep -Fq 'Completed context task 111' <<< "$anchor_output"
grep -Fq '"motifs":["MA0861.2"]' \
    "$run_root/packages/chrom-X/task-111/input_manifest.json"
if SLURM_ARRAY_TASK_ID=10 JASPAR_CONTEXT_TASK_OFFSET=100 \
    "$repository_root/scripts/run_motif_context_slurm_task.sh" \
        --run-root "$run_root" --scan-package "$scan_package" \
        --task-file "$run_root/plan/context_tasks.tsv" \
        --source "$source_tree" --duckdb "$temporary/duckdb" \
        --gtf "$temporary/annotation-wrong.gtf" >/dev/null 2>&1; then
    echo "E: Same-size changed GTF was accepted by the reuse gate." >&2
    exit 1
fi

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

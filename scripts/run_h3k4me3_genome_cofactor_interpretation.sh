#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
Usage: run_h3k4me3_genome_cofactor_interpretation.sh --analysis-run DIR
       --enrichment-package DIR --h3-package DIR --output-dir DIR
       --source DIR --source-commit HASH --duckdb FILE --rscript FILE

Build the compact all-JASPAR TP73/H3K4me3 interpretation tables and figures
after a finalized schema-9 cofactor analysis. The summary is published
atomically by its Python builder; plots are rendered into a staging directory
and promoted together. A completed summary or figure set is reused.

Options:
  --analysis-run DIR        Run containing final/h3k4me3_cofactor_analysis
  --enrichment-package DIR  Final TP73 cofactor-enrichment package
  --h3-package DIR          Final genome H3K4me3 signal package
  --output-dir DIR          Durable interpretation directory
  --source DIR              Repository root
  --source-commit HASH      Exact source commit required at execution
  --duckdb FILE             DuckDB executable
  --rscript FILE            Rscript executable
  -h, --help                Show this help
EOF
}

analysis_run=
enrichment_package=
h3_package=
output_dir=
source=
source_commit=
duckdb=
rscript=
while [[ $# -gt 0 ]]; do
    case "$1" in
        --analysis-run) analysis_run=${2:?}; shift 2 ;;
        --enrichment-package) enrichment_package=${2:?}; shift 2 ;;
        --h3-package) h3_package=${2:?}; shift 2 ;;
        --output-dir) output_dir=${2:?}; shift 2 ;;
        --source) source=${2:?}; shift 2 ;;
        --source-commit) source_commit=${2:?}; shift 2 ;;
        --duckdb) duckdb=${2:?}; shift 2 ;;
        --rscript) rscript=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "E: Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

for value in analysis_run enrichment_package h3_package output_dir source \
             source_commit duckdb rscript; do
    [[ -n ${!value} ]] || { usage >&2; exit 2; }
done
for path in "$analysis_run" "$enrichment_package" "$h3_package" \
            "$output_dir" "$source" "$duckdb" "$rscript"; do
    case "$path" in
        /data/sm718/*) ;;
        *) echo "E: Production paths must be below /data/sm718: $path" >&2; exit 2 ;;
    esac
done
[[ -x $duckdb ]] || { echo "E: DuckDB is unavailable: $duckdb" >&2; exit 1; }
[[ -x $rscript ]] || { echo "E: Rscript is unavailable: $rscript" >&2; exit 1; }
[[ $(git -C "$source" rev-parse --verify HEAD) == "$source_commit" ]] || {
    echo "E: Source commit differs from --source-commit." >&2
    exit 1
}
if ! git -C "$source" diff --quiet --ignore-submodules -- ||
   ! git -C "$source" diff --cached --quiet --ignore-submodules --; then
    echo "E: Interpretation requires a tracked-clean source tree." >&2
    exit 1
fi

analysis="$analysis_run/final/h3k4me3_cofactor_analysis"
enrichment="$enrichment_package/tables/jaspar2026/cofactor_primary_occupancy/part-000000.parquet"
intensity="$analysis/tables/intensity_effect.parquet"
context="$analysis/tables/context_stratified_intensity_effect.parquet"
gene_relation="$analysis/tables/gene_relation_stratified_intensity_effect.parquet"
gene_relation_occupancy="$analysis/tables/gene_relation_stratified_tp73_occupancy.parquet"
gradient="$analysis/tables/score_gradient.parquet"
interaction="$analysis/tables/tp73_interaction.parquet"
h3_summary="$h3_package/h3k4me3_change_summary_by_chromosome.parquet"
for path in "$analysis/manifest.json" "$enrichment" "$intensity" "$context" \
            "$gene_relation" "$gene_relation_occupancy" "$gradient" \
            "$interaction" "$h3_summary"; do
    [[ -s $path ]] || { echo "E: Finalized input is missing: $path" >&2; exit 1; }
done

if [[ -e $output_dir ]]; then
    [[ -s $output_dir/manifest.json ]] || {
        echo "E: Existing interpretation is incomplete: $output_dir" >&2
        exit 1
    }
    echo "I: Reusing compact interpretation: $output_dir" >&2
else
    "$source/scripts/summarize_h3k4me3_genome_cofactors.py" \
        --enrichment "$enrichment" --intensity "$intensity" \
        --context "$context" --gene-relation "$gene_relation" \
        --gene-relation-occupancy "$gene_relation_occupancy" \
        --score-gradient "$gradient" --interaction "$interaction" \
        --h3-summary "$h3_summary" --output-dir "$output_dir" \
        --duckdb "$duckdb"
fi

figures="$output_dir/figures"
expected=(
    h3k4me3_tp73_cofactor_effect_autosomes.png
    h3k4me3_tp73_context_correlations_autosomes.png
    h3k4me3_tp73_gene_relation_autosomes_promoter.png
    h3k4me3_tp73_gene_relation_autosomes_downstream.png
    h3k4me3_tp73_gene_relation_autosomes_gene_body.png
    h3k4me3_tp73_gene_relation_autosomes_intergenic.png
    context_correlations_autosomes.tsv
    gene_relation_correlations_autosomes.tsv
)
if [[ -d $figures ]]; then
    for name in "${expected[@]}"; do
        [[ -s $figures/$name ]] || {
            echo "E: Existing figure set is incomplete: $figures/$name" >&2
            exit 1
        }
    done
    echo "I: Reusing completed interpretation figures: $figures" >&2
    exit 0
fi

staging="$output_dir/.figures.staging-${SLURM_JOB_ID:-manual}-$$"
mkdir "$staging"
"$rscript" "$source/scripts/plot_h3k4me3_genome_cofactor_summary.R" \
    --joint "$output_dir/joint_primary_motif.tsv" \
    --context "$output_dir/context_primary_effect.tsv" \
    --gene-relation "$output_dir/gene_relation_primary_effect.tsv" \
    --output-effect "$staging/h3k4me3_tp73_cofactor_effect_autosomes.png" \
    --output-context "$staging/h3k4me3_tp73_context_correlations_autosomes.png" \
    --context-table "$staging/context_correlations_autosomes.tsv" \
    --output-gene-relation-prefix \
        "$staging/h3k4me3_tp73_gene_relation_autosomes" \
    --gene-relation-table "$staging/gene_relation_correlations_autosomes.tsv"
for name in "${expected[@]}"; do
    [[ -s $staging/$name ]] || {
        echo "E: Plotting omitted expected output: $staging/$name" >&2
        exit 1
    }
done
mv "$staging" "$figures"
echo "I: Published compact interpretation and figures: $output_dir" >&2

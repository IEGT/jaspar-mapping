#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-script-help.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

expected="$temporary/expected-scripts.txt"
actual="$temporary/actual-scripts.txt"

cat > "$expected" <<'EOF'
scripts/analyze_cutandrun_containment.py
scripts/analyze_dense_cutandrun_coverage.py
scripts/analyze_h3k4me3_cofactor_change.R
scripts/analyze_tp73_cofactor_enrichment.R
scripts/analyze_tp73_cofactor_score_surface.R
scripts/bedtools_map_serial.sh
scripts/benchmark_sparse_layout.py
scripts/build_density_capped_thresholds.py
scripts/build_fasta_index.py
scripts/build_h3k4me3_anchor_signal.py
scripts/build_motif_context.py
scripts/build_motif_score_thresholds.py
scripts/build_sparse_context_maxima.py
scripts/build_tp73_anchor_evidence.py
scripts/build_tp73_cutandrun_calibration_duckdb.sh
scripts/compare_tp73_cofactor_summaries.R
scripts/dense_tsv_to_parquet.sh
scripts/evaluate_tp73_cofactor_thresholds.R
scripts/evaluate_tp73_multifactor_model.R
scripts/evaluate_tp73_pair_stratified_model.R
scripts/export_bigwig_chrom_bedgraph.R
scripts/export_bigwig_chrom_bedgraph.py
scripts/export_dense_bed.py
scripts/fetch_region_to_embl.py
scripts/finalize_motif_context_run.py
scripts/fix_missing_bidirect.sh
scripts/fix_missing_bidirect_gz.sh
scripts/genelists.sh
scripts/inspect_chr1_dense_dry_run.sh
scripts/manage_genome_scan.py
scripts/manage_h3k4me3_cofactor_analysis.py
scripts/manage_h3k4me3_genome_signal.py
scripts/manage_motif_density_calibration.py
scripts/manage_motif_threshold_calibration.py
scripts/manage_tp73_cofactor_enrichment.py
scripts/manage_tp73_genome_context_maxima.py
scripts/manage_tp73_genome_evidence.py
scripts/plot_h3k4me3_genome_cofactor_summary.R
scripts/plot_h3k4me3_metaprofile.R
scripts/plot_informative_threshold_distribution.py
scripts/plot_tp73_score_distributions.R
scripts/query_genome_scan.py
scripts/run_chr1_2026_motif_panel.sh
scripts/run_chr1_patz1_tp73_dry_run.sh
scripts/run_genome_scan_slurm_chromosome.sh
scripts/run_genome_scan_slurm_finalize.sh
scripts/run_genome_scan_slurm_task.sh
scripts/run_h3k4me3_chr1_production.py
scripts/run_h3k4me3_cofactor_analysis_finalize.sh
scripts/run_h3k4me3_genome_cofactor_interpretation.sh
scripts/run_h3k4me3_genome_signal_finalize.sh
scripts/run_motif_context_slurm_task.sh
scripts/run_motif_density_calibration_finalize.sh
scripts/run_motif_density_calibration_slurm.sh
scripts/run_motif_threshold_anchor_setup.sh
scripts/run_motif_threshold_calibration_finalize.sh
scripts/run_motif_threshold_calibration_slurm_task.sh
scripts/run_negative_threshold_sensitivity_finalize.sh
scripts/run_negative_threshold_sensitivity_slurm_task.sh
scripts/run_tp73_chromosome_production.py
scripts/run_tp73_cofactor_enrichment_finalize.sh
scripts/run_tp73_cofactor_enrichment_setup.sh
scripts/run_tp73_cofactor_enrichment_slurm_batch.sh
scripts/run_tp73_cofactor_enrichment_slurm_task.sh
scripts/run_tp73_genome_context_maxima_batch.py
scripts/run_tp73_genome_context_maxima_finalize.sh
scripts/run_tp73_genome_evidence_finalize.sh
scripts/run_what_is_missing.sh
scripts/shift_bed.awk
scripts/stage_fasta_region.py
scripts/stage_motif_context_inputs.py
scripts/submit_genome_scan_slurm.sh
scripts/submit_h3k4me3_chr1_production_slurm.sh
scripts/submit_h3k4me3_cofactor_analysis_slurm.sh
scripts/submit_h3k4me3_genome_signal_slurm.sh
scripts/submit_motif_context_slurm.sh
scripts/submit_motif_density_calibration_slurm.sh
scripts/submit_motif_threshold_calibration_slurm.sh
scripts/submit_negative_threshold_sensitivity_slurm.sh
scripts/submit_tp73_cofactor_enrichment_slurm.sh
scripts/submit_tp73_genome_context_maxima_slurm.sh
scripts/submit_tp73_genome_evidence_slurm.sh
scripts/summarize_h3k4me3_genome_cofactors.py
scripts/summarize_tp73_cutandrun_threshold.R
scripts/summarize_tp73_patz1_cutandrun_threshold.R
EOF

(
    cd "$repository_root"
    find scripts -maxdepth 1 -type f -print | LC_ALL=C sort
) > "$actual"

if ! diff -u "$expected" "$actual"; then
    echo "E: Update the script-help inventory for every new or renamed CLI." >&2
    exit 1
fi

assert_help() {
    local relative_path=$1
    local script="$repository_root/$relative_path"
    local flag
    local output

    [[ -x $script ]] || {
        echo "E: User-facing script is not executable: $relative_path" >&2
        exit 1
    }

    if [[ $relative_path == *.R ]] && ! command -v Rscript >/dev/null 2>&1; then
        grep -q -- '--help' "$script" || {
            echo "E: Rscript is unavailable and $relative_path lacks static --help handling." >&2
            exit 1
        }
        echo "I: Rscript unavailable; statically checked $relative_path." >&2
        return
    fi

    for flag in -h --help; do
        if ! output=$(cd "$temporary" && "$script" "$flag" 2>&1); then
            echo "E: $relative_path $flag failed." >&2
            printf '%s\n' "$output" >&2
            exit 1
        fi
        if ! grep -Eiq '^usage:' <<< "$output"; then
            echo "E: $relative_path $flag did not print a Usage section." >&2
            printf '%s\n' "$output" >&2
            exit 1
        fi
    done
}

while IFS= read -r script; do
    assert_help "$script"
done < "$expected"

assert_help "OverlapTfPromoters/localMaxSkmelTADN.pl"

scanner_build=$("$repository_root/pssm_scan" --version-json)
python3 - "$scanner_build" <<'PY'
import json
import sys

build = json.loads(sys.argv[1])
required = {
    "program", "version", "source_commit", "source_dirty", "compiler",
    "cplusplus", "build_flags", "lto_enabled", "ndebug", "parquet_enabled",
    "arrow_version", "parquet_version",
}
assert required <= build.keys()
assert build["program"] == "pssm_scan"
assert isinstance(build["source_dirty"], bool)
assert "-O3" in build["build_flags"]
PY

panel_config=$("$repository_root/scripts/run_chr1_2026_motif_panel.sh" --print-config)
[[ $(grep -c $'^Motif\t' <<< "$panel_config") -eq 5 ]] || {
    echo "E: Chromosome 1 panel does not declare exactly five motifs." >&2
    exit 1
}
for motif in TP73:MA0861.2:16 E2F1:MA0024.3:12 SP1:MA0079.5:9 \
    PATZ1:MA1961.2:11 POU2F2:MA0507.3:13
do
    grep -Fq $'Motif\t'"$motif" <<< "$panel_config" || {
        echo "E: Chromosome 1 panel is missing $motif." >&2
        exit 1
    }
done

panel_help=$("$repository_root/scripts/run_chr1_2026_motif_panel.sh" --help)
grep -Fq -- '--source-commit SHA' <<< "$panel_help" || {
    echo "E: Chromosome 1 panel help omits its Git-free source mode." >&2
    exit 1
}

source_commit=$(git -C "$repository_root" rev-parse HEAD)
verified_commit=$(
    PATH=/nonexistent /bin/bash \
        "$repository_root/scripts/run_h3k4me3_genome_cofactor_interpretation.sh" \
        --source "$repository_root" --source-commit "$source_commit" \
        --verify-source-only
)
[[ $verified_commit == "$source_commit" ]] || {
    echo "E: Git-free H3K4me3 source verification returned another commit." >&2
    exit 1
}

shifted=$(printf '1\t100\t200\tplus\t0\t+\n1\t100\t200\tminus\t0\t-\n' |
    "$repository_root/scripts/shift_bed.awk")
expected_shifted=$'1\t600\t700\tplus\t0\t+\n1\t0\t0\tminus\t0\t-'
if [[ $shifted != "$expected_shifted" ]]; then
    echo "E: shift_bed.awk changed its interval transformation." >&2
    diff -u \
        <(printf '%s\n' "$expected_shifted") \
        <(printf '%s\n' "$shifted") >&2 || true
    exit 1
fi

unexpected=$(find "$temporary" -mindepth 1 -maxdepth 1 \
    ! -name expected-scripts.txt ! -name actual-scripts.txt -print -quit)
if [[ -n $unexpected ]]; then
    echo "E: A --help invocation created an unexpected file: $unexpected" >&2
    exit 1
fi

echo "Script help tests passed."

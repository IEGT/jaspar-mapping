#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
scanner_bed="$repository_root/pssm_scan"
scanner_parquet="$repository_root/pssm_scan_parquet"
fixture="$repository_root/test_files/synthetic_dense"
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-sparse-parity.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

python3 "$repository_root/scripts/build_fasta_index.py" \
    "$fixture/genome.fna" --output "$temporary/genome.fna.fai" >/dev/null

run_case() {
    local label=$1
    local score_mode=$2
    local pseudocount=$3
    local threshold=$4
    shift 4
    local bed_root="$temporary/$label/bed"
    local parquet_root="$temporary/$label/parquet/task_data/task_id=$label"
    local common=(
        --genome "$fixture/genome.fna"
        --fasta-index "$temporary/genome.fna.fai"
        --pssm "$fixture/JASPAR2026_synthetic.jaspar"
        --motif MA9999.1 --chr 1 --strand both --coordinate-mode bed
        --score-mode "$score_mode" --pseudocount "$pseudocount"
        --threshold "$threshold" --skip-N
        --motif-set-id synthetic_jaspar2026 --genome-id synthetic_acgtn_v1
        "$@"
    )

    "$scanner_bed" "${common[@]}" --outdir "$bed_root" >/dev/null
    "$scanner_parquet" "${common[@]}" --sparse-parquet \
        --outdir "$parquet_root" >/dev/null

    find "$bed_root" -name '*.bed' -type f -exec awk 'BEGIN { OFS="\t" }
        FNR > 1 { printf "%s\t%s\t%s\t%s\t%.3f\t%.6f\n",
            $1, $2, $3, $6, $5, $9 }' {} + \
        | LC_ALL=C sort > "$temporary/$label/bed.tsv"

    parquet_glob="$parquet_root/tables/jaspar2026/motif_hit/**/*.parquet"
    duckdb :memory: -csv -noheader -c "
        SELECT CAST(chrom AS VARCHAR), start, \"end\",
               CASE strand WHEN 'plus' THEN '+' WHEN 'minus' THEN '-' END,
               printf('%.3f', score), printf('%.6f', pwm_relative_score)
        FROM read_parquet('$parquet_glob', hive_partitioning=true)
        ORDER BY 1, 2, 3, 4, 5, 6;
    " | tr ',' '\t' > "$temporary/$label/parquet.tsv"

    if ! diff -u "$temporary/$label/bed.tsv" "$temporary/$label/parquet.tsv"; then
        echo "E: BED/Parquet parity failed for $label." >&2
        exit 1
    fi
}

# Inclusive equality at a known score of 1 exercises the float threshold edge.
run_case risk_equal log2_relative_risk 0 1

# Restricting coordinates and smoothing exercises the alternate score mode.
run_case odds_range log_odds 1 -10 --from 1 --to 4

# Both writers must publish valid, empty orientation files above every score.
run_case zero_hits log2_relative_risk 0 100
zero_files=$(find "$temporary/zero_hits/parquet" -name '*.parquet' -type f | wc -l)
[[ $zero_files -eq 2 ]] || {
    echo "E: Zero-hit sparse scan did not publish both orientation files." >&2
    exit 1
}

echo "Sparse BED/Parquet parity tests passed."

#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-script-help.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

expected="$temporary/expected-scripts.txt"
actual="$temporary/actual-scripts.txt"

cat > "$expected" <<'EOF'
scripts/bedtools_map_serial.sh
scripts/build_fasta_index.py
scripts/dense_tsv_to_parquet.sh
scripts/export_dense_bed.py
scripts/fetch_region_to_embl.py
scripts/fix_missing_bidirect.sh
scripts/fix_missing_bidirect_gz.sh
scripts/genelists.sh
scripts/inspect_chr1_dense_dry_run.sh
scripts/plot_tp73_score_distributions.R
scripts/run_chr1_patz1_tp73_dry_run.sh
scripts/run_what_is_missing.sh
scripts/shift_bed.awk
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

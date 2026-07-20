#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-context-inputs.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM
package="$temporary/scan_package"
staged="$temporary/staged"
mkdir -p "$package/tables/jaspar2026/scan_file_inventory"

write_hit_fixture() {
    local path=$1
    mkdir -p "$(dirname "$path")"
    duckdb :memory: -c "
        COPY (SELECT 10::BIGINT AS start, 20::BIGINT AS \"end\",
                     1.0::DOUBLE AS score)
        TO '$path' (FORMAT PARQUET);
    " >/dev/null
}

tp73_plus="task_data/task_id=batch-a/tables/hits/motif_id=MA0861.2/chrom=1/strand=plus/data.parquet"
tp73_minus="task_data/task_id=batch-a/tables/hits/motif_id=MA0861.2/chrom=1/strand=minus/data.parquet"
patz1_plus="task_data/task_id=batch-b/tables/hits/motif_id=MA1961.2/chrom=1/strand=plus/data.parquet"
patz1_minus="task_data/task_id=batch-b/tables/hits/motif_id=MA1961.2/chrom=1/strand=minus/data.parquet"
for relative in "$tp73_plus" "$tp73_minus" "$patz1_plus" "$patz1_minus"; do
    write_hit_fixture "$package/$relative"
done
tp73_plus_bytes=$(wc -c < "$package/$tp73_plus")
tp73_minus_bytes=$(wc -c < "$package/$tp73_minus")
patz1_plus_bytes=$(wc -c < "$package/$patz1_plus")
patz1_minus_bytes=$(wc -c < "$package/$patz1_minus")

(
    cd "$package"
    duckdb context.duckdb -c "
    CREATE TABLE inventory_source (
        task_id VARCHAR, output_relative_path VARCHAR, motif_id VARCHAR,
        chrom VARCHAR, strand VARCHAR, emitted_hits BIGINT, bytes BIGINT,
        sha256 VARCHAR
    );
    INSERT INTO inventory_source VALUES
        ('batch-a', 'tables/hits/motif_id=MA0861.2/chrom=1/strand=plus/data.parquet',
         'MA0861.2', '1', '+', 1, $tp73_plus_bytes, 'a'),
        ('batch-a', 'tables/hits/motif_id=MA0861.2/chrom=1/strand=minus/data.parquet',
         'MA0861.2', '1', '-', 1, $tp73_minus_bytes, 'b'),
        ('batch-b', 'tables/hits/motif_id=MA1961.2/chrom=1/strand=plus/data.parquet',
         'MA1961.2', '1', '+', 1, $patz1_plus_bytes, 'c'),
        ('batch-b', 'tables/hits/motif_id=MA1961.2/chrom=1/strand=minus/data.parquet',
         'MA1961.2', '1', '-', 1, $patz1_minus_bytes, 'd');
    COPY inventory_source TO
        'tables/jaspar2026/scan_file_inventory/data.parquet'
        (FORMAT PARQUET);
    DROP TABLE inventory_source;
    CREATE VIEW scan_file_inventory AS SELECT * FROM read_parquet(
        'tables/jaspar2026/scan_file_inventory/**/*.parquet'
    );
    " >/dev/null
)
printf '{"database":"context.duckdb"}\n' > "$package/manifest.json"

"$repository_root/scripts/stage_motif_context_inputs.py" \
    --package "$package" --output "$staged" \
    --motif MA0861.2 --motif MA1961.2 --chrom 1 >/dev/null

[[ $(find "$staged" -type f -name '*.parquet' | wc -l) -eq 4 ]]
[[ -f "$staged/input_manifest.json" ]]
[[ "$package/$tp73_plus" -ef "$staged/$tp73_plus" ]]
observed_motifs=$(duckdb -csv -noheader :memory: -c "
    SELECT string_agg(DISTINCT motif_id, ',' ORDER BY motif_id)
    FROM read_parquet('$staged/**/*.parquet', hive_partitioning=1);
")
[[ $observed_motifs == MA0861.2,MA1961.2 ]]

# A Hive-looking wrapper would override the genuine inner motif_id/chromosome
# labels when DuckDB reads the staged glob, so fail before creating it.
unsafe="$temporary/motif_id=MA1961.2/chrom=1"
if "$repository_root/scripts/stage_motif_context_inputs.py" \
    --package "$package" --output "$unsafe" \
    --motif MA0861.2 --motif MA1961.2 --chrom 1 >/dev/null 2>&1; then
    echo "E: Hive-conflicting context input wrapper was accepted" >&2
    exit 1
fi
[[ ! -e $unsafe ]]

# A completed identical selection is reusable; a changed selection is not.
"$repository_root/scripts/stage_motif_context_inputs.py" \
    --package "$package" --output "$staged" \
    --motif MA0861.2 --motif MA1961.2 --chrom 1 >/dev/null
if "$repository_root/scripts/stage_motif_context_inputs.py" \
    --package "$package" --output "$staged" \
    --motif MA0861.2 --chrom 1 >/dev/null 2>&1; then
    echo "E: mismatched existing context input was accepted" >&2
    exit 1
fi

echo "Motif-context input staging tests passed."

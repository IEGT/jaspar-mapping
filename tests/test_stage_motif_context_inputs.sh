#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-context-inputs.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM
package="$temporary/scan_package"
staged="$temporary/staged"
mkdir -p "$package/task_data/task_id=batch-a/tables/hits"
mkdir -p "$package/task_data/task_id=batch-b/tables/hits"

for strand in plus minus; do
    printf 'parquet\n' > "$package/task_data/task_id=batch-a/tables/hits/tp73-$strand.parquet"
    printf 'parquet\n' > "$package/task_data/task_id=batch-b/tables/hits/patz1-$strand.parquet"
done

duckdb "$package/context.duckdb" -c "
CREATE TABLE scan_file_inventory (
    task_id VARCHAR, output_relative_path VARCHAR, motif_id VARCHAR,
    chrom VARCHAR, strand VARCHAR, emitted_hits BIGINT, bytes BIGINT,
    sha256 VARCHAR
);
INSERT INTO scan_file_inventory VALUES
    ('batch-a', 'tables/hits/tp73-plus.parquet', 'MA0861.2', '1', '+', 2, 8, 'a'),
    ('batch-a', 'tables/hits/tp73-minus.parquet', 'MA0861.2', '1', '-', 2, 8, 'b'),
    ('batch-b', 'tables/hits/patz1-plus.parquet', 'MA1961.2', '1', '+', 2, 8, 'c'),
    ('batch-b', 'tables/hits/patz1-minus.parquet', 'MA1961.2', '1', '-', 2, 8, 'd');
" >/dev/null
printf '{"database":"context.duckdb"}\n' > "$package/manifest.json"

"$repository_root/scripts/stage_motif_context_inputs.py" \
    --package "$package" --output "$staged" \
    --motif MA0861.2 --motif MA1961.2 --chrom 1 >/dev/null

[[ $(find "$staged" -type f -name '*.parquet' | wc -l) -eq 4 ]]
[[ -f "$staged/input_manifest.json" ]]
[[ "$package/task_data/task_id=batch-a/tables/hits/tp73-plus.parquet" \
   -ef "$staged/task_data/task_id=batch-a/tables/hits/tp73-plus.parquet" ]]

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

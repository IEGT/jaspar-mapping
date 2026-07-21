#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-anchor-evidence.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping TP73 anchor-evidence test." >&2
    exit 0
}

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT * FROM (VALUES
        (10::BIGINT, 16::BIGINT, 1.0::FLOAT),
        (20::BIGINT, 26::BIGINT, 2.0::FLOAT),
        (30::BIGINT, 36::BIGINT, -0.5::FLOAT),
        (40::BIGINT, 46::BIGINT, 0.0::FLOAT)
    ) AS v(start, "end", score)
) TO '$temporary/plus.parquet' (FORMAT PARQUET);
COPY (
    SELECT * FROM (VALUES
        (10::BIGINT, 16::BIGINT, 3.0::FLOAT),
        (25::BIGINT, 31::BIGINT, 4.0::FLOAT)
    ) AS v(start, "end", score)
) TO '$temporary/minus.parquet' (FORMAT PARQUET);
SQL

printf '%s\n' \
    'track type=bedGraph' \
    $'chr1\t5\t12\t1' \
    $'chr1\t12\t17\t2' \
    $'chr1\t20\t32\t1' \
    $'chr2\t1\t100\t1' > "$temporary/anti.bedGraph"
printf '%s\n' \
    $'1\t19\t27\t1' \
    $'1\t39\t47\t1' > "$temporary/control.bedGraph"

"$repository_root/scripts/build_tp73_anchor_evidence.py" \
    --anchor-plus "$temporary/plus.parquet" \
    --anchor-minus "$temporary/minus.parquet" \
    --coverage supported_anti="$temporary/anti.bedGraph" \
    --coverage supported_control="$temporary/control.bedGraph" \
    --output "$temporary/anchors.parquet" --chrom 1 \
    --minimum-anchor-score 0 --duckdb "$duckdb"

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW anchors AS SELECT * FROM read_parquet('$temporary/anchors.parquet');
SELECT CASE WHEN (SELECT count(*) FROM anchors) <> 4
    THEN error('anchor score floor or orientation collapse is incorrect') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM anchors
    WHERE anchor_start = 10 AND anchor_end = 16 AND anchor_score = 3
      AND supported_anti AND NOT supported_control
) THEN error('adjacent coverage components did not immerse the anchor') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM anchors
    WHERE anchor_start = 20 AND NOT supported_anti AND supported_control
) THEN error('strict boundary equality was not enforced') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM anchors
    WHERE anchor_start = 25 AND supported_anti AND NOT supported_control
) THEN error('internal anchor was not immersed') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM anchors
    WHERE anchor_start = 40 AND NOT supported_anti AND supported_control
) THEN error('control immersion label is incorrect') END;
SQL

[[ -s $temporary/anchors.parquet.run_config.tsv &&
   -s $temporary/anchors.parquet.coverage_manifest.tsv ]] || {
    echo "E: Anchor evidence provenance sidecars are missing." >&2
    exit 1
}

echo "TP73 anchor-evidence tests passed."

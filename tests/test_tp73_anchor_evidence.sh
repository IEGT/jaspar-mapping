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

python3 - "$repository_root" <<'PY'
import importlib.util
import pathlib
import sys

path = pathlib.Path(sys.argv[1]) / "scripts" / "export_bigwig_chrom_bedgraph.py"
spec = importlib.util.spec_from_file_location("bigwig_export", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
assert module.normalize_chrom("chrM") == "MT"
assert module.normalize_chrom("MT") == "MT"
assert module.normalize_chrom("chr1") == "1"
assert module.chrom_matches("chr25", "MT", True)
assert not module.chrom_matches("chr25", "MT", False)
PY

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
      AND supported_anti AND depth_anti = 2
      AND NOT supported_control AND depth_control = 0
) THEN error('adjacent coverage components did not immerse the anchor') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM anchors
    WHERE anchor_start = 20 AND NOT supported_anti AND depth_anti = 0
      AND supported_control AND depth_control = 1
) THEN error('strict boundary equality was not enforced') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM anchors
    WHERE anchor_start = 25 AND supported_anti AND depth_anti = 1
      AND NOT supported_control AND depth_control = 0
) THEN error('internal anchor was not immersed') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM anchors
    WHERE anchor_start = 40 AND NOT supported_anti AND depth_anti = 0
      AND supported_control AND depth_control = 1
) THEN error('control immersion label is incorrect') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM anchors
    WHERE supported_anti <> (depth_anti > 0)
       OR supported_control <> (depth_control > 0)
) THEN error('support and effective depth disagree') END;
SQL

[[ -s $temporary/anchors.parquet.run_config.tsv &&
   -s $temporary/anchors.parquet.coverage_manifest.tsv ]] || {
    echo "E: Anchor evidence provenance sidecars are missing." >&2
    exit 1
}
grep -Fq $'effective_depth_rule\t' \
    "$temporary/anchors.parquet.run_config.tsv"
grep -Fq $'support_column\tdepth_column\t' \
    "$temporary/anchors.parquet.coverage_manifest.tsv"

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT * FROM (VALUES
        (10::BIGINT, 16::BIGINT, 'MA0861.2', '+', 3.0::FLOAT,
         'local_peak'),
        (10::BIGINT, 16::BIGINT, 'MA0861.2', '-', 2.0::FLOAT,
         'local_peak'),
        (30::BIGINT, 36::BIGINT, 'MA0861.2', '+', -0.5::FLOAT,
         'local_peak')
    ) AS v(start, "end", motif_id, strand, score,
           anchor_selection_class)
) TO '$temporary/context-anchor.parquet' (FORMAT PARQUET);
SQL
"$repository_root/scripts/build_tp73_anchor_evidence.py" \
    --anchor-source "$temporary/context-anchor.parquet" \
    --coverage supported_anti="$temporary/anti.bedGraph" \
    --output "$temporary/selected-anchors.parquet" --chrom 1 \
    --minimum-anchor-score -1 --duckdb "$duckdb"
"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW selected AS
SELECT * FROM read_parquet('$temporary/selected-anchors.parquet');
SELECT CASE WHEN (SELECT count(*) FROM selected) <> 2
    THEN error('schema-7 one-chromosome source was not collapsed by span') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM selected
    WHERE anchor_start = 10 AND anchor_score = 3 AND supported_anti
) THEN error('schema-7 tied orientations were not collapsed') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM selected
    WHERE anchor_start = 30 AND anchor_score = -0.5
) THEN error('negative local-peak anchor was lost') END;
SQL
grep -Fq $'schema7_local_peak_context_anchor\t' \
    "$temporary/selected-anchors.parquet.run_config.tsv"

: > "$temporary/empty.bedGraph"
"$repository_root/scripts/build_tp73_anchor_evidence.py" \
    --anchor-source "$temporary/context-anchor.parquet" \
    --coverage supported_empty="$temporary/empty.bedGraph" \
    --output "$temporary/empty-coverage-anchors.parquet" --chrom MT \
    --minimum-anchor-score -1 --allow-empty-coverage --duckdb "$duckdb"
"$duckdb" -batch :memory: >/dev/null <<SQL
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM read_parquet('$temporary/empty-coverage-anchors.parquet')
    WHERE supported_empty OR depth_empty <> 0
) THEN error('empty chromosome coverage did not become zero support') END;
SQL
grep -Fq $'zero_support\t' \
    "$temporary/empty-coverage-anchors.parquet.run_config.tsv"

if "$repository_root/scripts/build_tp73_anchor_evidence.py" \
    --anchor-source "$temporary/context-anchor.parquet" \
    --coverage supported_empty="$temporary/empty.bedGraph" \
    --output "$temporary/strict-empty.parquet" --chrom MT \
    --minimum-anchor-score -1 --duckdb "$duckdb" \
    > "$temporary/strict-empty.out" 2> "$temporary/strict-empty.err"; then
    echo "E: strict empty-coverage mode unexpectedly succeeded." >&2
    exit 1
fi
grep -Fq 'coverage has no positive rows' "$temporary/strict-empty.err"

printf '%s\n' $'chr25\t5\t17\t3' > "$temporary/mitochondrial.bedGraph"
"$repository_root/scripts/build_tp73_anchor_evidence.py" \
    --anchor-source "$temporary/context-anchor.parquet" \
    --coverage supported_mt="$temporary/mitochondrial.bedGraph" \
    --output "$temporary/mitochondrial-anchors.parquet" --chrom MT \
    --minimum-anchor-score -1 --mitochondrial-chromosome --duckdb "$duckdb"
"$duckdb" -batch :memory: >/dev/null <<SQL
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM read_parquet('$temporary/mitochondrial-anchors.parquet')
    WHERE anchor_start = 10 AND supported_mt AND depth_mt = 3
) THEN error('chr25 coverage did not match the run-scoped MT alias') END;
SQL

echo "TP73 anchor-evidence tests passed."

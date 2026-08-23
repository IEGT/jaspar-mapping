#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-fixed-threshold.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb -batch :memory: >/dev/null <<SQL
COPY (
  SELECT * FROM (VALUES
    ('source_v1', 'M1', 5.0::DOUBLE, -1.0::DOUBLE, 'calibrated'),
    ('source_v1', 'M2', 0.0::DOUBLE, -1.0::DOUBLE, 'calibrated'),
    ('another_set', 'M1', 9.0::DOUBLE, -1.0::DOUBLE, 'other')
  ) AS v(threshold_set_id, motif_id, recommended_threshold,
         source_minimum_score, calibration_status)
) TO '$temporary/source.parquet' (FORMAT PARQUET);
SQL

cat > "$temporary/motifs.tsv" <<'EOF'
motif_id	motif_name
M2	POU2F2
M1	POU2F1
EOF

"$repository_root/scripts/build_fixed_threshold_subset.py" \
    --source-registry "$temporary/source.parquet" \
    --source-threshold-set-id source_v1 \
    --motif-list "$temporary/motifs.tsv" --fixed-threshold 0 \
    --threshold-set-id pou_common_zero_v1 \
    --output-dir "$temporary/output"

duckdb -batch :memory: >/dev/null <<SQL
CREATE TEMP TABLE result AS
SELECT * FROM read_parquet('$temporary/output/motif_score_threshold.parquet');
SELECT CASE WHEN (SELECT count(*) FROM result) <> 2
  THEN error('fixed threshold subset has another cardinality') END;
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM result
  WHERE threshold_set_id <> 'pou_common_zero_v1'
     OR recommended_threshold <> 0 OR source_minimum_score <> -1
     OR calibration_status <> 'fixed_threshold_sensitivity'
     OR source_threshold_set_id <> 'source_v1'
) THEN error('fixed threshold subset changed its output contract') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM result
  WHERE motif_id='M1' AND motif_order=1 AND source_recommended_threshold=5
) THEN error('source threshold provenance or motif order was lost') END;
SQL

python3 - "$temporary/output/manifest.json" <<'PY'
import json
from pathlib import Path
import sys

manifest = json.loads(Path(sys.argv[1]).read_text())
assert manifest["threshold_set_id"] == "pou_common_zero_v1"
assert manifest["source_threshold_set_id"] == "source_v1"
assert manifest["fixed_threshold"] == 0
assert manifest["motif_count"] == 2
assert len(manifest["builder_sha256"]) == 64
assert len(manifest["output_sha256"]) == 64
PY

cat > "$temporary/missing.tsv" <<'EOF'
motif_id
MISSING
EOF
if "$repository_root/scripts/build_fixed_threshold_subset.py" \
    --source-registry "$temporary/source.parquet" \
    --source-threshold-set-id source_v1 \
    --motif-list "$temporary/missing.tsv" --fixed-threshold 0 \
    --threshold-set-id missing_v1 --output-dir "$temporary/missing-output" \
    >/dev/null 2>&1; then
    echo "E: Missing motif unexpectedly passed fixed-threshold validation." >&2
    exit 1
fi
[[ ! -e $temporary/missing-output ]] || {
    echo "E: Failed fixed-threshold build left a published output." >&2
    exit 1
}

echo "Fixed threshold subset tests passed."

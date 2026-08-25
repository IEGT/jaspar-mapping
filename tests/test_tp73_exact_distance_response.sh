#!/usr/bin/env bash

set -euo pipefail

if ! command -v duckdb >/dev/null 2>&1; then
    echo "E: duckdb CLI is required for this test." >&2
    exit 1
fi

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/tp73-exact-distance.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb -light-mode -batch -bail :memory: >/dev/null <<SQL
CREATE TABLE anchor_fixture AS
SELECT
    block_index,
    local_index,
    (500 + block_index * 5000 + local_index * 400)::BIGINT AS anchor_start,
    (516 + block_index * 5000 + local_index * 400)::BIGINT AS anchor_end,
    CASE
      WHEN block_index = 0 THEN local_index IN (0, 1, 2, 5)
      WHEN block_index = 1 THEN local_index IN (0, 1, 5)
      WHEN block_index = 2 THEN local_index IN (0, 1, 2, 5, 6)
      ELSE local_index IN (0, 1, 5, 6)
    END AS ta_anti
FROM range(4) block(block_index)
CROSS JOIN range(9) local(local_index);

COPY (
  SELECT
      '1'::VARCHAR AS chrom,
      anchor_start,
      anchor_end,
      3.0::DOUBLE AS anchor_score,
      ta_anti AS supported_tp73_saos2_TA,
      NOT ta_anti AS supported_negative_control_saos2_TA,
      NOT ta_anti AS supported_tp73_saos2_DN,
      ta_anti AS supported_negative_control_saos2_DN,
      ta_anti AS supported_tp73_skmel29_2_TA,
      NOT ta_anti AS supported_negative_control_skmel29_2_TA,
      NOT ta_anti AS supported_tp73_skmel29_2_DN,
      ta_anti AS supported_negative_control_skmel29_2_DN
  FROM anchor_fixture
) TO '$temporary/anchors.parquet' (FORMAT PARQUET);

COPY (
  SELECT anchor_start AS start, anchor_end AS "end", 3.0::DOUBLE AS score
  FROM anchor_fixture
  WHERE local_index % 2 = 0 OR local_index = 8
) TO '$temporary/anchor-plus.parquet' (FORMAT PARQUET);

COPY (
  SELECT anchor_start AS start, anchor_end AS "end", 3.0::DOUBLE AS score
  FROM anchor_fixture
  WHERE local_index % 2 = 1 OR local_index = 8
) TO '$temporary/anchor-minus.parquet' (FORMAT PARQUET);

COPY (
  SELECT anchor_start + 14 AS start, anchor_start + 24 AS "end",
         CASE WHEN local_index = 4 THEN 0.0 ELSE 3.0 END::DOUBLE AS score
  FROM anchor_fixture
  WHERE local_index IN (0, 2, 4)
) TO '$temporary/cofactor-plus.parquet' (FORMAT PARQUET);

COPY (
  SELECT anchor_start - 8 AS start, anchor_start + 2 AS "end",
         3.0::DOUBLE AS score
  FROM anchor_fixture
  WHERE local_index IN (1, 3)
  UNION ALL
  SELECT anchor_start + 14, anchor_start + 24, 1.0::DOUBLE
  FROM anchor_fixture
  WHERE local_index = 0
  UNION ALL
  SELECT anchor_start - 7, anchor_start + 2, 3.0::DOUBLE
  FROM anchor_fixture
  WHERE local_index = 5
) TO '$temporary/cofactor-minus.parquet' (FORMAT PARQUET);

COPY (
  SELECT * FROM read_parquet('$temporary/anchors.parquet')
  UNION ALL
  (SELECT * FROM read_parquet('$temporary/anchors.parquet') LIMIT 1)
) TO '$temporary/duplicate-anchors.parquet' (FORMAT PARQUET);
SQL

"$repository_root/scripts/build_tp73_exact_distance_counts.py" \
    --anchors "$temporary/anchors.parquet" \
    --anchor-plus "$temporary/anchor-plus.parquet" \
    --anchor-minus "$temporary/anchor-minus.parquet" \
    --cofactor-plus "$temporary/cofactor-plus.parquet" \
    --cofactor-minus "$temporary/cofactor-minus.parquet" \
    --motif-id MA9999.1 --motif-name SYNTHETIC --chrom 1 \
    --source-score-floor -1 --positive-threshold 2 \
    --inventory-output "$temporary/inventory.parquet" \
    --block-output "$temporary/block.parquet" \
    --block-size 5000 --memory-limit 1GB --max-temp-size 1GB \
    --temp-directory "$temporary/spill"

if "$repository_root/scripts/build_tp73_exact_distance_counts.py" \
    --anchors "$temporary/duplicate-anchors.parquet" \
    --anchor-plus "$temporary/anchor-plus.parquet" \
    --anchor-minus "$temporary/anchor-minus.parquet" \
    --cofactor-plus "$temporary/cofactor-plus.parquet" \
    --cofactor-minus "$temporary/cofactor-minus.parquet" \
    --motif-id MA9999.1 --motif-name SYNTHETIC --chrom 1 \
    --source-score-floor -1 --positive-threshold 2 \
    --inventory-output "$temporary/bad-inventory.parquet" \
    --block-output "$temporary/bad-block.parquet" \
    --block-size 5000 --memory-limit 1GB --max-temp-size 1GB \
    --temp-directory "$temporary/bad-spill" >/dev/null 2>&1
then
    echo "E: a duplicate TP73 anchor did not fail the SQL contract." >&2
    exit 1
fi

duckdb -light-mode -batch -bail :memory: >/dev/null <<SQL
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM read_parquet('$temporary/inventory.parquet')
  WHERE distance_frame = 'tp73_oriented' AND relative_orientation = 'all'
    AND center_offset_twice = 22 AND center_offset_bp = 11
    AND transaction_count = 32 AND anchors_source_present = 20
    AND itemset_count = 16 AND anchors_intermediate = 4
    AND anchors_negative = 12 AND physical_locus_pair_count = 20
) THEN error('the exact +11 bp oriented inventory is incorrect') END;

SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM read_parquet('$temporary/inventory.parquet')
  WHERE distance_frame = 'tp73_oriented' AND relative_orientation = 'same'
    AND center_offset_twice = 21 AND center_offset_bp = 10.5
    AND transaction_count = 32 AND itemset_count = 4
) THEN error('the half-base center offset was not preserved exactly') END;

SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM read_parquet('$temporary/inventory.parquet')
  WHERE distance_frame = 'tp73_oriented' AND relative_orientation = 'all'
    AND center_offset_twice = 24 AND anchors_source_present = 0
    AND anchors_negative = transaction_count
) THEN error('the exact-distance inventory is not zero-complete') END;

SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM read_parquet('$temporary/block.parquet')
  WHERE component_type = 'baseline' AND tp73_score_stratum = '[2,5)'
) THEN error('TP73 score strata were lost from block components') END;
SQL

"$repository_root/scripts/summarize_tp73_exact_distance_response.py" \
    --inventory "$temporary/inventory.parquet" \
    --block-components "$temporary/block.parquet" \
    --output-dir "$temporary/final" \
    --minimum-itemset-count 4 --minimum-blocks 2 \
    --minimum-periodicity-offsets 5 --memory-limit 1GB

duckdb -light-mode -batch -bail "$temporary/final/tp73_exact_distance.duckdb" \
    >/dev/null <<SQL
SELECT CASE WHEN abs((
  SELECT adjusted_odds_ratio
  FROM tp73_exact_distance_response
  WHERE motif_id = 'MA9999.1' AND distance_frame = 'tp73_oriented'
    AND relative_orientation = 'same' AND center_offset_twice = 22
    AND isoform = 'TA'
) - (5.0 / 3.0)) > 1e-9
THEN error('the TA exact-distance common odds ratio is incorrect') END;

SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM tp73_exact_distance_response
  WHERE motif_id = 'MA9999.1' AND distance_frame = 'tp73_oriented'
    AND relative_orientation = 'same' AND center_offset_twice = 22
    AND isoform = 'DN' AND evaluation_status = 'ok'
    AND adjusted_odds_ratio < 1
) THEN error('the DN exact-distance response is incorrect') END;

SELECT CASE WHEN (
  SELECT data_type FROM information_schema.columns
  WHERE table_name = 'tp73_exact_distance_response'
    AND column_name = 'center_offset_twice'
) <> 'BIGINT' THEN error('doubled-base coordinates lost their integer type') END;

SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM information_schema.columns
  WHERE table_name = 'tp73_exact_distance_peak_match'
    AND column_name = 'group_shift_hypothesis_bp'
) THEN error('empty derived tables did not retain their contract schema') END;
SQL

python3 - "$repository_root" <<'PY'
import importlib.util
import math
from pathlib import Path
import sys

module_path = Path(sys.argv[1]) / "scripts/summarize_tp73_exact_distance_response.py"
spec = importlib.util.spec_from_file_location("exact_distance", module_path)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)

responses = []
for isoform, phase in (("TA", 2.0), ("DN", 3.0)):
    for center in range(-30, 31):
        log_odds = 1.2 * math.cos(2.0 * math.pi * (center - phase) / 10.5)
        responses.append({
            "motif_id": "MA9999.1",
            "motif_name": "SYNTHETIC",
            "distance_frame": "tp73_oriented",
            "relative_orientation": "same",
            "isoform": isoform,
            "center_offset_bp": float(center),
            "adjusted_log_odds": log_odds,
            "adjusted_odds_ratio": math.exp(log_odds),
            "jackknife_se": 0.2,
            "q_value_bh_declared_offsets": 0.01,
            "itemset_count": 100,
            "evaluation_status": "ok",
        })

periodicity, contrast = module.periodicity_rows(responses, 10.5, 20)
assert len(periodicity) == 2
assert len(contrast) == 1
assert all(abs(row["harmonic_amplitude_log_odds"] - 1.2) < 1e-6
           for row in periodicity)

peaks = module.detect_peaks(responses, 1.0, 0.01, 3.0, 1)
for isoform in ("TA", "DN"):
    for direction in ("enriched", "depleted"):
        assert sum(row["isoform"] == isoform
                   and row["peak_direction"] == direction for row in peaks) <= 1

def peak(isoform, index, center):
    return {
        "motif_id": "MA9999.1",
        "motif_name": "SYNTHETIC",
        "distance_frame": "tp73_oriented",
        "relative_orientation": "same",
        "peak_direction": "enriched",
        "isoform": isoform,
        "peak_index": index,
        "peak_center_bp": center,
        "prominence_log_odds": 1.0,
    }

doublet = [
    peak("TA", 1, -42.0), peak("TA", 2, -31.0),
    peak("DN", 1, -31.0), peak("DN", 2, -20.0),
]
matches, pair_matches = module.match_peaks(doublet, 10.5, 16.0)
assert len(matches) == 2
assert all(row["group_shift_hypothesis_bp"] == 11.0 for row in matches)
assert all(row["dn_peak_is_closer_to_anchor"] for row in matches)
assert len(pair_matches) == 1
assert pair_matches[0]["ta_peak_spacing_bp"] == 11.0
assert pair_matches[0]["dn_peak_spacing_bp"] == 11.0
assert pair_matches[0]["dn_minus_ta_midpoint_shift_bp"] == 11.0
PY

echo "TP73 exact-distance response tests passed."

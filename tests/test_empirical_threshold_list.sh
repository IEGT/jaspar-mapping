#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)

python3 - "$repository_root" <<'PY'
import csv
import hashlib
import json
import pathlib
import sys

root = pathlib.Path(sys.argv[1])
matrix = root / "JASPAR2026_CORE_non-redundant_pfms_jaspar.txt"
thresholds = root / "thresholds" / "jaspar2026_grch38_chr1_tp73_context_v1.tsv"
metadata_path = thresholds.with_suffix(".json")

matrix_ids = []
with matrix.open(encoding="utf-8") as handle:
    for line in handle:
        if line.startswith(">"):
            matrix_ids.append(line[1:].split(maxsplit=1)[0])
assert len(matrix_ids) == len(set(matrix_ids)) == 2633

with thresholds.open(encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    assert reader.fieldnames == ["motif_id", "recommended_threshold"]
    rows = list(reader)

ids = [row["motif_id"] for row in rows]
assert len(ids) == len(set(ids)) == 2632
assert ids == sorted(ids)
assert set(ids) == set(matrix_ids) - {"MA0861.2"}

populated = {
    row["motif_id"]: int(row["recommended_threshold"])
    for row in rows
    if row["recommended_threshold"] != "NA"
}
assert len(populated) == 2615
assert len(rows) - len(populated) == 17
assert sum(row["recommended_threshold"] == "NA" for row in rows) == 17
assert min(populated.values()) == 0
assert max(populated.values()) == 14

expected_pilot = {
    "MA0024.3": 1,
    "MA0079.5": 3,
    "MA0138.3": 0,
    "MA0507.3": 0,
    "MA0740.2": 4,
    "MA0769.3": 6,
    "MA0790.2": 5,
    "MA0814.3": 5,
    "MA1961.2": 2,
}
assert {motif: populated[motif] for motif in expected_pilot} == expected_pilot

metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
digest = hashlib.sha256(thresholds.read_bytes()).hexdigest()
assert metadata["artifact"] == thresholds.name
assert metadata["columns"] == ["motif_id", "recommended_threshold"]
assert metadata["rows"] == len(rows)
assert metadata["populated_recommendations"] == len(populated)
assert metadata["null_recommendations"] == len(rows) - len(populated)
assert metadata["target_motif_id"] == "MA0861.2"
assert metadata["target_motif_included"] is False
assert metadata["threshold_tsv_sha256"] == digest
PY

echo "Empirical motif-threshold list tests passed."

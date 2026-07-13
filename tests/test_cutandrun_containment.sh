#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
fixture="$repository_root/test_files/synthetic_cutandrun"
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-cutandrun-containment.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

"$repository_root/scripts/analyze_cutandrun_containment.py" \
    --motifs "$fixture/tp73_motifs.bed" \
    --coverage-bed "$fixture/tp73_fragments.bed" \
    --output-dir "$temporary" \
    --chrom chr1 \
    --sample-id synthetic_tp73 \
    --score-mode synthetic_score >/dev/null

diff -u "$fixture/expected_motif_evidence.tsv" \
    "$temporary/motif_evidence.tsv"
diff -u "$fixture/expected_coverage_component_evidence.tsv" \
    "$temporary/coverage_component_evidence.tsv"
diff -u "$fixture/expected_threshold_curve.tsv" \
    "$temporary/threshold_curve.tsv"
diff -u "$fixture/expected_summary.json" "$temporary/summary.json"

"$repository_root/scripts/analyze_cutandrun_containment.py" \
    --motifs "$fixture/tp73_motifs.bed" \
    --coverage-bed "$fixture/tp73_coverage.bedGraph" \
    --output-dir "$temporary/bedgraph" \
    --chrom 1 \
    --sample-id synthetic_tp73_bedgraph \
    --score-mode synthetic_score >/dev/null

awk -F '\t' '
    NR == 1 { next }
    $5 == "TP73_A" { passed += ($8 == 3 && $9 == 3) }
    $5 == "TP73_B" { passed += ($8 == 2 && $9 == 2) }
    $5 == "TP73_C" { passed += ($8 == 4 && $9 == 4) }
    $5 == "TP73_D" { passed += ($8 == 5 && $9 == 0) }
    $5 == "TP73_E" { passed += ($8 == 6 && $9 == 0) }
    $5 == "TP73_F" { passed += ($8 == 7 && $9 == 7) }
    $5 == "TP73_G" { passed += ($8 == 9.5 && $9 == 9.5) }
    $5 == "TP73_H" { passed += ($8 == 0 && $9 == 0) }
    END { exit passed == 8 ? 0 : 1 }
' "$temporary/bedgraph/motif_evidence.tsv"
grep -Fq '"coverage_format": "bedgraph"' "$temporary/bedgraph/summary.json"
grep -Fq '"coverage_depth_semantics": "column_4_signal"' \
    "$temporary/bedgraph/summary.json"

if "$repository_root/scripts/analyze_cutandrun_containment.py" \
    --motifs "$fixture/tp73_motifs.bed" \
    --coverage-bed "$fixture/tp73_fragments.bed" \
    --output-dir "$temporary" >/dev/null 2>&1
then
    echo "E: Analyzer replaced existing outputs instead of refusing." >&2
    exit 1
fi

echo "Synthetic CUT&RUN containment test passed."

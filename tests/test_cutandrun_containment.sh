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

if "$repository_root/scripts/analyze_cutandrun_containment.py" \
    --motifs "$fixture/tp73_motifs.bed" \
    --coverage-bed "$fixture/tp73_fragments.bed" \
    --output-dir "$temporary" >/dev/null 2>&1
then
    echo "E: Analyzer replaced existing outputs instead of refusing." >&2
    exit 1
fi

echo "Synthetic CUT&RUN containment test passed."

#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-density-thresholds.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

cat > "$temporary/informative.tsv" <<'EOF'
motif_id	recommended_threshold
MOTIF_A	5
MOTIF_B	0
MOTIF_C	NA
EOF

cat > "$temporary/sensitivity.tsv" <<'EOF'
motif_id	best_negative_threshold	sensitivity_conclusion
MOTIF_B	-3	negative_threshold_higher_auc
EOF

cat > "$temporary/overrides.tsv" <<'EOF'
motif_id	informative_threshold	informative_source
TP73_TEST	-5	tp73_direct_calibration
EOF

make_distribution() {
    local motif_id=$1
    local counts=$2
    local output=$3
    {
        printf 'MotifID\tMotifName\tChromosome\tStrand\tScoreMode\tPseudocount\tBinScheme\tBinWidth\tValidWindows\tSkippedWindows\tMinScore\tMaxScore\tMeanScore\tScoreBinStart\tScoreBinEnd\tBinCount\tOrientationAggregation\n'
        local threshold count
        while IFS=: read -r threshold count; do
            printf '%s\t%s\t1\tboth\tlog2_relative_risk\t1.000000\tfixed\t1.000000\t1000\t0\t-5\t2\t0\t%s\t%s\t%s\tmax_score_per_alignment_span\n' \
                "$motif_id" "$motif_id" "$threshold" "$((threshold + 1))" "$count"
        done <<< "$counts"
    } > "$output"
}

# A starts from -1 because min(5,-1)=-1, then rises to 0 to meet <=5 loci.
make_distribution MOTIF_A $'-2:980\n-1:15\n0:3\n1:2' "$temporary/a.tsv"
# B keeps its informative -3 because exactly five loci survive there.
make_distribution MOTIF_B $'-4:995\n-3:2\n-2:1\n-1:1\n0:1' "$temporary/b.tsv"
# C has no informative threshold and uses the -1 fallback without density raising.
make_distribution MOTIF_C $'-2:998\n-1:2' "$temporary/c.tsv"
# The TP73 override begins at -5 but is density-limited to -4.
make_distribution TP73_TEST $'-5:995\n-4:2\n-3:2\n-2:1' "$temporary/tp73.tsv"

printf 'motif_id\tpath\tsha256\n' > "$temporary/manifest.tsv"
for entry in MOTIF_A:a MOTIF_B:b MOTIF_C:c TP73_TEST:tp73; do
    motif_id=${entry%%:*}
    stem=${entry#*:}
    path="$temporary/$stem.tsv"
    checksum=$(shasum -a 256 "$path" | awk '{print $1}')
    printf '%s\t%s\t%s\n' "$motif_id" "$path" "$checksum" \
        >> "$temporary/manifest.tsv"
done

"$repository_root/scripts/build_density_capped_thresholds.py" \
    --informative-thresholds "$temporary/informative.tsv" \
    --negative-sensitivity "$temporary/sensitivity.tsv" \
    --override-thresholds "$temporary/overrides.tsv" \
    --distribution-manifest "$temporary/manifest.tsv" \
    --output-tsv "$temporary/thresholds.tsv" \
    --output-json "$temporary/thresholds.json" \
    --threshold-set-id synthetic_density200_v1 \
    --genome-id synthetic_genome --motif-set-id synthetic_motifs \
    --minimum-spacing-bp 200 >/dev/null

python3 - "$temporary/thresholds.tsv" "$temporary/thresholds.json" <<'PY'
import csv
import json
import sys

with open(sys.argv[1], encoding="ascii") as stream:
    rows = {row["motif_id"]: row for row in csv.DictReader(stream, delimiter="\t")}
assert set(rows) == {"MOTIF_A", "MOTIF_B", "MOTIF_C", "TP73_TEST"}
assert rows["MOTIF_A"]["candidate_minimum_score"] == "-1"
assert rows["MOTIF_A"]["final_minimum_score"] == "0"
assert rows["MOTIF_A"]["density_limited"] == "true"
assert rows["MOTIF_A"]["loci_at_final_threshold"] == "5"
assert rows["MOTIF_B"]["informative_threshold"] == "-3"
assert rows["MOTIF_B"]["final_minimum_score"] == "-3"
assert rows["MOTIF_C"]["informative_threshold"] == "NA"
assert rows["MOTIF_C"]["final_minimum_score"] == "-1"
assert rows["TP73_TEST"]["candidate_minimum_score"] == "-5"
assert rows["TP73_TEST"]["final_minimum_score"] == "-4"
assert float(rows["TP73_TEST"]["mean_spacing_bp_at_final_threshold"]) == 200.0

metadata = json.load(open(sys.argv[2], encoding="ascii"))
assert metadata["motifs"] == 4
assert metadata["density_limited_motifs"] == 2
assert metadata["density_minimum_spacing_bp"] == 200.0
PY

echo "Density-capped threshold tests passed."

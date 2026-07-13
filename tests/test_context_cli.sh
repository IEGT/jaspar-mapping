#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-context-cli.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

cat > "$temporary/anchors.bed" <<'EOF'
Chromosome	From	To	Name	Score	Strand
1	100	116	TP73	5	+
1	200	216	TP73	6	-
EOF

cat > "$temporary/TP73_MA0861.2_both_1.bed" <<'EOF'
Chromosome	From	To	Name	Score	Strand
1	100	116	TP73	100	-
1	118	134	TP73	7	+
1	180	196	TP73	8	+
1	200	216	TP73	100	+
EOF

cat > "$temporary/SP1_MA0079.5_both_1.bed" <<'EOF'
Chromosome	From	To	Name	Score	Strand
1	80	89	SP1	4	-
1	115	124	SP1	3	+
1	220	229	SP1	2	-
EOF

"$repository_root/context" --flank 20 --batch-size 1 --temp-dir "$temporary" \
    "$temporary/anchors.bed" \
    "$temporary/TP73_MA0861.2_both_1.bed" \
    "$temporary/SP1_MA0079.5_both_1.bed" > "$temporary/context.tsv"

awk -F '\t' '
    NR == 1 {
        if ($7 != "ContextFlankBp" || $10 != "TP73_MA0861.2_Shift" ||
            $11 != "TP73_MA0861.2_GenomicShift" ||
            $14 != "TP73_MA0861.2_NumInWindow") exit 1
    }
    NR == 2 {
        if ($7 != 20 || $8 != "motif_center" || $9 != "exclude_same_named_locus" ||
            $10 != 18 || $11 != 18 || $12 != 7 || $13 != 1 || $14 != 1) exit 1
    }
    NR == 3 {
        if ($10 != 20 || $11 != -20 || $12 != 8 || $13 != 0 || $14 != 1) exit 1
    }
    END { if (NR != 3) exit 1 }
' "$temporary/context.tsv" || {
    echo "E: context CLI output did not preserve tandem distances/orientations." >&2
    cat "$temporary/context.tsv" >&2
    exit 1
}

"$repository_root/context" --help | grep -q '^Usage:'

cat > "$temporary/unsorted.bed" <<'EOF'
1	200	216	SP1	1	+
1	100	116	SP1	1	+
EOF
if "$repository_root/context" --temp-dir "$temporary" \
    "$temporary/anchors.bed" "$temporary/unsorted.bed" >/dev/null 2>&1; then
    echo "E: context accepted an unsorted secondary BED." >&2
    exit 1
fi

echo "Context CLI tests passed."

#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-h3k4me3.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping H3K4me3 cofactor-change test." >&2
    exit 0
}
command -v Rscript >/dev/null 2>&1 || {
    echo "I: Rscript unavailable; skipping H3K4me3 cofactor-change test." >&2
    exit 0
}
Rscript -e 'library(data.table)' >/dev/null 2>&1 || {
    echo "I: data.table unavailable; skipping H3K4me3 cofactor-change test." >&2
    exit 0
}

python3 - "$repository_root" "$temporary" <<'PY'
import errno
import importlib.util
import os
import pathlib
import sys

repository, temporary = map(pathlib.Path, sys.argv[1:])
spec = importlib.util.spec_from_file_location(
    "h3_signal", repository / "scripts" / "build_h3k4me3_anchor_signal.py"
)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
source = temporary / "cross-device-source.txt"
target = temporary / "cross-device-target.txt"
source.write_text("complete output\n")
real_replace = module.os.replace
raised = False

def forced_cross_device(left, right):
    global raised
    if not raised:
        raised = True
        raise OSError(errno.EXDEV, os.strerror(errno.EXDEV))
    return real_replace(left, right)

module.os.replace = forced_cross_device
module.promote_file(source, target)
assert target.read_text() == "complete output\n"
assert not list(temporary.glob(".cross-device-target.txt.tmp-*"))
PY

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT (100 + i * 500)::BIGINT AS start,
           (110 + i * 500)::BIGINT AS "end",
           (1 + (i % 7))::FLOAT AS score
    FROM range(16) AS r(i)
) TO '$temporary/plus.parquet' (FORMAT PARQUET);
COPY (
    SELECT (100 + i * 500)::BIGINT AS start,
           (110 + i * 500)::BIGINT AS "end",
           (1.5 + (i % 7))::FLOAT AS score
    FROM range(0, 16, 2) AS r(i)
) TO '$temporary/minus.parquet' (FORMAT PARQUET);
COPY (
    SELECT start, "end", 'MA0861.2'::VARCHAR AS motif_id,
           '+'::VARCHAR AS strand, score,
           'local_peak'::VARCHAR AS anchor_selection_class
    FROM read_parquet('$temporary/plus.parquet')
) TO '$temporary/context-anchor.parquet' (FORMAT PARQUET);
SQL

write_track() {
    local condition=$1
    local channel=$2
    local output=$3
    awk -v condition="$condition" -v channel="$channel" 'BEGIN {
        OFS="\t"
        if (channel == "negative_control") {
            print "1", 9000, 9010, 1
            exit
        }
        for (i = 0; i < 16; ++i) {
            start = 100 + i * 500
            end = start + 10
            if (channel == "tp73") {
                supported = (condition == "TA" && i % 2 == 0) ||
                            (condition == "DN" && i % 2 == 1) ||
                            (condition == "GFP" && i % 4 == 0)
                if (supported) print "1", start - 1, end + 1, 1 + (i % 3)
            } else {
                value = 1
                if (channel == "h3k4me3") {
                    value = 2
                    cofactor_positive = i % 4 < 2
                    if (condition == "TA" && cofactor_positive) value = 8
                    if (condition == "DN" && cofactor_positive) value = 1
                }
                print "1", start - 15, end + 15, value
            }
        }
    }' > "$output"
}

manifest="$temporary/tracks.tsv"
printf 'series_id\tcell_line\tcondition\tchannel\treplicate\tfilename\tanalysis_included\texclusion_reason\treference_build\tsignal_scale\n' > "$manifest"
for series in series_a series_b; do
    for condition in GFP TA DN; do
        for channel in h3k4me3 input tp73 negative_control; do
            filename="${series}_${condition}_${channel}.bedGraph"
            write_track "$condition" "$channel" "$temporary/$filename"
            printf '%s\t%s\t%s\t%s\tR1\t%s\ttrue\t\tGRCh38\tsynthetic_depth\n' \
                "$series" "$series" "$condition" "$channel" "$filename" \
                >> "$manifest"
        done
    done
done
for condition in GFP TA DN; do
    for channel in h3k4me3 input tp73 negative_control; do
        printf 'excluded_series\texcluded\t%s\t%s\tR1\tnot_present_%s_%s.bigWig\tfalse\tunresolved_test_series\tGRCh38\tsynthetic_depth\n' \
            "$condition" "$channel" "$condition" "$channel" >> "$manifest"
    done
done

awk 'BEGIN { FS = OFS = "\t" }
     NR == 2 { $2 = "inconsistent-cell-line" }
     { print }' "$manifest" > "$temporary/inconsistent-tracks.tsv"
if "$repository_root/scripts/build_h3k4me3_anchor_signal.py" \
    --anchor-plus "$temporary/plus.parquet" \
    --anchor-minus "$temporary/minus.parquet" \
    --track-manifest "$temporary/inconsistent-tracks.tsv" \
    --track-root "$temporary" --profile-only \
    --profile-output "$temporary/invalid-profile.tsv" \
    --chrom 1 --chrom-length 10000 --duckdb "$duckdb" \
    >"$temporary/inconsistent.stdout" 2>"$temporary/inconsistent.stderr"; then
    echo "E: inconsistent included-series metadata was accepted" >&2
    exit 1
fi
grep -Fq "included series series_a has inconsistent cell_line" \
    "$temporary/inconsistent.stderr"

"$repository_root/scripts/build_h3k4me3_anchor_signal.py" \
    --anchor-source "$temporary/context-anchor.parquet" \
    --track-manifest "$manifest" --track-root "$temporary" \
    --tp73-output "$temporary/tp73.parquet" \
    --signal-output "$temporary/signal.parquet" \
    --profile-output "$temporary/profile.tsv" \
    --window central_20:0:20 --chrom 1 --chrom-length 10000 \
    --profile-flank 40 --profile-bin-size 10 --profile-sample-size 16 \
    --threads 2 --memory-limit 1GB --duckdb "$duckdb"

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW signal AS SELECT * FROM read_parquet('$temporary/signal.parquet');
CREATE VIEW tp73 AS SELECT * FROM read_parquet('$temporary/tp73.parquet');
CREATE VIEW profile AS SELECT * FROM read_csv_auto(
    '$temporary/profile.tsv', delim='\t', header=true);
SELECT CASE WHEN (SELECT count(*) FROM signal) <> 192
    THEN error('signal output does not contain 16 anchors x 2 series x 3 conditions x 2 windows') END;
SELECT CASE WHEN EXISTS (SELECT 1 FROM signal WHERE series_id = 'excluded_series')
    THEN error('explicitly excluded series entered the signal output') END;
SELECT CASE WHEN (SELECT count(*) FROM profile) <> 32
    THEN error('GFP profile does not contain two channels x two series x eight bins') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM signal
    WHERE anchor_start = 100 AND series_id = 'series_a' AND condition = 'TA'
      AND window_name = 'central_20' AND h3k4me3_area = 320
      AND input_area = 40 AND effective_window_bp = 40
) THEN error('windowed integrated signal is incorrect') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM tp73
    WHERE anchor_start = 100
      AND supported_tp73_series_a_TA
      AND NOT supported_negative_control_series_a_TA
      AND depth_tp73_series_a_TA = 1
) THEN error('strict TP73 support was not retained') END;
SQL

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT '1'::VARCHAR AS chrom, (100 + i * 500)::BIGINT AS anchor_start,
           (110 + i * 500)::BIGINT AS anchor_end, motif_id,
           CASE WHEN i % 4 < 2 THEN 5.0 ELSE -2.0 END::FLOAT AS context_score,
           -2.0::DOUBLE AS source_score_floor,
           150::BIGINT AS context_flank_bp,
           'signed_interval_edge_distance'::VARCHAR AS context_distance_metric
    FROM range(16) AS r(i)
    CROSS JOIN (VALUES ('M_EFFECT'), ('M_OVERLAP')) AS m(motif_id)
) TO '$temporary/maxima.parquet' (FORMAT PARQUET);
SQL
printf 'motif_id\tpositive_threshold\tfactor_name\tpositive_threshold_source\tselection_semantics\nM_EFFECT\t4\tSynthetic effect\ttest_fixture\tprespecified\n' \
    > "$temporary/thresholds.tsv"
printf 'M_OVERLAP\t-1\tSynthetic overlap guard\ttest_fixture\tprespecified\n' \
    >> "$temporary/thresholds.tsv"

Rscript "$repository_root/scripts/analyze_h3k4me3_cofactor_change.R" \
    --signal "$temporary/signal.parquet" \
    --tp73-evidence "$temporary/tp73.parquet" \
    --cofactor-maxima "$temporary/maxima.parquet" \
    --thresholds "$temporary/thresholds.tsv" \
    --window central_20 --output-prefix "$temporary/result" \
    --series series_a --series series_b \
    --negative-references "-1,0" --pseudocount 1 --block-size 500 \
    --spline-df 1 --minimum-class-fraction 0.01 --minimum-class-count 2 \
    --minimum-interaction-cell-count 2 --duckdb "$duckdb"

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW intensity AS SELECT * FROM read_csv_auto(
    '$temporary/result_intensity_effect.tsv', delim='\t', header=true, nullstr='NA');
CREATE VIEW interaction AS SELECT * FROM read_csv_auto(
    '$temporary/result_tp73_interaction.tsv', delim='\t', header=true, nullstr='NA');
CREATE VIEW summary AS SELECT * FROM read_csv_auto(
    '$temporary/result_series_summary.tsv', delim='\t', header=true, nullstr='NA');
CREATE VIEW state_summary AS SELECT * FROM read_csv_auto(
    '$temporary/result_binding_state_summary.tsv', delim='\t', header=true, nullstr='NA');
CREATE VIEW occurrence_summary AS SELECT * FROM read_csv_auto(
    '$temporary/result_occurrence_summary.tsv', delim='\t', header=true, nullstr='NA');
CREATE VIEW config AS SELECT * FROM read_csv_auto(
    '$temporary/result_run_config.tsv', delim='\t', header=true, nullstr='NA');
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM intensity
    WHERE motif_id = 'M_EFFECT' AND isoform = 'TA'
      AND negative_reference_threshold = -1 AND evaluation_status = 'ok'
      AND try_cast(estimate AS DOUBLE) > 0
) THEN error('known positive TA-minus-GFP cofactor effect was not recovered') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM summary
    WHERE motif_id = 'M_EFFECT' AND isoform = 'TA'
      AND direction_consistency = 'all_series_positive'
) THEN error('two-series directional consistency is wrong') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM interaction
    WHERE contrast = 'cofactor_by_tp73_confirmation_interaction'
      AND evaluation_status = 'ok'
) THEN error('supported TP73 interaction was not estimated') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM intensity
    WHERE motif_id = 'M_OVERLAP' AND negative_reference_threshold = 0
      AND evaluation_status = 'overlapping_score_classes'
) THEN error('overlapping score classes were not rejected') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM state_summary
    WHERE motif_id = 'M_OVERLAP' AND negative_reference_threshold = 0
      AND comparison_status = 'overlapping_score_classes'
      AND cofactor_present IS NULL
) THEN error('binding-state output does not flag an invalid comparison') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM occurrence_summary
    WHERE motif_id = 'M_OVERLAP' AND negative_reference_threshold = 0
      AND comparison_status = 'overlapping_score_classes'
      AND cofactor_present IS NULL
) THEN error('occurrence output does not flag an invalid comparison') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM config
    WHERE all_zero_anchor_policy = 'retained'
      AND tp73_interaction_interpretation = 'secondary_descriptive_post_treatment'
) THEN error('analysis provenance omits the agreed estimand safeguards') END;
SQL

grep -Fq $'excluded_series\texcluded\t' \
    "$temporary/signal.parquet.track_manifest.tsv"
grep -Fq $'excluded_series\texcluded\tGFP\th3k4me3\tR1\tfalse\tfalse\tunresolved_test_series' \
    "$temporary/signal.parquet.track_manifest.tsv"
grep -Fq '"anchor_source_mode": "schema7_local_peak_context_anchor"' \
    "$temporary/signal.parquet.run_config.json"

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT * REPLACE (0.0::DOUBLE AS source_score_floor)
    FROM read_parquet('$temporary/maxima.parquet')
    WHERE motif_id = 'M_EFFECT'
) TO '$temporary/invalid-floor-maxima.parquet' (FORMAT PARQUET);
SQL
printf 'motif_id\tpositive_threshold\tfactor_name\nM_EFFECT\t4\tSynthetic effect\n' \
    > "$temporary/invalid-floor-thresholds.tsv"
if Rscript "$repository_root/scripts/analyze_h3k4me3_cofactor_change.R" \
    --signal "$temporary/signal.parquet" \
    --tp73-evidence "$temporary/tp73.parquet" \
    --cofactor-maxima "$temporary/invalid-floor-maxima.parquet" \
    --thresholds "$temporary/invalid-floor-thresholds.tsv" \
    --window central_20 --output-prefix "$temporary/invalid-floor-result" \
    --series series_a --series series_b --negative-references "-1,0" \
    --minimum-class-count 2 --minimum-interaction-cell-count 2 \
    --duckdb "$duckdb" \
    >"$temporary/invalid-floor.stdout" 2>"$temporary/invalid-floor.stderr"; then
    echo "E: an unobservable negative class was accepted" >&2
    exit 1
fi
grep -Fq "cofactor scan floor is above a requested negative reference" \
    "$temporary/invalid-floor.stderr"

echo "I: H3K4me3 cofactor-change synthetic test passed."

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
    --tp73-output "$temporary/evidence-only.parquet" --evidence-only \
    --chrom 1 --chrom-length 10000 --threads 2 --memory-limit 1GB \
    --duckdb "$duckdb"
"$duckdb" -batch :memory: >/dev/null <<SQL
SELECT CASE WHEN (SELECT count(*) FROM read_parquet(
    '$temporary/evidence-only.parquet')) <> 16
THEN error('evidence-only mode changed the anchor cardinality') END;
SQL
[[ ! -e $temporary/evidence-only-signal.parquet ]]
grep -Fq '"mode": "evidence_only"' \
    "$temporary/evidence-only.parquet.run_config.json"
awk -F '\t' 'NR > 1 && $4 == "h3k4me3" && $7 != "false" { exit 1 }' \
    "$temporary/evidence-only.parquet.track_manifest.tsv"

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

# Reusing completed TP73 evidence must stage only H3K4me3/input tracks and
# retain a complete zero-filled factorial grid when one chromosome track is
# empty.
: > "$temporary/series_b_DN_h3k4me3_empty.bedGraph"
awk 'BEGIN { FS = OFS = "\t" }
     NR == 1 || !($1 == "series_b" && $3 == "DN" && $4 == "h3k4me3") { print; next }
     { $6 = "series_b_DN_h3k4me3_empty.bedGraph"; print }' \
    "$manifest" > "$temporary/signal-only-tracks.tsv"
python3 - "$temporary/signal-only-tracks.tsv" "$temporary" \
    "$temporary/signal-track-files.tsv" <<'PY'
import csv
import hashlib
import pathlib
import sys

manifest, root, output = map(pathlib.Path, sys.argv[1:])
with manifest.open(newline="") as stream:
    rows = [row for row in csv.DictReader(stream, delimiter="\t")
            if row["analysis_included"] == "true"
            and row["channel"] in {"h3k4me3", "input"}]
fields = [
    "series_id", "condition", "channel", "filename", "resolved_source",
    "bytes", "mtime_ns", "sha256",
]
with output.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for row in rows:
        source = (root / row["filename"]).resolve()
        stat = source.stat()
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        writer.writerow({
            **{field: row[field] for field in (
                "series_id", "condition", "channel", "filename"
            )},
            "resolved_source": source,
            "bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "sha256": digest,
        })
PY
"$repository_root/scripts/build_h3k4me3_anchor_signal.py" \
    --tp73-evidence-input "$temporary/evidence-only.parquet" \
    --track-manifest "$temporary/signal-only-tracks.tsv" \
    --track-file-inventory "$temporary/signal-track-files.tsv" \
    --track-root "$temporary" \
    --signal-output "$temporary/signal-only.parquet" \
    --window central_20:0:20 --chrom 1 --chrom-length 10000 \
    --threads 2 --memory-limit 1GB --duckdb "$duckdb"
"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW reused AS SELECT * FROM read_parquet(
    '$temporary/signal-only.parquet');
SELECT CASE WHEN (SELECT count(*) FROM reused) <> 192
    THEN error('signal-only output changed the factorial row count') END;
SELECT CASE WHEN EXISTS (
    SELECT 1 FROM reused
    WHERE series_id = 'series_b' AND condition = 'DN'
      AND h3k4me3_area <> 0
) THEN error('empty H3K4me3 chromosome track was not zero-filled') END;
SELECT CASE WHEN (
    SELECT count(DISTINCT (chrom, anchor_start, anchor_end)) FROM reused
) <> 16 THEN error('signal-only mode changed the evidence anchor set') END;
SQL
grep -Fq '"mode": "signal_only"' \
    "$temporary/signal-only.parquet.run_config.json"
grep -Fq '"anchor_source_mode": "prebuilt_tp73_cutandrun_anchor_evidence"' \
    "$temporary/signal-only.parquet.run_config.json"
grep -Fq '"track_file_inventory_sha256":' \
    "$temporary/signal-only.parquet.run_config.json"
awk -F '\t' 'NR > 1 && $4 ~ /^(tp73|negative_control)$/ && $7 != "false" { exit 1 }' \
    "$temporary/signal-only.parquet.track_manifest.tsv"
python3 - "$temporary/signal-track-files.tsv" \
    "$temporary/signal-only.parquet.track_manifest.tsv" <<'PY'
import csv
import pathlib
import sys

expected_path, observed_path = map(pathlib.Path, sys.argv[1:])
with expected_path.open(newline="") as stream:
    expected = {
        (row["series_id"], row["condition"], row["channel"]): row["sha256"]
        for row in csv.DictReader(stream, delimiter="\t")
    }
with observed_path.open(newline="") as stream:
    observed = {
        (row["series_id"], row["condition"], row["channel"]): row["sha256"]
        for row in csv.DictReader(stream, delimiter="\t")
        if row["used_in_run"] == "true"
    }
assert observed == expected
PY
python3 - "$temporary/signal-track-files.tsv" <<'PY'
import csv
import os
import pathlib
import sys

with pathlib.Path(sys.argv[1]).open(newline="") as stream:
    source = pathlib.Path(next(csv.DictReader(stream, delimiter="\t"))["resolved_source"])
stat = source.stat()
os.utime(source, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
PY
if "$repository_root/scripts/build_h3k4me3_anchor_signal.py" \
    --tp73-evidence-input "$temporary/evidence-only.parquet" \
    --track-manifest "$temporary/signal-only-tracks.tsv" \
    --track-file-inventory "$temporary/signal-track-files.tsv" \
    --track-root "$temporary" \
    --signal-output "$temporary/stale-inventory-signal.parquet" \
    --window central_20:0:20 --chrom 1 --chrom-length 10000 \
    --threads 2 --memory-limit 1GB --duckdb "$duckdb" \
    >"$temporary/stale-inventory.stdout" \
    2>"$temporary/stale-inventory.stderr"; then
    echo "E: Changed source track was accepted by the pinned inventory." >&2
    exit 1
fi
grep -Fq "track differs from pinned file inventory" \
    "$temporary/stale-inventory.stderr"

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
    --analysis-role synthetic_held_out_validation \
    --spline-df 1 --minimum-class-fraction 0.01 --minimum-class-count 2 \
    --minimum-interaction-cell-count 2 --duckdb "$duckdb"

"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW intensity AS SELECT * FROM read_csv_auto(
    '$temporary/result_intensity_effect.tsv', delim='\t', header=true, nullstr='NA');
CREATE VIEW interaction AS SELECT * FROM read_csv_auto(
    '$temporary/result_tp73_interaction.tsv', delim='\t', header=true, nullstr='NA');
CREATE VIEW isoform_contrast AS SELECT * FROM read_csv_auto(
    '$temporary/result_isoform_contrast.tsv', delim='\t', header=true, nullstr='NA');
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
    SELECT 1 FROM intensity
    WHERE motif_id = 'M_EFFECT' AND isoform = 'TA'
      AND negative_reference_threshold = -1 AND evaluation_status = 'ok'
      AND adjusted_mean_change_positive IS NOT NULL
      AND adjusted_mean_change_negative IS NOT NULL
) THEN error('adjusted positive/negative margins were not recovered') END;
SELECT CASE WHEN NOT EXISTS (
    SELECT 1 FROM isoform_contrast
    WHERE motif_id = 'M_EFFECT' AND contrast = 'TA_minus_DN'
      AND negative_reference_threshold = -1 AND evaluation_status = 'ok'
) THEN error('paired anchor-level TA-minus-DN contrast was not recovered') END;
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
    WHERE all_zero_anchor_policy =
      'excluded_from_intensity_and_score_gradient_retained_for_occurrence'
      AND tp73_interaction_interpretation = 'secondary_descriptive_post_treatment'
      AND chromosomes = '1'
      AND analysis_role = 'synthetic_held_out_validation'
) THEN error('analysis provenance omits the agreed estimand safeguards') END;
SQL

grep -Fq $'excluded_series\texcluded\t' \
    "$temporary/signal.parquet.track_manifest.tsv"
grep -Fq $'excluded_series\texcluded\tGFP\th3k4me3\tR1\tfalse\tfalse\tunresolved_test_series' \
    "$temporary/signal.parquet.track_manifest.tsv"
grep -Fq '"anchor_source_mode": "schema_local_peak_context_anchor"' \
    "$temporary/signal.parquet.run_config.json"

# The genome evaluator consumes the compact change table and collapses
# orientation-duplicated schema-9 annotation only when physical-span fields
# agree.
"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
  WITH selected AS (
    SELECT * FROM read_parquet('$temporary/signal.parquet')
    WHERE window_name = 'central_20'
  ), pivoted AS (
    SELECT chrom, anchor_start, anchor_end, anchor_score, series_id, cell_line,
           replicate, window_name, min(inner_bp)::BIGINT AS inner_bp,
           max(outer_bp)::BIGINT AS outer_bp,
           max(segment_count)::INTEGER AS segment_count,
           max(effective_window_bp)::BIGINT AS effective_window_bp,
           max(h3k4me3_area) FILTER (condition = 'GFP') AS gfp_h3k4me3_area,
           max(input_area) FILTER (condition = 'GFP') AS gfp_input_area,
           max(h3k4me3_max) FILTER (condition = 'GFP') AS gfp_h3k4me3_max,
           max(input_max) FILTER (condition = 'GFP') AS gfp_input_max,
           max(h3k4me3_area) FILTER (condition = 'TA') AS ta_h3k4me3_area,
           max(input_area) FILTER (condition = 'TA') AS ta_input_area,
           max(h3k4me3_max) FILTER (condition = 'TA') AS ta_h3k4me3_max,
           max(input_max) FILTER (condition = 'TA') AS ta_input_max,
           max(h3k4me3_area) FILTER (condition = 'DN') AS dn_h3k4me3_area,
           max(input_area) FILTER (condition = 'DN') AS dn_input_area,
           max(h3k4me3_max) FILTER (condition = 'DN') AS dn_h3k4me3_max,
           max(input_max) FILTER (condition = 'DN') AS dn_input_max
    FROM selected
    GROUP BY chrom, anchor_start, anchor_end, anchor_score, series_id,
             cell_line, replicate, window_name
  ), normalized AS (
    SELECT *, log2((gfp_h3k4me3_area + 1) / (gfp_input_area + 1))
               AS gfp_log2_h3k4me3_input_ratio,
           log2((ta_h3k4me3_area + 1) / (ta_input_area + 1))
               AS ta_log2_h3k4me3_input_ratio,
           log2((dn_h3k4me3_area + 1) / (dn_input_area + 1))
               AS dn_log2_h3k4me3_input_ratio
    FROM pivoted
  )
  SELECT *,
         ta_log2_h3k4me3_input_ratio - gfp_log2_h3k4me3_input_ratio
             AS delta_ta_vs_gfp,
         dn_log2_h3k4me3_input_ratio - gfp_log2_h3k4me3_input_ratio
             AS delta_dn_vs_gfp,
         true AS has_any_h3k4me3_signal, true AS has_any_input_signal,
         false AS uninformative_all_h3k4me3_zero
  FROM normalized
) TO '$temporary/change.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);

COPY (
  SELECT '1'::VARCHAR AS chrom, (100 + i * 500)::BIGINT AS start,
         (110 + i * 500)::BIGINT AS "end", strand,
         CASE floor(i / 4)
           WHEN 0 THEN 'promoter_only'
           WHEN 1 THEN 'downstream_only'
           WHEN 2 THEN 'intron'
           ELSE 'strict_intergenic'
         END AS primary_genomic_context,
         floor(i / 4) = 3 AS strict_intergenic,
         floor(i / 4) = 0 AS overlaps_any_promoter,
         floor(i / 4) = 1 AS overlaps_any_downstream_region,
         floor(i / 4) = 2 AS in_any_transcript,
         CASE floor(i / 4)
           WHEN 0 THEN 'promoter'
           WHEN 1 THEN 'downstream'
           WHEN 2 THEN 'gene_body'
           ELSE 'intergenic'
         END AS gene_relation_class,
         false AS in_any_exon, false AS in_any_cds,
         (1000 + i * 100)::BIGINT AS nearest_tss_distance_bp,
         (1000 + i * 100)::BIGINT AS nearest_tss_genomic_distance_bp,
         i = 0 AS nearest_tss_has_mixed_strands,
         'downstream'::VARCHAR AS nearest_tss_relation,
         (2000 + i * 100)::BIGINT AS nearest_cds_distance_bp,
         (2000 + i * 100)::BIGINT AS nearest_cds_genomic_distance_bp,
         false AS nearest_cds_has_mixed_strands,
         'downstream'::VARCHAR AS nearest_cds_relation
  FROM range(16) AS r(i)
  CROSS JOIN (VALUES ('+'), ('-')) AS s(strand)
) TO '$temporary/schema9-annotation.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);

COPY (SELECT * FROM read_parquet('$temporary/change.parquet')
      WHERE anchor_start % 1000 = 100)
TO '$temporary/change-a.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * FROM read_parquet('$temporary/change.parquet')
      WHERE anchor_start % 1000 = 600)
TO '$temporary/change-b.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * FROM read_parquet('$temporary/schema9-annotation.parquet')
      WHERE start % 1000 = 100)
TO '$temporary/annotation-a.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * FROM read_parquet('$temporary/schema9-annotation.parquet')
      WHERE start % 1000 = 600)
TO '$temporary/annotation-b.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * FROM read_parquet('$temporary/maxima.parquet')
      WHERE motif_id = 'M_EFFECT')
TO '$temporary/maxima-effect.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * FROM read_parquet('$temporary/maxima.parquet')
      WHERE motif_id = 'M_OVERLAP')
TO '$temporary/maxima-overlap.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (
  SELECT * REPLACE (
    NOT supported_tp73_series_a_GFP
      AS supported_negative_control_series_a_GFP,
    (NOT supported_tp73_series_a_GFP)::INTEGER
      AS depth_negative_control_series_a_GFP,
    NOT supported_tp73_series_a_TA
      AS supported_negative_control_series_a_TA,
    (NOT supported_tp73_series_a_TA)::INTEGER
      AS depth_negative_control_series_a_TA,
    NOT supported_tp73_series_a_DN
      AS supported_negative_control_series_a_DN,
    (NOT supported_tp73_series_a_DN)::INTEGER
      AS depth_negative_control_series_a_DN,
    NOT supported_tp73_series_b_GFP
      AS supported_negative_control_series_b_GFP,
    (NOT supported_tp73_series_b_GFP)::INTEGER
      AS depth_negative_control_series_b_GFP,
    NOT supported_tp73_series_b_TA
      AS supported_negative_control_series_b_TA,
    (NOT supported_tp73_series_b_TA)::INTEGER
      AS depth_negative_control_series_b_TA,
    NOT supported_tp73_series_b_DN
      AS supported_negative_control_series_b_DN,
    (NOT supported_tp73_series_b_DN)::INTEGER
      AS depth_negative_control_series_b_DN
  )
  FROM read_parquet('$temporary/tp73.parquet')
) TO '$temporary/tp73-matched.parquet' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL

Rscript "$repository_root/scripts/analyze_h3k4me3_cofactor_change.R" \
    --change "$temporary/change-a.parquet" \
    --change "$temporary/change-b.parquet" \
    --tp73-evidence "$temporary/tp73-matched.parquet" \
    --cofactor-maxima "$temporary/maxima-effect.parquet" \
    --cofactor-maxima "$temporary/maxima-overlap.parquet" \
    --annotation "$temporary/annotation-a.parquet" \
    --annotation "$temporary/annotation-b.parquet" \
    --thresholds "$temporary/thresholds.tsv" \
    --output-prefix "$temporary/context-result" \
    --series series_a --series series_b --negative-references "-1,0" \
    --block-size 500 --spline-df 1 --minimum-class-fraction 0.01 \
    --minimum-class-count 2 --minimum-interaction-cell-count 2 \
    --duckdb "$duckdb"
"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW context_effect AS SELECT * FROM read_csv_auto(
  '$temporary/context-result_context_stratified_intensity_effect.tsv',
  delim='\t', header=true, nullstr='NA');
CREATE VIEW gene_relation_effect AS SELECT * FROM read_csv_auto(
  '$temporary/context-result_gene_relation_stratified_intensity_effect.tsv',
  delim='\t', header=true, nullstr='NA');
CREATE VIEW gene_relation_occupancy AS SELECT * FROM read_csv_auto(
  '$temporary/context-result_gene_relation_stratified_tp73_occupancy.tsv',
  delim='\t', header=true, nullstr='NA');
CREATE VIEW score_gradient AS SELECT * FROM read_csv_auto(
  '$temporary/context-result_score_gradient.tsv', delim='\t',
  header=true, nullstr='NA');
CREATE VIEW context_config AS SELECT * FROM read_csv_auto(
  '$temporary/context-result_run_config.tsv', delim='\t', header=true,
  nullstr='NA');
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM context_effect WHERE genomic_context_class = 'strict_intergenic'
) THEN error('strict-intergenic context stratum was not emitted') END;
SELECT CASE WHEN (SELECT count(DISTINCT gene_relation_class)
                  FROM gene_relation_effect) <> 4
  OR (SELECT count(*) FROM gene_relation_effect) <> 64
THEN error('four-way gene-relation strata were not emitted') END;
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM gene_relation_effect
  WHERE gene_relation_class NOT IN
    ('promoter', 'downstream', 'gene_body', 'intergenic')
) THEN error('unexpected gene-relation class was emitted') END;
SELECT CASE WHEN (SELECT count(*) FROM gene_relation_occupancy) <> 16
  OR (SELECT count(DISTINCT gene_relation_class)
      FROM gene_relation_occupancy) <> 4
THEN error('relation-specific TP73 occupancy rows were not emitted') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM gene_relation_occupancy
  WHERE motif_id = 'M_EFFECT' AND negative_reference_threshold = -1
    AND evaluation_status = 'ok' AND adjusted_odds_ratio > 0
) THEN error('matched relation-specific TP73 occupancy was not estimated') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM score_gradient
  WHERE estimate_unit = 'one_SD_increase_in_clamped_cofactor_score'
) THEN error('continuous cofactor-score sensitivity was not emitted') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM context_config
  WHERE input_mode = 'precomputed_change'
    AND schema_version = 6
    AND annotation_schema =
      'schema9_tp73_context_anchor_collapsed_to_physical_span'
    AND gene_relation_precedence =
      'promoter_then_downstream_then_gene_body_then_intergenic'
    AND strict_intergenic_definition =
      'no_transcript_promoter_or_downstream_region_overlap'
    AND gene_relation_tp73_occupancy_output
    AND gene_relation_tp73_occupancy_h3_zero_policy =
      'retained_independent_of_h3k4me3_signal'
) THEN error('genome inference provenance is incomplete') END;
SQL

# The baseline-adjusted variant is an explicit sensitivity run. It uses the
# same outcome and anchor set, but records the extra GFP term in provenance so
# it cannot be mistaken for the prespecified primary model.
Rscript "$repository_root/scripts/analyze_h3k4me3_cofactor_change.R" \
    --change "$temporary/change-a.parquet" \
    --change "$temporary/change-b.parquet" \
    --tp73-evidence "$temporary/tp73-matched.parquet" \
    --cofactor-maxima "$temporary/maxima-effect.parquet" \
    --cofactor-maxima "$temporary/maxima-overlap.parquet" \
    --annotation "$temporary/annotation-a.parquet" \
    --annotation "$temporary/annotation-b.parquet" \
    --thresholds "$temporary/thresholds.tsv" \
    --output-prefix "$temporary/baseline-adjusted-result" \
    --series series_a --series series_b --negative-references "-1,0" \
    --block-size 500 --spline-df 1 --minimum-class-fraction 0.01 \
    --minimum-class-count 2 --minimum-interaction-cell-count 2 \
    --adjust-gfp-baseline --duckdb "$duckdb"
"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW adjusted_effect AS SELECT * FROM read_csv_auto(
  '$temporary/baseline-adjusted-result_intensity_effect.tsv', delim='\t',
  header=true, nullstr='NA');
CREATE VIEW adjusted_config AS SELECT * FROM read_csv_auto(
  '$temporary/baseline-adjusted-result_run_config.tsv', delim='\t',
  header=true, nullstr='NA');
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM adjusted_effect
  WHERE motif_id='M_EFFECT' AND evaluation_status='ok'
) THEN error('GFP-baseline-adjusted sensitivity was not fitted') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM adjusted_config
  WHERE adjustment_variant='gfp_baseline_adjusted_sensitivity'
    AND gfp_baseline_adjustment LIKE
      'spline_or_linear_term_for_normalized_mark_GFP%'
    AND primary_adjustment LIKE '%normalized_GFP_mark_term%'
) THEN error('GFP-baseline sensitivity provenance is incomplete') END;
SQL

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
    SELECT * REPLACE (0.0::DOUBLE AS source_score_floor)
    FROM read_parquet('$temporary/maxima.parquet')
    WHERE motif_id = 'M_EFFECT'
) TO '$temporary/invalid-floor-maxima.parquet' (FORMAT PARQUET);
SQL
printf 'motif_id\tpositive_threshold\tfactor_name\nM_EFFECT\t4\tSynthetic effect\n' \
    > "$temporary/invalid-floor-thresholds.tsv"
Rscript "$repository_root/scripts/analyze_h3k4me3_cofactor_change.R" \
    --signal "$temporary/signal.parquet" \
    --tp73-evidence "$temporary/tp73.parquet" \
    --cofactor-maxima "$temporary/invalid-floor-maxima.parquet" \
    --thresholds "$temporary/invalid-floor-thresholds.tsv" \
    --window central_20 --output-prefix "$temporary/invalid-floor-result" \
    --series series_a --series series_b --negative-references "-1,0" \
    --minimum-class-count 2 --minimum-interaction-cell-count 2 \
    --duckdb "$duckdb" \
    >"$temporary/invalid-floor.stdout" 2>"$temporary/invalid-floor.stderr"
"$duckdb" -batch :memory: >/dev/null <<SQL
CREATE VIEW floor_result AS SELECT * FROM read_csv_auto(
  '$temporary/invalid-floor-result_intensity_effect.tsv', delim='\t',
  header=true, nullstr='NA');
CREATE VIEW floor_gradient AS SELECT * FROM read_csv_auto(
  '$temporary/invalid-floor-result_score_gradient.tsv', delim='\t',
  header=true, nullstr='NA');
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM floor_result
  WHERE negative_reference_threshold = -1
    AND evaluation_status = 'negative_reference_below_source_floor'
    AND NOT negative_reference_observable
) THEN error('censored negative reference was not retained as an explicit status') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM floor_gradient
  WHERE score_clamp_reference = -1
    AND evaluation_status = 'negative_reference_below_source_floor'
) THEN error('censored score gradient was not retained as an explicit status') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM floor_result
  WHERE negative_reference_threshold = 0 AND negative_reference_observable
) THEN error('observable threshold was lost beside a censored threshold') END;
SQL

echo "I: H3K4me3 cofactor-change synthetic test passed."

#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-h3-analysis.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping H3K4me3 analysis-manager test." >&2
    exit 0
}
command -v Rscript >/dev/null 2>&1 || {
    echo "I: Rscript unavailable; skipping H3K4me3 analysis-manager test." >&2
    exit 0
}
Rscript -e 'library(data.table)' >/dev/null 2>&1 || {
    echo "I: data.table unavailable; skipping H3K4me3 analysis-manager test." >&2
    exit 0
}
duckdb=$(command -v "$duckdb")

sha256_file() {
    shasum -a 256 "$1" | awk '{print $1}'
}

h3="$temporary/h3"
evidence="$temporary/evidence"
annotation="$temporary/annotation"
context="$temporary/context"
runtime="$temporary/runtime"
run_root="$temporary/run"
scratch="$temporary/scratch"
mkdir -p "$h3/tables" "$evidence/tables" "$annotation" "$context" \
    "$runtime/duckdb/bin" "$scratch"
ln -s "$duckdb" "$runtime/duckdb/bin/duckdb"

printf 'dataset\tsequence_order\tchrom\tchrom_length\tanalysis_partition\tprimary_inference\trelative_path\tbytes\tsha256\n' \
    > "$h3/chromosome_file_inventory.tsv"
printf 'task_index\tchrom\tdataset\tpackage_relative_path\trun_relative_path\tabsolute_path\tbytes\tmtime_ns\tis_parquet\n' \
    > "$annotation/files.tsv"

for chrom in {1..22}; do
    h3_relative="tables/change/chrom-$chrom.parquet"
    h3_path="$h3/$h3_relative"
    mkdir -p "$(dirname "$h3_path")"
    "$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
  WITH anchors AS (
    SELECT '$chrom'::VARCHAR AS chrom,
           (100 + i * 1000000)::BIGINT AS anchor_start,
           (116 + i * 1000000)::BIGINT AS anchor_end,
           (-0.5 + (i % 7))::FLOAT AS anchor_score, i
    FROM range(12) AS r(i)
  ), series AS (
    SELECT * FROM (VALUES
      ('saos2', 'SaOS-2'), ('skmel29_2', 'SK-Mel-29')
    ) AS s(series_id, cell_line)
  ), values AS (
    SELECT a.*, s.series_id, s.cell_line,
           2.0::DOUBLE AS gfp_h3k4me3_area,
           1.0::DOUBLE AS gfp_input_area,
           CASE WHEN i % 3 = 0
                THEN 8.0 + (($chrom + i) % 5) * 0.1
                ELSE 2.0 + (($chrom + i) % 5) * 0.1 END::DOUBLE
             AS ta_h3k4me3_area,
           1.0::DOUBLE AS ta_input_area,
           CASE WHEN i % 3 = 0 THEN 1.0 ELSE 2.0 END::DOUBLE
             AS dn_h3k4me3_area,
           1.0::DOUBLE AS dn_input_area
    FROM anchors a CROSS JOIN series s
  ), normalized AS (
    SELECT *, log2((gfp_h3k4me3_area + 1) / (gfp_input_area + 1))
                 AS gfp_log2_h3k4me3_input_ratio,
              log2((ta_h3k4me3_area + 1) / (ta_input_area + 1))
                 AS ta_log2_h3k4me3_input_ratio,
              log2((dn_h3k4me3_area + 1) / (dn_input_area + 1))
                 AS dn_log2_h3k4me3_input_ratio
    FROM values
  )
  SELECT * EXCLUDE (i), 'R1'::VARCHAR AS replicate,
         'flank_150_1000'::VARCHAR AS window_name,
         150::BIGINT AS inner_bp, 1000::BIGINT AS outer_bp,
         2::INTEGER AS segment_count, 1700::BIGINT AS effective_window_bp,
         gfp_h3k4me3_area AS gfp_h3k4me3_max,
         gfp_input_area AS gfp_input_max,
         ta_h3k4me3_area AS ta_h3k4me3_max,
         ta_input_area AS ta_input_max,
         dn_h3k4me3_area AS dn_h3k4me3_max,
         dn_input_area AS dn_input_max,
         ta_log2_h3k4me3_input_ratio - gfp_log2_h3k4me3_input_ratio
           AS delta_ta_vs_gfp,
         dn_log2_h3k4me3_input_ratio - gfp_log2_h3k4me3_input_ratio
           AS delta_dn_vs_gfp,
         true AS has_any_h3k4me3_signal, true AS has_any_input_signal,
         false AS uninformative_all_h3k4me3_zero
  FROM normalized
  ORDER BY anchor_start, series_id
) TO '$h3_path' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL
    bytes=$(wc -c < "$h3_path" | tr -d ' ')
    digest=$(sha256_file "$h3_path")
    printf 'h3k4me3_anchor_change\t%s\t%s\t%s\tautosome\ttrue\t%s\t%s\t%s\n' \
        "$((chrom - 1))" "$chrom" "$((chrom * 10000000 + 20000000))" \
        "$h3_relative" "$bytes" "$digest" \
        >> "$h3/chromosome_file_inventory.tsv"

    package="$temporary/annotation-packages/chrom-$chrom"
    annotation_path="$package/tables/jaspar2026/tp73_context_anchor/genome_id=GRCh38/chrom=$chrom/data.parquet"
    mkdir -p "$(dirname "$annotation_path")"
    "$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
  SELECT (100 + i * 1000000)::BIGINT AS start,
         (116 + i * 1000000)::BIGINT AS "end",
         CASE i % 4 WHEN 0 THEN 'promoter_only'
                    WHEN 1 THEN 'downstream_only'
                    WHEN 2 THEN 'intron'
                    ELSE 'strict_intergenic' END
           ::VARCHAR AS primary_genomic_context,
         (i % 4 = 3) AS strict_intergenic,
         (i % 4 = 0) AS overlaps_any_promoter,
         (i % 4 = 1) AS overlaps_any_downstream_region,
         (i % 4 = 2) AS in_any_transcript,
         CASE i % 4 WHEN 0 THEN 'promoter'
                    WHEN 1 THEN 'downstream'
                    WHEN 2 THEN 'gene_body'
                    ELSE 'intergenic' END::VARCHAR AS gene_relation_class,
         false AS in_any_exon,
         false AS in_any_cds,
         (1000 + i * 100)::BIGINT AS nearest_tss_distance_bp,
         (1000 + i * 100)::BIGINT AS nearest_tss_genomic_distance_bp,
         false AS nearest_tss_has_mixed_strands,
         'downstream'::VARCHAR AS nearest_tss_relation,
         (2000 + i * 100)::BIGINT AS nearest_cds_distance_bp,
         (2000 + i * 100)::BIGINT AS nearest_cds_genomic_distance_bp,
         false AS nearest_cds_has_mixed_strands,
         'downstream'::VARCHAR AS nearest_cds_relation
  FROM range(12) AS r(i)
) TO '$annotation_path' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL
    bytes=$(wc -c < "$annotation_path" | tr -d ' ')
    mtime=$(stat -f %m "$annotation_path")000000000
    printf '%s\t%s\ttp73_context_anchor\ttables/jaspar2026/tp73_context_anchor/data.parquet\tpackages/chrom-%s/tables/jaspar2026/tp73_context_anchor/data.parquet\t%s\t%s\t%s\ttrue\n' \
        "$((chrom - 1))" "$chrom" "$chrom" "$annotation_path" "$bytes" \
        "$mtime" >> "$annotation/files.tsv"
done

printf '{"state":"complete","primary_inference_partition":"autosome"}\n' \
    > "$h3/manifest.json"

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
  WITH chromosomes AS (
    SELECT range::VARCHAR AS chrom FROM range(1, 23)
  ), anchors AS (
    SELECT chrom, (100 + i * 1000000)::BIGINT AS anchor_start,
           (116 + i * 1000000)::BIGINT AS anchor_end,
           (-0.5 + (i % 7))::FLOAT AS anchor_score, i
    FROM chromosomes CROSS JOIN range(12) AS r(i)
  )
  SELECT * EXCLUDE (i),
         (i % 2 = 0) AS supported_tp73_saos2_GFP,
         (i % 2 = 0)::INTEGER AS depth_tp73_saos2_GFP,
         false AS supported_negative_control_saos2_GFP,
         0::INTEGER AS depth_negative_control_saos2_GFP,
         (i % 2 = 0) AS supported_tp73_saos2_TA,
         (i % 2 = 0)::INTEGER AS depth_tp73_saos2_TA,
         false AS supported_negative_control_saos2_TA,
         0::INTEGER AS depth_negative_control_saos2_TA,
         (i % 2 = 1) AS supported_tp73_saos2_DN,
         (i % 2 = 1)::INTEGER AS depth_tp73_saos2_DN,
         false AS supported_negative_control_saos2_DN,
         0::INTEGER AS depth_negative_control_saos2_DN,
         (i % 2 = 0) AS supported_tp73_skmel29_2_GFP,
         (i % 2 = 0)::INTEGER AS depth_tp73_skmel29_2_GFP,
         false AS supported_negative_control_skmel29_2_GFP,
         0::INTEGER AS depth_negative_control_skmel29_2_GFP,
         (i % 2 = 0) AS supported_tp73_skmel29_2_TA,
         (i % 2 = 0)::INTEGER AS depth_tp73_skmel29_2_TA,
         false AS supported_negative_control_skmel29_2_TA,
         0::INTEGER AS depth_negative_control_skmel29_2_TA,
         (i % 2 = 1) AS supported_tp73_skmel29_2_DN,
         (i % 2 = 1)::INTEGER AS depth_tp73_skmel29_2_DN,
         false AS supported_negative_control_skmel29_2_DN,
         0::INTEGER AS depth_negative_control_skmel29_2_DN
  FROM anchors
) TO '$evidence/tables/tp73_anchor_evidence_autosome.parquet'
  (FORMAT PARQUET, COMPRESSION ZSTD);
SQL
printf '{"state":"complete","primary_inference_partition":"autosome"}\n' \
    > "$evidence/manifest.json"

"$duckdb" -batch "$annotation/context.duckdb" >/dev/null <<SQL
CREATE TABLE context_run AS SELECT 9::INTEGER AS context_schema_version,
  'anchor_annotation'::VARCHAR AS task_kind;
CREATE TABLE context_file_inventory AS SELECT * FROM read_csv_auto(
  '$annotation/files.tsv', delim='\t', header=true);
SQL
printf '{"state":"complete","context_schema_version":9,"task_kind":"anchor_annotation","database":"context.duckdb"}\n' \
    > "$annotation/manifest.json"

mkdir -p "$context/source-tasks"
printf 'task_index\tmotif_id\tmotif_name\tanalysis_partition\tchromosome_count\tanchor_count\tscan_minimum_score\tapplied_context_threshold\trelative_path\tabsolute_path\tbytes\tsha256\n' \
    > "$context/inventory.tsv"
for task in 0 1; do
    motif="M$((task + 1))"
    maxima="$context/source-tasks/$motif.parquet"
    "$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
  WITH chromosomes AS (
    SELECT range::VARCHAR AS chrom FROM range(1, 23)
  )
  SELECT chrom, (100 + i * 1000000)::BIGINT AS anchor_start,
         (116 + i * 1000000)::BIGINT AS anchor_end,
         '$motif'::VARCHAR AS motif_id,
         CASE WHEN ($task = 0 AND i % 3 = 0)
                    OR ($task = 1 AND i % 4 = 0)
              THEN 5.0 ELSE -2.0 END::DOUBLE
           AS context_score,
         -2.0::DOUBLE AS source_score_floor,
         150::BIGINT AS context_flank_bp,
         'signed_interval_edge_distance'::VARCHAR AS context_distance_metric
  FROM chromosomes CROSS JOIN range(12) AS r(i)
) TO '$maxima' (FORMAT PARQUET, COMPRESSION ZSTD);
SQL
    bytes=$(wc -c < "$maxima" | tr -d ' ')
    digest=$(sha256_file "$maxima")
    printf '%s\t%s\tSynthetic %s\tautosome\t22\t264\t-2\t4\tsource-tasks/%s.parquet\t%s\t%s\t%s\n' \
        "$task" "$motif" "$motif" "$motif" "$maxima" "$bytes" "$digest" \
        >> "$context/inventory.tsv"
done
"$duckdb" -batch "$context/tp73_genome_context_maxima.duckdb" >/dev/null <<SQL
CREATE TABLE context_maxima_file_inventory AS SELECT * FROM read_csv_auto(
  '$context/inventory.tsv', delim='\t', header=true);
CREATE TABLE motif_score_threshold AS
SELECT motif_id, ('synthetic_' || motif_id)::VARCHAR AS threshold_set_id,
       4.0::DOUBLE AS recommended_threshold,
       'synthetic_reachable'::VARCHAR AS calibration_status
FROM context_maxima_file_inventory;
SQL
printf '{"state":"complete","database":"tp73_genome_context_maxima.duckdb"}\n' \
    > "$context/manifest.json"

batch_count=$(
    "$repository_root/scripts/manage_h3k4me3_cofactor_analysis.py" prepare \
        --run-root "$run_root" --h3-package "$h3" \
        --evidence-package "$evidence" --context-maxima-package "$context" \
        --annotation-catalog "$annotation" --runtime-prefix "$runtime" \
        --source "$repository_root" --scratch-root "$scratch" \
        --run-id synthetic_h3_cofactor --motifs-per-batch 2 \
        --block-size 1000 --spline-df 1 --minimum-class-fraction 0.001 \
        --minimum-class-count 2 --minimum-interaction-cell-count 2 \
        --minimum-free-run-gb 0 --minimum-free-scratch-gb 0
)
[[ $batch_count -eq 1 ]]
"$repository_root/scripts/manage_h3k4me3_cofactor_analysis.py" preflight \
    --run-root "$run_root" >/dev/null
"$repository_root/scripts/manage_h3k4me3_cofactor_analysis.py" preflight \
    --run-root "$run_root" >/dev/null
"$repository_root/scripts/manage_h3k4me3_cofactor_analysis.py" run-batch \
    --run-root "$run_root" --batch-index 0 --rscript Rscript
"$repository_root/scripts/manage_h3k4me3_cofactor_analysis.py" run-batch \
    --run-root "$run_root" --batch-index 0 --rscript Rscript

status=$(
    "$repository_root/scripts/manage_h3k4me3_cofactor_analysis.py" status \
        --run-root "$run_root"
)
grep -Fq $'complete_motifs\t2' <<< "$status"
grep -Fq $'pending_motifs\t0' <<< "$status"

mkdir -p "$temporary/finalizer-tmp"
finalize=(
    "$repository_root/scripts/manage_h3k4me3_cofactor_analysis.py" finalize
    --run-root "$run_root" --duckdb "$duckdb" --threads 1
    --memory-limit 1GB --temp-directory "$temporary/finalizer-tmp"
)
"${finalize[@]}"
"${finalize[@]}"

final="$run_root/final/h3k4me3_cofactor_analysis"
for provenance_file in run_config.json fixed_inputs.tsv cofactor_tasks.tsv \
        preflight.json; do
    [[ -f "$final/provenance/$provenance_file" ]] || {
        echo "E: Missing finalized provenance file: $provenance_file" >&2
        exit 1
    }
done
(
    cd "$final"
    "$duckdb" -batch h3k4me3_cofactor_analysis.duckdb >/dev/null <<'SQL'
SELECT CASE WHEN (SELECT count(*) FROM intensity_effect) <> 16
  THEN error('combined intensity output is incomplete') END;
SELECT CASE WHEN (SELECT count(*) FROM tp73_interaction) <> 48
  THEN error('combined interaction output is incomplete') END;
SELECT CASE WHEN (SELECT count(*) FROM score_gradient) <> 16
  THEN error('combined score-gradient output is incomplete') END;
SELECT CASE WHEN (SELECT count(*)
                  FROM gene_relation_stratified_intensity_effect) <> 64
  THEN error('combined four-way gene-relation output is incomplete') END;
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM intensity_effect
  WHERE p_value IS NOT NULL AND q_value_bh_all_motifs IS NULL
) THEN error('global all-motif BH values are missing') END;
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM gene_relation_stratified_intensity_effect
  WHERE p_value IS NOT NULL AND q_value_bh_all_motifs IS NULL
) THEN error('gene-relation global BH values are missing') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM intensity_effect
  WHERE q_value_bh_task IS NOT NULL
) THEN error('task-local diagnostic q values were not retained') END;
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM run_config
  WHERE input_reference_semantics <> 'package_provenance_selector'
     OR execution_inputs_staged_to_scratch <> true
     OR change NOT LIKE 'provenance/fixed_inputs.tsv#%'
     OR cofactor_maxima NOT LIKE 'provenance/cofactor_tasks.tsv#motif_id=%'
) THEN error('portable evaluator provenance was not retained') END;
SQL
)

echo "H3K4me3 cofactor-analysis manager tests passed."

#!/usr/bin/env bash

set -euo pipefail

repository_root=$(cd "$(dirname "$0")/.." && pwd)
temporary=$(mktemp -d "${TMPDIR:-/tmp}/jaspar-h3-summary.XXXXXX")
trap 'rm -rf "$temporary"' EXIT HUP INT TERM

duckdb=${DUCKDB:-duckdb}
command -v "$duckdb" >/dev/null 2>&1 || {
    echo "I: DuckDB unavailable; skipping H3K4me3 genome-summary test." >&2
    exit 0
}

"$duckdb" -batch :memory: >/dev/null <<SQL
COPY (
  SELECT * FROM (VALUES
    ('1','saos2','SaOS-2',10,8,2,-0.2,-1.0),
    ('2','saos2','SaOS-2',20,15,5,0.1,-2.0),
    ('X','saos2','SaOS-2',100,100,0,9.0,9.0)
  ) AS t(chrom,series_id,cell_line,anchors,anchors_with_h3k4me3_signal,
         uninformative_all_h3k4me3_zero,mean_delta_ta_vs_gfp,
         mean_delta_dn_vs_gfp)
) TO '$temporary/h3_summary.parquet' (FORMAT PARQUET);

COPY (
  SELECT * FROM (VALUES
    ('M1','Factor1',1,-1,'ok',true,2.0,1.5,2.5,0.10,0.01,
     'anti_p73_enriched',100,200),
    ('M2','Factor2',2,-1,'ok',true,0.5,0.4,0.6,0.10,0.02,
     'anti_p73_depleted',80,220),
    ('M3','Factor3',3,1,'negative_reference_below_source_floor',false,
     NULL,NULL,NULL,NULL,NULL,'not_estimable',40,0)
  ) AS t(motif_id,factor_name,positive_threshold,source_score_floor,
         evaluation_status,negative_reference_observable,adjusted_odds_ratio,
         confidence_interval_95_lower,confidence_interval_95_upper,
         block_clustered_standard_error,q_value_bh_all_jaspar,
         association_direction,anchors_positive,
         anchors_negative_reference)
) TO '$temporary/enrichment.parquet' (FORMAT PARQUET);

COPY (
  WITH motifs AS (
    SELECT * FROM (VALUES ('M1','Factor1',-1),('M2','Factor2',-1),
                               ('M3','Factor3',1))
      AS m(motif_id,factor_name,source_score_floor)
  ), designs AS (
    SELECT * FROM (VALUES ('saos2','TA'),('saos2','DN'),
                               ('skmel29_2','TA'),('skmel29_2','DN'))
      AS d(series_id,isoform)
  ), refs AS (SELECT * FROM (VALUES (-1),(0)) AS r(reference))
  SELECT motif_id,factor_name,1 AS positive_threshold,source_score_floor,
         reference AS negative_reference_threshold,series_id,isoform,
         CASE WHEN motif_id='M3' AND reference=-1
              THEN 'negative_reference_below_source_floor' ELSE 'ok' END
           AS evaluation_status,
         CASE WHEN motif_id='M1' THEN 0.2 ELSE -0.2 END::DOUBLE AS estimate,
         0.05::DOUBLE AS standard_error,
         -0.3::DOUBLE AS confidence_interval_95_lower,
         0.3::DOUBLE AS confidence_interval_95_upper,
         0.01::DOUBLE AS q_value_bh_all_motifs,
         100::BIGINT AS anchors_positive,200::BIGINT AS anchors_negative
  FROM motifs CROSS JOIN designs CROSS JOIN refs
) TO '$temporary/intensity.parquet' (FORMAT PARQUET);

COPY (
  SELECT e.motif_id,e.factor_name,1 AS positive_threshold,
         e.source_score_floor,-1 AS negative_reference_threshold,
         d.series_id,d.isoform,'strict_intergenic' AS genomic_context_class,
         'ok' AS evaluation_status,0.1::DOUBLE AS estimate,
         0.02::DOUBLE AS q_value_bh_all_motifs
  FROM (VALUES ('M1','Factor1',-1),('M2','Factor2',-1))
    AS e(motif_id,factor_name,source_score_floor)
  CROSS JOIN (VALUES ('saos2','TA'),('skmel29_2','TA'))
    AS d(series_id,isoform)
) TO '$temporary/context.parquet' (FORMAT PARQUET);

COPY (
  SELECT e.motif_id,e.factor_name,1 AS positive_threshold,
         e.source_score_floor,-1 AS negative_reference_threshold,
         d.series_id,d.isoform,r.gene_relation_class,
         'ok' AS evaluation_status,0.1::DOUBLE AS estimate,
         0.02::DOUBLE AS q_value_bh_all_motifs
  FROM (VALUES ('M1','Factor1',-1),('M2','Factor2',-1))
    AS e(motif_id,factor_name,source_score_floor)
  CROSS JOIN (VALUES ('saos2','TA'),('saos2','DN'),
                     ('skmel29_2','TA'),('skmel29_2','DN'))
    AS d(series_id,isoform)
  CROSS JOIN (VALUES ('promoter'),('downstream'),('gene_body'),('intergenic'))
    AS r(gene_relation_class)
) TO '$temporary/gene_relation.parquet' (FORMAT PARQUET);

COPY (
  SELECT e.motif_id,-1 AS negative_reference_threshold,
         r.gene_relation_class,
         CASE r.gene_relation_class
           WHEN 'promoter' THEN CASE e.motif_id WHEN 'M1' THEN 4.0 ELSE 0.25 END
           WHEN 'downstream' THEN CASE e.motif_id WHEN 'M1' THEN 0.25 ELSE 4.0 END
           WHEN 'gene_body' THEN CASE e.motif_id WHEN 'M1' THEN 1.5 ELSE 0.75 END
           ELSE CASE e.motif_id WHEN 'M1' THEN 0.8 ELSE 1.2 END
         END::DOUBLE AS adjusted_odds_ratio,
         0.1::DOUBLE AS confidence_interval_95_lower,
         5.0::DOUBLE AS confidence_interval_95_upper,
         0.01::DOUBLE AS p_value,0.02::DOUBLE AS q_value_bh_all_motifs,
         CASE WHEN adjusted_odds_ratio > 1 THEN 'anti_p73_enriched'
              ELSE 'anti_p73_depleted' END AS association_direction,
         'ok'::VARCHAR AS evaluation_status,'synthetic matched fit'::VARCHAR
           AS evaluation_note,
         100::BIGINT AS anchors_total,40::BIGINT AS anchors_positive,
         60::BIGINT AS anchors_negative_reference,
         300::BIGINT AS discordant_observations,
         6::BIGINT AS matched_samples,20::BIGINT AS genomic_blocks
  FROM (VALUES ('M1'),('M2')) AS e(motif_id)
  CROSS JOIN (VALUES ('promoter'),('downstream'),('gene_body'),('intergenic'))
    AS r(gene_relation_class)
) TO '$temporary/gene_relation_occupancy.parquet' (FORMAT PARQUET);

COPY (
  SELECT e.motif_id,e.factor_name,1 AS positive_threshold,
         e.source_score_floor,-1 AS score_clamp_reference,
         d.series_id,d.isoform,'ok' AS evaluation_status,
         0.1::DOUBLE AS estimate,0.02::DOUBLE AS q_value_bh_all_motifs
  FROM (VALUES ('M1','Factor1',-1),('M2','Factor2',-1))
    AS e(motif_id,factor_name,source_score_floor)
  CROSS JOIN (VALUES ('saos2','TA'),('skmel29_2','TA'))
    AS d(series_id,isoform)
) TO '$temporary/gradient.parquet' (FORMAT PARQUET);

COPY (
  SELECT e.motif_id,e.factor_name,1 AS positive_threshold,
         e.source_score_floor,-1 AS negative_reference_threshold,
         d.series_id,d.isoform,
         'cofactor_by_tp73_confirmation_interaction' AS contrast,
         'ok' AS evaluation_status,0.1::DOUBLE AS estimate,
         0.02::DOUBLE AS q_value_bh_all_motifs
  FROM (VALUES ('M1','Factor1',-1),('M2','Factor2',-1))
    AS e(motif_id,factor_name,source_score_floor)
  CROSS JOIN (VALUES ('saos2','TA'),('skmel29_2','TA'))
    AS d(series_id,isoform)
) TO '$temporary/interaction.parquet' (FORMAT PARQUET);
SQL

summary="$temporary/summary"
"$repository_root/scripts/summarize_h3k4me3_genome_cofactors.py" \
    --enrichment "$temporary/enrichment.parquet" \
    --intensity "$temporary/intensity.parquet" \
    --context "$temporary/context.parquet" \
    --gene-relation "$temporary/gene_relation.parquet" \
    --gene-relation-occupancy "$temporary/gene_relation_occupancy.parquet" \
    --score-gradient "$temporary/gradient.parquet" \
    --interaction "$temporary/interaction.parquet" \
    --h3-summary "$temporary/h3_summary.parquet" \
    --output-dir "$summary" --duckdb "$duckdb"

"$duckdb" -batch :memory: >/dev/null <<SQL
SELECT CASE WHEN (SELECT count(*) FROM read_parquet(
  '$summary/motif_coverage.parquet')) <> 3
  THEN error('coverage registry lost a censored motif') END;
SELECT CASE WHEN (SELECT count(*) FROM read_parquet(
  '$summary/joint_primary_motif.parquet')) <> 2
  THEN error('strict joint table has the wrong estimable set') END;
SELECT CASE WHEN (SELECT count(*) FROM read_parquet(
  '$summary/correlation_robustness.parquet')) <> 8
  THEN error('conditioned/unconditioned correlation diagnostics are incomplete') END;
SELECT CASE WHEN EXISTS (SELECT 1 FROM read_parquet(
  '$summary/correlation_robustness.parquet')
  WHERE motifs <> 2 OR pearson_correlation IS NULL
    OR spearman_rank_correlation IS NULL)
  THEN error('correlation diagnostics lost their estimable motif set') END;
SELECT CASE WHEN (SELECT count(*) FROM read_parquet(
  '$summary/h3_design_diagnostics.parquet')) <> 4
  THEN error('H3 design diagnostics are incomplete') END;
SELECT CASE WHEN EXISTS (SELECT 1 FROM read_parquet(
  '$summary/joint_primary_motif.parquet')
  WHERE tp73_standard_error IS NULL OR ta_saos_se IS NULL)
  THEN error('joint output dropped per-axis uncertainty') END;
SELECT CASE WHEN (SELECT count(*) FROM read_parquet(
  '$summary/joint_reference_zero_motif.parquet')) <> 3
  THEN error('score-zero sensitivity did not recover its estimable motif') END;
SELECT CASE WHEN EXISTS (SELECT 1 FROM read_parquet(
  '$summary/joint_primary_motif.parquet') WHERE motif_id='M3')
  THEN error('censored motif entered strict inference') END;
SELECT CASE WHEN EXISTS (SELECT 1 FROM read_parquet(
  '$summary/h3_autosome_baseline.parquet')
  WHERE anchors<>30 OR abs(mean_delta_ta_vs_gfp) > 1e-12
    OR abs(mean_delta_dn_vs_gfp + 5.0/3.0) > 1e-12)
  THEN error('autosome baseline is wrong or includes chromosome X') END;
SELECT CASE WHEN (SELECT count(*) FROM read_parquet(
  '$summary/gene_relation_primary_effect.parquet')) <> 32
  THEN error('four-way gene-relation summary is incomplete') END;
SELECT CASE WHEN NOT EXISTS (
  SELECT 1 FROM read_parquet('$summary/gene_relation_primary_effect.parquet')
  WHERE motif_id='M1' AND gene_relation_class='promoter'
    AND tp73_adjusted_odds_ratio=4.0
    AND global_tp73_adjusted_odds_ratio=2.0
) THEN error('gene-relation summary reused the global TP73 odds ratio') END;
SQL

[[ -s $summary/manifest.json && -s $summary/summary_metrics.tsv ]]
grep -Fq '"schema_version": 4' "$summary/manifest.json"
grep -Fq '"ranking_status": "provisional_pending_baseline_covariates_and_empirical_null"' \
    "$summary/manifest.json"
if command -v Rscript >/dev/null 2>&1 && Rscript -e \
    'quit(status=ifelse(requireNamespace("data.table",quietly=TRUE) && requireNamespace("ggplot2",quietly=TRUE),0,1))'; then
    "$repository_root/scripts/plot_h3k4me3_genome_cofactor_summary.R" \
        --joint "$summary/joint_primary_motif.tsv" \
        --context "$summary/context_primary_effect.tsv" \
        --gene-relation "$summary/gene_relation_primary_effect.tsv" \
        --output-effect "$temporary/effect.png" \
        --output-context "$temporary/context.png" \
        --context-table "$temporary/context.tsv" \
        --output-gene-relation-prefix "$temporary/gene-relation" \
        --gene-relation-table "$temporary/gene-relation.tsv" \
        --label-motifs M1
    for relation in promoter downstream gene_body intergenic; do
        [[ -s $temporary/gene-relation-$relation.png ]]
    done
    [[ -s $temporary/effect.png && -s $temporary/context.png \
       && -s $temporary/gene-relation.tsv ]]
fi
if "$repository_root/scripts/summarize_h3k4me3_genome_cofactors.py" \
    --enrichment "$temporary/enrichment.parquet" \
    --intensity "$temporary/intensity.parquet" \
    --context "$temporary/context.parquet" \
    --gene-relation "$temporary/gene_relation.parquet" \
    --gene-relation-occupancy "$temporary/gene_relation_occupancy.parquet" \
    --score-gradient "$temporary/gradient.parquet" \
    --interaction "$temporary/interaction.parquet" \
    --h3-summary "$temporary/h3_summary.parquet" \
    --output-dir "$summary" --duckdb "$duckdb" >/dev/null 2>&1; then
    echo "E: Extractor overwrote an existing immutable summary." >&2
    exit 1
fi

echo "H3K4me3 genome cofactor-summary tests passed."

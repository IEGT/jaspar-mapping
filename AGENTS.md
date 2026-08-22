# AGENTS.md

## Promoter Architecture And ML Query Surface

For work on promoter architecture, expression-linked ML features, OpenClaw
integration, or agent-facing query surfaces, read
[`docs/promoter_architecture_ml_schema.md`](docs/promoter_architecture_ml_schema.md)
before designing or changing code.

For motif-score calibration against raw CUT&RUN fragments, also read
[`docs/cutandrun_motif_score_calibration.md`](docs/cutandrun_motif_score_calibration.md).
For joint TP73-score/cofactor-score CUT&RUN response surfaces, fixed negative
references, or empirical cofactor score bands, also read
[`docs/tp73_cofactor_score_response.md`](docs/tp73_cofactor_score_response.md).
For all-JASPAR enrichment/depletion around TP73 anchors, depth-stratified
CUT&RUN summaries, or their Slurm production run, also read
[`docs/jaspar2026_chr1_all_cofactor_enrichment_run.md`](docs/jaspar2026_chr1_all_cofactor_enrichment_run.md).
For species-aware vertebrate filtering and the isoform-by-exclusive-distance
TP73 cofactor ranking, also read
[`docs/tp73_distance_species_aware_cofactor_enrichment.md`](docs/tp73_distance_species_aware_cofactor_enrichment.md).
Join all-motif rankings to the normalized JASPAR 2026 catalog described in
[`docs/jaspar_2026_provenance.md`](docs/jaspar_2026_provenance.md). Preserve the
all-taxon sequence-signature result, but use `tax_group = 'vertebrates'` for the
primary human cofactor interpretation and report `includes_homo_sapiens` as a
separate source-provenance sensitivity. A matrix source species is not a claim
about the only species in which its sequence pattern can function.
For the whole-genome extension, autosome/sex/mitochondrial partitioning,
zero-complete context counts, or scan-floor censoring of negative references,
also read
[`docs/jaspar2026_tp73_whole_genome_followup.md`](docs/jaspar2026_tp73_whole_genome_followup.md).
For fixed-threshold genomic transfer of the prespecified cofactor panel and
H3K4me3 change effects, also read
[`docs/jaspar2026_chr2_heldout_validation_20260810.md`](docs/jaspar2026_chr2_heldout_validation_20260810.md).
The completed schema-7 local-peak result is reported in
[`docs/jaspar2026_chr1_localpeak_enrichment_20260810.md`](docs/jaspar2026_chr1_localpeak_enrichment_20260810.md).

For GFP-referenced H3K4me3 change at TP73 anchors, read
[`docs/h3k4me3_cofactor_change.md`](docs/h3k4me3_cofactor_change.md). The
completed chromosome-1 result is reported in
[`docs/h3k4me3_chr1_production_20260810.md`](docs/h3k4me3_chr1_production_20260810.md).
For the schema-9 annotation rebuild and restart-safe whole-genome H3K4me3
signal and all-JASPAR cofactor analysis, also read
[`docs/h3k4me3_whole_genome_production.md`](docs/h3k4me3_whole_genome_production.md).
Gene-relation figures must join `gene_relation_stratified_intensity_effect`
to `gene_relation_stratified_tp73_occupancy` by motif, negative reference, and
relation class; never substitute the global TP73 occupancy estimate on their
x-axis.
The completed all-autosome joint TP73 enrichment/H3K4me3 interpretation is
reported in
[`docs/h3k4me3_tp73_cofactor_results_autosomes_20260812.md`](docs/h3k4me3_tp73_cofactor_results_autosomes_20260812.md).

For TP73-local motif spacing, tandem TP73 sites, or transcript/intron context,
also read [`docs/tp73_motif_context.md`](docs/tp73_motif_context.md).
For motif-specific convenient thresholds, score-floor calibration, or applying
thresholds to context queries, also read
[`docs/motif_score_thresholds.md`](docs/motif_score_thresholds.md).

For whole-genome JASPAR 2026 scans, Slurm batching/requeue behavior, production
inventories, or cross-species identity, also read
[`docs/jaspar2026_genome_scan_plan.md`](docs/jaspar2026_genome_scan_plan.md).
For the informative-or-`-1` per-motif policy, chromosome-1 physical-locus
density calibration, or the 200 bp production ceiling, also read
[`docs/jaspar2026_informative_density200_scan.md`](docs/jaspar2026_informative_density200_scan.md).
For scanner execution, scratch staging, finalization, checksum verification,
build provenance, or exact-file DuckDB queries, also read
[`docs/scanner_maintenance.md`](docs/scanner_maintenance.md).

For TP73 context analyses, never use a density-capped threshold as the source
retention floor. Select maxima from the low-floor scan, and keep source-floor,
score-zero, and empirical operating-threshold counts as separate fields.

Treat [`sql/schema.sql`](sql/schema.sql) and [`sql/queries.sql`](sql/queries.sql)
as the draft query contract:

- Use BED 0-based half-open coordinates for exported motif and promoter tables.
- Treat `transcription_start_site` and `promoter` as physical interval
  dimensions. Recover transcript/gene ownership through `transcript_tss` and
  `promoter_gene`; do not duplicate or collapse shared TSSs into a gene-grain
  table.
- Use Q20 with `tp73_anchor_nearest_tss` when all tied nearest TSSs and their
  associated genes matter. Use Q21 with `tp73_anchor_promoter` for canonical
  many-to-many promoter membership, and always bind `promoter_definition_id`.
  Membership means positive half-open interval overlap; mere abutment is not
  membership.
- Use Q22 with `tp73_anchor_nearest_cds` when all tied nearest physical CDS
  segments and transcript/gene owners matter. Use Q23 only for the compact,
  deterministic anchor annotation covariates used by statistical models.
- Use Q24 with `tp73_anchor_downstream_region` for versioned many-to-many
  transcript-end downstream membership. The compact `gene_relation_class`
  applies promoter, downstream, gene-body, then intergenic precedence without
  replacing the independent relationship tables.
- Treat TSSs, promoters, and `coding_sequence_segment` rows as physical
  dimensions. Recover their transcript/gene ownership and CDS phase through
  the corresponding bridge tables; exon and intron rows remain
  transcript-specific because isoforms can assign different boundaries.
- Keep large generated data out of Git; package it as versioned Parquet plus a
  rebuildable DuckDB query index.
- Prefer the materialized `promoter_card` for low-latency agent lookups, bind
  its `score_mode` and `pseudocount` explicitly, and use named parameterized
  queries instead of ad hoc joins from agents.
- Keep ML features in long/tidy tables; pivot to wide matrices only at the ML
  export boundary.
- Keep sequence-derived promoter architecture separate from CUT&RUN-derived
  evidence unless an analysis explicitly declares an "architecture + binding"
  feature block.

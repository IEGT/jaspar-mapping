-- Optional comparative-genomics layer. Load after genome_scan_schema.sql in a
-- package that contains the corresponding Parquet files.

CREATE OR REPLACE VIEW syntenic_interval AS
SELECT
    compara_release,
    compara_division,
    method_link_type,
    method_link_species_set_id,
    ensembl_synteny_region_id,
    source_genome_id,
    source_chrom,
    source_start,
    source_end,
    target_genome_id,
    target_chrom,
    target_start,
    target_end,
    target_strand,
    mapping_rank,
    mapping_quality,
    aligned_base_fraction,
    reciprocal,
    synteny_source_id,
    synteny_source_sha256
FROM read_parquet(
    'tables/jaspar2026/syntenic_interval/**/*.parquet',
    hive_partitioning = true,
    union_by_name = true
);

-- Ensembl synteny regions are intentionally too coarse for projecting a short
-- motif. Pairwise LASTZ_NET blocks provide the exact mapping substrate. The
-- importer assigns alignment_block_id when a bulk MAF record does not expose
-- an Ensembl genomic_align_block_id.
CREATE OR REPLACE VIEW comparative_alignment_block AS
SELECT
    compara_release,
    compara_division,
    method_link_type,
    method_link_species_set_id,
    alignment_block_id,
    ensembl_genomic_align_block_id,
    source_genome_id,
    source_chrom,
    source_start,
    source_end,
    source_strand,
    source_cigar,
    target_genome_id,
    target_chrom,
    target_start,
    target_end,
    target_strand,
    target_cigar,
    alignment_score,
    alignment_source_id,
    alignment_source_sha256
FROM read_parquet(
    'tables/jaspar2026/comparative_alignment_block/**/*.parquet',
    hive_partitioning = true,
    union_by_name = true
);

CREATE OR REPLACE VIEW orthologous_motif_hit AS
SELECT
    orthology_group_id,
    compara_release,
    method_link_species_set_id,
    ensembl_synteny_region_id,
    alignment_block_id,
    projection_status,
    source_genome_id,
    source_chrom,
    source_start,
    source_end,
    source_motif_set_id,
    source_motif_id,
    source_strand,
    source_score,
    target_genome_id,
    target_chrom,
    target_start,
    target_end,
    target_motif_set_id,
    target_motif_id,
    target_strand,
    target_score,
    target_hit_present,
    motif_sequence_conserved,
    score_delta,
    mapping_rank,
    reciprocal,
    aligned_base_fraction,
    derivation_run_id
FROM read_parquet(
    'tables/jaspar2026/orthologous_motif_hit/**/*.parquet',
    hive_partitioning = true,
    union_by_name = true
);

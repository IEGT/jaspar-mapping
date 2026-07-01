# TF Mapping Schema Notes

## Project context

These notes extend the source tree for Möller et al. 2025, which introduced
the current JASPAR-mapping and CUT&RUN/expression analysis workflow for TAp73
and DNp73 cancer programs (PMID: 41594604; PMCID: PMC12839168; DOI:
10.3390/biom16010063).

The updated data model should remain layered:

1. Raw motif mapping: genomewide JASPAR hits with coordinates, score, strand,
   motif identity, and matched sequence.
2. Motif architecture: assay-independent annotations derived from the motif
   hit and motif model, including half-site/spacer structure and
   oligomerisation compatibility. This layer should be split into a compact
   motif-level sidecar and separate per-hit derived architecture files.
3. Experimental evidence: CUT&RUN support, coverage, peak overlap, condition,
   antibody, cell line, and isoform context mapped onto motif or architecture
   hits. This layer should move toward long-form artifacts instead of
   preserving the current wide matrix as the primary exchange format.
4. Regulatory and gene context: promoters, UTRs, gene bodies, nearby genes,
   gene sets, and expression-data associations.

## Filtering semantics

There are two distinct filtering concepts that should stay separate in the
schema and in file names:

1. Sequence-score filtering: the default scanner mode normalizes JASPAR counts to a
   summed log2 relative-risk score, `log2(frequency / background)`, using a
   uniform A/C/G/T background. This is the active published/default scoring
   path. The scanner also supports an explicit `--score-mode log_odds` option.
2. Evidence or contrast filtering: downstream CUT&RUN, expression, promoter,
   or enrichment analyses can filter on log ratios, log odds, peak support,
   or other contrast-specific metrics. These should not be conflated with the
   motif sequence score.

The raw mapping layer should therefore carry explicit score provenance:

- `sequence_score`
- `sequence_score_mode`, initially `log2_relative_risk`
- `sequence_score_threshold`
- `background_model`, initially `uniform_acgt_0.25`
- `jaspar_version`

For score-distribution sidecars, use a separate histogram table rather than
emitting every subthreshold genomic window. The first intended sampling run
uses chromosome 1 on the forward strand only, with fields such as `motif_id`,
`motif_name`, `chromosome`, `strand`, `score_mode`, `bin_scheme`, `bin_width`,
`score_bin_start`, `score_bin_end`, `bin_count`, `valid_windows`, and
`skipped_windows`. The default adaptive log-ratio binning uses widths of 0.2
for scores >= -10, 1 for scores >= -50, 5 for scores >= -250, 10 for scores
>= -1000, 100 for scores >= -10000, and 500 below that. In
`log2_relative_risk` mode, `skipped_windows` includes windows that hit
zero-probability motif entries and therefore receive the scanner's sentinel
score rather than a finite negative score.

The experimental and regulatory layers should carry separate contrast
provenance:

- `effect_metric`, for example `log2_ratio`, `log_odds`, `peak_overlap`, or
  `coverage_delta`
- `contrast`, for example `TAp73_vs_DNp73`
- `effect_threshold`
- `direction`, for example `TAp73_enriched`, `DNp73_enriched`, or
  `bidirectional`
- `support_metric`, such as peak count, read coverage, or expression change

Log-odds motif scoring is available only as an explicit scanner option; changing
the meaning of existing score columns silently would make old and new runs hard
to compare.

## Oligomerisation metadata proposal

The current analysis mostly treats motif identity, score, count, and strand as
the data carried forward into the table layer. For display, accessibility, and
future OpenClaw use, it would be helpful to add a small amount of redundant
motif architecture metadata so a table row can be understood without repeatedly
dereferencing JASPAR or external annotations.

The clearest example is the TP53 family. Wilson et al. 2025 describe TP53,
TP63, and TP73 binding sites as two RRRCWWGYYY half-sites, separated by a
0-13 bp spacer, and bound by the TF as a tetramer, specifically a dimer of
dimers. The same paper emphasizes that nucleosome binding depends not only on
motif composition but also on accessibility and helical orientation.

Suggested motif-level fields:

- `motif_id`: JASPAR accession, for example `MA0861.2`
- `motif_name`: JASPAR motif name, for example `TP73`
- `tf_family`: broad family label, for example `TP53_family`
- `binding_unit_model`: controlled value such as `monomer`, `homodimer`,
  `heterodimer`, `dimer_of_dimers`, `tetramer`, or `unknown`
- `oligomeric_state`: redundant integer or short label, for example `4` or
  `tetramer`
- `dna_binding_unit_count`: number of DNA-contacting motif units
- `half_site_count`: number of half-sites represented by the motif model
- `half_site_pattern`: optional consensus, for example `RRRCWWGYYY`
- `spacer_min_bp` and `spacer_max_bp`: allowed half-site spacer range
- `architecture_source`: citation or PMCID supporting the architecture
- `architecture_confidence`: controlled value such as `curated`,
  `database_inferred`, or `unknown`

Suggested per-hit redundant fields:

- `motif_id`
- `motif_name`
- `binding_unit_model`
- `oligomeric_state`
- `architecture_class`, for example `full_site`, `half_site`, or
  `degenerate_full_site`
- `spacer_bp`, when a motif can be decomposed into half-sites
- `half_site_1_score` and `half_site_2_score`, where applicable
- `oligomer_compatible`, a boolean indicating whether the matched sequence has
  the expected architecture for the motif family
- `nucleosome_context_model`, initially `not_assessed`, reserving space for
  later values such as `edge_accessible`, `dyad_centered`, or
  `orientation_sensitive`

For the current wide tables, these fields can be introduced as companion
metadata rather than as many repeated columns. A compact sidecar table keyed by
`motif_id` would avoid unnecessary bloat, while generated long-format exports
for display can repeat the fields per hit for accessibility.

References:

- Möller et al. 2025, Integrative Multimodal Profiling of TAp73 and DNp73
  Reveals Isoform-Specific Transcriptomic Coregulator Landscapes in Cancer
  Programs. https://pubmed.ncbi.nlm.nih.gov/41594604/
- Seelan et al. 2002, The human p73 promoter: characterization and
  identification of functional E2F binding sites.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC1531693/
- Wilson et al. 2025, Nucleosome binding by TP53, TP63, and TP73 is
  determined by the composition, accessibility, and helical orientation of
  their binding sites. https://pmc.ncbi.nlm.nih.gov/articles/PMC11960462/

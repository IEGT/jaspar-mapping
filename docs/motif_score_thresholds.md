# Motif score-threshold registry

## Purpose

`motif_score_threshold` stores versioned, evidence-derived proposals for turning
a continuous motif score into an inclusive binary feature. It does **not**
replace either the raw score or `scan_threshold_policy`:

- `scan_threshold_policy.minimum_score` is the storage floor used while
  scanning. It determines which alignments remain recoverable without a rescan.
- `motif_score_threshold.recommended_threshold` is a downstream interpretation
  chosen for one target, assay, geometry, and validation design.
- ML analyses should ordinarily retain the continuous score and use the
  convenient threshold as a comparison feature, visualization aid, or simple
  query default.

A score is not an assay-independent binding constant. Every consumer must bind
`threshold_set_id`, target motif, genome/motif set, and full scoring
configuration. A threshold calibrated for TP73 CUT&RUN context must not silently
be reused for another transcription factor, promoter outcome, species, score
mode, or pseudocount.

## Grain and key

One row describes one motif in one named threshold set. Its logical key is:

```text
threshold_set_id, genome_id, motif_set_id, motif_id, threshold_role,
target_motif_id, score_mode, pseudocount, background_model_id,
pseudocount_scheme, calibration_stratum_id
```

The first fill uses:

```text
threshold_role = tp73_context_binary_feature
target_motif_id = MA0861.2
score_mode = log2_relative_risk
pseudocount = 1
context_distance_metric = signed_interval_edge_distance
context_flank_bp = 150
context_min_interval_distance_bp = NULL
context_max_interval_distance_bp = 150
context_relation_filter = any
calibration_stratum_id = all_tp73_anchors
source_minimum_score = -1
```

TP73's own direct CUT&RUN threshold is a different role (for example,
`anchor_cutrun_support`) and must be imported from its direct calibration rather
than inferred by treating TP73 as its own neighboring cofactor.

## Columns

Identity and scoring columns pin the threshold set, genome, motif set, motif,
target, score mode, pseudocount/background semantics, and the lower storage
floor from which candidates could be recovered.

Evidence columns pin the calibration scope, CUT&RUN dataset, outcome definition,
context geometry, fold definition, anchor/observation counts, and source metric
URI/checksum. `schema_version` versions the row contract; `source_dirty` makes
exploratory local runs honest when they were produced before their
implementation commit existed.
The stratum ID plus JSON definition permits later singleton/tandem, genomic
region, or other biologically declared strata without turning them into hidden
filename conventions. Signed interval bounds and the relation filter support
overlap, adjacency, and wider distance-band calibrations.

Selection columns record the candidate grid, tested range/count, optimized
metric, tie rule, and the fraction used to define a near-optimal range. The
initial rule maximizes held-out macro ROC AUC gain over the TP73-score baseline;
ties choose the lower integer threshold. `useful_threshold_min/max` includes
tested thresholds reaching at least 90% of the selected positive gain.
By default, a candidate is eligible only when both the retained and absent
classes contain at least 1% of anchors; the bound and eligible-candidate count
are stored rather than hidden in the selector.

Result columns contain the recommended inclusive threshold, useful range,
retained-anchor fraction, selected metric/gain, adjusted odds ratio, effect
direction, and descriptive sample/fold consistency counts. The sample-fold
cells are correlated summaries, not independent evidence.

## Threshold-qualified occurrence counts

Once a threshold set has been selected, the sparse context builder can
materialize `tp73_motif_threshold_count`. Its grain is one physical TP73
alignment locus, neighboring motif accession, and fully identified threshold
set. The initial feature schema is version 1.
`n_neighbor_loci_above_threshold` is the number of distinct neighboring
alignment spans whose best orientation-specific score meets the motif's
inclusive threshold and recorded interval-distance/relation bounds.

The count deliberately collapses `+` and `-` records at an identical genomic
span. Those records are alternative motif orientations, not two physical
places. Every requested anchor/motif combination is emitted, including an
explicit zero when no locus qualifies. This makes the field suitable for a
long-form ML feature or a later wide pivot without conflating missing rows with
zero occurrences.

The term "locus" is intentional. Passing a sequence-score threshold does not
show that the corresponding protein bound the DNA. CUT&RUN remains a separate
experimental label or feature block.

For a populated registry, rerun `scripts/build_sparse_context_maxima.py` with
`--threshold-parquet` and `--threshold-set-id`. The tool validates that every
requested motif has exactly one populated threshold row, that score and source
floors are compatible, and that the registered interval geometry fits inside
the supplied capture flank. The output retains `context_score` as the
unthresholded maximum alongside the threshold-qualified count.

## Null and status semantics

- `exploratory_positive_gain`: a positive held-out gain was observed and the
  recommended threshold is populated.
- `no_positive_gain`: candidates were evaluated but none improved the selected
  metric; the recommended threshold remains `NULL`.
- `insufficient_class_support`: candidates were evaluated but none met the
  recorded retained/absent support bound; the recommendation remains `NULL`.
- `no_finite_metric`: supported candidates exist but the selected metric is not
  finite; the recommendation remains `NULL`.
- `pending`: the motif was expected but no metric rows were supplied; the
  recommended threshold remains `NULL`.

Consumers must not convert either null state to zero. A separate, explicitly
named fallback policy may be used for display, but it is not calibration.

## Initial chromosome-1 fill

The first local set evaluates the strongest interval-defined score within 150
bp of each TP73 anchor for:

| Motif | JASPAR accession |
| --- | --- |
| E2F1 | `MA0024.3` |
| SP1 | `MA0079.5` |
| REST | `MA0138.3` |
| POU2F2 | `MA0507.3` |
| KLF14 | `MA0740.2` |
| TCF7 | `MA0769.3` |
| POU4F1 | `MA0790.2` |
| TFAP2C | `MA0814.3` |
| PATZ1 | `MA1961.2` |

All source scans retain scores at least -1. The KLF14 chromosome-1 scan is
generated with the same JASPAR 2026, GRCh38, score, pseudocount, coordinate,
strand, and sparse-Parquet configuration as the existing eight-motif panel.
`scripts/build_sparse_context_maxima.py` creates the rectangular long-form
maximum table; `scripts/evaluate_tp73_cofactor_thresholds.R --thresholds auto`
tests every observed non-negative integer floor independently for each motif;
and `scripts/build_motif_score_thresholds.py` writes the registry directly to
Parquet.

The calibration is exploratory because threshold selection and reporting use
the same five contiguous chromosome-1 folds. Before claiming a biological
cutoff, add anchor-clustered uncertainty and validate the selected
transformations on held-out chromosomes or independent CUT&RUN data.

The populated values, retained fractions, and interpretation limits are
reported in
[`tp73_context_convenient_thresholds_chr1_20260721.md`](tp73_context_convenient_thresholds_chr1_20260721.md).

## Query contract

`sql/schema.sql` exposes the complete `motif_score_threshold` registry and the
populated convenience subset `motif_convenient_threshold`. Q16 in
`sql/queries.sql` inspects one explicitly named set, including pending/no-gain
rows. Q17 applies only populated thresholds to `motif_context_pair` while
preserving the underlying continuous records and enforcing the calibrated
interval-distance and relation bounds.
Q18 reads the materialized, zero-complete `tp73_motif_threshold_count` surface
for one cofactor and returns one row per physical TP73 anchor.

For a genome-wide all-JASPAR fill, process motif partitions independently and
append rows to a new immutable threshold-set version. Never mutate the first
chromosome-1 set when the candidate grid, context geometry, evidence labels, or
validation design changes.

The production chromosome-1 all-motif execution contract, reproducible
CUT&RUN anchor reconstruction, Slurm layout, status reporting, and finalization
gate are specified in
[`jaspar2026_chr1_all_motif_threshold_run.md`](jaspar2026_chr1_all_motif_threshold_run.md).

On Haumea, the scalable unit should be one neighboring motif against the fixed
TP73 anchor/label table. Resolve its exact plus/minus chromosome files through
`scan_file_inventory`, derive maxima and metric rows in a requeueable task, and
reuse one pinned fold manifest and TP73-only baseline across tasks. Do not make
each task discover files with a package-wide glob. Finalization should verify
one metric/registry row per expected motif and preserve explicit pending rows
for failed or not-yet-run tasks.

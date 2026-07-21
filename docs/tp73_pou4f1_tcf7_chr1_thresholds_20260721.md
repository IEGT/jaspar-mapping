# POU4F1 and TCF7 score thresholds around TP73 sites on chromosome 1

## Question and conclusion

This analysis asks whether the very frequent POU4F1 (`MA0790.2`) and TCF7
(`MA0769.3`) motif matches are useful predictors of the strict TP73 CUT&RUN
immersion label, and whether score floors above zero improve their usefulness.
POU2F2 (`MA0507.3`) is included as the established inhibitory comparison.

The low-score matches are indeed too frequent to be useful as binary features:
at score 0, POU4F1 occurs within 150 bp of 95.41% of TP73 anchors and TCF7
within 99.58%. Raising the floor restores discrimination. The exploratory
operating points are:

| Motif | Suggested binary floor | Anchors retained | Macro ROC AUC gain | Median adjusted odds ratio |
| --- | ---: | ---: | ---: | ---: |
| POU2F2 `MA0507.3` | 0 | 40.03% | +0.03078 | 0.736 |
| POU4F1 `MA0790.2` | 5 | 70.71% | +0.04278 | 0.653 |
| TCF7 `MA0769.3` | 6 | 69.51% | +0.02285 | 0.765 |

POU4F1 therefore does reproduce the earlier inhibitory association and, by
this ranking metric, adds more discrimination than the POU2F2 indicator.
TCF7 is also informative after thresholding, but its effect is smaller.
These are associations with the technical CUT&RUN support label, not proof of
causal inhibition or direct protein binding.

## Inputs and geometry

The analysis uses 310,782 chromosome-1 TP73 anchors and six matched anti-p73
versus control CUT&RUN comparisons. Only discordant anti/control observations
enter the binary outcome, yielding 532,116 `(anchor, sample)` observations.

Cofactor maxima were recomputed from the chromosome-1 sparse scans retained at
score -1. They use schema-v4 geometry: signed interval-edge distance no greater
than 150 bp. Thus abutting BED intervals have distance 0 and overlaps have a
negative distance. The derived lossless center prefilter is 179 bp for these
motifs. `scripts/build_sparse_context_maxima.py` writes the long Parquet table
and a companion `.run_config.tsv` containing every exact source path.

The full chromosome scan contains the following orientation-specific records:

| Motif | Width | Score >= -1 | Score >= 0 | Records at suggested floor |
| --- | ---: | ---: | ---: | ---: |
| POU2F2 `MA0507.3` | 13 bp | 661,894 | 528,420 | 528,420 at 0 |
| POU4F1 `MA0790.2` | 12 bp | 15,326,568 | 11,585,518 | 2,210,504 at 5 |
| TCF7 `MA0769.3` | 7 bp | 12,925,512 | 9,930,039 | 1,244,166 at 6 |

These record counts preserve orientation and are not counts of distinct
orientation-collapsed loci.

## Validation design

The baseline logistic model contains a sample intercept and a four-degree
natural spline of the TP73 score. Each augmented model adds one binary
indicator saying that the best interval-defined cofactor score reaches the
tested floor. Five equal-width contiguous chromosome regions are held out in
turn; all sample copies of a TP73 anchor remain in the same fold.

The TP73-only macro ROC AUC is 0.52092. The complete integer threshold sweep
from 0 through 12 gives these main observations:

- POU2F2 is most useful at score 0 and loses information monotonically as the
  floor rises.
- POU4F1 has its highest ROC AUC gain at 5. A floor of 3 gives the best average
  precision and log-loss gains, so the broad useful range is approximately
  3-6 rather than one biologically exact cutoff.
- TCF7 has its highest ROC AUC and average-precision gains at 6, with 5-6 the
  useful range. Its seven-base matrix produces discrete score classes; floors
  8 and 9 consequently select exactly the same anchors.
- The raw odds ratio is below one in all six samples at every tested floor for
  all three motifs. At the suggested floors, the ROC AUC improves in 30/30
  sample-fold cells for POU2F2, 29/30 for POU4F1, and 28/30 for TCF7. These
  cells are correlated descriptive summaries, not independent replicates.
- At their suggested floors, the three binary features are only weakly
  correlated across anchors (phi 0.032-0.094). A joint conditional model is
  still required before calling any effect independent of the other motifs.

The thresholds were selected and displayed using the same cross-validation
folds. There is no outer selection split, anchor-clustered confidence interval,
or held-out chromosome in this exploratory run.

## Relation to the schema-v4 review

Claude's review of the interval geometry, derived prefilter, canonical
same-motif pairing, and context-scoped pair counting is consistent with the
implementation and tests. Two qualifications matter downstream:

1. The source TP73 scan retains scores down to -5, but the current context
   anchor selection uses local peaks at score at least -1. Those are distinct
   provenance decisions; the context package does not currently expose every
   retained -5 TP73 scan record as an anchor.
2. `cutandrun_score_calibration --feature-parquet` still defines a cofactor
   neighborhood by center distance. On these anchors, its POU2F2 score-at-least
   -1 indicator selected 134,025 anchors, whereas schema-v4 interval geometry
   selected 143,479. The interval calculation adds 9,454 anchors and changes
   the selected maximum for 4,761 anchors shared by both sets. The present
   analysis therefore reuses only the TP73/CUT&RUN labels from that feature
   file and recomputes every cofactor score with schema-v4 geometry.

## Recommendation

Keep the general sparse scan floor at -1. It is the storage floor, not the
eventual biological decision threshold, and it permits this sensitivity
analysis without rescanning the genome. For simple binary summaries use score
5 for POU4F1 and score 6 for TCF7, while retaining score 0 for POU2F2. For the
planned multivariable model, prefer the actual maximum score (with an explicit
missing/below--1 state) over prematurely reducing it to one binary feature.

Before a biological claim, repeat the model jointly with the other cofactors,
add anchor-clustered uncertainty, and validate the selected transformations on
held-out chromosomes or independent CUT&RUN data.

## Reproducible outputs

The local run is under:

`dry_runs/tp73_tcf7_pou4f1_threshold_sweep_20260721/`

The coherent accession-keyed output prefix is
`tp73_cofactor_thresholds_accessions`. Its `_run_config.tsv`,
`_fold_manifest.tsv`, `_threshold_metrics.tsv`, `_sample_fold_metrics.tsv`,
`_fold_coefficients.tsv`, and `_raw_association.tsv` tables are generated by
`scripts/evaluate_tp73_cofactor_thresholds.R`. The figure is
`_threshold_sweep.png`. Both new scripts have self-documenting `--help` output,
and their synthetic tests require no genome or JASPAR download.

# TP73 Pair Architecture And Cofactors: Chromosome-1 CUT&RUN Model, 2026-07-18

## Summary

Distinct nearby TP73 motif-model alignments carry reproducible information
about TP73-specific CUT&RUN immersion, but the effect is modest compared with
the local cofactor scores. In the primary interleaved chromosome-block split:

- TP73 score alone reached mean within-sample ROC AUC `0.5174`;
- TP73 score plus TP73 pair architecture reached `0.5272`;
- TP73 plus PATZ1, TFAP2C, POU2F2, and E2F1 reached `0.5951`; and
- the additive TP73, pair-architecture, and four-cofactor model reached
  `0.5997`.

Adding pair architecture to TP73 score improved AUC in 27 of 30 sample/fold
cells. Adding it to TP73 plus cofactors improved AUC in 23 of 30 cells. An
explicit pair-class-by-cofactor interaction model did not improve held-out
performance: its AUC was `0.5992`, a mean change of `-0.00047`, and it won in
only 10 of 30 cells. The present evidence therefore supports retaining pair
architecture as an additive feature, not yet fitting separate cofactor rules
for every pair class.

In a stricter sensitivity analysis using five contiguous chromosome-1 spans,
the additive full-model AUC fell to `0.5867`, but pair architecture still added
`0.0084` AUC over TP73 alone and `0.0041` over TP73 plus cofactors. The absolute
primary AUC should therefore be treated as optimistic; the smaller incremental
architecture result is more stable than the absolute performance estimate.

These classes describe sequence architecture. They do not directly observe a
TP73 oligomer or quaternary structure.

## Inputs And Target

The analysis used chromosome 1, JASPAR 2026 `MA0861.2`,
`log2_relative_risk`, pseudocount 1, and TP73 scores at least `-1`. The
calibrator selects one orientation per genomic start by maximum TP73 score.
This leaves 417,723 candidate starts. Of these, 337,382 had at least one
discordant anti-p73/control observation and contributed 715,678 sample-anchor
observations to the model.

The six biological comparisons are TA and DN in Saos-2 and two SK-Mel-29
replicates. For each candidate and sample:

- outcome 1 means the full 16-bp TP73 alignment was strictly immersed in an
  anti-p73 coverage component but not its matched control component;
- outcome 0 means strict immersion in control but not anti-p73; and
- candidates supported by both or neither were excluded from model fitting.

Strict immersion requires coverage to extend beyond both sides of the scored
alignment. It is a technical occupancy label, not proof of direct binding.
Sample identity is included as a nuisance intercept. Feature slopes are shared
across all six samples in this first analysis; the model does not yet test
sample- or isoform-specific slopes.

Local cofactor features are the highest orientation-collapsed scores within
150 bp for:

| Factor | JASPAR motif |
|---|---|
| PATZ1 | `MA1961.2` |
| TFAP2C | `MA0814.3` |
| POU2F2 | `MA0507.3` |
| E2F1 | `MA0024.3` |

A cofactor is represented by a retained-hit indicator at score `-1` plus a
nonlinear continuous score excess above that floor.

## Pair Definition

A tandem partner is a distinct, non-overlapping `MA0861.2` alignment whose
edge-to-edge gap from the anchor is at most 20 bp. Shifted overlapping
alignments are not tandem partners. Orientation records at the same partner
span are collapsed into one locus and marked ambiguous if both strands are
reported.

This produces five anchor classes:

- `singleton`;
- `tandem_same_orientation`;
- `tandem_opposite_orientation`;
- `tandem_orientation_ambiguous`; and
- `tandem_mixed_orientation` for multiple resolved partner orientations.

This is an external tandem-alignment layer. It is distinct from the two
half-sites within one p53-family response element, which still require the
separate `hit_architecture` layer. A tandem class is compatible with one or
more TP73 complexes but cannot establish their protein stoichiometry.

## Direct Association

The table reports mean support across the six samples before discarding
concordant observations. `Anti-only fraction` is calculated among discordant
anti-p73/control observations.

| Pair class | Anchors per sample | Anti support | Control support | Difference | Anti-only fraction |
|---|---:|---:|---:|---:|---:|
| Singleton | 380,082 | 21.223% | 14.747% | +6.476 pp | 0.5706 |
| Tandem, opposite | 10,281 | 24.938% | 15.611% | +9.326 pp | 0.5998 |
| Tandem, same | 13,487 | 24.758% | 15.164% | +9.594 pp | 0.6102 |
| Tandem, ambiguous | 13,326 | 27.993% | 15.332% | +12.661 pp | 0.6539 |
| Tandem, mixed | 547 | 30.225% | 16.789% | +13.437 pp | 0.6331 |

The ordering is suggestive, but not uniform across samples. The anti-minus-
control difference is positive in 3 of 6 samples for singletons, 4 of 6 for
same- and opposite-orientation tandems, and 5 of 6 for ambiguous and mixed
tandems. Sample-specific differences range from negative to strongly positive,
so the mean must not be read as a universal effect.

The six sample ranges and the sample/fold win counts below are descriptive.
Nearby loci, repeated sample observations of one locus, and cross-validation
fits share information, so neither samples nor the 30 sample/fold cells are
independent replicates. This exploratory run reports no confidence interval or
P value for the small architecture increment.

The added TP73 score interval `[-1,0)` preserves the same ordering. Mean
anti-minus-control differences are +6.425 percentage points for singletons,
+9.683 for same-orientation tandems, +8.412 for opposite-orientation tandems,
+9.938 for ambiguous tandems, and +10.091 for mixed tandems. Thus the weaker
TP73 candidates still contain pair-architecture information and are useful for
learning, even if a later storage policy is stricter.

## Primary Held-Out Prediction

Five folds assign 5 Mb genomic blocks cyclically by coordinate.
All six sample copies of an anchor remain in one fold. Primary metrics are
macro-means over the 30 sample/fold cells. The cells provide balanced
diagnostic summaries, not 30 independent tests.

| Model | ROC AUC | Average precision | Log loss |
|---|---:|---:|---:|
| Sample only | 0.5000 | 0.5730 | 0.5289 |
| TP73 score | 0.5174 | 0.5889 | 0.5283 |
| Pair architecture only | 0.5145 | 0.5835 | 0.5279 |
| TP73 + pair architecture | 0.5272 | 0.5988 | 0.5273 |
| Cofactors only | 0.5916 | 0.6237 | 0.5218 |
| TP73 + cofactors | 0.5951 | 0.6277 | 0.5211 |
| TP73 + pair + cofactors | **0.5997** | **0.6331** | **0.5201** |
| TP73 + pair x cofactors | 0.5992 | 0.6322 | 0.5201 |

Relative to TP73 alone, adding pair architecture gains `0.0098` AUC and
`0.0099` average precision. Relative to TP73 plus the four cofactors, adding
pair architecture gains `0.0046` AUC and `0.0054` average precision. Adding the
four cofactors to TP73 plus pair architecture gains `0.0724` AUC and improves
AUC and average precision in all 30 sample/fold cells.

The best additive model reaches AUC `0.5962` for singletons, `0.5981` for
same-orientation tandems, `0.5929` for opposite-orientation tandems, and
`0.6145` for orientation-ambiguous tandems. The mixed class reaches only
`0.5277` on 1,161 discordant observations and is too small for a stable
class-specific conclusion.

### Contiguous-span sensitivity

A second run divides the observed chromosome-1 coordinate range into five
equal-width contiguous spans instead of distributing 5 Mb blocks cyclically.
All copies of an anchor still remain together. This makes each fold a long
spatial extrapolation and lowers the absolute cofactor-model performance.

| Model | ROC AUC | Average precision | Log loss |
|---|---:|---:|---:|
| Sample only | 0.5000 | 0.5626 | 0.5375 |
| TP73 score | 0.5180 | 0.5782 | 0.5370 |
| Pair architecture only | 0.5126 | 0.5714 | 0.5364 |
| TP73 + pair architecture | 0.5264 | 0.5856 | 0.5359 |
| Cofactors only | 0.5784 | 0.6059 | 0.5295 |
| TP73 + cofactors | 0.5826 | 0.6098 | 0.5289 |
| TP73 + pair + cofactors | **0.5867** | **0.6138** | **0.5278** |
| TP73 + pair x cofactors | 0.5863 | 0.6132 | 0.5278 |

Pair architecture improves AUC over TP73 score in 28 of 30 diagnostic cells
and over TP73 plus cofactors in 23 of 30. The interaction model improves over
the additive model in only 14 of 30 cells and loses `0.00038` mean AUC. This
sensitivity supports the additive feature choice, while the AUC reduction of
about `0.0130` for the full model shows why contiguous or chromosome-held-out
validation is needed for final performance claims.

## Cofactor Combinations

The additive model's held-out predictions remain compatible with the earlier
cofactor findings. The common PATZ1 + TFAP2C + E2F1 combination has the highest
mean predicted anti-p73-only fraction in singletons (`0.605`, n = 271,663),
same-orientation tandems (`0.632`, n = 10,159), opposite-orientation tandems
(`0.624`, n = 8,904), and ambiguous tandems (`0.664`, n = 10,190).

The exploratory heat map displays combinations with at least 100 observations.
Among those displayed cells, combinations dominated by POU2F2 are lower:
POU2F2 alone is `0.486` in singletons (n = 5,038), and TFAP2C + POU2F2 is
`0.540`, `0.524`, and `0.581` in same (n = 1,459), opposite (n = 893), and
ambiguous (n = 1,920) tandems.

These are predictions for observed candidate populations, not controlled
causal contrasts. Cofactor presence is correlated with the continuous scores,
sequence composition, repeat context, and one another. In particular, TFAP2C
is retained near almost every anchor at this floor. The failure of the explicit
interaction model to improve held-out metrics means there is not yet evidence
that separate pair-class-specific cofactor slopes are warranted. The additive
model was selected and displayed using the same cross-validation predictions;
there is no outer model-selection split, so the combination display is
exploratory rather than an independent estimate of the chosen model.

## Reproduce

Build pair features directly from sparse TP73 Parquet, without rescanning the
genome:

```sh
scripts/build_motif_context.py \
  --motif-hits 'TP73_SPARSE/task_data/*/tables/jaspar2026/motif_hit/**/*.parquet' \
  --output RUN/context --anchor-motif MA0861.2 \
  --motif-set-id jaspar2026_core_nonredundant \
  --genome-id homo_sapiens_grch38_ensembl113_primary \
  --anchor-minimum-score -5 --partner-minimum-score -1 \
  --score-mode log2_relative_risk --pseudocount 1 --chrom 1 \
  --capture-flank 500 --context-flank 150 --tandem-flank 20 \
  --threads 4 --memory-limit 4GB --max-temp-size 20GB
```

Generate one exact feature Parquet for each cofactor with
`cutandrun_score_calibration --feature-parquet`, as described in
[`cutandrun_motif_score_calibration.md`](cutandrun_motif_score_calibration.md),
then run:

```sh
scripts/evaluate_tp73_pair_stratified_model.R \
  --pair-parquet RUN/context/tables/jaspar2026/tp73_pair_feature/chrom=1/data_0.parquet \
  --patz1-parquet RUN/patz1_features.parquet \
  --tfap2c-parquet RUN/tfap2c_features.parquet \
  --pou2f2-parquet RUN/pou2f2_features.parquet \
  --e2f1-parquet RUN/e2f1_features.parquet \
  --output-prefix RUN/tp73_pair_stratified --context-floor -1
```

Repeat the evaluator with `--fold-mode contiguous` and a different output
prefix for the contiguous-span sensitivity. The default `interleaved` mode
retains the primary cyclic 5 Mb block scheme.

The evaluator writes exact formulas, run configuration, sample/fold metrics,
the coordinate extent and anchor count of every fold, pair-class metrics,
cofactor-combination summaries, and `_overview.png`. The completed ignored
local run is
`dry_runs/tp73_pair_stratified_chr1_cutrun_real_20260718/`.

## Limits And Next Tests

- Chromosome 1 and the 20-bp tandem rule were not selected on independent data.
- Interleaved and contiguous chromosome-1 folds do not replace held-out
  chromosomes or independent CUT&RUN.
- Formal uncertainty must cluster repeated observations by anchor or genomic
  block; the present point estimates and correlated cell win counts are not an
  inferential test.
- Shared feature slopes may hide TA/DN- or sample-specific relationships. Test
  sample-by-feature interactions or stratified fits before biological claims.
- The additive-versus-interaction choice needs nested validation when its
  selected predictions are used for a final model assessment.
- Orientation ambiguity is a property of motif scoring, not direct evidence of
  an ambiguous protein arrangement.
- Repeat/Alu class, low complexity, GC, mappability, accessibility, promoter,
  TSS, intron, and CUT&RUN depth remain potential covariates.
- The next genome-wide model should keep the shared additive architecture,
  validate on held-out chromosomes, and add internal half-site architecture
  before claiming a relationship to TP73 quaternary structure.

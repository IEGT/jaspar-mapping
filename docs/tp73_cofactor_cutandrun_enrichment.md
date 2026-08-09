# TP73 Cofactor Enrichment Across CUT&RUN Depth

## Scope

This analysis asks whether a threshold-qualified cofactor motif within 150 bp
of a putative TP73 site is associated with strict anti-p73 CUT&RUN immersion,
beyond the same association in the matched control track. It also describes
how that contrast changes across TP73 scores and stronger CUT&RUN depth tiers.

The outcome remains technical. Strict immersion means that a merged
positive-coverage component extends beyond both ends of the motif-model
alignment. It is not proof of direct protein binding, recruitment, or causal
cofactor activity.

## Cofactor Score Classes

For a cofactor operating point `T` and negative reference `N`, the classes are:

```text
positive:      best context score >= T
negative:      best context score < N, or no hit retained at the -1 source floor
intermediate:  N <= best context score < T
```

The inequalities are deliberate. A score exactly equal to `N` is not in the
strict negative class. Intermediate anchors are excluded from the corresponding
positive-versus-negative comparison.

Two cumulative negative references are reported:

- `< -1` is the primary sensitivity contrast. It holds the negative set fixed
  while the positive operating point varies and does not call the `[-1,0)`
  scores negative.
- `< 0` is the secondary compatibility contrast with the original score-zero
  analysis. When `T = 0`, it partitions every anchor into `< 0` and `>= 0`
  without an intermediate class.

The checked-in
[`tp73_chr1_cofactor_operating_points_v1.tsv`](../thresholds/tp73_chr1_cofactor_operating_points_v1.tsv)
contains the nine positive thresholds discussed so far. They were selected by
the historical threshold-complement analysis and are used here as fixed,
prespecified operating points. They are not relabelled as thresholds selected
under the new fixed-reference design.

## Descriptive Estimand

For antibody track `a`, motif class `M`, and one depth-tier event `Y`, report:

```text
RR_a = P(Y = 1 | M positive) / P(Y = 1 | M negative)
specificity ratio = RR_anti-p73 / RR_matched-control
```

Thus a log2 specificity ratio above zero is anti-p73 enrichment and a value
below zero is depletion relative to the matched control. The companion
difference-in-differences is:

```text
[P(Y | positive, anti) - P(Y | negative, anti)]
  - [P(Y | positive, control) - P(Y | negative, control)]
```

This conditioning matches the primary model; it does not use the transposed
`P(motif | depth)` quantity. Raw event probabilities are retained. Relative
risks use Jeffreys `Beta(0.5,0.5)` smoothed probabilities so sparse secondary
cells do not emit infinite log ratios.

The TP73 score strata are half-open:

```text
[-5,-1), [-1,0), [0,1), [1,2), [2,5), [5,10), [10,15), [15,+Inf)
```

An additional `all` row is the prespecified unstratified result. The current
310,782-anchor chromosome-1 evidence file starts at TP73 score 0, so its two
negative TP73 strata are empty. They become informative only after rebuilding
the anchor evidence from the wider TP73 production store.

## CUT&RUN Depth Tiers

`effective_max_depth` is zero unless strict immersion holds. The first tier is
therefore strict immersion (`depth > 0`). Five stronger tiers use the 50th,
75th, 90th, 95th, and 99th percentiles among positive effective depths,
calculated separately within every anti-p73 and control track. The anti and
control events consequently compare within-track ranks rather than raw depths
from libraries of different size.

CUT&RUN depth is discrete and often equals one. Several quantiles can resolve
to the same event set. The depth-tier manifest records the cutoff, achieved
event fraction, and whether it duplicates the previous tier; duplicate rows are
kept rather than pretending the requested ranks were separable.

## Primary Matched Test

There is one primary test per motif, at its supplied positive threshold and the
`< -1` negative reference. For each anchor and biological sample pair, retain
only discordant strict-immersion outcomes:

```text
outcome = 1  anti-p73 immersed, matched control not immersed
outcome = 0  matched control immersed, anti-p73 not immersed
```

The logistic model contains sample fixed effects, a natural spline of TP73
score, and the cofactor positive indicator. Conditioning on the discordant
anti/control pair makes the cofactor coefficient the
`cofactor presence x antibody` interaction: a cofactor main effect shared by
both tracks cancels. This is the model counterpart of the descriptive
specificity ratio, not a pooled accessibility association.

Uncertainty uses a sandwich covariance clustered by 5 Mb genomic blocks. The
block is much wider than the 150 bp context and expected short-range CUT&RUN
autocorrelation. The result reports the number of effective blocks. Benjamini-
Hochberg adjustment covers only the motifs supplied in that invocation; the
nine-motif result is a panel-level exploratory adjustment, not an all-JASPAR
q-value.

## Current Limits

The present anchor evidence does not carry GC, mappability, repeat/Alu class,
accessibility, or GFP coverage. The matched control reduces generic coverage
confounding but does not make those covariates unnecessary. They are recorded
as unavailable in every primary row. The current thresholds were selected and
tested on chromosome 1, no null motifs are included, and no independent
chromosome or data set is held out. These limits preclude a biological or
causal claim from this first panel.

Near-ubiquitous motifs are flagged when either anchor class is below the
declared 1% minimum. Their descriptive rows remain available, but a later model
should use count or continuous score rather than forcing an unstable presence
indicator.

## Run

The local nine-motif chromosome-1 command is:

```sh
scripts/analyze_tp73_cofactor_enrichment.R \
  --anchor-evidence \
    dry_runs/tp73_multifactor_chr1_cutrun_real_20260717/patz1_features.parquet \
  --cofactor-maxima \
    dry_runs/jaspar2026_chr1_convenient_thresholds_20260721/all_cofactor_maxima_v4.parquet \
  --thresholds thresholds/tp73_chr1_cofactor_operating_points_v1.tsv \
  --output-prefix RUN/tp73_cofactor_enrichment
```

The evaluator writes:

- `_class_counts.tsv`: independently inspectable positive, negative, and
  intermediate counts for both cumulative references and every TP73 stratum;
- `_depth_tier_manifest.tsv`: per-track quantile cutoffs and achieved support;
- `_descriptive.tsv`: sample-level probabilities, relative risks, specificity
  ratios, and difference-in-differences;
- `_macro_summary.tsv`: compact sample-macro descriptive summaries;
- `_primary_occupancy.tsv`: one matched occupancy interaction per motif with
  block-clustered uncertainty and panel-level BH adjustment; and
- `_run_config.tsv`: paths, boundaries, geometry, provenance, and declared
  inference limits.

The synthetic test needs no genome or JASPAR download:

```sh
bash tests/test_tp73_cofactor_enrichment.sh
```

For all JASPAR motifs, use one motif per restart-safe cluster task, combine the
primary rows, and apply BH once across the complete prespecified motif family.
Depth and TP73-score strata remain secondary until a separate validation set
is available.

The production orchestration, exact input inventory, restart contract, and
combined DuckDB/Parquet query surface are documented in
[`jaspar2026_chr1_all_cofactor_enrichment_run.md`](jaspar2026_chr1_all_cofactor_enrichment_run.md).

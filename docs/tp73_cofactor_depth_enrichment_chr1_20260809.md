# TP73 Cofactor And CUT&RUN Depth On Chromosome 1

## Result

The first fixed-reference enrichment run evaluates nine cofactor motifs around
310,782 chromosome-1 TP73 anchors. Cofactor maxima use schema-v4 interval-edge
geometry within 150 bp. Positive thresholds are the prespecified historical
operating points; the primary negative reference is strict score `< -1`.

The primary model conditions on discordant anti-p73/matched-control immersion,
adjusts for sample and a TP73-score spline, and clusters uncertainty over 5 Mb
genomic blocks. Its odds ratio estimates the cofactor-by-antibody interaction,
not cofactor prevalence pooled across both antibodies.

| Factor | JASPAR | Positive | Positive anchors | Negative `< -1` | Matched OR (95% block CI) |
| --- | --- | ---: | ---: | ---: | ---: |
| E2F1 | `MA0024.3` | `>= 1` | 41.01% | 36.35% | 1.287 (1.223-1.354) |
| SP1 | `MA0079.5` | `>= 3` | 47.20% | 31.68% | 1.434 (1.358-1.515) |
| REST | `MA0138.3` | `>= 0` | 34.19% | 56.42% | 1.252 (1.208-1.298) |
| POU2F2 | `MA0507.3` | `>= 0` | 40.03% | 53.83% | 0.719 (0.685-0.756) |
| KLF14 | `MA0740.2` | `>= 4` | 43.76% | 15.55% | 1.597 (1.494-1.707) |
| TCF7 | `MA0769.3` | `>= 6` | 69.51% | 0.23% | not fitted: negative class too small |
| POU4F1 | `MA0790.2` | `>= 5` | 70.71% | 3.09% | 0.359 (0.306-0.422) |
| TFAP2C | `MA0814.3` | `>= 5` | 61.38% | 2.56% | 1.664 (1.558-1.777) |
| PATZ1 | `MA1961.2` | `>= 2` | 56.47% | 25.79% | 1.545 (1.459-1.637) |

Eight motifs pass the prespecified 1% class-support gate. Their block-clustered
confidence intervals exclude one in this exploratory panel. The p-values and
panel-level BH values are retained in the output, but threshold selection on
the same chromosome, absent confounders, and absent null motifs make them
unsuitable as confirmatory evidence.

## Two Negative References

The cumulative `< 0` reference increases negative support as intended:

| Factor | Negative `< -1` | Negative `< 0` | Overall log2 specificity, `< -1` | Overall log2 specificity, `< 0` |
| --- | ---: | ---: | ---: | ---: |
| E2F1 | 36.35% | 47.63% | +0.362 | +0.312 |
| SP1 | 31.68% | 35.02% | +0.560 | +0.545 |
| REST | 56.42% | 65.81% | +0.326 | +0.301 |
| POU2F2 | 53.83% | 59.97% | -0.497 | -0.466 |
| KLF14 | 15.55% | 22.26% | +0.794 | +0.714 |
| TCF7 | 0.23% | 0.42% | -1.390 | -1.257 |
| POU4F1 | 3.09% | 4.59% | -1.148 | -1.056 |
| TFAP2C | 2.56% | 5.53% | +1.060 | +0.968 |
| PATZ1 | 25.79% | 31.29% | +0.712 | +0.655 |

Every direction is preserved. The `< 0` rows therefore provide the requested
historical comparison without redefining the intermediate `[-1,T)` anchors as
negative in the primary `< -1` analysis. TCF7 remains too close to ubiquitous
for a defensible binary primary model under either reference. Its negative
class rises only from 721 to 1,311 anchors.

## CUT&RUN Depth

The descriptive specificity contrast becomes larger in absolute magnitude for
stronger CUT&RUN events. The table below gives the mean sample-level log2
specificity ratio for the primary `< -1` reference. Positive values indicate
anti-p73 enrichment over matched control; negative values indicate depletion.

| Factor | Immersed | Positive-depth q90 | q95 | q99 |
| --- | ---: | ---: | ---: | ---: |
| E2F1 | +0.362 | +1.019 | +1.237 | +1.748 |
| SP1 | +0.560 | +1.673 | +2.007 | +2.660 |
| REST | +0.326 | +0.753 | +0.864 | +1.081 |
| POU2F2 | -0.497 | -1.298 | -1.608 | -2.410 |
| KLF14 | +0.794 | +2.155 | +2.539 | +2.801 |
| TCF7 | -1.390 | -2.960 | -3.584 | -4.187 |
| POU4F1 | -1.148 | -2.887 | -3.566 | -4.679 |
| TFAP2C | +1.060 | +2.676 | +2.757 | +1.927 |
| PATZ1 | +0.712 | +1.912 | +2.238 | +2.848 |

All six sample pairs agree in direction for every unstratified entry shown.
TFAP2C remains enriched at q99 but is the one non-monotonic tail. The q99 cells
are sparse and descriptive; they should not be interpreted as separate
hypothesis tests.

Depth ranks are normalized independently within every track. Many median and
75th-percentile cutoffs equal depth one and duplicate strict immersion. The
90th-99th percentile rows provide most of the actual intensity separation, and
the manifest reports each achieved tail fraction rather than assuming that a
nominal percentile split was possible.

## TP73 Score Strata

The current anchor table contains TP73 scores from 0 upward. Across `[0,1)`,
`[1,2)`, `[2,5)`, `[5,10)`, and `[10,15)`, all nine motifs retain their overall
direction in all six samples except occasional single-sample reversals for
KLF14, PATZ1, SP1, TCF7, and TFAP2C. The mean magnitude is broadly stable rather than
being confined to one TP73 score interval. The `[15,+Inf)` tail is less stable,
and TCF7 remains underpowered.

The planned TP73 source floor reaches -5, but this older CUT&RUN anchor file
does not. Negative TP73-score strata must therefore be revisited after the
anchor evidence is rebuilt from the production scan; the evaluator already
retains their empty rows so that omission is visible.

## Interpretation

This result supports the qualitative control pattern seen earlier:

- TFAP2C, KLF14, and PATZ1 have the strongest positive associations; SP1 is
  also strongly positive.
- POU4F1 has a stronger negative association than POU2F2 but a much smaller
  reference class. POU2F2 remains the more balanced binary discriminator.
- TCF7 is directionally negative but should move to a count or continuous-score
  representation instead of a binary presence call.
- Increasing anti-p73 CUT&RUN depth generally sharpens rather than erases the
  factor directions, which is the behavior expected of useful predictors.

These remain sequence-context associations with technical immersion. GC,
mappability, repeat/Alu class, accessibility, and GFP evidence are absent from
this anchor package. The next inferential analysis must add those covariates,
null motifs, and independent chromosome or data-set validation.

## Reproduce

The command and output contract are documented in
[`tp73_cofactor_cutandrun_enrichment.md`](tp73_cofactor_cutandrun_enrichment.md).
The ignored local result package is:

```text
dry_runs/tp73_cofactor_enrichment_chr1_20260809/
```

It contains six compact TSV outputs totaling about 1.7 MiB. The run is
computed directly from Parquet; no intermediate BED or TSV motif table is
materialized.

# TP73 and Cofactor Score Response

## Question and estimand

The fixed-threshold TP73/CUT&RUN cofactor analysis asks whether a cofactor
motif above its prespecified operating point is associated with strict
anti-p73 CUT&RUN immersion more strongly than the matched negative-control
track. Its primary adjusted odds ratio is intentionally one number per motif.

`scripts/analyze_tp73_cofactor_score_surface.R` adds two secondary questions:

1. Does the cofactor association change as the TP73 motif score increases?
2. At the same TP73 score, do stronger cofactor motifs have a different
   association from weaker cofactor motifs?

The outcome and matched design remain unchanged. For every anchor and sample,
only anti-p73/control-discordant strict-immersion pairs enter the logistic
model; `outcome = 1` denotes anti-p73 immersion. Therefore score-band
coefficients remain cofactor-by-antibody associations rather than generic
sequence accessibility effects.

## Score classes

The strict negative reference is fixed at cofactor score `< -1`, including an
absent sparse hit only when `source_score_floor <= -1`. Scores at or above `-1`
are divided into non-overlapping empirical bands with default lower boundaries
at their 50th, 75th, 90th, 95th, and 99th percentiles:

```text
negative, [-1,q50), [q50,q75), [q75,q90), [q90,q95), [q95,q99), [q99,+Inf)
```

Tied quantiles are collapsed. These are ordered score bands, not repeated
nested `score >= T` groups, so one physical anchor enters exactly one class.
The band containing the historical operating threshold is marked, but the
operating threshold does not determine the score bands.

The evaluator refuses the model when the negative reference is censored, when
the negative class covers less than 1% of anchors, or when any modeled band has
fewer than 100 anchors by default.

## Model and outputs

The fitted model is:

```text
anti_p73_wins ~ sample + chromosome
                + ns(TP73_score, df=4) * cofactor_score_band
```

Chromosome is omitted for a single-chromosome input. Sandwich uncertainty is
clustered by chromosome-qualified 5 Mb blocks. The output surfaces evaluate
the fitted contrasts at fixed TP73 scores plus the observed 5th, 25th, 50th,
75th, and 95th percentiles.

The evaluator writes:

- `_motif_status.tsv`: observability, support, model status, and the omnibus
  TP73-score-by-cofactor-band interaction;
- `_score_band.tsv`: exact score bounds, counts, and operating-point location;
- `_score_surface.tsv`: every band and TP73 grid point, with odds ratios versus
  the `< -1` reference, the immediately weaker band, and the weakest positive
  band;
- `_tp73_score_contrast.tsv`: the ratio between the cofactor odds ratio at the
  observed 95th and 5th TP73-score percentiles; and
- `_run_config.tsv`: inputs, included/excluded samples, geometry, model, and
  scientific status.

An odds ratio versus the negative class answers whether that score band is an
enriched or depleted sequence context. The strongest-versus-weakest-positive
contrast asks whether increasing the cofactor score adds information without
changing the negative reference. The TP73 high-versus-low ratio of odds ratios
asks whether increasing TP73 score modifies that cofactor association.

## Chromosome-1 pilot

The local pilot reuses 310,782 chromosome-1 anchors and nine rectangular
cofactor-maxima tables. It explicitly excludes both `skmel29_1` samples. The
older TP73 anchor set starts just above score zero, so this pilot cannot assess
TP73 scores in `[-1,0)` and is not a substitute for the whole-autosome run.
Threshold selection and evaluation also both used chromosome 1.

Eight motifs pass the 1% negative-reference guard. TCF7 `MA0769.3` does not:
only 721 anchors (0.23%) are in its `< -1` class, so it is reported as
`underpowered_score_reference` rather than fitted.

At the median TP73 score (2.53), the strongest empirical cofactor band compared
with the weakest positive band gives:

| Motif | Strongest band | OR | 95% CI | Interpretation |
|---|---:|---:|---:|---|
| E2F1 `MA0024.3` | `>= 10.41` | 1.963 | 1.479-2.605 | stronger enrichment |
| SP1 `MA0079.5` | `>= 15.06` | 4.054 | 2.990-5.496 | stronger enrichment |
| REST `MA0138.3` | `>= 13.26` | 1.215 | 0.999-1.477 | borderline increase |
| POU2F2 `MA0507.3` | `>= 16.37` | 0.934 | 0.697-1.250 | no clear score gain |
| KLF14 `MA0740.2` | `>= 14.52` | 4.424 | 3.293-5.943 | stronger enrichment |
| POU4F1 `MA0790.2` | `>= 14.74` | 0.521 | 0.426-0.637 | stronger depletion |
| TFAP2C `MA0814.3` | `>= 14.43` | 2.311 | 1.916-2.788 | stronger enrichment |
| PATZ1 `MA1961.2` | `>= 17.42` | 2.507 | 1.835-3.425 | stronger enrichment |

The omnibus TP73-by-cofactor score interaction is detectable for E2F1, SP1,
REST, POU2F2, KLF14, TFAP2C, and PATZ1 in this pilot; POU4F1 is not significant
(`p = 0.073`). This does not imply that stronger TP73 scores uniformly amplify
the cofactor association. For the strongest cofactor band, none of the direct
95th-versus-5th TP73-score contrasts excludes one after `skmel29_1` is removed.
The surface can therefore be non-flat without being monotonic from the lower to
upper TP73-score tail.

These p-values are exploratory and have not been adjusted across motifs,
bands, or displayed contrasts. The strongest 1% bands also have the smallest
counts and widest high-TP73 confidence intervals.

## Production use

The completed whole-autosome fixed-threshold package remains the authoritative
primary result. Its compact tables are under:

```text
/data/sm718/jaspar_mapping_runs/
  jaspar2026_grch38_tp73_cofactor_enrichment_autosomes_v1/
  final/cofactor_enrichment/
```

The score surface can be run immediately only for motifs whose context maxima
have `source_score_floor <= -1`. Motifs retained above that floor are censored,
not negative. Extending the score surface to every JASPAR motif therefore
requires rebuilding those motifs' 150 bp TP73 context maxima to a `-1` floor;
it does not require rerunning TP73 CUT&RUN evidence or the genome-wide
fixed-threshold analysis.

Example:

```sh
scripts/analyze_tp73_cofactor_score_surface.R \
  --anchor-evidence ANCHOR_EVIDENCE.parquet \
  --cofactor-maxima COFACTOR_MAXIMA.parquet \
  --thresholds MOTIF_THRESHOLDS.tsv \
  --output-prefix RESULTS/tp73_cofactor_score_surface \
  --negative-reference -1
```

The synthetic test uses known TP73/cofactor score interaction and a censored
motif and requires no genome download:

```sh
bash tests/test_tp73_cofactor_score_surface.sh
```

# TP73 Chromosome-1 CUT&RUN Threshold, 2026-07-13

## Conclusion

For `MA0861.2` genome-wide storage, use:

- **`log2_relative_risk`, pseudocount 1: minimum score `0.0`**
- **`log2_relative_risk`, pseudocount 0: minimum score `0.0`**
- **`log_odds`, pseudocount 1: minimum score about `38`**

The log-odds cutoff is `37.8` at the stored 0.2-bin resolution. It retains
almost the same number of chromosome-1 starts as risk score 0 and performs
slightly better. Keep `score_mode` and `pseudocount` with the threshold because
the numerical scales are not interchangeable.

Score 0 is an operational storage threshold, not a sharp biological boundary.
The association rises gradually. A less stringent, explicitly data-derived
lower bound is risk score `-2.0` with pseudocount 1 (`-2.2` without smoothing).
Rounding that lower bound to 0 removes another 44% of chromosome-1 candidates
while retaining essentially the same median anti-p73-versus-IgG specificity.

## Inputs

The analysis used:

- JASPAR 2026 TP73 motif `MA0861.2`, length 16;
- all 248,956,407 possible motif-model starts on GRCh38 chromosome 1;
- the maximum of the plus- and minus-orientation scores at each start;
- 230,478,786 starts with valid sequence scores after the existing `N` policy;
- nine anti-p73 bedGraph tracks: TA, DN, and GFP conditions for `saos2`,
  `skmel29_1`, and `skmel29_2`;
- the nine corresponding IgG BigWigs, exported for chromosome 1 without
  changing their positive depth values; and
- risk scores with pseudocounts 0 and 1 plus log odds with pseudocount 1.

The GFP tracks also use the p73 antibody. They are therefore condition
baselines, not negative-antibody controls. Threshold selection used the six TA
and DN anti-p73 tracks against their matched IgG tracks. Sample names are
preserved as supplied; this analysis does not infer that `_1` and `_2` are
interchangeable biological replicates.

## Coverage Rule

Positive bedGraph runs that overlap or directly adjoin are one coverage
component. A scored span `[start,start+16)` is supported only if:

~~~text
component_start < start AND component_end > start + 16
~~~

The effective depth is the maximum bedGraph depth within the 16-bp span when
that strict immersion rule holds, and zero otherwise. This implements the
requested technical rule while accepting that adjacent runs may jointly
immerse a span.

## Practical Criterion

For each threshold and TA/DN sample:

1. motif precision was normalized by that track's unthresholded support
   prevalence;
2. mean effective maximum depth was normalized by that track's unthresholded
   mean; and
3. each normalized anti-p73 value was divided by its matched normalized IgG
   value.

The declared lower-bound rule requires **every one of the six TA/DN
comparisons** to show at least a 5% anti-p73 advantage for both immersion and
depth. This 5% is a practical effect-size requirement, not a p-value. Treating
the hundreds of millions of overlapping motif spans as independent
observations would produce misleadingly small standard errors.

| Score definition | Lower bound | Starts retained | Operational cutoff | Starts retained |
|---|---:|---:|---:|---:|
| Risk, pseudocount 0 | -2.2 | 556,569 | 0.0 | 303,456 |
| Risk, pseudocount 1 | -2.0 | 550,839 | 0.0 | 310,782 |
| Log odds, pseudocount 1 | 33.4 | 668,335 | 37.8 | 315,829 |

At add-one risk score 0, 0.1348% of sequence-valid starts are retained, or one
retained start per approximately 801 total chromosome-1 alignment starts.
Across the six TA/DN comparisons:

- the least and median normalized immersion advantages over IgG are 1.053 and
  1.232; and
- the least and median normalized depth advantages over IgG are 1.084 and
  1.370.

At the equal-storage log-odds cutoff 37.8, those least/median ratios improve to
1.060/1.257 for immersion and 1.098/1.422 for depth.

## Score Comparison

Median binned ROC AUC for predicting strict immersion within each track was:

| Score definition | Anti-p73 | IgG |
|---|---:|---:|
| Risk, pseudocount 0 | 0.5134 | 0.5067 |
| Risk, pseudocount 1 | 0.5168 | 0.5093 |
| Log odds, pseudocount 1 | 0.5212 | 0.5116 |

The effect is real but modest because the label is positive coverage at every
chromosome start, not a curated set of isolated p73 peaks. Log odds ranks the
full score range best. Add-one smoothing improves global discrimination over
the unsmoothed risk score, while the two risk variants are nearly equivalent
around score 0. The earlier impression that unsmoothed risk had better recall
does not hold across the full chromosome and matched controls.

## Reproduce And Inspect

Build the streaming calibrator:

~~~sh
make cutandrun_score_calibration
./cutandrun_score_calibration --help
~~~

Convert a matched BigWig control chromosome without expanding the whole
genome:

~~~sh
scripts/export_bigwig_chrom_bedgraph.R \
  --input neg_skmel29_2_TA_R1.bigWig \
  --output RUN/control_bedgraph/neg_skmel29_2_TA_chr1.bedGraph \
  --chrom 1
~~~

Pass each plus/minus dense Parquet pair and repeat `--coverage ID=FILE` for the
tracks analyzed together. The executable writes only compact score histograms,
threshold curves, summaries, and provenance:

~~~sh
./cutandrun_score_calibration \
  --plus-parquet PLUS.parquet \
  --minus-parquet MINUS.parquet \
  --coverage skmel29_2_TA=tp73_skmel29_2_TA_R1.clipped.bedGraph \
  --coverage skmel29_2_DN=tp73_skmel29_2_DN_R1.clipped.bedGraph \
  --output-dir RUN/risk_p1 \
  --motif-id MA0861.2 --motif-length 16 --chrom 1 \
  --score-mode log2_relative_risk --pseudocount 1 --bin-width 0.2
~~~

After the three anti-p73 and three matched `igg_*` result directories exist:

~~~sh
scripts/summarize_tp73_cutandrun_threshold.R --run-root RUN
scripts/build_tp73_cutandrun_calibration_duckdb.sh --run-root RUN
duckdb -readonly RUN/tp73_calibration.duckdb
~~~

The local completed run is under
`dry_runs/tp73_chr1_cutrun_real_20260713/` and is intentionally ignored by
Git. Its `tp73_calibration.duckdb` contains `threshold_curve`,
`score_histogram`, `calibration_summary`, `run_config`,
`threshold_aggregate`, and `threshold_recommendation`.

The three generated figures are:

- `tp73_threshold_support.png`;
- `tp73_threshold_depth.png`; and
- `tp73_threshold_storage_tradeoff.png`.

Before using the cutoff as a publication claim, repeat the calibration on
held-out chromosomes or a held-out biological dataset. Chromosome 1 selected
this operational threshold; it must not also be presented as independent
validation.

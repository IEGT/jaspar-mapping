# TP73 Plus PATZ1 Chromosome-1 CUT&RUN Analysis, 2026-07-14

## Conclusion

PATZ1 adds information to TP73 sequence score on chromosome 1. Keep the
existing TP73 `log2_relative_risk >= 0` storage threshold as the primary rule,
but store the strongest local PATZ1 score as a separate architecture feature.

For each TP73 alignment span scoring at least 0, requiring the strongest PATZ1
alignment within 150 bp center-to-center to score at least 0:

- retains 206,513 of 310,782 TP73 starts (66.45%);
- raises the median normalized anti-p73-versus-IgG strict-immersion signal by
  8.39% relative to a TP73-only selection with exactly the same expected
  storage count; and
- raises the corresponding median maximum-depth signal by 15.89%.

All six TA/DN comparisons improve in that matched-size comparison. The weakest
improvements are only 0.81% for immersion and 1.81% for depth, so this is not
yet a reason to discard every TP73 >= 0 occurrence lacking PATZ1 support.

## Definition Of "Together"

This analysis follows the historical rule: for every inspected TP73 alignment,
select the single highest-scoring PATZ1 alignment in the local window. It uses:

- TP73 `MA0861.2`, length 16;
- PATZ1 `MA1961.2`, length 11;
- add-one `log2_relative_risk` for both motifs;
- the maximum score of the two motif orientations at each genomic start;
- the highest PATZ1 score whose alignment center is within +/-150 bp of the
  TP73 alignment center; and
- the leftmost PATZ1 alignment, then `+`, as the deterministic exact-score tie
  rule.

The PATZ1 score is therefore a maximum over as many as approximately 300 local
alignment starts. PATZ1 score 0 in this table does not have the same null
interpretation as score 0 for one predeclared alignment.

The selected PATZ1 alignment may overlap the TP73 alignment. This preserves the
historical maximum-score definition but means that the result alone cannot
distinguish independent PATZ1/TP73 motif architecture from two PSSMs extracting
information from overlapping sequence. A non-overlap sensitivity analysis is
the next focused check before interpreting this as cooperative binding syntax.

## Real CUT&RUN Comparison

The labels and controls are identical to the TP73-only calibration documented
in [`tp73_chr1_cutandrun_threshold_20260713.md`](tp73_chr1_cutandrun_threshold_20260713.md):

- nine anti-p73 tracks and nine matched IgG tracks;
- threshold decisions from the six TA and DN pairs;
- strict immersion of the complete TP73 alignment span; and
- maximum effective depth set to zero unless strict immersion holds.

At PATZ1 cutoff 0, the absolute anti-p73/IgG results across the six pairs are:

| Statistic | Least-enriched sample | Median |
|---|---:|---:|
| Normalized strict immersion | 1.0685 | 1.3789 |
| Normalized maximum depth | 1.1361 | 1.6311 |

For the matched-storage comparison, the TP73-only expectation is interpolated
within one tied 0.2-score bin. This avoids favoring the joint rule because the
nearest complete TP73 bin contains 212,422 rather than 206,513 starts.

| Joint / count-matched TP73-only | Least-improved sample | Median |
|---|---:|---:|
| Strict immersion | 1.0081 | 1.0839 |
| Maximum depth | 1.0181 | 1.1589 |

The strongest exploratory PATZ1 gate for which both measures improve in every
sample is PATZ1 >= 3.4. It retains 139,507 starts and has median matched-storage
improvements of 12.77% for immersion and 32.02% for depth. This cutoff was
selected on chromosome 1 and must not be treated as independently validated.

## Lower-Score Rescue

A piecewise rule can retain all TP73 >= 0 starts and selectively recover the
`-2 <= TP73 < 0` tier:

```text
TP73 >= 0
OR (-2 <= TP73 < 0 AND best local PATZ1 >= p)
```

At `p = 0`, this retains 477,289 starts, 53.58% more than TP73 >= 0 alone. At
the same expected storage count it improves median immersion by 3.38% and
median depth by 4.30% over TP73 score alone; the least improvements are 0.53%
and 1.00%. PATZ1 >= 4.6 gives a smaller 403,079-start set with median gains of
5.10% and 7.23%.

This demonstrates useful rescue information, but it does not overturn the
TP73 threshold of 0: rescue increases storage, and chromosome 1 is the same
data used to choose the rule.

## Geometry Check

For TP73 >= 0 and best PATZ1 >= 0:

- 41.86% of winning PATZ1 alignments have the same reported orientation as
  TP73;
- 42.32% of supported alignments have the same orientation;
- mean signed PATZ1-minus-TP73 center distance is -0.49 bp; and
- mean absolute center distance is 75.38 bp.

The nearly unchanged orientation fraction and centered signed distance do not
suggest an orientation preference in this first aggregate. Exact spacing and
overlap require the retained pair-level architecture layer rather than these
compact score histograms.

## Reproduce And Inspect

The completed ignored run is:

```text
dry_runs/tp73_patz1_chr1_cutrun_real_20260714/
```

The streaming executable accepts the two PATZ1 dense Parquet files through
`--context-plus-parquet` and `--context-minus-parquet`, plus
`--context-motif-length 11 --context-flank 150 --minimum-anchor-score -10`.
It writes `joint_score_histogram.tsv` and `joint_run_config.tsv` without a
chromosome-wide joined intermediate.

Reduce anti-p73 and IgG outputs with:

```sh
scripts/summarize_tp73_patz1_cutandrun_threshold.R \
  --run-root dry_runs/tp73_patz1_chr1_cutrun_real_20260714
```

Key outputs are:

- `tp73_patz1_gate_curve.tsv` and `tp73_patz1_gate_sample_curve.tsv`;
- `tp73_patz1_rescue_curve.tsv` and `tp73_patz1_rescue_sample_curve.tsv`;
- `tp73_patz1_selected.tsv`;
- `tp73_patz1_absolute_effect.png`; and
- `tp73_patz1_matched_storage.png`.

The original one-dimensional TP73 outputs from both joint passes are
byte-for-byte identical to the prior TP73-only outputs. This verifies that the
second motif stream did not alter TP73 scores or CUT&RUN containment.

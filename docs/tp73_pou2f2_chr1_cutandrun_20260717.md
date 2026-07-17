# TP73 With POU2F2: Chromosome-1 CUT&RUN Analysis, 2026-07-17

## Conclusion

The strongest nearby POU2F2 motif score is a strong negative predictor of TP73
CUT&RUN evidence on chromosome 1. This confirms the anticipated direction as a
sequence-context association: TP73 alignments near stronger POU2F2 alignments
are less likely to be strictly immersed in anti-p73 CUT&RUN coverage and have
lower effective depth than an equally sized set selected by TP73 score alone.

This does not yet establish direct antagonism by POU2F2. The POU2F2 motif may
also mark an AT-rich sequence or chromatin context that is unfavorable to TP73.
The POU2F2 score should enter the machine-learning layer as a continuous,
potentially negative feature rather than become a hard exclusion rule.

## Model And Analysis

The context model is JASPAR 2026
[`MA0507.3`](https://jaspar.elixir.no/matrix/MA0507.3/), a 13-bp human POU2F2
CAP-SELEX profile. The analysis otherwise matches the prior cofactors:

- TP73 `MA0861.2` anchors and POU2F2 use add-one `log2_relative_risk`;
- both motif orientations are scored densely across chromosome 1;
- each TP73 alignment receives the highest orientation-collapsed POU2F2 score
  within +/-150 bp center-to-center;
- nine anti-p73 and nine matched IgG bedGraph tracks are streamed directly;
- threshold decisions use the six TA and DN comparisons; and
- strict immersion and maximum effective depth are evaluated for the complete
  16-bp TP73 alignment span.

Every POU2F2 gate is compared with a TP73-only rule at exactly the same expected
storage count. Interpolation is confined to one tied 0.2-score bin. This is
essential here: POU2F2 is not called detrimental merely because a joint gate
selects fewer starts.

## Results

There are 310,782 chromosome-1 TP73 starts scoring at least 0. The table reports
percentage change from the exactly count-matched TP73-only expectation. Within
each cell, values are the minimum, median, and maximum across the six TA/DN
comparisons. Thus even the maximum being negative means that all six agree on
the detrimental direction.

| POU2F2 policy | Starts | Retained | Strict immersion change, min / median / max | Effective-depth change, min / median / max |
|---|---:|---:|---:|---:|
| POU2F2 >= 0 | 115,658 | 37.22% | -48.45% / -24.71% / -5.29% | -53.31% / -39.81% / -11.26% |
| Approximately 50%, POU2F2 >= -2.2 | 155,490 | 50.03% | -43.52% / -21.74% / -4.26% | -48.29% / -36.41% / -9.65% |
| Exploratory strongest detrimental gate, POU2F2 >= 7.8 | 22,421 | 7.21% | -63.22% / -34.59% / -11.28% | -69.93% / -46.43% / -15.39% |

The score-zero result is the most useful predeclared comparison. In absolute
anti-p73-versus-IgG terms, it has median normalized strict immersion 0.9895 and
median normalized depth 0.9695. Individual absolute ratios range from 0.9134
to 1.3672 for immersion and 0.7742 to 1.0369 for depth. A subset can therefore
remain enriched over IgG in one sample while still performing substantially
worse than TP73 score alone at the same selected count.

The threshold 7.8 was chosen from this chromosome and these six comparisons.
It is descriptive and must not be used as a validated biological or storage
cutoff.

## Geometry

At POU2F2 >= 0:

- 48.08% of winning POU2F2 alignments have the same reported orientation as
  TP73;
- the median same-orientation fraction among supported starts is 48.31%;
- mean signed POU2F2-minus-TP73 center distance is -0.39 bp; and
- mean absolute center distance is 75.11 bp.

This first aggregate provides no evidence for a preferred orientation or a
specific offset. The near-midpoint mean absolute distance is also compatible
with a broad neighborhood effect. Pair-level distance distributions are needed
before interpreting the result as a particular antagonistic architecture.

## Machine-Learning Interpretation

Retain the continuous highest local POU2F2 score for every TP73 anchor,
including negative scores. Also derive a `POU2F2_score_ge_0` indicator, distance,
orientation, overlap, and local-hit count. The model can then learn a negative
main effect and interactions with positive cofactors without deleting the
underlying TP73 observation.

Focused sensitivity checks should precede a mechanistic claim:

1. exclude POU2F2 alignments overlapping the TP73 alignment;
2. compare 20, 50, 100, and 150 bp context radii;
3. match or adjust for GC content, dinucleotide composition, promoter class,
   and accessibility;
4. compare other POU-family motifs and sequence-matched decoy motifs; and
5. validate on held-out chromosomes or an independent experiment.

## Reproduce And Inspect

The ignored run root is:

```text
dry_runs/tp73_pou2f2_chr1_cutrun_real_20260717/
```

Its dense plus and minus Parquet files each contain all 248,956,410 possible
13-bp chromosome-1 alignment starts. The anti-p73 and IgG provenance is pinned
in `joint_risk_p1/joint_run_config.tsv` and
`igg_joint_risk_p1/joint_run_config.tsv`. The ordinary TP73 score histograms,
threshold curves, and summaries from both passes are byte-for-byte identical
to the prior E2F1 context runs, confirming that adding POU2F2 did not alter the
anchor calculation.

Regenerate the summaries and plots with:

```sh
scripts/summarize_tp73_patz1_cutandrun_threshold.R \
  --run-root dry_runs/tp73_pou2f2_chr1_cutrun_real_20260717 \
  --context-label 'POU2F2 MA0507.3'
```

Important outputs are:

- `tp73_pou2f2_ma0507_3_gate_sample_curve.tsv`, with all six paired results;
- `tp73_pou2f2_ma0507_3_gate_curve.tsv`, with min/median/max summaries;
- `tp73_pou2f2_ma0507_3_selected.tsv`, with score-zero, retention, and
  direction-consistent exploratory rows; and
- `tp73_pou2f2_ma0507_3_matched_storage.png`, with the complete replicate
  envelope.

The positive and negative half-retention comparison is under:

```text
dry_runs/tp73_cofactor_comparison_20260717/
```

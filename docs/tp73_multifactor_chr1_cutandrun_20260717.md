# TP73 With PATZ1, TFAP2C, And POU2F2: Chromosome-1 Model, 2026-07-17

## Conclusion

The three local motif scores together are a **moderate but consistent**
predictor of TP73-specific CUT&RUN immersion among chromosome-1 TP73 candidates.
In held-out genomic blocks, the model containing TP73 plus all three context
scores reached mean within-sample ROC AUC `0.6039`, compared with `0.5210` for
TP73 score alone. Mean average precision increased from `0.5935` to `0.6394`.

POU2F2 was the strongest single context feature, but the three-score model was
better than any single cofactor model. The AUC and average-precision gains over
TP73 alone were positive in all 30 sample-by-fold comparisons.

This is exploratory evidence of complementary sequence context, not yet a
validated predictor or a mechanistic interaction claim.

## Prediction Target

The source set contains the 310,782 chromosome-1 TP73 `MA0861.2` alignment
starts with add-one `log2_relative_risk >= 0`. For each start and each of six
TA/DN samples, strict immersion of the complete 16-bp TP73 span was compared
between anti-p73 CUT&RUN and its matched control:

- outcome `1`: immersed in anti-p73 but not in the control;
- outcome `0`: immersed in the control but not in anti-p73; and
- excluded: immersed in both or in neither.

This gives 532,116 discordant observations: 330,382 anti-p73-only and 201,734
control-only. Sample identity is included as a nuisance term, but primary
metrics are calculated within each sample and held-out fold so differences in
sample prevalence cannot create motif discrimination.

## Features And Validation

For every TP73 candidate, the feature tables retain the exact highest
orientation-collapsed motif score within +/-150 bp center-to-center for:

- PATZ1 `MA1961.2`;
- TFAP2C `MA0814.3`; and
- POU2F2 `MA0507.3`.

All scores use JASPAR 2026, add-one `log2_relative_risk`, and both orientations.
Models are additive logistic regressions with four-degree-of-freedom natural
splines. Contiguous 5-Mb chromosome blocks are assigned round-robin to five
folds; all sample copies of one genomic coordinate remain in the same fold.
The primary values below are macro-averages over the 30 sample/fold cells.

| Model | ROC AUC | Average precision | Log loss |
|---|---:|---:|---:|
| Sample identity only | 0.5000 | 0.5754 | 0.5309 |
| TP73 | 0.5210 | 0.5935 | 0.5302 |
| Three cofactors, without TP73 score | 0.6003 | 0.6353 | 0.5214 |
| TP73 + PATZ1 | 0.5785 | 0.6173 | 0.5253 |
| TP73 + TFAP2C | 0.5760 | 0.6140 | 0.5259 |
| TP73 + POU2F2 | 0.5841 | 0.6311 | 0.5234 |
| TP73 + PATZ1 + TFAP2C + POU2F2 | **0.6039** | **0.6394** | **0.5206** |

Relative to TP73 alone, the complete model improves mean ROC AUC by `0.0830`,
average precision by `0.0459`, and log loss by `0.0096`. AUC and average
precision improve in all 30 paired sample/fold cells; log loss improves in 23.
The mean AUC gain is positive in all six samples and ranges from `0.0172` in
`skmel29_1_TA` to `0.1656` in `skmel29_2_DN`.

The three cofactors without continuous TP73 score already reach AUC `0.6003`.
This does **not** imply that TP73 sequence is unimportant: every observation
was admitted by the TP73 score-zero rule. It says that, after this entry filter,
variation in the three local context scores carries considerably more ranking
information than variation in TP73 score itself.

## Threshold Sensitivity: 0 Versus -1

Two threshold effects were tested separately. First, the original TP73
`score >= 0` anchors were held fixed while cofactor information was censored as
it would be in a sparse table retaining context matches at either `0` or `-1`.
For each cofactor, the model receives a retained-hit indicator and a nonlinear
function of score excess above the chosen floor. It receives no ordering among
scores below that floor. Labels, folds, and all TP73 scores are unchanged.

| Model | AUC at context floor 0 | AUC at context floor -1 | AUC change | AP change | Log-loss change |
|---|---:|---:|---:|---:|---:|
| TP73 + PATZ1 | 0.576150 | 0.577151 | +0.001002 | +0.000205 | -0.000131 |
| TP73 + TFAP2C | 0.575978 | 0.575966 | -0.000012 | -0.000017 | +0.000005 |
| TP73 + POU2F2 | 0.558439 | 0.563123 | +0.004684 | +0.001915 | -0.000394 |
| Three cofactors, without TP73 score | 0.591737 | 0.592753 | +0.001016 | +0.000614 | -0.000150 |
| TP73 + all three | 0.595544 | 0.596472 | +0.000928 | +0.000563 | -0.000152 |

Lowering the context floor to `-1` improves PATZ1 AUC in all five held-out
genomic folds and 26 of 30 sample/fold cells. POU2F2 improves in all five folds
and 29 of 30 cells, with the largest mean gain. TFAP2C is effectively unchanged:
it improves in only one fold and 12 of 30 cells. Among the 250,369 TP73 anchors
that contribute at least one discordant observation, the lower floor retains an
additional 14,217 PATZ1, 6,920 TFAP2C, and 14,469 POU2F2 local maxima. Thus the
`-1` extension adds useful information for PATZ1 and especially POU2F2, but not
detectably for TFAP2C in this analysis.

Second, the TP73 entry threshold itself was lowered to `-1`. The source cohort
grows from 310,782 to 417,723 orientation-collapsed candidates. After exclusion
of anchors with no discordant anti-p73/control observation in any sample, the
model contains 337,382 anchors and 715,678 observations. Performance in the
newly admitted TP73 `[-1,0)` stratum is:

| Model | ROC AUC | Average precision | Log loss |
|---|---:|---:|---:|
| TP73 | 0.496445 | 0.565326 | 0.522810 |
| TP73 + PATZ1 | 0.567284 | 0.600536 | 0.518785 |
| TP73 + TFAP2C | 0.560967 | 0.591800 | 0.519557 |
| TP73 + POU2F2 | 0.537976 | 0.579169 | 0.520884 |
| Three cofactors, without TP73 score | **0.584343** | **0.612028** | 0.517312 |
| TP73 + all three | 0.584275 | 0.611962 | **0.516941** |

TP73 score alone does not rank binding in this narrow low-score interval, but
the context motifs remain informative. Relative to TP73 alone, the complete
model improves AUC in all 30 sample/fold cells. Training on the expanded cohort
does not materially change predictions for the original `TP73 >= 0` stratum:
complete-model AUC is `0.596311` instead of `0.596472`, and log loss is
`0.522903` instead of `0.522896`.

These results favor `-1` over `0` as a sparse retention threshold when these
cofactors will be used for TP73 prediction. The gain is modest on the original
anchors, and the exact-score analysis remains an informative upper comparison
(complete-model AUC `0.6039`), so `-1` should be viewed as a storage/recall
compromise rather than proof that all lower scores are irrelevant.

## Chromosome-1 Storage And Timing Benchmark

A matched local scan on 2026-07-17 used TP73 `MA0861.2`, both orientations,
GRCh38 chromosome 1, JASPAR 2026, add-one `log2_relative_risk`, BED coordinates,
and skipped windows containing `N`. Only `--threshold` differed. The source
commit was `bb884295e73db50f8cef25176823822ef982659d` with the documented local
changes still uncommitted.

| Threshold | Strand-specific rows | Scanner wall time | CPU time | Raw BED | gzip BED | Lean Parquet |
|---:|---:|---:|---:|---:|---:|---:|
| -1 | 550,860 | 46.08 s warm (51.38 s cold) | 37.02 s warm | 39.517 MiB | 7.547 MiB | 3.784 MiB |
| 0 | 416,103 | 46.40 s warm | 36.75 s warm | 29.757 MiB | 5.720 MiB | 2.828 MiB |
| Change at -1 | +32.385% | -0.32 s (-0.69%) warm | +0.27 s (+0.73%) | +32.798% | +31.950% | +33.812% |

The initially requested order was `-1` followed by `0`. Because the first run
had a cold filesystem cache, `-1` was repeated into `/private/tmp`; its two BED
files are byte-identical to the cold run. The warm measurements show no material
threshold-dependent runtime difference: every chromosome window must be scored
in either case, and output volume is too small to dominate. The storage increase
tracks the 32.4% increase in retained rows.

For this initial benchmark, `gzip BED` is a preserved gzip copy of the native
BED files. `Lean Parquet` was materialized directly from the existing dense
Parquet source, partitioned by strand, with only `start` and float32 `score`
stored in each row; it was a lower-bound storage estimate. The scanner now has
direct `--sparse-parquet` output with the complete sparse physical schema.

The direct writer was then measured with the same inputs and scan settings. It
streams bounded Arrow record batches to two ZSTD Parquet files, one per strand,
without first writing BED or TSV:

| Threshold | Rows | Direct wall time | Direct CPU time | Direct Parquet |
|---:|---:|---:|---:|---:|
| -1 | 550,860 | 20.31 s | 19.74 s | 7.731 MiB |
| 0 | 416,103 | 19.49 s | 19.25 s | 6.085 MiB |
| Change at -1 | +32.385% | +4.21% | +2.55% | +27.052% |

Relative to the native outputs, direct Parquet is 80.44% smaller than raw BED
at threshold `-1` and 79.55% smaller at threshold `0`. It is 2.43% and 6.38%
larger than the corresponding gzip BED files. Unlike gzip BED, however, these
files are typed and directly queryable, and their generation avoids text
formatting plus a second compression pass. Direct elapsed time was 55.92%
lower than the warm native BED run at `-1` and 58.00% lower at `0`.

The direct files are about twice the size of `Lean Parquet` because that
lower-bound estimate deliberately kept only `start` and `score`. The production
physical schema additionally stores `end` and `pwm_relative_score`, uses int64
coordinates, and reserves a nullable `matched_seq` field. Column metadata shows
that the all-null `matched_seq` field is negligible; the useful endpoint and
relative-score columns account for nearly all of the difference.

The ignored local benchmark is in
`dry_runs/tp73_chr1_threshold_benchmark_20260717/`; its
`benchmark_summary.tsv` records exact timings and byte counts. The direct
benchmark and its corresponding summary are in
`dry_runs/tp73_chr1_direct_sparse_parquet_benchmark_20260717/`.

## Reproduce And Inspect

The optional exact export is produced by:

```sh
./cutandrun_score_calibration \
  --plus-parquet TP73_PLUS.parquet --minus-parquet TP73_MINUS.parquet \
  --context-plus-parquet COFACTOR_PLUS.parquet \
  --context-minus-parquet COFACTOR_MINUS.parquet \
  --coverage SAMPLE=TRACK.bedGraph \
  --feature-parquet RUN/cofactor_features.parquet \
  --output-dir RUN/cofactor_calibration \
  --motif-id MA0861.2 --motif-length 16 \
  --context-motif-id COFACTOR_ID --context-motif-length COFACTOR_LENGTH \
  --context-flank 150 --minimum-anchor-score 0 \
  --score-mode log2_relative_risk --pseudocount 1 \
  --context-pseudocount 1 --chrom 1
```

After making PATZ1, TFAP2C, and POU2F2 feature Parquets:

```sh
scripts/evaluate_tp73_multifactor_model.R \
  --patz1-parquet RUN/patz1_features.parquet \
  --tfap2c-parquet RUN/tfap2c_features.parquet \
  --pou2f2-parquet RUN/pou2f2_features.parquet \
  --output-prefix RUN/tp73_multifactor \
  --context-floor -1
```

Omit `--context-floor` to retain the exact local-score analysis. With a floor,
the evaluator also writes `_context_availability.tsv`, `_context_basis.tsv`, and
separate metrics for the original `TP73 >= 0` and added `TP73 [-1,0)` strata.

The completed ignored local run is
`dry_runs/tp73_multifactor_chr1_cutrun_real_20260717/`. Its compact outputs
include `tp73_multifactor_model_metrics.tsv`,
`tp73_multifactor_sample_fold_metrics.tsv`, and
`tp73_multifactor_comparison.png`.

## Limits

- Chromosome 1, the motifs, the score-zero entry rule, and the 150-bp radius
  were not selected on independent data.
- Five block folds reduce local sequence leakage but do not replace validation
  on held-out chromosomes or an independent CUT&RUN experiment.
- Strict immersion is a technical operational label, not proof of direct
  binding or absence.
- The present model learns nonlinear main effects but no motif interactions,
  distances, orientations, promoter class, accessibility, or GC composition.
- The strong sample-specific differences make pooled AUC (`0.7805` for the
  complete model) a secondary and potentially misleading summary; the
  within-sample macro result is the primary one.

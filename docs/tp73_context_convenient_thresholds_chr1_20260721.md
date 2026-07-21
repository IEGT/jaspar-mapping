# TP73-context convenient motif thresholds on chromosome 1

## Scope

This exploratory threshold set asks whether the strongest occurrence of each
neighboring motif within 150 bp of a putative TP73 site improves prediction of
strict anti-p73 CUT&RUN immersion. It does not define a universal binding
threshold and does not change the scan storage floor of -1.

The calibration uses JASPAR 2026 motifs scored as log2 relative risks with a
pseudocount of 1. Distances are signed interval-edge distances; overlapping
motifs therefore have negative distances and abutting motifs have distance 0.
For every TP73 anchor, the maximum cofactor score at interval distance at most
150 bp is converted to an inclusive binary feature `score >= threshold`.

The analysis contains 310,782 chromosome-1 TP73 anchors and 532,116 discordant
anchor/sample observations (strict immersion in an anti-p73 sample only versus
its matched control only). Held-out predictions use five contiguous,
equal-width chromosome-1 folds. Candidate thresholds are the observed
non-negative integers. A candidate is eligible only when both the retained and
absent classes contain at least 1% of anchors.

The convenient threshold maximizes the held-out macro ROC AUC gain over the
TP73-score baseline, with ties resolved toward the lower threshold. The useful
range contains tested thresholds reaching at least 90% of the selected positive
gain.

## Results

| Motif | JASPAR | Threshold | Useful range | Anchors retained | Delta ROC AUC | Adjusted OR | Direction |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| E2F1 | `MA0024.3` | 1 | 0-1 | 41.01% | 0.017766 | 1.2166 | positive association |
| SP1 | `MA0079.5` | 3 | 0-5 | 47.20% | 0.030285 | 1.3303 | positive association |
| REST | `MA0138.3` | 0 | 0 | 34.19% | 0.019490 | 1.2288 | positive association |
| POU2F2 | `MA0507.3` | 0 | 0 | 40.03% | 0.030778 | 0.7364 | negative association |
| KLF14 | `MA0740.2` | 4 | 2-5 | 43.76% | 0.031968 | 1.3588 | positive association |
| TCF7 | `MA0769.3` | 6 | 5-6 | 69.51% | 0.022853 | 0.7646 | negative association |
| POU4F1 | `MA0790.2` | 5 | 3-6 | 70.71% | 0.042777 | 0.6529 | negative association |
| TFAP2C | `MA0814.3` | 5 | 5-7 | 61.38% | 0.037259 | 1.4371 | positive association |
| PATZ1 | `MA1961.2` | 2 | 0-4 | 56.47% | 0.035599 | 1.4192 | positive association |

All nine selected features improved ROC AUC in at least 25 of the 30
sample-fold summaries, and their raw association directions agreed across all
six anti-p73 samples. Those sample-fold summaries are correlated diagnostics,
not independent replicates. The odds ratios describe association with the
technical CUT&RUN immersion outcome; the direction labels must not be read as
causal facilitation or inhibition by the corresponding protein.

The focused values obtained earlier are recovered exactly: POU2F2 0, TCF7 6,
and POU4F1 5. REST receives a threshold of 0. KLF14 (`MA0740.2`) receives 4;
its useful 2-5 range indicates that the precise integer is not sharply
identified.

## Data and provenance

KLF14 was added with a direct sparse-Parquet scan of both chromosome-1 strands
at the -1 storage floor. The scan emitted 1,005,324 plus-strand and 1,010,980
minus-strand records in 9.70 seconds of recorded scan time. The resulting
nine-motif context table contains 2,797,038 rows: one maximum for each of nine
motifs at each TP73 anchor.

The populated registry is:

```text
dry_runs/jaspar2026_chr1_convenient_thresholds_20260721/
  tables/jaspar2026/motif_score_threshold/part-000000.parquet
```

Its threshold-set ID is `tp73_chr1_cutrun_context_roc_auc_v1`. Every row pins
the motif and target accessions, score configuration, context geometry,
candidate-selection rule, retained fraction, evidence dataset, source metric
checksum, source commit, and dirty-worktree state. `sql/schema.sql` exposes the
complete registry as `motif_score_threshold` and positive recommendations as
`motif_convenient_threshold`.

## Interpretation and next validation

These values are convenient defaults for simple TP73-context queries, not
filters to apply destructively to the motif-hit store. Continuous scores remain
the preferred ML input, and all hits down to -1 remain available for alternate
thresholds, distance bands, pair architectures, and genomic strata.

Selection and reported performance currently use the same chromosome-1 folds.
The next threshold-set version should freeze these transformations and test
them on held-out chromosomes or independent anti-p73 CUT&RUN data, with
anchor-clustered uncertainty. Thresholds should also be re-estimated separately
where tandem TP73 architecture, promoter position, or distance band is treated
as a biological stratum.

# TP73 With E2F1 And TFAP2C: Chromosome-1 CUT&RUN Analysis, 2026-07-14

## Conclusion

The strongest nearby TFAP2C motif score is a robust predictor of TP73
CUT&RUN evidence beyond the TP73 sequence score alone. The effect is comparable
to PATZ1 and is reproduced by all three current JASPAR 2026 TFAP2C matrices.

E2F1 is a weaker candidate. Its local score improves the median result, most
clearly for maximum CUT&RUN depth, but the approximately half-retention rule
does not improve strict immersion in every one of the six TA/DN comparisons.

These are sequence-context associations with anti-p73 CUT&RUN. They do not show
that TFAP2C or E2F1 physically occupies the selected alignment, nor that either
factor causes TP73 binding. Keep the TP73 `log2_relative_risk >= 0` storage rule
and store local cofactor scores as separate architecture features.

## Motif Models

The E2F1 model is JASPAR [`MA0024.3`](https://jaspar.elixir.no/matrix/MA0024.3/),
a 12-bp HT-SELEX model.

JASPAR 2026 contains three distinct current TFAP2C profiles, so all three were
run before inspecting their relative performance:

- [`MA0814.3`](https://jaspar.elixir.no/matrix/MA0814.3/), 9 bp, ChIP-seq/ReMap,
  is the primary model;
- [`MA0524.3`](https://jaspar.elixir.no/matrix/MA0524.3/), 11 bp, is the
  HT-SELEX DBD sensitivity model; and
- [`MA0815.1`](https://jaspar.elixir.no/matrix/MA0815.1/), 13 bp, is an
  HT-SELEX motif variant.

The accessions remain separate. Their raw scores and score cutoffs are not
interchangeable.

## Common Analysis

For E2F1 and every TFAP2C model, the analysis uses:

- dense chromosome-1 scores on both orientations;
- add-one `log2_relative_risk` for TP73 and the context motif;
- TP73 `MA0861.2`, with every start scoring at least -10 entering the compact
  joint histogram;
- the highest orientation-collapsed context score whose alignment center is
  within +/-150 bp of the TP73 alignment center;
- nine anti-p73 and nine matched IgG tracks;
- the six TA and DN pairs for threshold comparisons; and
- strict immersion of the complete 16-bp TP73 alignment span plus maximum
  effective CUT&RUN depth.

Each joint rule is compared with a TP73-only rule at exactly the same expected
number of stored starts. The comparator interpolates within one tied 0.2-score
bin, as documented for the PATZ1 analysis. This controls the major advantage
that any additional threshold obtains merely by retaining fewer TP73 starts.

## Count-Matched Comparison

The fairest cross-model comparison chooses the context-score threshold nearest
to retaining half of the 310,782 chromosome-1 TP73 starts scoring at least 0.
The table reports improvement over the exactly count-matched TP73-only rule.

| Context model | Context cutoff | Starts | Strict immersion, weakest | Strict immersion, median | Maximum depth, weakest | Maximum depth, median |
|---|---:|---:|---:|---:|---:|---:|
| PATZ1 `MA1961.2` | 2.8 | 155,643 | +0.91% | +11.45% | +2.98% | +28.28% |
| E2F1 `MA0024.3` | 0.0 | 154,605 | -0.12% | +3.45% | +0.24% | +15.38% |
| TFAP2C `MA0814.3` primary | 6.0 | 153,187 | +1.61% | +10.65% | +3.33% | +25.36% |
| TFAP2C `MA0524.3` sensitivity | 4.8 | 156,212 | +1.40% | +11.35% | +3.67% | +26.59% |
| TFAP2C `MA0815.1` sensitivity | -3.0 | 155,723 | +1.42% | +11.48% | +3.72% | +27.33% |

The agreement among TFAP2C matrices is more informative than their small
ordering differences. In particular, the negative cutoff for `MA0815.1` is a
feature of that model's score distribution, not evidence that negative scores
have a universal TFAP2C interpretation.

For the primary `MA0814.3` model, the exploratory strongest all-sample gate is
TFAP2C >= 6.8. It retains 133,256 starts and has median count-matched gains of
12.41% for strict immersion and 29.61% for depth; the weakest gains are 1.93%
and 4.00%. This cutoff was selected from the same chromosome and samples and
therefore is descriptive, not a validated operational threshold.

The reducer's exploratory strongest E2F1 gate that improves both measures in
all six comparisons is E2F1 >= -0.2. It retains 163,053 starts and yields
median gains of 3.59% for strict immersion and 15.11% for depth. Its weakest
gains are only 0.25% and 0.79%, respectively.

## Interpretation

TFAP2C should enter the promoter-architecture layer as an accession-specific
feature. The primary feature should use `MA0814.3`; the other two scores should
remain available for sensitivity analysis and later model selection. E2F1 is
also worth retaining, but its present evidence is not as uniform.

The result has four important limits:

1. The highest score is selected from roughly 300 nearby alignment starts, so
   it is a local maximum rather than the score of one predeclared position.
2. The winning context alignment may overlap the TP73 alignment. A non-overlap
   sensitivity run is required before calling the signal cooperative motif
   syntax.
3. GC content, promoter accessibility, and other sequence architecture can
   raise both motif scores and CUT&RUN coverage. Matched decoy motifs and
   chromosome-held-out modeling are needed to estimate specificity.
4. The same chromosome and six comparisons were used to explore cutoffs and
   report effects. Independent chromosomes or experiments must validate any
   fixed threshold.

## Reproduce And Inspect

The ignored run roots are:

```text
dry_runs/tp73_e2f1_chr1_cutrun_real_20260714/
dry_runs/tp73_tfap2c_ma0814_3_chr1_cutrun_real_20260714/
dry_runs/tp73_tfap2c_ma0524_3_chr1_cutrun_real_20260714/
dry_runs/tp73_tfap2c_ma0815_1_chr1_cutrun_real_20260714/
```

Each `joint_risk_p1/joint_run_config.tsv` and
`igg_joint_risk_p1/joint_run_config.tsv` records the exact dense Parquet,
coverage, motif, flank, pseudocount, bin-width, and tie-rule inputs.

Reduce an individual pair of anti-p73 and IgG runs with, for example:

```sh
scripts/summarize_tp73_patz1_cutandrun_threshold.R \
  --run-root dry_runs/tp73_tfap2c_ma0814_3_chr1_cutrun_real_20260714 \
  --context-label 'TFAP2C MA0814.3'
```

The reducer's historical filename remains PATZ1-specific, but its labels and
calculations are cofactor-generic. Recreate the combined table and plot with:

```sh
scripts/compare_tp73_cofactor_summaries.R \
  --input 'PATZ1 MA1961.2=dry_runs/tp73_patz1_chr1_cutrun_real_20260714/tp73_patz1_selected.tsv' \
  --input 'E2F1 MA0024.3=dry_runs/tp73_e2f1_chr1_cutrun_real_20260714/tp73_e2f1_selected.tsv' \
  --input 'TFAP2C MA0814.3 (primary)=dry_runs/tp73_tfap2c_ma0814_3_chr1_cutrun_real_20260714/tp73_tfap2c_ma0814_3_selected.tsv' \
  --input 'TFAP2C MA0524.3 (sensitivity)=dry_runs/tp73_tfap2c_ma0524_3_chr1_cutrun_real_20260714/tp73_tfap2c_ma0524_3_selected.tsv' \
  --input 'TFAP2C MA0815.1 (sensitivity)=dry_runs/tp73_tfap2c_ma0815_1_chr1_cutrun_real_20260714/tp73_tfap2c_ma0815_1_selected.tsv' \
  --output-prefix dry_runs/tp73_cofactor_comparison_20260714/tp73_cofactor_comparison
```

This writes `tp73_cofactor_comparison.tsv` with the score-zero,
approximately-half-retention, and exploratory policies, plus
`tp73_cofactor_comparison.png` for the fair half-retention comparison.

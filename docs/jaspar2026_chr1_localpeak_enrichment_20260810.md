# JASPAR 2026 Cofactors at TP73 Local Peaks (2026-08-10)

## Result package

The schema-7 chromosome-1 analysis is published at:

```text
/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_localpeak_enrichment_v3/final/cofactor_enrichment_localpeak_v3
```

It combines 2,632 non-TP73 JASPAR 2026 motifs with 305,492 physical TP73
local-peak alignment spans selected at score `>= -1`. Cofactor-positive classes
use the fixed historical operating-point registry; the primary negative class
is absent or has a context maximum `< -1`. Six matched TP73/control comparisons
cover TA, DN, and GFP in SaOS-2 and SK-Mel-29. `skmel29_1` is excluded by the
source-controlled evidence manifest.

The final package contains one DuckDB database and seven Zstandard-compressed
Parquet tables. All 2,632 tasks completed. Seventeen motifs without a registry
recommendation retain a null recommendation and use threshold zero only as an
explicit descriptive fallback.

## Primary occupancy screen

The primary model conditions on discordant TP73/control strict-immersion
outcomes, adjusts for sample and TP73 score, and clusters uncertainty by 5 Mb
genomic block. Of 2,632 motifs, 1,571 pass the declared class-support gate:

| Direction | Estimable motifs | All-JASPAR BH q <= 0.05 |
| --- | ---: | ---: |
| TP73 enriched | 782 | 723 |
| TP73 depleted | 789 | 744 |
| Not estimable | 1,061 | 0 |

The large significant fraction reflects the large anchor population and broad
sequence/context structure. It must not be read as 1,467 independently
validated cofactors. This is a chromosome-1 hypothesis screen using operating
points also derived on chromosome 1; GC, mappability, repeat class, chromatin
accessibility, and independent chromosomes remain validation work.

## Prespecified cofactors

The nine motifs selected before the all-JASPAR result retain the directions
seen in the historical analysis. Odds ratios compare the cofactor-positive and
strict-negative classes in the matched TP73/control model.

| Factor | Threshold | Positive anchors | Odds ratio (95% CI) | BH q | Direction |
| --- | ---: | ---: | ---: | ---: | --- |
| E2F1 | 1 | 41.1% | 1.309 (1.239, 1.384) | 3.83e-21 | enriched |
| KLF14 | 4 | 42.4% | 1.601 (1.484, 1.727) | 3.97e-33 | enriched |
| PATZ1 | 2 | 55.9% | 1.553 (1.449, 1.665) | 2.66e-34 | enriched |
| POU2F2 | 0 | 40.4% | 0.731 (0.693, 0.772) | 1.01e-28 | depleted |
| POU4F1 | 5 | 72.6% | 0.366 (0.313, 0.429) | 6.55e-35 | depleted |
| REST | 0 | 33.6% | 1.270 (1.214, 1.330) | 4.45e-24 | enriched |
| SP1 | 3 | 46.7% | 1.455 (1.364, 1.552) | 3.83e-29 | enriched |
| TCF7 | 6 | 71.0% | not estimable | NA | class support below gate |
| TFAP2C | 5 | 60.9% | 1.691 (1.538, 1.860) | 1.53e-26 | enriched |

TCF7 remains useful in descriptive tables, but its near-ubiquity makes the
binary operating-point model underpowered under the predeclared support rule.

## Increasing CUT&RUN depth

The descriptive specificity ratio compares the positive-versus-negative risk
ratio in TP73 CUT&RUN with the corresponding risk ratio in the matched control.
For the prespecified panel, stronger depth usually magnifies the strict-
immersion direction:

| Factor | Strict immersion log2 ratio | 99th-percentile depth log2 ratio |
| --- | ---: | ---: |
| E2F1 | +0.335 | +1.506 |
| KLF14 | +0.631 | +2.611 |
| PATZ1 | +0.597 | +2.536 |
| POU2F2 | -0.399 | -1.810 |
| POU4F1 | -0.976 | -3.805 |
| REST | +0.291 | +0.774 |
| SP1 | +0.472 | +2.357 |
| TCF7 | -1.474 | -3.492 |
| TFAP2C | +0.877 | +0.861 |

TFAP2C is non-monotone: it reaches `+2.463` at the 90th-percentile tier and
then weakens. The 50th-percentile event is identical to strict immersion in
these discrete depth tracks and is retained explicitly in the database.
Across all motifs with descriptive estimates, the median log2 specificity
ratio moves from `-0.258` at strict immersion to `-1.190` at the 99th-
percentile tier; the number with a negative ratio rises from 1,554 to 1,811 of
2,559. These depth summaries are descriptive and do not have the primary
model's inferential status.

The retained TP73 anchors begin at score `-1`, so the predeclared `[-5,-1)`
stratum is correctly empty. All strata from `[-1,0)` upward are populated. At
strict immersion, the across-motif median log2 specificity ratio becomes less
negative as TP73 score rises (`-0.216` in `[-1,0)` and `-0.091` at `>= 15`),
but this is an aggregate descriptive trend rather than a motif-level threshold
test.

## Historical comparison

The earlier all-JASPAR package used the historical score-zero TP73 anchor
population. It had 1,614 estimable motifs: 771 enriched and 843 depleted. In
the schema-7 local-peak rerun:

- 1,562 motifs are estimable in both analyses;
- 1,536 of those retain direction and 26 change direction;
- the two adjusted log-odds estimates correlate at `r = 0.9942`; and
- the prespecified nine-motif panel retains every estimable direction.

The notable significant reversal is TP63 (`MA0525.2`): odds ratio `1.122` in
the historical analysis and `0.923` in the local-peak analysis. It should be
treated as a sensitivity finding, especially because TP53-family motifs can
overlap the TP73 anchor definition.

## Provenance correction

The first consolidation inherited the historical default run ID even though
its inputs and output path were the new local-peak run. Commit `3e685a5`
separated immutable plan identity from publication identity. Slurm job
`5702503` then rebuilt only the combined package from the exact completed task
outputs, publishing:

```text
run_id       jaspar2026_chr1_tp73_localpeak_enrichment_v3
plan_run_id  jaspar2026_chr1_tp73_cofactor_enrichment_v2
```

The predecessor was retained. `EXCEPT ALL` comparisons report zero differing
rows for every one of the seven tables; the repair changes provenance, not
scientific values. The corrected manifest records metric source commit
`4ddec92` and finalization commit `3e685a5`.

See
[`jaspar2026_chr1_all_cofactor_enrichment_run.md`](jaspar2026_chr1_all_cofactor_enrichment_run.md)
for the execution/query contract and
[`h3k4me3_chr1_production_20260810.md`](h3k4me3_chr1_production_20260810.md)
for the separate GFP-referenced H3K4me3 change analysis.

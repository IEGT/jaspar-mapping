# Chromosome-1 H3K4me3 Change Induction (2026-08-10)

## Result package

The first schema-7 production analysis completed as Slurm job `5689195` in
26 minutes 3 seconds with peak RSS 7.3 GiB. It was executed from source commit
`24f2c3b` and published atomically at:

```text
/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_h3k4me3_production_v4/final
```

The package contains 305,492 distinct chromosome-1 TP73 local-peak alignment
spans selected at score `>= -1`. The source contained 305,528 strand-aware
records; tied orientations at an identical physical span were collapsed by
their maximum score. The evidence Parquet has exactly one row per span.

The signal table is rectangular: 305,492 anchors x two accepted experimental
series x three conditions x five windows = 9,164,760 rows. `saos2` and
`skmel29_2` are analyzed separately; `skmel29_1` is excluded by the input
manifest and contributes zero rows. Nine prespecified cofactor motifs produce
2,749,428 context-maximum rows. The complete package has 17 checksummed files
and occupies about 173 MiB.

## Estimand

For anchor `i`, series `r`, and condition `c`, the windowed mark is:

```text
H(i,r,c) = log2((H3K4me3_area(i,r,c) + 1) /
                 (input_area(i,r,c) + 1))
```

The primary outcome is `H(TA) - H(GFP)`; `H(DN) - H(GFP)` is the predeclared
parallel analysis. The primary window is the two flanks 150-1,000 bp from the
TP73 alignment, chosen using GFP-only profiles. A positive cofactor effect
means that this condition-minus-GFP change is larger in cofactor-positive than
in cofactor-negative anchors after adjustment for TP73 motif score. It does
not mean that the absolute mark increased.

Positive cofactors meet their prespecified score threshold. The primary
negative reference is a context maximum below `-1`; intermediate anchors are
excluded. Standard errors are clustered over 5 Mb chromosome-1 blocks.

## Global changes

Mean H3K4me3/input levels confirm that the two systems cannot be pooled as
replicates:

| Series | GFP | TA | DN | TA - GFP | DN - GFP |
| --- | ---: | ---: | ---: | ---: | ---: |
| SaOS-2 | -4.406 | -4.413 | -5.951 | -0.007 | -1.545 |
| SK-Mel-29 | -1.377 | -1.148 | -1.722 | +0.229 | -0.345 |

The numbers of anchors strictly immersed by TP73 CUT&RUN but not by the
matched negative-control component are:

| Series | GFP | TA | DN |
| --- | ---: | ---: | ---: |
| SaOS-2 | 29,076 | 23,842 | 6,297 |
| SK-Mel-29 | 12,219 | 71,149 | 82,952 |

These are technical occupancy labels, not independent biological replicates.

## TA change

The table gives the adjusted cofactor-positive minus cofactor-negative effect
on `TA - GFP`, with block-clustered 95% confidence intervals.

| Cofactor | SaOS-2 effect (95% CI) | SK-Mel-29 effect (95% CI) |
| --- | ---: | ---: |
| E2F1 | -0.053 (-0.114, +0.009) | +0.071 (+0.029, +0.114) |
| KLF14 | -0.040 (-0.121, +0.040) | +0.193 (+0.117, +0.269) |
| PATZ1 | -0.040 (-0.104, +0.025) | +0.185 (+0.115, +0.255) |
| POU2F2 | +0.042 (-0.006, +0.090) | -0.136 (-0.184, -0.088) |
| POU4F1 | +0.247 (+0.120, +0.374) | -0.537 (-0.709, -0.365) |
| REST | +0.011 (-0.040, +0.062) | +0.144 (+0.098, +0.190) |
| SP1 | -0.043 (-0.103, +0.016) | +0.175 (+0.110, +0.240) |
| TCF7 | not estimable | not estimable |
| TFAP2C | -0.091 (-0.245, +0.062) | +0.224 (+0.077, +0.371) |

REST is the only estimable motif with the same positive direction in both
systems, and its SaOS-2 interval includes zero. The other seven estimable
motifs reverse direction. This reproduces the pilot's central observation:
the strong SK-Mel-29 associations do not behave as replicated TA effects in
SaOS-2. Replacing every score-`>= 0` span with schema-7 local peaks at `>= -1`
does not change any of the pilot's per-series TA directions.

## DN change

The corresponding `DN - GFP` contrasts are directionally consistent across
the two systems for every estimable motif:

| Cofactor | SaOS-2 effect (95% CI) | SK-Mel-29 effect (95% CI) |
| --- | ---: | ---: |
| E2F1 | +0.204 (+0.134, +0.274) | +0.157 (+0.110, +0.204) |
| KLF14 | +0.386 (+0.294, +0.478) | +0.332 (+0.265, +0.400) |
| PATZ1 | +0.345 (+0.260, +0.430) | +0.329 (+0.266, +0.392) |
| POU2F2 | -0.232 (-0.283, -0.182) | -0.242 (-0.286, -0.199) |
| POU4F1 | -0.576 (-0.690, -0.463) | -0.716 (-0.815, -0.618) |
| REST | +0.147 (+0.087, +0.208) | +0.174 (+0.134, +0.214) |
| SP1 | +0.306 (+0.239, +0.374) | +0.309 (+0.256, +0.362) |
| TCF7 | not estimable | not estimable |
| TFAP2C | +0.527 (+0.404, +0.650) | +0.444 (+0.305, +0.582) |

Overall DN H3K4me3 is lower than GFP in both systems. Positive values here
therefore describe a smaller loss at cofactor-positive anchors, while negative
values describe a larger loss. They must not be called DN-induced H3K4me3.
TCF7 has only 449 primary negative anchors and fails the declared support gate.

## TP73-confirmation interaction

TP73 confirmation is post-treatment and may mediate H3K4me3 change, so its
interaction is secondary and descriptive. In the TA comparison, KLF14, PATZ1,
REST, and SP1 interactions are positive in both systems, while POU4F1 is
negative in both. Only selected SK-Mel-29 contrasts survive within-family BH
correction; SaOS-2 intervals remain broad. DN interactions are generally weak:
none survives its within-family correction in both systems.

The package also retains gained/lost H3K4me3 occurrence summaries and
TP73-binding-state strata. These are supporting views rather than alternative
primary outcomes.

## Interpretation and limits

The subsequent chromosome-2 holdout materially sharpens this interpretation.
SK-Mel-29 TA effects transfer for all eight estimable motifs, while no SaOS-2
TA effect excludes zero on chromosome 2. The chromosome-1 SaOS-2 POU4F1 effect
and the apparent cross-system REST direction do not replicate as supported
effects. See
[`jaspar2026_chr2_heldout_validation_20260810.md`](jaspar2026_chr2_heldout_validation_20260810.md)
for the fixed-threshold comparison.

- The production anchor definition confirms rather than resolves the
  system-specific TA result seen in the pilot.
- DN cofactor contrasts are remarkably consistent, but they stratify the
  magnitude of a global DN-associated loss; they do not demonstrate induction.
- There is one track per condition. The two accepted cell systems are the only
  biological replication level, and no pooled replicate p-value is reported.
- Motif presence is sequence-derived. The analysis associates motif context
  with a processed H3K4me3 change surface; it does not establish cofactor
  binding or causality.
- The exact upstream BigWig normalization remains pending provenance. The
  manifest records source file hashes, and absolute signal interpretation must
  remain conditional on that missing processing history.

See [`h3k4me3_cofactor_change.md`](h3k4me3_cofactor_change.md) for the full
method and restart contract, and
[`h3k4me3_gfp_window_calibration_chr1_20260809.md`](h3k4me3_gfp_window_calibration_chr1_20260809.md)
for the GFP-only window choice.

# Chromosome-1 GFP-Referenced H3K4me3 Cofactor Pilot (2026-08-09)

## Scope

This is a structural and statistical pilot, not a publishable biological
result. It uses 310,782 orientation-collapsed TP73 motif-alignment spans with
log-risk-ratio score at least 0 on chromosome 1, the separately analyzed
`saos2` and `skmel29_2` experimental
series, and the nine prespecified motifs in
[`config/h3k4me3_chr1_pilot_cofactors_v1.tsv`](../config/h3k4me3_chr1_pilot_cofactors_v1.tsv).
`skmel29_1` was explicitly excluded and no file from that series was opened.

The primary window was predeclared from GFP-only profiles as the two flanks
150-1000 bp from each TP73 motif-alignment centre. For each condition, the
outcome is `log2((H3K4me3 area + 1)/(input area + 1))`; TA-minus-GFP and
DN-minus-GFP changes are modeled separately. Positive cofactors meet their
motif-specific inclusive threshold. The primary negative class is absent or
has a maximum context score below -1; intermediate anchors are excluded.

## Extraction

The persistent products are 9,323,460 window rows (142 MiB) and 310,782 strict
TP73/control evidence rows (2.9 MiB). The signal table is exactly rectangular:
two series x three conditions x five windows for every anchor. Source paths,
sizes, modification times, and SHA-256 values are recorded in the generated
track manifest.

Mean mark levels show why condition differencing and separate series models are
essential:

| Series | GFP | TA | DN | mean TA-GFP | mean DN-GFP |
|---|---:|---:|---:|---:|---:|
| SaOS-2 | -4.424 | -4.408 | -5.929 | +0.016 | -1.505 |
| SK-Mel-29 | -1.388 | -1.121 | -1.688 | +0.267 | -0.300 |

Strict TP73-positive/negative-control-negative anchor counts were respectively
29,877/27,231/7,892 in SaOS-2 GFP/TA/DN and
12,500/78,406/90,782 in SK-Mel-29. These labels describe strict immersion in
the processed tracks, not biological replication.

## Primary total association

The table reports the adjusted cofactor-positive minus cofactor-negative effect
on TA-minus-GFP H3K4me3 change, in log2-ratio units. `q` is BH-adjusted within
series, isoform, and negative-reference family.

| Cofactor | Threshold | SaOS-2 effect (q) | SK-Mel-29 effect (q) |
|---|---:|---:|---:|
| E2F1 | 1 | -0.060 (0.652) | +0.072 (0.0028) |
| KLF14 | 4 | -0.036 (0.652) | +0.263 (9.4e-10) |
| PATZ1 | 2 | -0.024 (0.652) | +0.229 (3.1e-8) |
| POU2F2 | 0 | +0.014 (0.652) | -0.155 (2.2e-8) |
| POU4F1 | 5 | +0.101 (0.652) | -0.589 (8.8e-13) |
| REST | 0 | +0.015 (0.652) | +0.164 (6.1e-9) |
| SP1 | 3 | -0.024 (0.652) | +0.218 (3.7e-9) |
| TCF7 | 6 | not estimable | not estimable |
| TFAP2C | 5 | -0.047 (0.652) | +0.179 (0.036) |

The headline is heterogeneity, not replication. Except for a small positive
REST contrast, TA total-association signs differ between the two systems. The
directions are stable across the five signal windows within each system for
most motifs, so this is not explained by the selected annulus alone. TCF7 has
only 721 primary negative anchors and fails the declared 0.5% class-support
gate. Its historical `< 0` comparison is estimable but rests on only 1,311
negative anchors and is not promoted to the primary analysis.

## TP73-confirmation interaction

TP73 confirmation is post-treatment and can mediate the response, so these
interactions are explicitly secondary and descriptive. A positive value means
that the cofactor contrast is larger among strict TP73-confirmed anchors than
among unconfirmed anchors.

| Cofactor | SaOS-2 interaction (q) | SK-Mel-29 interaction (q) |
|---|---:|---:|
| E2F1 | +0.055 (0.827) | +0.022 (0.709) |
| KLF14 | +0.286 (0.537) | +0.248 (0.0072) |
| PATZ1 | +0.105 (0.827) | +0.128 (0.081) |
| POU2F2 | -0.022 (0.881) | -0.088 (0.095) |
| POU4F1 | -0.317 (0.627) | -0.200 (0.081) |
| REST | +0.114 (0.632) | +0.055 (0.278) |
| SP1 | +0.227 (0.380) | +0.155 (0.0072) |
| TCF7 | not estimable | not estimable |
| TFAP2C | -0.059 (0.881) | -0.052 (0.786) |

Across all five windows, PATZ1 and SP1 interactions remain positive in both
systems, while POU2F2 remains negative in both. KLF14 is positive in all
SK-Mel-29 windows but has one near-zero sign reversal in SaOS-2. These are
useful hypotheses. They are not independent replications, and only the
SK-Mel-29 KLF14 and SP1 primary-window interactions survive their within-family
BH correction.

## Robustness and limits

- All five stored windows were evaluated with the final script. The historical
  `< 0` negative reference preserves the primary direction in 32 of 36
  estimable series/isoform/motif comparisons; the four changes are near-zero
  contrasts rather than stable reversals.
- Standard errors are clustered by 5 Mb chromosome-1 blocks. There is one R1
  track per condition and only two accepted biological systems; no pooled
  replicate P-value is reported.
- All-zero H3K4me3/input anchors are retained. The primary model does not
  condition on TP73 confirmation.
- The exact BigWig normalization, clipping, and channel semantics are not in
  the repository. Absolute effect interpretation and publication must wait for
  that provenance. The current manifest deliberately records
  `processed_bigwig_scale_pending_provenance`.
- Motif presence is sequence-derived. The analysis establishes association
  with a processed H3K4me3 change surface, not cofactor binding or causality.

The durable local outputs are under
`dry_runs/h3k4me3_cofactor_change_chr1_20260809/`; they are reproducible from
the source-controlled manifests and commands in
[`h3k4me3_cofactor_change.md`](h3k4me3_cofactor_change.md).

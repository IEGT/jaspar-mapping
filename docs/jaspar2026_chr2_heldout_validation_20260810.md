# Chromosome-2 Genomic Holdout Validation (2026-08-10)

## Design

Chromosome 2 was held out from the chromosome-1 threshold selection. The nine
prespecified cofactor score thresholds, the TP73 local-peak floor (`>= -1`),
the 150 bp context radius, the `flank_150_1000` H3K4me3 window, and every model
setting were carried over unchanged. This is therefore a genomic holdout under
the same CUT&RUN experiments. It tests transfer to an untouched sequence
region, but it is not a new biological replicate.

The corrected result package is:

```text
/data/sm718/jaspar_mapping_runs/jaspar2026_chr2_tp73_heldout_validation_v2/final
```

It was computed from source commit `d95c6ff` by Slurm job `5705248`. The
package was submitted with four CPUs, 32 GiB, and the `requeue` partition. It
completed in 21 minutes 50 seconds with peak RSS 7.4 GiB and occupies about
165 MiB on `/data`.

## Cardinality

The chromosome-2 package contains 311,385 physical TP73 local-peak alignment
spans. The evidence table has exactly one row per span. Its two accepted
experimental series, three conditions, and five windows yield 9,341,550
H3K4me3/input rows. Nine motifs yield 2,802,465 rectangular context-maximum
rows. `skmel29_1` remains excluded and contributes no rows.

Strict TP73 immersion without simultaneous negative-control immersion occurs
at the following numbers of anchors:

| Series | GFP | TA | DN |
| --- | ---: | ---: | ---: |
| SaOS-2 | 24,732 | 17,921 | 4,008 |
| SK-Mel-29 | 12,464 | 56,883 | 71,029 |

These are technical occupancy labels from the existing tracks, not independent
binding measurements.

## Occupancy Transfer

Positive classes use the chromosome-1-selected thresholds. Negative references
have no context hit or a maximum score strictly below `-1`; intermediate scores
are excluded. Odds ratios come from the same matched TP73/control model used
on chromosome 1.

| Factor | Threshold | Chr2 positive | Chr2 negative | Chr1 OR | Chr2 OR (95% CI) |
| --- | ---: | ---: | ---: | ---: | ---: |
| E2F1 | 1 | 38.6% | 39.41% | 1.309 | 1.284 (1.246, 1.323) |
| KLF14 | 4 | 38.1% | 19.77% | 1.601 | 1.511 (1.431, 1.596) |
| PATZ1 | 2 | 50.1% | 31.85% | 1.553 | 1.470 (1.399, 1.546) |
| POU2F2 | 0 | 44.8% | 48.73% | 0.731 | 0.757 (0.731, 0.784) |
| POU4F1 | 5 | 76.7% | 1.98% | 0.366 | 0.350 (0.319, 0.384) |
| REST | 0 | 31.1% | 59.89% | 1.270 | 1.237 (1.191, 1.284) |
| SP1 | 3 | 41.6% | 37.03% | 1.455 | 1.392 (1.333, 1.454) |
| TCF7 | 6 | 74.3% | 0.13% | not estimable | not estimable |
| TFAP2C | 5 | 55.1% | 3.31% | 1.691 | 1.476 (1.340, 1.624) |

All eight motifs estimable on both chromosomes retain their direction. The
chromosome-specific adjusted log-odds estimates correlate at `r = 0.9963`; the
mean absolute log-odds difference is 0.0526. Thus E2F1, KLF14, PATZ1, REST,
SP1, and TFAP2C transfer as enriched sequence contexts, while POU2F2 and
POU4F1 transfer as depleted contexts. All eight chromosome-2 panel-adjusted
tests pass `q < 0.05`.

Chromosome-1 q-values were corrected over all 2,632 motifs, whereas the
chromosome-2 q-values cover this prespecified nine-motif panel. The transfer
comparison therefore uses coefficients, confidence intervals, and direction;
it does not compare q-values across those different testing families.

TCF7 does not become estimable: only 390 chromosome-2 anchors (0.13%) are in
its primary negative class, compared with 449 on chromosome 1. Its binary
operating point remains unsuitable for this model.

The descriptive depth behavior also transfers. The table compares strict
immersion with the 99th-percentile positive-depth event on chromosome 2:

| Factor | Strict log2 specificity ratio | Q99 log2 specificity ratio |
| --- | ---: | ---: |
| E2F1 | +0.338 | +1.548 |
| KLF14 | +0.605 | +2.903 |
| PATZ1 | +0.559 | +3.088 |
| POU2F2 | -0.380 | -1.725 |
| POU4F1 | -1.117 | -3.575 |
| REST | +0.278 | +0.833 |
| SP1 | +0.455 | +2.107 |
| TCF7 | -1.389 | -4.567 |
| TFAP2C | +0.693 | +2.129 |

Every strict and Q99 direction agrees with chromosome 1. Unlike the
chromosome-1 Q99 summary, TFAP2C remains strongly positive at Q99 on
chromosome 2. Depth tiers are descriptive repeated thresholds on the same
tracks and are not independent tests.

## Global H3K4me3 Change

Input-normalized chromosome-2 H3K4me3 changes preserve the chromosome-1
system-level pattern:

| Series | Chr1 TA-GFP | Chr2 TA-GFP | Chr1 DN-GFP | Chr2 DN-GFP |
| --- | ---: | ---: | ---: | ---: |
| SaOS-2 | -0.007 | -0.037 | -1.545 | -1.644 |
| SK-Mel-29 | +0.229 | +0.259 | -0.345 | -0.370 |

The absolute input-normalized levels differ by chromosome, but the condition
changes are close. This supports treating the cofactor contrasts as modifiers
of a series-specific change rather than as absolute H3K4me3 abundance.

## H3K4me3 Cofactor Effects

The primary coefficient is cofactor-positive minus cofactor-negative in the
condition-minus-GFP input-normalized H3K4me3 change. The summary below compares
the eight estimable motifs between chromosomes.

| Isoform | Series | Direction concordance | Effect correlation | Chr1 CIs excluding 0 | Chr2 CIs excluding 0 |
| --- | --- | ---: | ---: | ---: | ---: |
| DN | SaOS-2 | 8/8 | 0.987 | 8/8 | 8/8 |
| DN | SK-Mel-29 | 8/8 | 0.997 | 8/8 | 8/8 |
| TA | SaOS-2 | 5/8 | -0.131 | 1/8 | 0/8 |
| TA | SK-Mel-29 | 8/8 | 0.994 | 8/8 | 8/8 |

The TA effects show the chromosome-1 and chromosome-2 estimates side by side:

| Factor | SaOS-2 chr1 | SaOS-2 chr2 | SK-Mel-29 chr1 | SK-Mel-29 chr2 |
| --- | ---: | ---: | ---: | ---: |
| E2F1 | -0.053 (-0.114, +0.009) | -0.031 (-0.081, +0.019) | +0.071 (+0.029, +0.114) | +0.095 (+0.038, +0.152) |
| KLF14 | -0.040 (-0.121, +0.040) | -0.014 (-0.078, +0.049) | +0.193 (+0.117, +0.269) | +0.308 (+0.234, +0.382) |
| PATZ1 | -0.040 (-0.104, +0.025) | -0.019 (-0.070, +0.032) | +0.185 (+0.115, +0.255) | +0.221 (+0.157, +0.285) |
| POU2F2 | +0.042 (-0.006, +0.090) | +0.003 (-0.038, +0.045) | -0.136 (-0.184, -0.088) | -0.150 (-0.206, -0.094) |
| POU4F1 | +0.247 (+0.120, +0.374) | -0.014 (-0.157, +0.129) | -0.537 (-0.709, -0.365) | -0.869 (-1.136, -0.602) |
| REST | +0.011 (-0.040, +0.062) | -0.016 (-0.063, +0.032) | +0.144 (+0.098, +0.190) | +0.157 (+0.096, +0.218) |
| SP1 | -0.043 (-0.103, +0.016) | -0.019 (-0.073, +0.035) | +0.175 (+0.110, +0.240) | +0.217 (+0.155, +0.279) |
| TFAP2C | -0.091 (-0.245, +0.062) | +0.027 (-0.087, +0.141) | +0.224 (+0.077, +0.371) | +0.274 (+0.156, +0.391) |

The SK-Mel-29 TA architecture transfers almost exactly: every direction and
every interval excluding zero is reproduced on chromosome 2. SaOS-2 has no
chromosome-2 interval excluding zero. The chromosome-1 positive POU4F1 result
does not transfer, and the small REST and TFAP2C sign changes remain compatible
with zero on both relevant SaOS-2 estimates. Consequently, no motif in this
panel has supported TA-associated H3K4me3 modulation across both cell systems
and both chromosomes.

DN effects transfer in direction, magnitude, and interval support for all
eight motifs in both systems. They still describe differential attenuation or
accentuation of a global DN-associated H3K4me3 loss; they are not evidence of
DN-induced H3K4me3.

## Secondary TP73 Interaction

The post-treatment cofactor-by-TP73-confirmation interaction is less stable,
as expected. For TA in SK-Mel-29, 7/8 directions transfer (`r = 0.849`); the
only reversal is TFAP2C, whose chromosome-1 estimate was close to zero. SaOS-2
TA interactions transfer in only 4/8 directions and no interval excludes zero
on either chromosome. These interactions remain descriptive and do not alter
the primary change conclusion.

## Interpretation

- The fixed chromosome-1 cofactor thresholds transfer exceptionally well for
  TP73 versus control occupancy on chromosome 2.
- Increasing CUT&RUN depth strengthens the same enriched/depleted directions,
  providing a useful held-out check on the descriptive depth response.
- The SK-Mel-29 TA H3K4me3 cofactor pattern and both DN patterns are genomic,
  not peculiar to chromosome 1.
- The SaOS-2 TA pattern is weak and does not provide a replicated cofactor
  effect. The earlier POU4F1 and REST observations should not be promoted to
  cross-system biological claims.
- This remains one set of experimental tracks. Chromosome 2 supplies disjoint
  genomic observations, not additional cell lines or biological replicates.
- GC, mappability, repeat/ALU class, accessibility, and the upstream BigWig
  normalization history remain unresolved covariates/provenance limits.

The next defensible expansion is to apply the fixed panel to additional
chromosomes and treat chromosome as the replication unit for genomic transfer.
The all-JASPAR chromosome-1 screen should remain hypothesis generation until
more motifs are evaluated outside chromosome 1.

## Publication Repair

The first completed package (`jaspar2026_chr2_tp73_heldout_validation_v1`) had
valid data and checksums, but three text sidecars retained their temporary
`attempts/...` paths after atomic directory promotion. It is retained as an
audit predecessor. Commit `d95c6ff` canonicalizes published output paths and
durable anchor/motif input locators before checksumming. The `v2` package is
the publication target; ephemeral node-local bedGraph paths remain explicitly
marked as non-persistent execution traces while their source BigWigs are
identified and checksummed in the track manifest.

The reproducibility comparison covered the three Parquet data tables, the GFP
metaprofile, five H3K4me3 result tables, and five occupancy result tables.
Thirteen of these fourteen scientific files are byte-identical between `v1`
and `v2`. The cofactor-maxima Parquet has a different physical serialization,
but both versions contain 2,802,465 rows and a bidirectional `EXCEPT ALL`
comparison finds zero differing rows. No corrected TSV or JSON sidecar retains
an `attempts/...` path, a scratch anchor path, or a scratch motif path.

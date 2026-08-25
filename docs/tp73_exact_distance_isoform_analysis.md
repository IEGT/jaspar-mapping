# Exact-distance TP73 isoform and cofactor analysis

## Question

This analysis asks whether the CUT&RUN association of a candidate motif changes
with its exact position relative to a TP73 motif, and whether that positional
response differs between full-length TAp73alpha and DNp73beta. In particular,
it is designed to retain:

- single preferred positions that move toward or away from TP73;
- paired response peaks separated by about one DNA helical turn;
- a paired pattern that is preserved between isoforms but translated as a
  whole; and
- TP73/cofactor relative orientation.

The first implementation is candidate-level. It processes one motif and one
chromosome at a time into compact sufficient statistics. It does not materialize
an anchor-by-motif-by-position table.

## Relationship to TF-COMB

[TF-COMB](https://github.com/loosolab/TF-COMB) is the methodological reference
for transaction support, direction-aware motif distance, smoothing, preferred
distances, orientation, and downstream pair-network reporting. This project
adapts those ideas to its existing TP73 anchor and matched CUT&RUN design; it
does not require TF-COMB at run time and does not copy its implementation.

There are two different transaction universes:

1. A future **unconditioned regulatory-region universe** can estimate global
   motif-pair support, confidence, lift, cosine, and background-corrected
   co-occurrence.
2. The present **TP73-conditioned universe** contains one transaction per
   physical TP73 anchor. It estimates the conditional support of a candidate
   motif and its association with anti-p73 versus negative-control support.

In the second universe every transaction contains TP73. Consequently,
TP73-to-candidate lift is always 1, reverse confidence is always 1, and cosine
is just the square root of support. Reporting those values would add no
information. The package therefore reports support and forward confidence,
labels the conditioning explicitly, and omits the degenerate metrics. This is
not a genome-wide motif-pair enrichment analysis.

## Geometry

All intervals are BED 0-based half-open intervals. Two complementary distances
are retained:

- `interval_distance_bp` is the signed edge distance already used for the
  overlap/adjacent/near/far bands. Abutting intervals have distance 0 and
  overlaps have negative distance.
- `center_offset_twice` is the signed motif-center offset in doubled-base
  coordinates. It preserves exact half-base offsets when one motif has odd
  length and the other has even length. `center_offset_bp` is its readable
  floating-point representation.

The genomic frame is positive toward increasing reference coordinates. The
TP73-oriented frame is positive in the direction of the best-supported TP73
motif strand. A TP73 span with equal best scores on both strands remains in the
genomic frame but is excluded from TP73-oriented curves. Cofactor loci found on
both strands at one physical span are counted once. A higher strand-specific
score can resolve their orientation; an exact tie remains ambiguous.

The primary exact-spacing reports use TP73-oriented motif centers. The signed
edge-distance bands remain the primary surface for direct comparison with the
earlier analysis. Neither coordinate is called a transcription factor's
"window" on DNA.

## Score classes

The source scan floor and positive operating threshold are independent:

- `strict_negative`: no physical locus retained at the source floor, normally
  score `< -1`;
- `intermediate`: a retained locus exists at the exact offset but its best
  score is below the declared positive threshold; and
- `positive`: the best score at that offset reaches the declared threshold.

Intermediate anchors are not silently used as negatives in the CUT&RUN odds
ratio. The exact-distance inventory remains zero-complete and reports all three
classes. This preserves the project's low-floor context contract even for a
motif whose convenient operating threshold is much higher.

`--source-score-floor` is mandatory. It records and applies the requested
floor, but it cannot restore hits omitted by an upstream scan. A production
manager must verify from the scan-package inventory that both strand files were
actually retained at that floor. A density-capped production threshold must
not be presented as a score `-1` source.

## CUT&RUN estimand

For each series and isoform, the binary outcome is the existing strict-immersion
comparison between anti-p73 and its negative control. Only discordant anchors
contribute to the matched binary contrast. The common odds ratio is
Mantel-Haenszel adjusted over:

- cell-line series;
- 5 Mb genomic block; and
- TP73 score stratum.

The exact-position feature compares `positive` with `strict_negative`; the
intermediate class is omitted. The two supported series are Saos-2 and
`skmel29_2`; `skmel29_1` remains excluded. TAp73 and DNp73 are estimated
separately, followed by a direct TA-versus-DN log-odds contrast. A leave-one-5
Mb-block-out jackknife supplies uncertainty. Series-specific estimates and a
direction-consistency flag remain visible.

Effect contrasts use `log(OR_TA) - log(OR_DN)`. Positional shifts use
`DN_position - TA_position`; negative and positive signs therefore describe
movement in the TP73-oriented coordinate frame, while
`dn_peak_is_closer_to_anchor` states the absolute-distance interpretation
directly.

This is an association with the technical CUT&RUN occupancy label, not proof
that a protein binds directly or physically interacts with TP73.

## Periodicity and peak matching

The prespecified periodic component has period 10.5 bp. For each isoform,
distance frame, and orientation, a weighted quadratic trend is compared with
the same trend plus sine and cosine terms. The output records harmonic
amplitude, phase, and weighted residual reduction.

Preferred positions are found on a lightly Gaussian-smoothed log-odds curve.
Peak calls are descriptive candidates. TA and DN peaks are matched only after
choosing a coherent group-level shift hypothesis. This prevents a doublet at
TA positions `-42,-31` and DN positions `-31,-20` from being reduced to one
nearest-neighbour match at `-31`: the two preserved 11 bp spacings and the
11 bp movement of their midpoint are retained.

The 10.5 bp term is a useful DNA-helical-rotation hypothesis and positive
control. It is not a universal gate for biological relevance. Peak-location
uncertainty and a chromosome/block-preserving spatial shift null are still
required before treating a peak or phase difference as inferential evidence.

## Command-line stages

Build one chromosome/motif statistic package:

```bash
scripts/build_tp73_exact_distance_counts.py \
  --anchors ANCHOR_EVIDENCE_CHROM.parquet \
  --anchor-plus TP73_PLUS_CHROM.parquet \
  --anchor-minus TP73_MINUS_CHROM.parquet \
  --cofactor-plus MOTIF_PLUS_CHROM.parquet \
  --cofactor-minus MOTIF_MINUS_CHROM.parquet \
  --motif-id MA0024.3 --motif-name E2F1 --chrom 1 \
  --source-score-floor -1 --positive-threshold 0 \
  --inventory-output E2F1.chr1.inventory.parquet \
  --block-output E2F1.chr1.blocks.parquet
```

Repeat the two input switches for every chromosome of the same motif, then
finalize them by repeating `--inventory` and `--block-components`:

```bash
scripts/summarize_tp73_exact_distance_response.py \
  --inventory E2F1.chr1.inventory.parquet \
  --inventory E2F1.chr2.inventory.parquet \
  --block-components E2F1.chr1.blocks.parquet \
  --block-components E2F1.chr2.blocks.parquet \
  --output-dir E2F1_exact_distance
```

Both commands refuse to replace existing outputs. The chromosome output is
restart-friendly, and the finalizer records input checksums in `manifest.json`.

## Output contract

The finalized directory contains Parquet tables and one portable DuckDB file:

- `distance_inventory`: zero-complete occurrence/support counts;
- `run_config`: score, geometry, block, motif, and chromosome provenance;
- `distance_response`: isoform-specific CUT&RUN odds-ratio curves;
- `isoform_contrast`: direct TA-versus-DN comparisons;
- `periodicity` and `periodicity_isoform_contrast`;
- `peak`, `peak_match`, and `peak_pair_match`;
- `report.md`; and
- `manifest.json`.

Every table has a stable schema even when it has no rows. The DuckDB tables use
the same names prefixed with `tp73_exact_distance_`.

## Validation and expansion

The synthetic integration test fixes the following invariants:

- exact 10.5 and 11 bp offsets;
- reverse-strand sign conversion;
- one count per physical cofactor span across strand files;
- unresolved TP73 orientation excluded only from the oriented frame;
- zero-complete source/intermediate/positive counts;
- TP73 score strata retained in sufficient statistics;
- stable schemas for empty derived outputs; and
- coherent matching of an 11 bp doublet translated by 11 bp.

The next production stage should run a prespecified candidate panel first,
including all POU motifs, PATZ1, SP1, E2F1, REST, TFAP2C, TCF7, TP53, and TP63.
An all-JASPAR exact-distance sweep should follow only after the spatial-null and
peak-location validation are fixed. Arbitrary cofactor-cofactor grammar needs a
separate A-B pair kernel and an unconditioned transaction universe; the present
TP73-conditioned result must not be relabelled as that network.

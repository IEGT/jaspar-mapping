# TP73 CUT&RUN association across named POU motifs

## Scope

This sensitivity analysis includes every JASPAR 2026 vertebrate matrix whose
official matrix name begins with `POU`, case-insensitively. The locked list is
[`config/tp73_pou_named_motifs_v1.tsv`](../config/tp73_pou_named_motifs_v1.tsv).
It contains 18 matrices: 15 uniquely named single-factor motifs, two different
POU6F1 matrices, and two POU-SOX2 heteromeric matrices. Seventeen are derived
from human source material; `Pou5f1::Sox2` (`MA0142.1`) is mouse-derived.

HNF1A, HNF1B, HMBOX1, and mouse Hnf1A are classified by JASPAR as POU-domain
factors but do not match the declared `POU*` name rule. They remain related
controls and are not silently mixed into this family panel.

The source result is the completed autosome 1-22, Saos-2 plus `skmel29_2`
analysis with `skmel29_1` excluded. The local complete extracts are kept under
the ignored directory
`dry_runs/tp73_distance_species_isoform_highlights_autosomes_20260823/` as
`pou_family_calibrated_by_band.tsv` and
`pou_family_calibrated_cross_band.tsv`.

## Calibrated-threshold result

Every one of the 18 POU-named matrices is depleted among TP73 CUT&RUN-supported
anchors for both TAp73 and DNp73 in every one of the six exclusive distance
bands. Every association has the same direction in Saos-2 and `skmel29_2`.
For every matrix and every distance band, depletion is stronger for DNp73 than
for TAp73 (`OR_TA / OR_DN > 1`).

| Matrix | Calibrated threshold | TA OR range | DN OR range | TA/DN OR range |
|---|---:|---:|---:|---:|
| POU5F1 `MA1115.2` | 6 | 0.388-0.763 | 0.307-0.583 | 1.260-1.406 |
| POU5F1B `MA0792.1` | 6 | 0.444-0.719 | 0.320-0.496 | 1.386-1.594 |
| POU3F1 `MA0786.2` | 6 | 0.493-0.716 | 0.339-0.464 | 1.452-1.633 |
| POU6F1 `MA0628.2` | 3 | 0.534-0.831 | 0.374-0.723 | 1.150-1.490 |
| POU3F3 `MA0788.1` | 5 | 0.547-0.723 | 0.383-0.497 | 1.428-1.597 |
| POU2F1 `MA0785.2` | 5 | 0.550-0.748 | 0.393-0.537 | 1.394-1.505 |
| POU4F1 `MA0790.2` | 5 | 0.556-0.747 | 0.388-0.567 | 1.316-1.634 |
| POU4F3 `MA0791.2` | 5 | 0.571-0.744 | 0.390-0.533 | 1.397-1.634 |
| Pou5f1::Sox2 `MA0142.1` | 5 | 0.567-0.771 | 0.413-0.587 | 1.312-1.378 |
| POU6F1 `MA1549.2` | 5 | 0.594-0.860 | 0.446-0.705 | 1.220-1.346 |
| POU3F2 `MA0787.1` | 4 | 0.612-0.758 | 0.430-0.546 | 1.381-1.503 |
| POU2F3 `MA0627.3` | 6 | 0.616-0.767 | 0.437-0.557 | 1.354-1.467 |
| POU3F4 `MA0789.1` | 4 | 0.668-0.772 | 0.474-0.549 | 1.342-1.437 |
| POU6F2 `MA0793.2` | 4 | 0.650-0.900 | 0.501-0.844 | 1.067-1.335 |
| POU2F2 `MA0507.3` | 0 | 0.720-0.831 | 0.501-0.681 | 1.221-1.438 |
| POU4F2 `MA0683.2` | 0 | 0.747-0.833 | 0.482-0.710 | 1.173-1.551 |
| POU1F1 `MA0784.3` | 0 | 0.725-0.826 | 0.508-0.664 | 1.244-1.426 |
| POU2F1::SOX2 `MA1962.1` | 0 | 0.749-0.815 | 0.550-0.670 | 1.217-1.360 |

The coherence is biologically interesting but does not make these 18
independent confirmations. POU-family matrices can have similar sequence
preferences, and calibrated thresholds range from 0 to 6. Motif-similarity
clustering and a shared-threshold analysis are therefore required before
assigning the pattern to individual POU proteins.

## Historical identity

The publication-tag source explicitly selected POU2F2 `MA0507.2`. Solitary
POU2F1 was not mislabeled as POU2F2; the old code separately mentioned the
POU2F1::SOX2 heteromeric matrix. POU2F2 therefore remains in the historical
prespecified panel, while schema 4 adds all 18 named POU matrices under the
separate `pou_family_panel` criterion.

## Common score-zero sensitivity

The calibrated analysis is primary. A second run fixes the positive threshold
to score 0 for all 18 matrices while retaining the same strict-negative rule
of no score-`>= -1` locus. This tests whether the family-wide depletion and the
stronger DNp73 effect survive a common score interpretation.

The source-controlled command first derives a small threshold package from the
production registry without changing or copying the scan data:

```bash
scripts/build_fixed_threshold_subset.py \
  --source-registry PRODUCTION_THRESHOLDS.parquet \
  --source-threshold-set-id SOURCE_THRESHOLD_SET \
  --motif-list config/tp73_pou_named_motifs_v1.tsv \
  --fixed-threshold 0 \
  --threshold-set-id jaspar2026_grch38_tp73_pou_common0_v1 \
  --output-dir /data/sm718/jaspar_mapping_runs/thresholds/pou_common0_v1
```

The distance-enrichment manager then sees only those 18 registry rows and runs
the compact count/block aggregation against the existing score-floor `-1`
Parquet. No genome or motif rescan is required. The common-zero result must be
reported beside, not in place of, the calibrated result because very frequent
motifs may lose a useful strict-negative comparison class at score 0.

The production sensitivity run completed across autosomes 1-22. All 108
motif-by-distance estimates remain depleted for both isoforms, all 108 retain
the same direction in Saos-2 and `skmel29_2`, and DNp73 depletion remains
stronger than TAp73 depletion in every cell. Thus the family-wide direction is
not an artifact of using motif-specific thresholds between 3 and 6.

| Matrix | Common-0 TA OR range | Common-0 DN OR range | Common-0 TA/DN OR range |
|---|---:|---:|---:|
| POU5F1 `MA1115.2` | 0.440-0.800 | 0.365-0.654 | 1.179-1.268 |
| POU5F1B `MA0792.1` | 0.505-0.759 | 0.408-0.588 | 1.224-1.322 |
| POU3F1 `MA0786.2` | 0.544-0.760 | 0.425-0.573 | 1.266-1.352 |
| POU6F1 `MA0628.2` | 0.567-0.811 | 0.443-0.667 | 1.216-1.390 |
| POU3F3 `MA0788.1` | 0.584-0.764 | 0.451-0.584 | 1.297-1.368 |
| POU2F1 `MA0785.2` | 0.598-0.777 | 0.472-0.612 | 1.250-1.309 |
| POU4F1 `MA0790.2` | 0.589-0.779 | 0.451-0.625 | 1.246-1.408 |
| POU4F3 `MA0791.2` | 0.602-0.781 | 0.462-0.626 | 1.249-1.404 |
| Pou5f1::Sox2 `MA0142.1` | 0.607-0.798 | 0.485-0.643 | 1.215-1.265 |
| POU6F1 `MA1549.2` | 0.641-0.890 | 0.531-0.810 | 1.098-1.221 |
| POU3F2 `MA0787.1` | 0.641-0.794 | 0.502-0.635 | 1.251-1.350 |
| POU2F3 `MA0627.3` | 0.651-0.789 | 0.517-0.618 | 1.244-1.306 |
| POU3F4 `MA0789.1` | 0.700-0.790 | 0.544-0.631 | 1.217-1.312 |
| POU6F2 `MA0793.2` | 0.687-0.880 | 0.557-0.786 | 1.119-1.250 |
| POU2F2 `MA0507.3` | 0.720-0.831 | 0.501-0.681 | 1.221-1.438 |
| POU4F2 `MA0683.2` | 0.747-0.833 | 0.482-0.710 | 1.173-1.551 |
| POU1F1 `MA0784.3` | 0.725-0.826 | 0.508-0.664 | 1.244-1.426 |
| POU2F1::SOX2 `MA1962.1` | 0.749-0.815 | 0.550-0.670 | 1.217-1.360 |

POU2F1 is more depleted than POU2F2 in most non-overlapping bands, especially
at gaps of 51-150 bp, whereas POU2F2 is more depleted in the overlap band.
This makes a recollection of POU2F1 as the stronger distal inhibitor
understandable, but it does not change the source audit: the publication code
selected POU2F2.

The score-zero package is
`/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_tp73_distance_pou_common0_v1`
and its schema-4 final was produced by source commit `e20ff49`. The calibrated
schema-4 comparison is under
`jaspar2026_grch38_tp73_distance_species_vertebrates_v1/final/`
`distance_enrichment_isoform_highlights_v4`. Both manifests record the full
finalizer commit and a clean source tree.

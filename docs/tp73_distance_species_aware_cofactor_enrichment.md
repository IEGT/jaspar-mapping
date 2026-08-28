# Species-aware TP73 cofactor enrichment by distance

## Purpose

This analysis ranks JASPAR motif signatures associated with successful TP73
CUT&RUN support separately for TAp73 and DNp73 and for six mutually exclusive
distance bands around each TP73 anchor. It replaces a distance-agnostic pooled
ranking with a compact calculation that does not materialize an
anchor-by-motif matrix.

The primary human-genome analysis is limited to JASPAR matrices whose official
`tax_group` is `vertebrates`. The complete JASPAR 2026 metadata source and its
checksum are described in [jaspar_2026_provenance.md](jaspar_2026_provenance.md).

## Taxonomic interpretation

Two related but different JASPAR attributes are retained:

- `tax_group` is the broad JASPAR classification used to define the primary
  motif family. There are 1,019 vertebrate matrices in the 2,633-matrix CORE
  non-redundant PFM, including TP73, and therefore 1,018 candidate cofactors.
- `jaspar_matrix_species` is a normalized one-to-many table of NCBI taxonomy
  IDs and scientific names. These are the organisms contributing to matrix
  construction. They are provenance, not proof that the sequence preference
  can bind only in those species.

`includes_homo_sapiens` provides a convenient exact-source sensitivity subset.
The final package contains top/bottom tables for both all selected vertebrate
matrices and this human-source subset. Mouse-derived vertebrate matrices are not
discarded from the primary analysis merely because the human TF has not been
used to construct that particular JASPAR matrix.

Matrices from plants, insects, fungi, nematodes, urochordates, and diatoms do
not enter the primary ranking. A separate run with `--tax-group all` can retain
them as an explicitly all-taxon sequence-signature sensitivity analysis.

## Cofactor classes

For each TP73 anchor and motif, plus- and minus-strand records with the same
physical interval are collapsed. The maximum score is selected independently
within each exclusive interval-distance band:

| Band | Signed interval-edge distance |
|---|---:|
| `overlap` | `< 0` |
| `adjacent_0_5` | `0` through `5` bp |
| `gap_6_20` | `6` through `20` bp |
| `gap_21_50` | `21` through `50` bp |
| `gap_51_100` | `51` through `100` bp |
| `gap_101_150` | `101` through `150` bp |

Distance zero means that one half-open motif interval ends exactly where the
other starts. Negative values are overlaps. Each anchor-band is classified as:

- **positive:** its best score is at least the motif-specific operating
  threshold;
- **intermediate:** a retained locus is present at score `>= -1` but below the
  operating threshold;
- **strict negative:** no retained score-`>= -1` locus occurs in that band.

Intermediate anchors are intentionally excluded from the positive-versus-
negative odds ratio. This preserves the low scan floor and prevents the
operating threshold from redefining absence.

## CUT&RUN estimand

The two valid series are `saos2` and `skmel29_2`; `skmel29_1` remains excluded.
TA and DN are estimated separately. Within each series and isoform, only TP73
anchors discordant between anti-p73 support and its matched negative-control
support contribute to the odds ratio. The response therefore compares whether
positive versus strictly negative cofactor contexts are preferentially found
on the anti-p73-supported side of the paired technical contrast.

The isoform comparison is a separate estimand, not a comparison of two
significance labels. For every motif and distance band it reports
`OR_TA / OR_DN`, equivalently `log(OR_TA) - log(OR_DN)`. Values above one mean
that the cofactor association odds ratio is larger for TAp73; values below one
mean that it is larger for DNp73. This interpretation still requires the two
individual odds ratios: for example, a larger value can mean stronger
enrichment or weaker depletion. Its uncertainty is calculated by deleting the
same 5 Mb genomic block from both isoform estimates before taking their
difference, thereby retaining the covariance induced by their shared anchors.

The common odds ratio is a Mantel-Haenszel estimate stratified by series and by
TP73 anchor-score bands `[-5,-1)`, `[-1,0)`, `[0,1)`, `[1,2)`, `[2,5)`,
`[5,10)`, `[10,15)`, and `[15,+Inf)`. Uncertainty is a delete-one-cluster
jackknife over 5 Mb genomic blocks. BH adjustment is performed separately for
each isoform and distance band, using every planned matrix in the declared
taxonomic family as the family size.

These are motif-context associations, not yet unique TF identities. Similar
matrices and TF-family redundancy require a later consolidation analysis.

## Motif frequency and historical compatibility

Schema 5 reports motif frequency beside every enrichment/depletion estimate.
Repeated motif loci within one TP73 anchor and exclusive distance band count as
one presence, matching the binary context interpretation used by the earlier
JASPAR-2022 analysis. The following quantities remain distinct:

- `all_tp73_anchor_vicinity_frequency` is the fraction of all TP73 anchors
  with a motif score at or above the declared positive threshold in the band;
- `anti_supported_positive_anchor_fraction_discordant` is the corresponding
  fraction on the anti-p73-supported side of the matched-discordant CUT&RUN
  comparison;
- `control_supported_positive_anchor_fraction_discordant` is the matched
  control-side fraction; and
- `anti_to_control_positive_anchor_log2_ratio_discordant` is their descriptive
  log2 ratio.

The historical scatter plot used all CUT&RUN-supported anchors. The modern
supported frequencies are intentionally restricted to anchors discordant
between anti-p73 and matched control, because that is the population entering
the primary matched analysis. Their frequency denominator includes positive,
intermediate, and strict-negative contexts, with intermediate contexts counted
as not positive. The primary odds ratio instead compares positive with strict
negative and excludes intermediate contexts. The historical-style log2
frequency ratio is therefore a compatibility view, not a replacement for the
stratified adjusted odds ratio or its block-jackknife uncertainty.

The finalizer-only schema extension reuses existing task sufficient statistics;
it does not rescan motifs or recompute chromosome-level geometry. Generate both
frequency plots from the finalized TSV with:

```bash
scripts/plot_tp73_distance_frequency_enrichment.R \
  --input FINAL/cofactor_distance_frequency_enrichment.tsv \
  --output-adjusted frequency_vs_adjusted_log2_odds.png \
  --output-frequency-ratio frequency_vs_descriptive_log2_ratio.png \
  --output-table frequency_plot_data.tsv
```

The first figure is primary: adjusted log2 odds on the x axis and
anti-p73-supported motif frequency on the y axis. The second reproduces the
older frequency-versus-log-ratio presentation as closely as the matched design
allows, with its descriptive scope stated in the subtitle.

## Restart-safe production

The kernel is
[`build_tp73_distance_cofactor_counts.py`](../scripts/build_tp73_distance_cofactor_counts.py).
It emits only chromosome-level class counts and block sufficient statistics.
The manager
[`manage_tp73_distance_cofactor_enrichment.py`](../scripts/manage_tp73_distance_cofactor_enrichment.py)
pins exact scan files, thresholds, source hashes, chromosome scope, and JASPAR
metadata.

Each Slurm task handles one motif. It stages one chromosome's anchor and two
strand files at a time on `/scratch`, then atomically publishes compact
chromosome checkpoints below the dedicated `/data/sm718` run root. Requeue and
manual repetition reuse validated checkpoints. Completed motif outputs and the
final package are also atomic and checksum-validated.

Schema-2 finalization can reuse schema-1 task sufficient statistics without
rerunning motif jobs. Its manifest distinguishes the task source commit from
the finalizer source commit and records exact scientific-source hashes. An
existing schema-1 final directory is never overwritten; use a new
`--final-name` for the schema-2 derivative. If the finalizer is newer than the
task source, the submission command supplies `--finalizer-source-commit`; the
compute node verifies file hashes and does not need a Git executable.

Schema 3 added presentation tables only. Schema 4 additionally appends every
vertebrate JASPAR matrix whose name begins with `POU` (case-insensitive) under
the separate `pou_family_panel` criterion. Both schemas reuse the same block
components and do not alter any scientific estimate.

## Colleague-facing highlight criteria

The complete side-by-side comparison retains `OR_TA`, `OR_DN`, and
`OR_TA / OR_DN`, their confidence intervals and q-values, and ranks for TA
enrichment/depletion, DN enrichment/depletion, and absolute isoform difference.
The compact highlight table then selects up to 20 rows per distance band and
direction for four explicitly labelled criteria:

1. `TA_association`: strongest significant enrichment and depletion for TAp73;
2. `DN_association`: the same selection for DNp73;
3. `isoform_difference`: strongest significant direct differences, separated
   into larger TA and larger DN odds ratios;
4. `opposite_direction`: the strict subset in which both individual 95%
   intervals lie on opposite sides of OR 1, both individual BH tests and the
   direct isoform test pass at 0.05.

POU2F2, SP1, PATZ1, REST, and E2F1 are appended under
`prespecified_panel` in every distance band. This does not alter any
data-driven rank. Their rows retain their actual ranks under all three complete
ranking systems, making absence from a top 20 distinguishable from absence
from the analysis.

Schema 4 also appends the complete named POU family under `pou_family_panel`.
This is deliberately distinct from the historical panel and includes both
POU6F1 matrices and POU-SOX2 heteromeric matrices. It does not silently include
HNF1A, HNF1B, or HMBOX1 merely because JASPAR classifies them as POU-domain
factors; those are related controls rather than names matching the declared
`POU*` rule.

The submission entry point is:

```bash
scripts/submit_tp73_distance_cofactor_enrichment_slurm.sh \
  --run-root /data/sm718/jaspar_mapping_runs/RUN_ID \
  --scan-package /data/sm718/jaspar_mapping_runs/SCAN/package \
  --anchor-evidence /data/sm718/jaspar_mapping_runs/EVIDENCE.parquet \
  --thresholds /data/sm718/jaspar_mapping_runs/thresholds.parquet \
  --threshold-set-id THRESHOLD_SET_ID \
  --jaspar-catalog /data/sm718/jaspar_mapping_runs/resources/jaspar/2026/CORE/catalog_core_nonredundant \
  --duckdb /home/sm718/micromamba/bin/duckdb \
  --source /data/sm718/jaspar_mapping_runs/source_checkouts/jaspar-mapping-COMMIT \
  --run-id RUN_ID
```

The exact paths must be resolved from their manifests before submission; the
placeholders above are not production defaults.

## Query surface

The final DuckDB database contains:

- `cofactor_distance_enrichment`: one row per motif, isoform, and distance
  band, including class support, series-specific odds ratios, block-jackknife
  uncertainty, taxonomic group, and flattened source-species label;
- `cofactor_distance_frequency_enrichment`: the compact plotting/query surface
  with all-anchor, anti-p73-supported, and control-supported motif frequencies,
  the descriptive anti/control log2 frequency ratio, and the primary adjusted
  log2 odds ratio;
- `cofactor_distance_isoform_contrast`: one row per motif and distance band,
  with the TA and DN estimates side by side, their odds-ratio ratio, a paired
  block-jackknife confidence interval, and BH-adjusted direct isoform test;
- `cofactor_distance_isoform_comparison` and its `_human_source` subset: the
  complete side-by-side table with TA, DN, and absolute isoform-difference
  ranks;
- `cofactor_distance_highlight_20` and its `_human_source` subset: the four
  labelled highlight selections plus the prespecified five-motif panel;
- `cofactor_distance_top_bottom_20`: vertebrate top and bottom 20 per isoform
  and distance band, carrying the other isoform and direct contrast beside each
  ranked result;
- `cofactor_distance_top_bottom_20_human_source`: the corresponding exact-human
  source sensitivity ranking;
- `cofactor_distance_class_count`: positive, intermediate, and strict-negative
  anchor counts;
- `cofactor_distance_block_component`: the sufficient statistics from which
  estimates and block uncertainty are reconstructed;
- `jaspar_matrix`, `jaspar_matrix_species`, and `jaspar_motif_set_matrix`: the
  normalized official metadata and exact PFM membership.

For example:

```sql
SELECT direction, rank, motif_id, motif_name, adjusted_odds_ratio,
       ta_adjusted_odds_ratio, dn_adjusted_odds_ratio,
       ta_vs_dn_odds_ratio_ratio, ta_vs_dn_q_value_bh_tax_group,
       source_species
FROM cofactor_distance_top_bottom_20
WHERE isoform = 'TA' AND distance_band = 'gap_6_20'
ORDER BY CASE direction WHEN 'enriched' THEN 1 ELSE 2 END, rank;
```

Sort `cofactor_distance_isoform_contrast` directly when the question is which
motifs differ most between isoforms. Do not infer an isoform difference merely
because one isoform's confidence interval excludes one and the other's does
not. A distance-agnostic TA/DN-pooled ranking may be retained as a descriptive
supplement, but it is not the primary cofactor result.

The presentation table can be queried without reconstructing any ranking:

```sql
SELECT selection_criterion, selection_direction, selection_rank,
       motif_name, motif_id, ta_adjusted_odds_ratio, dn_adjusted_odds_ratio,
       ta_vs_dn_odds_ratio_ratio, ta_q_value_bh_tax_group,
       dn_q_value_bh_tax_group, q_value_bh_tax_group
FROM cofactor_distance_highlight_20_human_source
WHERE distance_band = 'gap_21_50'
ORDER BY selection_criterion, selection_direction,
         selection_rank NULLS LAST, motif_name;
```

Use `cofactor_distance_top_bottom_20_human_source` for the exact-human source
sensitivity view, or join `jaspar_matrix_species` by `matrix_id = motif_id` when
individual taxonomy IDs rather than a display label are needed.

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

The common odds ratio is a Mantel-Haenszel estimate stratified by series and by
TP73 anchor-score bands `[-5,-1)`, `[-1,0)`, `[0,1)`, `[1,2)`, `[2,5)`,
`[5,10)`, `[10,15)`, and `[15,+Inf)`. Uncertainty is a delete-one-cluster
jackknife over 5 Mb genomic blocks. BH adjustment is performed separately for
each isoform and distance band, using every planned matrix in the declared
taxonomic family as the family size.

These are motif-context associations, not yet unique TF identities. Similar
matrices and TF-family redundancy require a later consolidation analysis.

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
- `cofactor_distance_top_bottom_20`: vertebrate top and bottom 20 per isoform
  and distance band;
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
       confidence_interval_95_lower, confidence_interval_95_upper,
       source_species
FROM cofactor_distance_top_bottom_20
WHERE isoform = 'TA' AND distance_band = 'gap_6_20'
ORDER BY CASE direction WHEN 'enriched' THEN 1 ELSE 2 END, rank;
```

Use `cofactor_distance_top_bottom_20_human_source` for the exact-human source
sensitivity view, or join `jaspar_matrix_species` by `matrix_id = motif_id` when
individual taxonomy IDs rather than a display label are needed.

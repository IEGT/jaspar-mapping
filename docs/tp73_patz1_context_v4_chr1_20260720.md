# TP73/PATZ1 chromosome-1 context-v4 pilot

This is an operational benchmark for the first production-shaped schema-v4
context build. It is not a biological result or a replacement for the
CUT&RUN-calibrated analyses.

## Inputs and policy

- Genome identity: `homo_sapiens_grch38_ensembl113_primary`
- Motif set: `jaspar2026_core_nonredundant`
- Chromosome: `1`
- Motifs: TP73 `MA0861.2` and PATZ1 `MA1961.2`, both orientations
- Score mode: `log2_relative_risk`, pseudocount `1`
- TP73 anchor floor: `-1`
- Anchor selection: physical-locus local peak within a 150 bp start radius
- Tandem eligibility: both TP73 orientation records score at least `0`, with a
  non-overlapping edge gap no larger than 20 bp
- Context and cofactor-pair interval distance: 150 bp
- Output tier: `selected`
- Annotation: Ensembl release 113 GTF, converted directly to BED coordinates

The four source files were exposed through a hard-link tree, so no motif data
were copied. The local TP73 source was retained only to `-1`; this is sufficient
for the selected anchors and the score-`0` tandem rule. The production Haumea
scan retains TP73 to `-5`, and will additionally preserve those deeper rejected
TP73 alignments as inspectable context neighbors.

## Execution

The build used two DuckDB threads, a 6 GB memory ceiling, and a 30 GB spill
ceiling. It completed in 470.36 seconds (7 minutes 50 seconds). Sampled DuckDB
resident memory remained below 1.9 GB. The staging tree reached approximately
10 GB while joins were spilling, then collapsed to a 1.3 GB finished package.

The package is locally available at:

```text
dry_runs/jaspar2026_chr1_requested_panel_20260720/
  context_v4_local_patz1_selected/
```

## Cardinalities

| Relation | Rows |
|---|---:|
| `tp73_anchor_locus` | 305,492 |
| `tp73_pair_feature` | 305,528 |
| `tp73_context_anchor` | 305,528 |
| `motif_context_pair` | 1,278,832 |
| `tp73_motif_context_summary` | 782,352 |
| `cofactor_motif_locus` | 1,906,567 |
| `cofactor_motif_pair` | 7,684,131 |
| `cofactor_locus_pair_feature` | 1,906,567 |
| `tp73_cofactor_pair_context` | 4,527,996 |
| `tp73_cofactor_pair_summary` | 1,563,537 |
| `motif_transcript_context` | 1,612,458 |
| `transcript` | 35,122 |
| `intron` | 167,554 |

There are 305,492 physical TP73 local-peak spans but 305,528 orientation
anchors because equal-scoring strand alternatives at a retained physical locus
remain explicit. Their `anchor_hit_id` and `(chrom,start,strand)` keys are both
unique.

PATZ1 occurs within the 150 bp interval context of 224,999 TP73 orientation
anchors. TP73 pair classes are:

| Pair class | Orientation anchors |
|---|---:|
| `singleton` | 297,935 |
| `tandem_same_orientation` | 3,299 |
| `tandem_opposite_orientation` | 2,806 |
| `tandem_orientation_ambiguous` | 1,412 |
| `tandem_mixed_orientation` | 76 |

Primary transcript-region summaries are 218,380 intronic, 77,118 intergenic,
9,145 transcribed non-intronic, and 885 crossing an intron boundary. These are
anchor-level conveniences; the per-transcript table retains isoform ambiguity.

## Acceptance checks

The completed package has zero violations for:

- local-peak status, score floor, and nonnegative regional prominence;
- unique anchor IDs and orientation keys;
- exhaustive TP73 partner-locus orientation counts and singleton nullability;
- symmetric tandem score floors, distinct spans, non-overlap, and the 20 bp
  edge-gap limit;
- the 150 bp signed-interval capture boundary;
- one annotation summary row per anchor; and
- textual chromosome types across every Hive-partitioned query view.

The canonical `motif_context_run_config` view and backward-compatible
`context_run_config` alias both expose schema version 4. Native, R, DuckDB,
context, exact-input-staging, script-help, scanner, and SIGUSR1 tests pass.

## Cluster consequence

The PATZ1 expansion shows why selected context should run as independent
motif/chromosome tasks. Dense motifs such as TCF7 and POU4F1 can create much
larger same-motif pair layers. The Haumea pilot therefore uses one cofactor per
array task, exact inventory files, durable output under `/data/sm718`, and
DuckDB spill under job-local `/scratch`. A failure or resource underestimate is
isolated to one motif rather than invalidating the panel.

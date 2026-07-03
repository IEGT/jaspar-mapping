# Promoter Architecture → TA-vs-DN Differential: Schema & Query Draft

Status: **draft**. Companion files: [`sql/schema.sql`](../sql/schema.sql) (DuckDB
DDL), [`sql/queries.sql`](../sql/queries.sql) (canned queries). This extends the
four-layer model in [`docs/tf_mapping_schema_notes.md`](tf_mapping_schema_notes.md).

## Goal

Make promoter architecture and gene expression jointly queryable with low
latency so that (a) an ML model can predict **isoform-specific differential
expression (TAp73 vs DNp73)** from sequence-derived promoter architecture, and
(b) OpenClaw / Claude / Codex can reach any slice of it through one stable
interface without bespoke plumbing.

## Design principles

The model is the four layers from the schema notes, stored **long/tidy** and
partitioned by chromosome, never as the ~8,000-column wide matrix. Layers 1–2
are deterministic from sequence and regenerable on demand; layer 3 (CUT&RUN) is
irreducible experimental data; layer 4 is annotation. Apache Parquet is the
canonical at-rest format (ZSTD, sorted within partition by `start`). A single
**DuckDB** file is the query surface: zero-copy `VIEW`s over Parquet for the raw
layers, and a small number of **materialized** `TABLE`s for the one interval
join and the feature aggregation we refuse to recompute per query. DuckDB is
embedded, single-file, and speaks SQL — which is exactly what makes it a
reasonable target for an agent (no server, no session state, predictable
latency, introspectable schema).

Coordinates are BED 0-based half-open everywhere; run the 2026 scan with
`--coordinate-mode bed` so layer 1 lands in that convention and no join has to
reconcile two coordinate systems.

## Physical layout

```
tables/jaspar2026/
  manifest.json            # source commit, JASPAR 2026 sha256, Ensembl release, row counts
  checksums.sha256
  motif_metadata/…          # tiny motif dimension: motif_id, name, length, source
  motif_hit/chrom=*/…        # layer 1  (large, regenerable)
  motif_score_dense/
    motif_id=*/score_mode=*/pseudocount=*/chrom=*/strand=*/…
                            # dense calibration blocks, score for every PSSM alignment start
  motif_cutandrun_score_bin_stats/
    motif_id=*/score_mode=*/pseudocount=*/chrom=*/…
                            # CUT&RUN enrichment by score bin
  motif_architecture/…       # layer 2a (tiny, motif-keyed dimension)
  hit_architecture/chrom=*/… # layer 2b (optional, per-hit, keyed to motif_hit)
  cutandrun/chrom=*/…        # layer 3  (irreducible experimental)
  gene/…  promoter/…  gene_set/…            # layer 4 dimensions
  expression/…  expression_differential/…   # layer 4 label source
  jaspar2026.duckdb        # query index over Parquet plus materialized hot tables
```

Only non-zero facts are stored for the genome-wide hit layer — `motif_hit`
holds actual hits, `cutandrun` one row per (locus, sample),
`promoter_arch_feature` one row per (gene, motif). Dense calibration runs are
the explicit exception: for selected motifs/chromosomes they store one score
per motif-model alignment start and strand in block vectors, with repeated
identity moved into partition paths.

## Layers and keys

Full column lists live in `sql/schema.sql`; the essentials:

- **`motif_hit`** — `chrom, start, end, motif_id, motif_name, strand, score, score_mode, pseudocount, pwm_relative_score, matched_seq`. Straight from `pssm_scan`. `score_mode` plus `pseudocount` define the motif scoring transform; `pwm_relative_score` (0–1) is the cross-motif-comparable score to threshold on.
- **`motif_score_dense`** — logical view over dense score blocks: `chrom, start, end, motif_id, motif_name, strand, score_mode, pseudocount, score`. Physical files store only `block_start` plus a `FLOAT[] scores` vector; `motif_id`, `score_mode`, `pseudocount`, `chrom`, and `strand` live in hive partitions.
- **`motif_cutandrun_score_bin_stats`** — dense-score calibration output by score bin and CUT&RUN sample: window counts, covered-window counts, baseline fraction, enrichment ratio, log2 enrichment, and signal summaries.
- **`motif_architecture`** — one row per `motif_id`: family, `binding_unit_model`, half-site pattern, spacer range, `architecture_confidence`. Curated for the TP53 family, `unknown` elsewhere.
- **`hit_architecture`** *(optional)* — per-hit decomposition (`spacer_bp`, half-site scores, `oligomer_compatible`), derived from `matched_seq`, kept separate so layer 1 stays stable.
- **`cutandrun`** — `chrom, start, end, cell_line, isoform, antibody, replicate, signal`. The sample facets are columns, so a TA-vs-DN contrast is a `GROUP BY`, not a filename parse.
- **`gene / promoter / gene_set`** — annotation dimensions; `promoter` carries the strand and TSS used for distance features.
- **`expression_differential`** — `gene_id, cell_line, log2fc_ta_vs_dn, padj, …`. **This is the ML label.**

Three materialized tables carry all the cost: `promoter_motif_hit` (the single
interval join of hits into promoter windows, with strand-aware `tss_distance`),
`promoter_arch_feature` (per-gene, per-motif architecture summary), and
`promoter_card` (one compact row per gene for agent lookup). The `ml_ta_vs_dn`
view joins the long feature table to the label.

## Dense Chr1 Calibration

Before choosing a genome-wide storage threshold, run a dense chromosome-1
calibration for TP73 (`MA0861.2`) on both strands:

- `--score-mode log2_relative_risk --pseudocount 1`
- `--score-mode log_odds --pseudocount 1`

Each run stores one score per motif-model alignment start and strand as dense
block vectors. Position is implicit (`start = block_start + offset`), and
`end` is derived from `motif_metadata.motif_length`. This interval is the DNA
span used for scoring the PSSM alignment; it should not be read as an asserted
physical footprint of the TF complex. The theoretical anchor remains score
`0`: for log2 relative risk it is the point where the scored sequence has no
net information gain over equal A/C/G/T background. The calibration then asks
whether CUT&RUN support diverges from the random-alignment baseline above, at,
or below that anchor, and whether log odds gives a cleaner or more symmetric
separation.

Use the same score bins for both dense runs:

| Bin | Interval |
| --- | --- |
| 1 | `(-Inf,-10000)` |
| 2 | `[-10000,-200)` |
| 3 | `[-200,-50)` |
| 4 | `[-50,-10)` |
| 5 | `[-10,-5)` |
| 6 | `[-5,-2)` |
| 7 | `[-2,-1)` |
| 8 | `[-1,-0.5)` |
| 9 | `[-0.5,0)` |
| 10 | `[0,0.5)` |
| 11 | `[0.5,1)` |
| 12 | `[1,2)` |
| 13 | `[2,5)` |
| 14 | `[5,10)` |
| 15 | `[10,+Inf)` |

The comparison statistic is per score bin, strand, and CUT&RUN sample:
`n_windows`, `n_covered_windows`, `overlap_fraction`, `baseline_fraction`,
`enrichment_ratio`, `log2_enrichment`, and optional mean/max signal. The
genome-wide `motif_hit` threshold should be set just below the point where the
coverage enrichment collapses toward baseline.

## The ML framing (TA vs DN)

**Unit of analysis:** the gene (optionally per `cell_line`, since the contrast
is computed per line). **Label:** `log2fc_ta_vs_dn` for regression, or
`label_ta_up` (sign) for classification — both in `ml_ta_vs_dn`. **Features:**
promoter architecture aggregated per gene — counts and strong-hit counts per
motif / TF family, max relative score, minimum distance to TSS, and architecture
flags such as `has_dimer_of_dimers`. The wide training frame is produced by
`PIVOT` at read time (`queries.sql` Q3), so the stored tables stay long and the
feature set can grow without a schema migration.

**Keep the model sequence-only by default.** Features come from layers 1–2
(architecture), which are deterministic from sequence. The CUT&RUN differential
(layer 3) is the *cause-adjacent* signal for TA-vs-DN expression, so feeding it
in as a feature largely leaks the label. Treat CUT&RUN as either a validation
target (Q4) or a clearly-separated second feature block for an
"architecture + binding" ablation — never silently mixed into the sequence-only
matrix.

**Leakage cautions that matter here:**

- **Positional leakage** — split train/test by **chromosome**, not random rows, so paralogous/neighbouring promoters don't straddle the split.
- **Paralog leakage** — the TP53 family (TP53/TP63/TP73) shares motif structure; hold out whole families/paralog groups when estimating generalization.
- **Label leakage via p73 itself** — decide deliberately whether the `TP73` motif (`MA0861.2`) and p73 CUT&RUN belong in the feature set; they are the manipulated factor.
- **Gene-set imbalance** — differential genes cluster in programs (EMT, p53 pathway); report per-set performance (Q6), not just pooled.

## How agents query it "reasonably"

The contract for OpenClaw / Claude / Codex is deliberately small:

1. **One package, one schema.** From the exported package root, open `tables/jaspar2026/jaspar2026.duckdb` and `.read sql/schema.sql`; every layer is then a named view/table over package-relative Parquet paths. This document plus `schema.sql` is the entire surface an agent needs to read to know what exists.
2. **Named, parameterized queries.** `sql/queries.sql` is the canned set (region browse, gene architecture, ML pull, CUT&RUN cross-check, enrichment, gene-set slice). Agents fill `$params`; they don't invent joins. Region and gene lookups prune by the chrom partition + `start` zonemap, so latency is a row-group scan, not a full pass.
3. **Hot promoter cards.** Ordinary agent questions start with `promoter_card`,
   a one-row summary keyed by gene. It contains the promoter span, motif-hit
   totals, strongest PWM-relative score, nearest hit to the TSS, TP53-family
   dimer-of-dimers flag, and count of represented TF families. This keeps
   "tell me about TP73" fast while still pointing deeper queries at the long
   architecture tables.
4. **Stable, versioned identity.** `manifest.json` records the source commit, the JASPAR 2026 SHA-256, the Ensembl release, and row counts, so an answer is reproducible and an agent can state exactly which package it queried.
5. **Read-only by default.** The `.duckdb` is a rebuildable index over Parquet; agents read it. Rebuilds come from the exporter, not from queries.

If a networked surface is later needed, the same file sits behind a thin DuckDB
HTTP endpoint or a small MCP tool exposing exactly the `queries.sql` set — but
the file-plus-SQL form already covers local agent use with no infrastructure.

## Minimal agent tool surface

OpenClaw, Claude, and Codex should expose a small read-only surface rather than
free-form access to raw genome-wide files:

- `resolve_gene(gene_name)` -> stable `gene_id`, chromosome, strand, TSS.
- `get_promoter_card(gene_name)` -> one compact row from `promoter_card`.
- `get_promoter_architecture(gene_name)` -> query Q2 over
  `promoter_arch_feature`.
- `query_region(chrom, start, end)` -> query Q1 over `motif_hit`.
- `query_dense_score_region(motif_id, chrom, strand, score_mode, pseudocount, start, end)` -> query Q7 over `motif_score_dense`.
- `get_score_calibration(motif_id, chrom, score_mode, pseudocount)` -> query Q8 over `motif_cutandrun_score_bin_stats`.
- `get_score_calibration_bins()` -> query Q9 over `score_calibration_bin`.
- `get_cutandrun_promoter_signal(gene_name, cell_line)` -> query Q4.
- `export_ml_matrix(cell_line, feature_set)` -> query Q3 or a pinned Parquet
  export recorded in `manifest.json`.

This is intentionally boring: each tool maps to a named SQL statement, all
parameters are typed, and no agent needs to discover table joins on the fly.

## Open questions to settle before building the exporter

- **Promoter window definition** — the GeneLists `*.promoter.bed` extents (upstream/downstream of TSS) fix every distance feature; pin and record them in the manifest.
- **Differential source** — how `log2fc_ta_vs_dn` is computed from E-MTAB-14704 (per cell line, which contrast, which shrinkage), and whether `skmel29` and `saos2` are modelled separately or pooled.
- **Motif-version consistency** — the run uses `MA0861.2`; the legacy datatable pipeline still references `MA0861.1` (see the earlier review). Reconcile before joining old and new outputs.
- **Ensembl release** — genome vs GTF release for the 2026 run (Makefile pulls 113; the paper used 112). Coordinates are stable on GRCh38, but layer-4 annotation is not.
- **`hit_architecture` scope** — build per-hit decomposition for the TP53 family only at first, or attempt it for all families with `architecture_confidence` flags?

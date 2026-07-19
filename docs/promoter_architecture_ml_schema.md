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
manifest.json              # finalization marker and immutable run identity
task_data/task_id=*/       # one atomically promoted, validated scan batch
  tables/jaspar2026/motif_hit/
    motif_set_id=*/genome_id=*/motif_id=*/score_mode=*/pseudocount=*/
    background_model_id=*/pseudocount_scheme=*/minimum_score=*/
    minimum_pwm_relative_score=*/maximum_pwm_relative_score=*/
    chrom=*/strand=*/n_policy=*/matched_sequence=*/…
tables/jaspar2026/
  motif_set/… genome/… sequence_region/…
  scan_run/… scan_threshold_policy/… scan_task/… scan_file_inventory/…
  motif_metadata/…          # tiny motif dimension: motif_id, name, length, source
  motif_score_dense/
    motif_set_id=*/genome_id=*/motif_id=*/score_mode=*/pseudocount=*/
    background_model_id=*/pseudocount_scheme=*/chrom=*/strand=*/…
                            # range and N policy are encoded in each part filename
                            # dense calibration blocks, score for every PSSM alignment start
  motif_cutandrun_score_bin_stats/
    motif_id=*/score_mode=*/pseudocount=*/chrom=*/…
                            # CUT&RUN enrichment by score bin
  motif_cutandrun_containment_curve/
    motif_id=*/score_mode=*/pseudocount=*/chrom=*/sample_id=*/…
                            # merged-coverage immersion ROC/PR calibration
  motif_architecture/…       # layer 2a (tiny, motif-keyed dimension)
  hit_architecture/chrom=*/… # layer 2b (optional, per-hit, keyed to motif_hit)
  motif_context_pair/genome_id=*/chrom=*/… # TP73 anchor × neighboring hit
  tp73_pair_feature/genome_id=*/chrom=*/… # orientation-collapsed partners
  tp73_context_anchor/genome_id=*/chrom=*/… # tandem/local-context summary
  cutandrun_coverage/chrom=*/… # layer 3 positive coverage intervals/signal
  cutandrun/chrom=*/…        # layer 3 aggregated locus signal
  gene/… promoter/… transcript/… intron/… gene_set/… # layer 4 dimensions
  motif_transcript_context/… # layer 4 per-anchor/per-transcript bridge
  expression/…  expression_differential/…   # layer 4 label source
  jaspar2026.duckdb        # query index over Parquet plus materialized hot tables
```

Only threshold-retained facts are stored for the genome-wide hit layer —
`motif_hit` holds selected sequence matches, `cutandrun` one row per (locus, sample),
`promoter_arch_feature` one row per (gene, motif). Dense calibration runs are
the explicit exception: for selected motifs/chromosomes they store one score
per motif-model alignment start and strand in block vectors, with repeated
identity moved into partition paths.

## Layers and keys

Full column lists live in `sql/schema.sql`; the essentials:

The executable whole-genome batching, validation, inventory, and cross-species
identity contract is specified in
[`jaspar2026_genome_scan_plan.md`](jaspar2026_genome_scan_plan.md). The raw-only
DuckDB surface is `sql/genome_scan_schema.sql`; optional synteny bridges are in
`sql/cross_species_schema.sql`.

- **`motif_hit`** — keyed by `genome_id, motif_set_id, chrom, start, end, motif_id, strand, score_mode, pseudocount`. Its logical configuration also exposes `background_model_id`, `pseudocount_scheme`, `minimum_score`, both PWM-relative bounds, `n_policy`, matched-sequence policy, and coordinate mode. The physical file still stores only `start`, `end`, float32 `score`, float32 `pwm_relative_score`, and nullable `matched_seq`; row-constant identity lives in Hive partitions and `motif_name` comes from `motif_metadata`.
- **`genome / sequence_region / scan_run / scan_threshold_policy / scan_task / scan_file_inventory`** — immutable production provenance and completeness tables. The file inventory includes expected and skipped-N windows, sentinel/threshold/PWM rejections, emitted rows, bytes, SHA-256, task state, Slurm IDs, scanner checksum, and source commit. A zero-hit Parquet file still has an inventory row.
- **`motif_score_dense`** — logical view over dense score blocks, also keyed by explicit `genome_id` and `motif_set_id`. Physical files store only `block_start` plus a `FLOAT[] scores` vector; run identity lives in Hive partitions, while requested range and N policy are encoded in collision-resistant part filenames.
- **`motif_cutandrun_score_bin_stats`** — dense-score calibration output by score bin and CUT&RUN sample: window counts, covered-window counts, baseline fraction, enrichment ratio, log2 enrichment, and signal summaries.
- **`motif_cutandrun_containment_curve`** — score-threshold calibration from
  strict immersion in merged positive coverage: coverage-component recall,
  motif precision, motif recall, false-positive rate, effective depth, ROC AUC,
  and average precision. See
  [`cutandrun_motif_score_calibration.md`](cutandrun_motif_score_calibration.md).
- **`motif_architecture`** — one row per `motif_id`: family, `binding_unit_model`, half-site pattern, spacer range, `architecture_confidence`. Curated for the TP53 family, `unknown` elsewhere.
- **`hit_architecture`** *(optional)* — per-hit decomposition (`spacer_bp`, half-site scores, `oligomer_compatible`), derived from `matched_seq`, kept separate so layer 1 stays stable.
- **`motif_context_pair`** — every retained TP73-anchor/neighbor occurrence
  within the broad capture radius, with genomic and anchor-oriented motif-center
  distance, relative orientation, provisional-context membership, and a strict
  same-motif/distinct-span tandem flag. See
  [`tp73_motif_context.md`](tp73_motif_context.md).
- **`tp73_pair_feature`** — one row per TP73 anchor with mutually exclusive
  singleton/same/opposite/mixed/ambiguous pair class, partner-locus counts,
  orientation-specific gaps and scores, and conservative two-site score
  summaries. Strand alternatives at one neighboring alignment span remain in
  the raw pair table but collapse to one orientation-ambiguous partner locus.
  These are sequence-architecture predictions, not observations of protein
  quaternary structure.
- **`tp73_context_pair_feature`** — the chromosome-wide raw context pairs
  decorated with their anchor's pair class and a 5-bp oriented-distance bin.
  This is the pair-stratified feature surface for per-site CUT&RUN models and
  is not restricted to promoters.
- **`tp73_context_anchor`** — one feature-ready row per TP73 occurrence with
  tandem summary, local motif counts, nearest signed TSS distance, and
  transcript/intron indicators. It never replaces the pair table.
- **`cutandrun_coverage`** — BED/bedGraph/bigWig-derived positive coverage
  components and signal with `sample_id` and sample facets. Overlapping and
  adjacent intervals are merged for motif-immersion labels; depth remains a
  separate signal.
- **`cutandrun`** — `chrom, start, end, cell_line, isoform, antibody, replicate, signal`. The sample facets are columns, so a TA-vs-DN contrast is a `GROUP BY`, not a filename parse.
- **`gene / promoter / transcript / intron / gene_set`** — annotation
  dimensions regenerated from the pinned GTF. Transcript-level intron state is
  kept in `motif_transcript_context` because isoforms can classify one locus
  differently.
- **`expression_differential`** — `gene_id, cell_line, log2fc_ta_vs_dn, padj, …`. **This is the ML label.**

Five materialized tables carry the reusable join/aggregation cost:
`promoter_motif_hit` (the single
interval join of hits into promoter windows, deduplicated to one row per gene
and unique hit while retaining the closest transcript's strand-aware
`tss_distance` and the number of overlapping transcripts),
`promoter_arch_feature` (per-gene, per-motif architecture summary),
`promoter_pair_feature` (TP73 anchors grouped by pair class),
`promoter_motif_pair_feature` (neighbor motifs grouped by TP73 pair class,
orientation, side, and 5-bp distance bin), and
`promoter_card` (one compact row per gene and scoring configuration for agent
lookup). `ml_ta_vs_dn`, `ml_ta_vs_dn_pair`, and
`ml_ta_vs_dn_motif_pair` join the corresponding sequence-only feature tables
to the label without mixing their grains.

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

For coverage data, additionally apply strict immersion in merged positive
coverage rather than ordinary overlap. Evaluate one genomic span using the
maximum score across orientations as the primary unstranded analysis, retain
orientation-specific results secondarily, and report ROC AUC plus the more
imbalance-sensitive precision-recall metrics. The storage threshold must be
chosen across replicates, not from a synthetic fixture or ROC AUC alone.

## The ML framing (TA vs DN)

**Unit of analysis:** the gene (optionally per `cell_line`, since the contrast
is computed per line). **Label:** `log2fc_ta_vs_dn` for regression, or
`label_ta_up` (sign) for classification — both in `ml_ta_vs_dn`. **Features:**
promoter architecture aggregated per gene — counts and strong-hit counts per
motif / TF family, max relative score, minimum distance to TSS, and architecture
flags such as `has_dimer_of_dimers`. Query Q3 returns a parameterized long
family-feature frame for one cell line and scoring configuration. The ML
client/exporter pivots that bound result to wide form, so DuckDB does not have
to discover dynamic columns while parameters are still unresolved and the
stored tables remain long without a schema migration.

Pair-aware models use two additional, opt-in feature blocks. The first counts
TP73 anchors by `pair_class`, including singletons. The second stratifies each
neighboring motif by the anchor pair class, relative orientation, anchor-facing
side, and a stable 5-bp oriented-distance bin. This supports interactions such
as `PATZ1 score x TP73 pair class` while preserving the original continuous
scores. Descriptive performance and calibration should also be reported
separately for singleton, same-orientation, opposite-orientation, and
ambiguous/mixed anchors; sparse strata should not automatically become
independent models.

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
   a summary keyed by gene, score mode, and pseudocount. It contains the
   promoter span, motif-hit totals, strongest PWM-relative score, nearest hit
   to the TSS, TP53-family dimer-of-dimers motif-model flag, count of represented
   TF families, and singleton/tandem TP73 anchor counts. Keeping configurations
   in separate rows prevents incompatible runs being combined while keeping
   "tell me about TP73" fast.
4. **Stable, versioned identity.** `manifest.json` records the source commit, the JASPAR 2026 SHA-256, the Ensembl release, and row counts, so an answer is reproducible and an agent can state exactly which package it queried.
5. **Read-only by default.** The `.duckdb` is a rebuildable index over Parquet; agents read it. Rebuilds come from the exporter, not from queries.

If a networked surface is later needed, the same file sits behind a thin DuckDB
HTTP endpoint or a small MCP tool exposing exactly the `queries.sql` set — but
the file-plus-SQL form already covers local agent use with no infrastructure.

## Minimal agent tool surface

OpenClaw, Claude, and Codex should expose a small read-only surface rather than
free-form access to raw genome-wide files:

- `resolve_gene(gene_name)` -> stable `gene_id`, chromosome, strand, TSS.
- `get_promoter_card(gene_name, score_mode, pseudocount)` -> one compact row
  from `promoter_card`.
- `get_promoter_architecture(gene_name, score_mode, pseudocount)` -> query Q2 over
  `promoter_arch_feature`.
- `query_region(chrom, start, end)` -> query Q1 over `motif_hit`.
- `query_dense_score_region(motif_id, chrom, strand, score_mode, pseudocount, start, end)` -> query Q7 over `motif_score_dense`.
- `get_score_calibration(motif_id, chrom, score_mode, pseudocount)` -> query Q8 over `motif_cutandrun_score_bin_stats`.
- `get_score_calibration_bins()` -> query Q9 over `score_calibration_bin`.
- `get_cutandrun_containment_curve(motif_id, chrom, score_mode, pseudocount, sample_id)` -> merged-coverage immersion threshold metrics.
- `get_tp73_tandem_partners(chrom, start, strand, score_mode, pseudocount)` -> query Q10.
- `get_tp73_gene_context(gene_name, score_mode, pseudocount)` -> query Q11.
- `get_tp73_promoter_pair_features(gene_name, score_mode, pseudocount)` ->
  query Q12.
- `get_tp73_neighbor_features(gene_name, neighbor_motif_id, score_mode, pseudocount)`
  -> query Q13.
- `get_tp73_site_neighbor_features(chrom, start, strand, neighbor_motif_id, score_mode, pseudocount)`
  -> query Q14, without a promoter restriction.
- `get_cutandrun_promoter_signal(gene_name, cell_line)` -> query Q4.
- `export_ml_matrix(cell_line, score_mode, pseudocount, feature_set)` -> query
  Q3, pivot its bound long result in the client, or read a pinned Parquet export
  recorded in `manifest.json`.

This is intentionally boring: each tool maps to a named SQL statement, all
parameters are typed, and no agent needs to discover table joins on the fly.

## Open questions to settle before building the exporter

- **Promoter window definition** — regenerate the historical GeneLists
  `*.promoter.bed` files with the corrected GTF-to-BED conversion before they
  seed this layer. Pin the chosen upstream/downstream extents and coordinate
  convention in the manifest because they fix every distance feature.
- **Differential source** — how `log2fc_ta_vs_dn` is computed from E-MTAB-14704 (per cell line, which contrast, which shrinkage), and whether `skmel29` and `saos2` are modelled separately or pooled.
- **Motif-version consistency** — the run uses `MA0861.2`; the legacy datatable pipeline still references `MA0861.1` (see the earlier review). Reconcile before joining old and new outputs.
- **Ensembl release** — genome vs GTF release for the 2026 run (Makefile pulls 113; the paper used 112). Coordinates are stable on GRCh38, but layer-4 annotation is not.
- **`hit_architecture` scope** — build per-hit decomposition for the TP53 family only at first, or attempt it for all families with `architecture_confidence` flags?

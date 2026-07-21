# TP73 Local Motif Context And Tandem Sites

## Purpose

The published analysis (PMID 41594604, PMCID PMC12839168) summarized other
motifs around predicted TP73 occurrences in a wide table. The updated analysis
keeps the same biological question but treats local motif syntax as a de novo,
long-form machine-learning feature layer. Every retained pair remains
inspectable; max-score, count, tandem, and gene-context features are derived
from that layer rather than replacing it.

The historical result remains available from tag `v0.0.0`. Current defaults are
new-analysis defaults and must not be used to reinterpret the published table.

The first pair-stratified chromosome-1 CUT&RUN benchmark is reported in
[`tp73_pair_stratified_chr1_cutandrun_20260718.md`](tp73_pair_stratified_chr1_cutandrun_20260718.md).
It treats these fields as sequence-architecture features and explicitly tests
additive and pair-by-cofactor predictors.

The first production-shaped schema-v4 package and its resource/cardinality
measurements are recorded in
[`tp73_patz1_context_v4_chr1_20260720.md`](tp73_patz1_context_v4_chr1_20260720.md).

The evaluator offers `--fold-mode interleaved` for cyclic genomic blocks and
`--fold-mode contiguous` for equal-width contiguous spans within each
chromosome. It writes a fold manifest with the actual coordinate extents.
Sample/fold metrics are correlated diagnostics; they must not be counted as
independent replicates for uncertainty or significance claims.

## Distance semantics

All coordinates are BED 0-based, half-open. The primary proximity measurement
for two motif-model alignment spans is:

```text
interval_distance = max(start1, start2) - min(end1, end2)
```

Directly abutting spans have distance `0`, one unoccupied base gives `1`, and
an overlap of `k` bases gives `-k`. `interval_relation` distinguishes identical
spans, either direction of containment, partial overlap, abutting, and disjoint
spans. The mutually exclusive bands are `overlap`, `adjacent_0_5`, `gap_6_20`,
`gap_21_50`, `gap_51_100`, and `gap_101_150`. Cumulative `within_5`,
`within_20`, `within_50`, `within_100`, and `within_150` flags include overlaps
while leaving their relation visible.

For a motif-model span `[start,end)`, its center is `(start + end) / 2`;
centers may lie at half-base coordinates. Center offsets encode direction but
do not decide context membership.

For a TP73 anchor `a` and neighboring motif occurrence `n`:

```text
genomic_distance = center(n) - center(a)
oriented_distance = genomic_distance                 if strand(a) = '+'
                    -genomic_distance                 if strand(a) = '-'
```

`relative_orientation` is `same` when the two reported motif orientations are
equal and `opposite` otherwise. Both distances are stored. This matters because
the oriented sign is useful for motif syntax, whereas the genomic sign remains
unambiguous if a near-palindromic TP53-family match has uncertain orientation.
The scored interval remains a PSSM alignment span, not an asserted physical TF
footprint.

## Context and pair limits

- `capture_flank_bp = 150`: retain genome-wide relationships with signed
  interval distance at most 150 bp. A wider value is an explicit,
  promoter-scoped analysis, not the genome-wide default.
- `context_flank_bp = 150`: interval-distance radius used by ordinary context
  summaries. It may be narrower than capture in a targeted package.
- `tandem_flank_bp = 20`: maximum edge-to-edge gap for a distinct,
  non-overlapping occurrence of the same TP73 motif to be a tandem partner.
- `cofactor_pair_flank_bp = 150`: maximum signed interval distance considered
  when describing distinct same-motif cofactor loci as a sequence pair.

The center prefilter used to accelerate capture is the deliberately
conservative bound `capture + max_anchor_span + max_neighbor_span`, with widths
derived from the retained input. Cofactor-pair capture analogously uses
`pair_flank + 2 * max_neighbor_span`. The exact interval predicate is applied
afterward. Schema-v4 provenance records the observed widths, derived
prefilters, and interval capture geometry.

An opposite-strand record with exactly the same alignment span is retained as
an orientation alternative but is not a tandem partner and does not inflate
the anchor's ordinary neighboring-hit/locus counts. A tandem pair must use
the same motif accession, have a different non-overlapping span, and have an
edge-to-edge gap of at most 20 bp. Both orientation-specific TP73 scores must
also be at least the configured tandem minimum; the production value is `0`.
Shifted overlapping matches remain in the pair table but are not called
tandem. Pair rows retain center distance, overlap, gap, and every partner; the
anchor table also
provides the partner count and the nearest partner's distance, side,
orientation, and score.

Anchor retention and tandem eligibility are deliberately separate. The raw
TP73 scan remains available to `-5`, but the conservative production context
uses `--anchor-minimum-score -1`. With `--anchor-selection-mode local_peak`, an
orientation is retained as an anchor only if it is the best score at its
physical span and no strictly stronger physical TP73 locus starts within
`--anchor-local-peak-flank` (150 bp). The rule applies to positive and negative
scores. Equal-scoring regional maxima are retained rather than broken
arbitrarily. A lower positive TP73 occurrence can therefore contribute to the
selected peak's tandem architecture without materializing a nearly identical
cofactor neighborhood as a second anchor.

`--tandem-minimum-score` is applied symmetrically to both TP73 members before
`is_tandem_tp73` is set; `--partner-minimum-score` remains a legacy command-line
alias. The floors, selection mode, local competitor score, prominence, and
selection class are recorded in the derived package. A rejected TP73 row is
still available in the raw scan and may still occur as a neighbor of a retained
anchor. Thus region-specific analyses can relax the derived-layer policy
without rescanning the entire genome.

## Pair-feature semantics

`motif_context_pair` deliberately remains orientation-specific. A
near-palindromic alignment may therefore have both `+` and `-` records at the
same neighboring span. For ML features, `tp73_pair_feature` collapses those
records to one partner locus keyed by its alignment span. Its orientation is
`ambiguous` when both strand alternatives are present; this avoids counting an
uncertain orientation as two physical partners while retaining both source
records for inspection.

Every TP73 anchor receives exactly one `pair_class`:

- `singleton`: no distinct, non-overlapping TP73 partner within the tandem gap;
- `tandem_same_orientation`: all resolved partner loci have the anchor's
  reported orientation;
- `tandem_opposite_orientation`: all resolved partner loci have the other
  orientation;
- `tandem_mixed_orientation`: distinct resolved partner loci include both; or
- `tandem_orientation_ambiguous`: at least one partner locus has both strand
  alternatives.

Counts, nearest gaps, best partner scores, the best minimum score of a pair,
and the best score sum remain separate columns. These fields describe sequence
architecture compatible with one or more bound complexes. They do not assert
that a particular oligomer or quaternary structure was experimentally present.
The two half-sites within one TP53-family motif alignment remain a separate
question represented by the optional `hit_architecture` layer.

The package contract requires one `tp73_pair_feature` row per `anchor_hit_id`,
the same/opposite/ambiguous partner counts to sum to the total partner-locus
count, and singleton partner-derived scores to remain `NULL` rather than zero.
`tests/test_motif_context.sh` enforces these invariants on a synthetic package.

`tp73_context_pair_feature` joins the resulting anchor pair class back onto
every raw neighboring-motif record. This chromosome-wide view is the direct
input for pair-stratified CUT&RUN predictors; promoter aggregation is an
additional downstream view, not a prerequisite.

The 150 bp value is not frozen. It is also the default range used by SpaMo, a
method that independently formalized strand-separated motif-spacing tests, but
that precedent is not evidence that 150 bp is optimal for TP73
([Whitington et al. 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3159476/),
[SpaMo documentation](https://meme-suite.org/meme/doc/spamo.html)).

## Cofactor pair architecture

Generic cofactor pairing is computed after collapsing identical-span strand
alternatives into one locus keyed by genome, chromosome, motif accession, start,
and end. Such a locus is `plus`, `minus`, or `ambiguous`. Canonical pairs use
two distinct loci of the same motif accession and are stored once in genomic
left-to-right order. Their arrangement is `codirectional_plus` (`>>`),
`codirectional_minus` (`<<`), `convergent` (`><`), `divergent` (`<>`), or
`ambiguous`; ambiguity takes precedence when either member has both strand
alternatives.

This is a TP73 context package, not a second chromosome-wide motif inventory.
A cofactor locus captured around at least one TP73 anchor seeds pair formation;
every same-motif partner within `--cofactor-pair-flank` is retained even when
that second member falls just outside the anchor's context. Pairs for which
neither member occurs in any TP73 context are omitted. This preserves the
one-member boundary case while avoiding a global pair table that cannot
contribute to a TP73 feature.

`cofactor_motif_locus` therefore contains context loci plus their retained pair
partners and marks the former with `in_any_tp73_context`.
`cofactor_motif_pair` retains exact member geometry and scores, while
`cofactor_locus_pair_feature` describes same-motif partner architecture only
for the context loci whose features can be attached to TP73. A canonical pair
is attached to an individual TP73 anchor when at least one member is within the
configured context distance. `tp73_cofactor_pair_context` records whether one
or both members qualify and attributes side/direction through the nearest
qualifying member. `tp73_cofactor_pair_summary` is the compact nonduplicating
feature surface. These are sequence-compatible pair architectures, not
observations of protein oligomerisation.

## Transcript and intron context

When `--gtf` is supplied, `build_motif_context.py` reads the pinned GTF directly
through DuckDB and writes Parquet; no BED or TSV annotation intermediate is
created. GTF 1-based inclusive coordinates are converted to BED coordinates.
Transcript bounds come from transcript/exon records, and introns are the gaps
between consecutive exons of each transcript.

The output deliberately preserves transcript ambiguity:

- `nearest_tss_distance_bp` is signed in the nearest transcript's direction;
  positive means downstream of transcription initiation.
- `motif_transcript_context` has one row per anchor and containing transcript.
- `motif_transcript_context_pair` derives cofactor direction in each
  transcript's frame, signed TSS distances for both spans, and whether either
  scored span crosses the TSS. It is a view and does not multiply the primary
  chromosome-wide ML rows.
- `fully_within_intron` requires the complete motif-model span to lie within an
  intron. `overlaps_intron` also detects boundary-crossing spans.
- `in_any_intron` is a convenient anchor-level summary, but does not erase the
  fact that the same locus can be intronic for one isoform and non-intronic for
  another.

## Building the Parquet package

For a finalized production scan, first select exact files through its immutable
inventory. This creates hard links only; it neither copies payload data nor
walks the package-wide Parquet tree:

```sh
scripts/stage_motif_context_inputs.py \
  --package RUN/package --output CONTEXT_RUN/input/MA1961.2/1 \
  --motif MA0861.2 --motif MA1961.2 --chrom 1
```

The resulting `input_manifest.json` pins the source package-manifest checksum,
selection, inventory paths, recorded file sizes, and checksums. The preserved
task/Hive paths can be passed to the builder as
`CONTEXT_RUN/input/MA1961.2/1/**/*.parquet`. Keep the input tree on the same
filesystem as the scan package so hard-link creation remains metadata-only.
Use neutral names for its parent directories rather than Hive labels such as
`motif_id=...` or `chrom=...`: DuckDB reads labels from every path component,
so an outer label would override the genuine identity in the linked scan paths.
The stager rejects such conflicting wrappers.

On Slurm, submit one cofactor/chromosome package per array task. The worker
uses exact inventory hard links on `/data`, a unique DuckDB spill directory on
node-local `/scratch`, and refuses to overwrite either an incompatible input
selection or an existing output:

```sh
scripts/submit_motif_context_slurm.sh \
  --run-root /data/sm718/jaspar_mapping_runs/jaspar2026_grch38_context_v4 \
  --scan-package /data/sm718/jaspar_mapping_runs/jaspar2026_grch38_sparse_v3/package \
  --gtf /data/sm718/resources/ensembl/113/Homo_sapiens.GRCh38.113.gtf.gz \
  --motif MA1961.2 --chrom 1 \
  --account cluster --partition requeue \
  --max-concurrent 1 --cpus 4 --memory 32G \
  --memory-limit 24GB --max-temp-size 100GB --dry-run
```

Inspect the rendered command before removing `--dry-run`. For a selected
panel, repeat `--motif`; for genome-wide summary production, use
`--output-tier summary` and repeat or comma-separate chromosomes. The immutable
task plan makes a changed selection fail rather than silently mix packages.

```sh
scripts/build_motif_context.py \
  --motif-hits 'RUN/task_data/*/tables/jaspar2026/motif_hit/**/*.parquet' \
  --gtf Homo_sapiens.GRCh38.113.gtf.gz \
  --output RUN/tp73_context \
  --anchor-motif MA0861.2 \
  --motif-set-id jaspar2026_core_nonredundant \
  --genome-id homo_sapiens_grch38_ensembl113_primary \
  --anchor-minimum-score -1 --tandem-minimum-score 0 \
  --anchor-selection-mode local_peak --anchor-local-peak-flank 150 \
  --score-mode log2_relative_risk --pseudocount 1 \
  --capture-flank 150 --context-flank 150 --tandem-flank 20 \
  --cofactor-pair-flank 150 --output-tier selected
```

`--output-tier selected` retains inspectable raw TP73-neighbor, cofactor-pair,
and pair-attribution rows. `--output-tier summary` writes schema-stable empty
raw surfaces while retaining `tp73_motif_context_summary` and
`tp73_cofactor_pair_summary`; use that mode for motif-at-a-time all-JASPAR
production.

The package contains:

- `tables/jaspar2026/motif_context_pair/genome_id=*/chrom=*/`: all retained local pairs;
- `tables/jaspar2026/tp73_anchor_locus/genome_id=*/chrom=*/`: one physical
  TP73 alignment span with explicit orientation state;
- `tables/jaspar2026/tp73_motif_context_summary/`: nonempty motif-, band-,
  side-, and orientation-stratified groups partitioned by neighboring motif;
- `tables/jaspar2026/cofactor_motif_locus.parquet` and
  `cofactor_motif_pair.parquet`: TP73-context loci plus necessary outside pair
  partners, and canonical pairs having at least one context member;
- `tables/jaspar2026/tp73_cofactor_pair_context.parquet` and
  `tp73_cofactor_pair_summary.parquet`: pair attribution and compact features;
- `tables/jaspar2026/tp73_pair_feature/genome_id=*/chrom=*/`: one orientation-collapsed
  pair-state row per TP73 anchor, including its retention class and regional
  score prominence;
- `tables/jaspar2026/tp73_context_anchor/genome_id=*/chrom=*/`: one ML-friendly row per
  TP73 occurrence;
- `tables/jaspar2026/motif_transcript_context.parquet`: per-transcript location;
- regenerated `tables/jaspar2026/transcript.parquet` and
  `tables/jaspar2026/intron.parquet` annotation dimensions;
- `tables/jaspar2026/context_run_config.parquet`: complete radius/source
  semantics; and
- `context.duckdb`: small read-only views over those Parquet files.

Signed interval capture, TP73 locus grain, distance bands, conservative
regional-peak anchor selection, symmetric tandem score eligibility, and generic
cofactor pair architecture are context schema version 4. Explicit
genome/motif-set identity and independent anchor/partner score floors appeared
in version 3;
the TP73 pair-class contract appeared in version 2. Older context packages can
be rebuilt from their motif-hit Parquet files; no genome rescan is required.

The builder reports source size and free output-disk space before starting.
Raw TP73-neighbor and same-motif cofactor-pair tables can be substantially
larger than the motif-hit input, so this is a warning rather than a reliable
output-size estimate. Process all-JASPAR summaries one neighboring motif at a
time; retain raw pair tables only for selected motifs. DuckDB is bounded by
`--memory-limit` and `--max-temp-size`. On a cluster, point
`--temp-directory` at a new, empty job-specific directory on node-local
`/scratch`; the final package remains on durable storage while join spill does
not load the shared filesystem. The builder never removes an externally
supplied spill directory, leaving that lifecycle to the scheduler.

Both canonical `+`/`-` strands and the direct sparse writer's Hive values
`plus`/`minus` are accepted and normalized to `+`/`-`. `motif_name` is optional
in the source Parquet; when omitted by the sparse writer, the context package
uses the motif accession as its display name.

Inspect from the package root:

```sh
duckdb context.duckdb -c \
  "SELECT * FROM tp73_context_anchor WHERE has_tandem_tp73 LIMIT 20"

duckdb context.duckdb -c \
  "SELECT relative_orientation, anchor_oriented_center_distance_bp, count(*) \
   FROM motif_context_pair WHERE is_tandem_tp73 GROUP BY ALL ORDER BY 2"

duckdb context.duckdb -c \
  "SELECT pair_class, count(*) AS n_anchors, \
          sum(n_tandem_tp73_partner_loci) AS n_partner_loci \
   FROM tp73_pair_feature GROUP BY pair_class ORDER BY pair_class"
```

## Assessing the 150 bp boundary

The genome-wide package deliberately stops at 150 bp. Count neighboring loci
by the stored exclusive bands, orientation, TP73 score stratum, and transcript
region, then compare them with chromosome-, GC-, mappability-, regulatory-, and
TP73-score-matched controls. Raw counts are strongly affected by motif
frequency, so estimate excess occurrence density as well.

Persistent excess in `gap_101_150` is evidence for a targeted follow-up, not an
automatic genome-wide expansion. Promoter-associated anchors or selected motif
accessions may be rebuilt with a wider explicit capture while retaining the
same interval geometry and provenance. Radius and feature stability must
ultimately be checked on held-out chromosomes and biological replicates rather
than tuned and reported on chromosome 1 alone.

The implemented pair-aware ML surface includes TP73 anchor counts by pair
class, partner-locus counts and scores, and neighboring-motif counts stratified
by anchor pair class, relative orientation, anchor-facing side, and stable 5-bp
oriented-distance bins. Models should retain the continuous TP73 and neighbor
scores and add interactions such as pair-class-by-PATZ1 or tandem-by-intron.
Report performance within pair classes, but use one regularized model initially
rather than fragmenting sparse classes into independent fits. CUT&RUN remains
an experimental-evidence block and is not silently mixed into the sequence-only
architecture features.

For threshold sensitivity against the strict CUT&RUN immersion label,
`scripts/build_sparse_context_maxima.py` derives interval-defined maxima for a
specified set of exact single-chromosome motif-hit partitions. It writes a
source-file run config beside its Parquet output and rejects multi-chromosome
anchors because the sparse partition files do not repeat chromosome in every
row. `scripts/evaluate_tp73_cofactor_thresholds.R` then compares score floors
with contiguous chromosome folds. The chromosome-1 POU2F2, POU4F1, and TCF7
result is documented in
[`tp73_pou4f1_tcf7_chr1_thresholds_20260721.md`](tp73_pou4f1_tcf7_chr1_thresholds_20260721.md).

Do not conflate the TP73 source scan floor with context anchor eligibility. The
production TP73 scan may retain scores down to -5 while a context package can
conservatively select local-peak anchors only from scores at least -1. Both
values must remain explicit provenance.

## Legacy wide exporter

`context` remains available for the historical R workflow, now with explicit
`--flank` (default 150), center distances, correct counts, strict input errors,
embedded provenance columns, exact-anchor exclusion, and bounded batching. It
opens one neighboring file at a time and uses a checked temporary spool rather
than requiring thousands of simultaneous file descriptors. `_Shift` is the
anchor-oriented center displacement; `_GenomicShift` retains genomic direction.
The main BED is expected to contain one fixed-width anchor motif and all inputs
must be chromosome/start sorted.
Only the strongest neighbor per motif is represented in that wide compatibility
table. New analyses should use `motif_context_pair` whenever all observations or
radius selection matter.

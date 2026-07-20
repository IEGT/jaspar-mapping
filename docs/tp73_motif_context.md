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

The evaluator offers `--fold-mode interleaved` for cyclic genomic blocks and
`--fold-mode contiguous` for equal-width contiguous spans within each
chromosome. It writes a fold manifest with the actual coordinate extents.
Sample/fold metrics are correlated diagnostics; they must not be counted as
independent replicates for uncertainty or significance claims.

## Distance semantics

All coordinates are BED 0-based, half-open. For a motif-model alignment span
`[start,end)`, its center is `(start + end) / 2`; centers may lie at half-base
coordinates when motif lengths differ.

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

## Three radii

- `capture_flank_bp = 500`: retain every local pair this far from TP73. This is
  the reusable observation layer and is deliberately wider than the current
  feature definition.
- `context_flank_bp = 150`: provisional radius for ordinary local-context
  features. It means 150 bp on each side, measured center to center.
- `tandem_flank_bp = 20`: maximum edge-to-edge gap for a distinct,
  non-overlapping occurrence of the same TP73 motif to be a tandem partner.

An opposite-strand record with exactly the same alignment span is retained as
an orientation alternative but is not a tandem partner and does not inflate
the anchor's ordinary neighboring-hit/locus counts. A tandem pair must use
the same motif accession, have a different non-overlapping span, and have an
edge-to-edge gap of at most 20 bp. Shifted overlapping matches remain in the
pair table but are not called tandem. Pair rows retain center distance, overlap,
gap, and every partner; the anchor table also
provides the partner count and the nearest partner's distance, side,
orientation, and score.

Anchor retention and tandem-partner eligibility are deliberately separate.
`--anchor-minimum-score` controls which TP73 rows become modeled anchors;
`--partner-minimum-score` controls only whether another TP73 row can set
`is_tandem_tp73`. Both values are stored on every derived pair/anchor row and in
`context_run_config`. A source scan retained to `-5` can therefore be reused to
compare partner floors `-5`, `-1`, and `0` without changing the raw observations
or rescanning sequence. A below-floor TP73 neighbor remains inspectable in the
raw pair table but is not called a tandem partner at that floor.

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
- `fully_within_intron` requires the complete motif-model span to lie within an
  intron. `overlaps_intron` also detects boundary-crossing spans.
- `in_any_intron` is a convenient anchor-level summary, but does not erase the
  fact that the same locus can be intronic for one isoform and non-intronic for
  another.

## Building the Parquet package

```sh
scripts/build_motif_context.py \
  --motif-hits 'RUN/task_data/*/tables/jaspar2026/motif_hit/**/*.parquet' \
  --gtf Homo_sapiens.GRCh38.113.gtf.gz \
  --output RUN/tp73_context \
  --anchor-motif MA0861.2 \
  --motif-set-id jaspar2026_core_nonredundant \
  --genome-id homo_sapiens_grch38_ensembl113_primary \
  --anchor-minimum-score -5 --partner-minimum-score -1 \
  --score-mode log2_relative_risk --pseudocount 1 \
  --capture-flank 500 --context-flank 150 --tandem-flank 20
```

The package contains:

- `tables/jaspar2026/motif_context_pair/genome_id=*/chrom=*/`: all retained local pairs;
- `tables/jaspar2026/tp73_pair_feature/genome_id=*/chrom=*/`: one orientation-collapsed
  pair-state row per TP73 anchor;
- `tables/jaspar2026/tp73_context_anchor/genome_id=*/chrom=*/`: one ML-friendly row per
  TP73 occurrence;
- `tables/jaspar2026/motif_transcript_context.parquet`: per-transcript location;
- regenerated `tables/jaspar2026/transcript.parquet` and
  `tables/jaspar2026/intron.parquet` annotation dimensions;
- `tables/jaspar2026/context_run_config.parquet`: complete radius/source
  semantics; and
- `context.duckdb`: small read-only views over those Parquet files.

Explicit genome/motif-set identity and independent anchor/partner score floors
are context schema version 3. The pair-class contract first appeared in version
2. A version-1 context
package still preserves its raw pair observations, but it must be rebuilt from
the existing motif-hit Parquet files to obtain `tp73_pair_feature`; no genome
rescan is required.

The builder reports source size and free output-disk space before starting.
The pair table can be substantially larger than the motif-hit input, so this is
a warning rather than a reliable output-size estimate; DuckDB is also bounded
by `--memory-limit` and `--max-temp-size`.

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

## Selecting the eventual context radius

Use the 500 bp capture table to estimate distance distributions without a new
genome scan:

1. Count all neighboring occurrences in 1 bp bins and display stable 5 bp bins,
   separately for motif, same/opposite orientation, TP73 score stratum, and
   transcript region.
2. Compare with chromosome-, GC-, mappability-, regulatory-region-, and
   TP73-score-matched control anchors. Estimate excess occurrence density rather
   than raw counts alone.
3. Find the smallest rounded radius containing approximately 95% of positive
   excess density and require the following 50 bp to be statistically
   indistinguishable from background. If signal remains at 500 bp, recapture at
   1,000 bp.
4. Verify radius and feature stability on held-out chromosomes and biological
   replicates. Do not tune and report performance on chromosome 1 alone.

The implemented pair-aware ML surface includes TP73 anchor counts by pair
class, partner-locus counts and scores, and neighboring-motif counts stratified
by anchor pair class, relative orientation, anchor-facing side, and stable 5-bp
oriented-distance bins. Models should retain the continuous TP73 and neighbor
scores and add interactions such as pair-class-by-PATZ1 or tandem-by-intron.
Report performance within pair classes, but use one regularized model initially
rather than fragmenting sparse classes into independent fits. CUT&RUN remains
an experimental-evidence block and is not silently mixed into the sequence-only
architecture features.

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

# A Small Promoter Cofactor Package For Glen And GENtle

## Scope

The exporter `scripts/export_promoter_collaboration.py` builds a relocatable,
read-only analytical package from **finalized** TP73 cofactor results. Its first
production target is the common-score-zero, Ensembl **extended-promoter**
autosomal analysis. That source contains 34,035 eligible TP73 anchors; retain
them all, including negative CUT&RUN observations. A selected cofactor panel
must not become a selected set of positive anchors.

The complete vertebrate motif-level results are small and travel with the
package. Detailed positions are restricted to a deterministic panel of up to
100 motifs by default. Exact-name matches for the prespecified candidates are
selected first, followed by a round-robin of the TA-enriched, DN-enriched,
TA-depleted, DN-depleted and absolute isoform-difference rankings in each of
the six exclusive distance bands. This is an exploratory presentation panel,
not a new statistical test or a new BH correction family. Original estimates,
confidence intervals, adjusted p-values and support flags are unchanged.

Source species flags are retained. The overview covers vertebrate motifs, not
only matrices originally measured in humans. No substitute motif is invented
for a requested name absent from the result catalog. In particular Glen's
SUB1, HMGB3 and YBX3 have no exact-name matrix in the current source, whereas
IRF9 has MA0653.1. The generated manifest lists availability explicitly.

The target is less than 1 GB; **the measured publication ceiling is 2 GB**,
including the DuckDB index and companion files. No size claim is made before
the production files are measured. An oversized attempt remains unpublished;
use a new run with a smaller panel rather than silently dropping rows. The
underlying complete genome scan and raw CUT&RUN tracks remain on Haumea.

## Contents And Limits

- `anchors.parquet`: physical TP73 anchor keys and the available summarized
  CUT&RUN support from the upstream evidence split, joined back to the pinned
  source evidence for matching maximum-depth measurements. The
  within-package numeric `anchor_id` is not stable across different cohorts;
  use assembly, chromosome, start and end for cross-package joins.
- `anchor_promoter`, `promoter`, `promoter_gene`: the existing, audited
  many-to-many extended regulatory-promoter membership, physical features and
  gene links. Native Ensembl ownership and TSS-derived links retain their
  distinct provenance. These are NOT substituted fixed 2 kb/500 bp promoters.
- `cofactor_distance_isoform_comparison`, `cofactor_distance_enrichment` and
  `cofactor_distance_frequency_enrichment`: all upstream motif overview rows,
  including non-estimable results. JASPAR matrix/species metadata also travel.
- `feature_MA*.parquet`: one strongest physical cofactor locus per included
  anchor, motif and exclusive distance band; counts at scores >= -1 and >= 0;
  best position, interval distance, genomic side and both retained strand
  scores. A strongest-strand tie is `.`. The other orientation being absent
  means its score was below the source floor, not zero.
- `cofactors.duckdb`, `schema.sql`, `overview.html`, a standalone query/export
  helper, and checksummed provenance. The index contains views, not duplicate
  analytical tables, and can be rebuilt from `schema.sql` in a new database.

Coordinates remain **GRCh38 BED 0-based half-open** throughout. Distance is
`max(anchor_start, hit_start) - min(anchor_end, hit_end)`: overlaps are negative
and abutting intervals have distance zero. Bands are overlap, 0-5, 6-20, 21-50,
51-100 and 101-150 bp. The exact interval predicate follows a conservatively
widened center-bin prefilter; it does not lose long motifs near the boundary.
For equal best scores choose smallest absolute interval distance, then lowest
start and end. Both strands at one physical span count as **one occurrence**.

Only TP73 anchors must be within the promoter cohort. Cofactor neighbours can
extend outside a promoter: clipping them would change the published estimand.
No cofactor-pair or TP73 quaternary-structure claim can be reconstructed from
this strongest-only subset. That requires the full context/pair package.

Missing sparse rows imply zero retained occurrences only for included motifs
and anchors. The `anchor_features()` query restores all six rows, zero counts
and NULL best scores. It rejects motifs outside the selected panel. A span
without an included anchor is **not** evidence of no motif genome-wide.
Counts at arbitrary other cutoffs cannot be reconstructed from maxima. Presence
at any cutoff >= -1 can be tested; an all-150 count is the sum of exclusive
counts, but an all-150 frequency must be recalculated, never summed.

There are no raw bigWigs, genome FASTA, H3K4me3 model effects or X/Y results in
this first package. It reports CUT&RUN association, not proven causal cofactors.
SAOS2 and SK-MEL-29_2 remain distinct experimental series; SK-MEL-29_1 is
excluded. Scores, source floor and analytical threshold are different concepts.

## Execution And Restart Safety

Requirements: Python 3.9+ and DuckDB CLI, with no Python data dependencies or
runtime downloads. `--help` documents every command and the submission wrapper.

1. `prepare` requires the completed upstream schema-6 distance finalizer and
   audited promoter membership. It validates hashes, pins the result/config,
   copies the small overview and anchor tables, and derives a panel. It reads
   the **exact** low-floor scan file inventory from that analysis. A missing,
   duplicated or floor-censored motif/chromosome/strand is an error.
   The explicit genome and score configuration is recovered from the pinned
   scan's Hive partition keys and retained in the collaborator manifest; mixed
   configurations are rejected.
   The original anchor-evidence Parquet is staged once to scratch to retrieve
   the selected anchors' depth values, verifying its source hash and exact
   agreement with the split support flags. Depth is the upstream maximum over
   the motif interval conditional on strict immersion in a covered component;
   unsupported anchors retain depth zero. This is not a raw read-depth track.
2. `motif --task-index N` stages the current chromosome's two exact source
   Parquet files to node-local scratch, checks their hashes, collapses physical
   loci and retains per-band maxima/counts. It does not rescan DNA. A durable
   checksum checkpoint is published for each chromosome.
3. `finalize` requires every selected motif/chromosome checkpoint and compares
   reconstructed presence counts to the upstream finalized frequencies. It
   writes relative-path views and publishes only after the complete byte cap
   passes. Requeue reuses validated checkpoints; locks serialize duplicate
   workers. Interrupted durable attempts are preserved and never counted as
   complete. Only disposable job scratch is cleaned.

`submit_promoter_collaboration_slurm.sh` queues this chain after the existing
promoter finalizer succeeds. It defaults to 4 concurrent motif tasks, 2 CPUs,
32 GB and 2 hours per motif; preparation/finalization use 8 GB and 1 hour.
All use account `cluster`, partition `requeue`. Durable files stay under a new
dedicated `/data/sm718` directory. Source reaches Haumea via Git, not source
file copying. Submission IDs are journalled immediately; an existing journal
blocks accidental duplicate chains. `afterok` can block forever after an
upstream failure, so inspect failed dependencies rather than assuming export
success. Interrupted attempt directories can use extra space beyond the final
package cap; they are not part of the delivered package.

## GENtle Access

The first integration uses GENtle's existing **assembly-bound region import**,
not its full-genome motif scan provider. The reduced package must never pretend
to be a complete scan with complete retention. No GENtle Rust change is needed
for this route.

```sh
python3 scripts/export_promoter_collaboration.py inspect --package /path/to/package
python3 scripts/export_promoter_collaboration.py region --package /path/to/package \
  --motif MA0653.1 --chrom 1 --start 1000000 --end 1100000 \
  --output /path/to/new_region_example
```

The example is a coordinate query, not a claim that IRF9 occurs there. The
command rejects unrepresented chromosomes/motifs, empty anchor selections,
existing destinations and results exceeding the default 2,000 total anchor-band/annotation-row
cap. It exports all six bands for each overlapping included anchor, preserving
full cofactor spans even when they extend outside the requested viewing range.

The new directory contains `regions.bed`, `evidence.json` and
`gentle_workflow.json`. From that directory execute
`gentle_cli workflow gentle_workflow.json`. This imports the spans, then exports
GENtle's canonical `regions.gentle.json`; GENtle computes its own content and
identity hashes. Raw motif scores and CUT&RUN evidence remain in the sidecar,
not in the BED score column: BED score is explicitly the placeholder 0.
The sidecar also carries promoter memberships and gene ownership, separately
from the anchor rows. Promoter spans are included as additional BED regions;
no observations are duplicated for promoters shared by several genes.
The workflow does not load a reference genome or automatically put the
statistical sidecar into GENtle's UI. A future native collaborator-package
adapter should expose these same bounded queries and explicit coverage limits.

The tests run the complete exporter on a tiny two-chromosome fixture, verify
negative anchors, overlapping promoter ownership, distance 150 vs 151, strand
ties, cross-chromosome keys, counts, floor/cap failures, requeue and relocation.
When a local `gentle_cli` is available, they also import the generated workflow
and verify a canonical GRCh38 region export. Run:

```sh
python3 tests/test_promoter_collaboration.py
```

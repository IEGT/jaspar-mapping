# Whole-genome TP73 CUT&RUN and cofactor follow-up

This document defines the production extension of the chromosome-1 TP73
analysis to the finalized JASPAR 2026 GRCh38 scan. It covers TP73/control
CUT&RUN evidence, convenient-threshold context counts, and all-motif
enrichment/depletion. The schema-9 annotation rebuild and GFP-referenced
H3K4me3 change extension are defined as subsequent, independently restartable
stages.

## Scientific partitions

The primary inference set is autosomes 1-22. X and Y are built and retained as
a separate sensitivity set. The mitochondrial sequence is retained only as a
bystander control because high mitochondrial copy number can distort CUT&RUN
coverage; it never contributes to the primary fit.

Mitochondrial labels `M`, `MT`, `25`, and the corresponding `chr` forms are
aliases only inside that explicitly identified partition. The immutable task
plan records both the scan chromosome and the schema-7 annotation chromosome,
so an atlas using `25` can safely consume an anchor package using `MT`. No such
aliasing is applied to autosomes or sex chromosomes.

## Stage 1: TP73/control evidence

`submit_tp73_genome_evidence_slurm.sh` plans one restartable task per finalized
scan sequence region represented in `scan_file_inventory`. Every admitted
region must have exactly one finalized file per catalog motif and strand;
FASTA-catalog contigs that were not part of the completed scan are not silently
promoted into evidence tasks. Each task resolves exactly one schema-7
`tp73_context_anchor` file and retains physical TP73 local peaks at score
`>= -1`. It exports only TP73 and matched negative-control tracks for SaOS-2
and SkMel-29 series 2 in GFP, TA, and DN conditions. SkMel-29 series 1 remains
excluded by design.

For every anchor and track, support means that one merged positive-coverage
component starts strictly before the motif span and ends strictly after it.
Adjacent positive intervals may join. The accompanying depth is the maximum
coverage over the motif span within that embracing component. A chromosome
absent from a track is represented by support `false` and depth `0`, not by a
missing chromosome or dropped anchor.

The final package contains separate Parquet tables for autosomes, sex
chromosomes, and the mitochondrial bystander control, plus the exact
chromosome-file inventory and a DuckDB catalog. Every support column must equal
`depth > 0`.

## Stage 2: spatial and threshold context

The existing schema-7 context-band packages remain the authoritative spatial
layer. Their `anchor_motif_band_feature` surface preserves the mutually
exclusive interval-distance bands `overlap`, `adjacent_0_5`, `gap_6_20`,
`gap_21_50`, `gap_51_100`, and `gap_101_150`, with one coupled strongest locus,
orientation, side, and same-motif pair architecture per band. Do not replace
that surface with one undifferentiated 150 bp maximum.

The complementary whole-autosome threshold run is prepared by
`submit_tp73_genome_context_maxima_slurm.sh`. For every non-TP73 motif it emits
one rectangular row per TP73 anchor over the complete 150 bp interval-distance
radius. The row stores the strongest score and the number of distinct physical
motif spans reaching the convenient threshold. Opposite-strand reports of one
span count once; anchors with no qualifying span remain explicit zeroes.

Geometry provenance distinguishes the maximum span observed in a chromosome's
hit payload, the catalog-declared motif span, and the effective maximum used by
the conservative prefilter. Thus an empty hit partition records an observed
span of zero without losing the motif length needed to prove complete capture.

The applied context threshold is

```text
max(coalesce(chromosome_1_context_recommendation, 0),
    finalized_genome_scan_retention_floor)
```

Both inputs, the applied result, fallback status, and whether the scan floor
raised the recommendation are retained. Preflight proves that the source
registry covers every finalized non-TP73 scan motif. A motif with empty hit
partitions still receives valid geometry from `motif_metadata.motif_length`.

Array jobs handle a bounded motif batch, stage the 22 chromosome anchor files
once, and process chromosomes sequentially. Each completed motif is published
atomically before the next motif starts. Requeue or timeout therefore reuses
already completed motifs. The compact final catalog does not duplicate the
large per-motif payloads; it provides exact paths and synonymous
`context_maxima_files(file_paths)` and
`tp73_motif_threshold_count_files(file_paths)` query macros.

## Stage 3: TP73 enrichment/depletion

The existing all-JASPAR analyzer is then run against the autosomal anchor
evidence and each motif's zero-complete context table. It reports:

- TP73 score strata `[-5,-1)`, `[-1,0)`, `[0,1)`, `[1,2)`, `[2,5)`,
  `[5,10)`, `[10,15)`, and `[15,+Inf)`;
- strict immersion and positive-depth quantile tiers for every TP73/control
  track pair;
- descriptive anti-p73 versus control risk-ratio differences; and
- one primary discordant matched anti-p73/control occupancy model per motif,
  with sample and chromosome fixed effects, a TP73-score spline, and standard
  errors clustered by chromosome-qualified 5 Mb blocks.

The requested negative references remain strict `< -1` and `< 0`; raising the
positive threshold never turns intermediate scores into negatives. A negative
reference is observable only when the motif's actual scan retention floor is
at or below that reference. When it is not, absent context rows are censored
rather than known negatives, and that contrast is emitted as
`negative_reference_below_source_floor` instead of being fitted. This is
essential for motifs whose 200 bp density cap raised the production scan floor.

Enrichment Slurm jobs stage the shared autosomal anchor table once for a small
motif batch. Each motif's compact results are still independently validated and
published, so requeued jobs skip durable work. The finalizer performs one
Benjamini-Hochberg correction across the complete planned non-TP73 motif family.

These are whole-autosome exploratory associations. The convenient operating
points were selected on chromosome 1, so more bases improve precision but do
not constitute independent threshold validation. Spatial interpretation should
use the context-band layer; the zero-complete 150 bp table answers the separate
threshold/count question.

## Production sequence

Use dedicated paths below `/data/sm718` and a tracked-clean checkout. Render
each command with `--dry-run` first.

```sh
scripts/submit_tp73_genome_evidence_slurm.sh \
  --run-root "$EVIDENCE_RUN" \
  --annotation-run "$ANNOTATION_RUN" \
  --scan-package "$SCAN_RUN/package" \
  --track-root "$CUTANDRUN_ROOT" \
  --runtime-prefix "$RUNTIME" \
  --max-concurrent 20 --cpus 4 --memory 32G --time 02:00:00 \
  --dry-run
```

After its finalizer succeeds:

```sh
scripts/submit_tp73_genome_context_maxima_slurm.sh \
  --run-root "$CONTEXT_RUN" \
  --scan-package "$SCAN_RUN/package" \
  --evidence-package "$EVIDENCE_RUN/final/genome_evidence" \
  --threshold-registry "$CHR1_CONTEXT_THRESHOLDS" \
  --runtime-prefix "$RUNTIME" \
  --motifs-per-batch 3 --max-concurrent 20 \
  --cpus 4 --memory 32G --time 02:00:00 \
  --dry-run
```

After the context finalizer succeeds, pass that context run as
`--source-threshold-run` and the autosomal evidence Parquet as
`--anchor-evidence` to `submit_tp73_cofactor_enrichment_slurm.sh`. Its default
two-motif batches request 2 CPUs, 32 GB, and one hour.

Stages 4 and 5 rebuild the shared nuclear annotation at schema 9, extract
H3K4me3/input signal in five fixed windows by reusing the finalized TP73 anchor
evidence, and fit the all-JASPAR GFP-referenced change models. Their exact
inputs, Slurm commands, restart contract, genomic-context covariates, and
global multiple-testing family are specified in
[`h3k4me3_whole_genome_production.md`](h3k4me3_whole_genome_production.md).

## Local contract tests

`tests/test_tp73_genome_evidence_manager.sh` builds a 25-region fixture,
including scan label `25` matched to annotation label `MT`.
`tests/test_tp73_genome_context_maxima.sh` builds 22 tiny autosomes and checks a
real nearby hit, completely empty motif partitions, threshold raising, restart
reuse, final catalog construction, and enrichment-plan handoff.
`tests/test_tp73_cofactor_enrichment.sh` verifies multi-chromosome adjustment
and proves that a negative reference below the retained scan floor is not
treated as observed evidence.

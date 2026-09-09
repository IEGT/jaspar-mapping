# Ensembl regulatory context at TP73 anchors

Status, 2026-09-09: importer, anchor membership, TP73 evidence selection and
H3K4me3 subset refitting implemented and tested locally. The full independent
coordinate audit passes. Restart-safe annotation has completed on Haumea for
all 22 autosomes. Both statistical Slurm managers now pin and validate the
regulatory cohort, including subset-aware checkpoints and finalization.

## Scientific contract

The whole-genome anchor cohort remains the scientific primary. Ensembl
regulatory annotations are additional, independently selectable contexts,
not replacements for the existing gene-relation annotation:

- Extended Ensembl promoters are primary **among promoter-focused screens**.
- Core-only promoters are a nested sensitivity; the existing strand-aware
  TSS upstream-2000/downstream-500 definition is a parallel comparator, not
  an intersection with Ensembl.
- Enhancers, open chromatin, other feature types and unannotated anchors
  remain available. No regulatory overlap is not synonymous with intergenic.
- Ensembl reference-tissue annotation is not demonstrated activity in SaOS-2
  or SK-Mel-29. Restriction changes the estimand; GFP-baseline dependence,
  GC/CpG, mappability and accessibility remain sensitivity questions. A
  GFP-referenced change does not make this selection automatically harmless.

Select anchors by positive half-open interval overlap, excluding abutment.
Leave their cofactor neighbourhoods and existing score thresholds unchanged.
A cofactor may lie outside the promoter containing its TP73 anchor. The
"both motifs in the same promoter" sensitivity is deliberately **not**
implemented here: the production maxima do not retain neighbour coordinates.
That sensitivity must revisit the retained low-floor positions, not infer
membership from the maximum score alone. Neither path needs a new DNA scan.

Regulatory memberships are overlapping Boolean properties, not an exclusive
classification. Core membership implies extended membership; enhancer and
promoter membership may coexist. Frequency is the number of unique positive
anchors divided by unique eligible anchors in the chosen subset, including
intermediate cofactor scores in the denominator. Frequencies from exclusive
distance bands must **not** be summed: their positive-anchor sets can overlap.
Per-anchor locus counts can be summed across bands; the all-150 maximum is
the maximum of band maxima.

## Verified Ensembl source

The [2025-12 track hub](https://regulation.ensembl.org/2025-12/trackhub/hub.txt)
lists the human GRCh38 assembly. Its
[track documentation](https://regulation.ensembl.org/2025-12/trackhub/docs/trackhub.html)
links both BigBed and GFF. The documentation distinguishes the TSS
identification window (90 upstream/10 downstream) from promoter-core
construction (490 upstream/10 downstream), and describes merging overlapping
cores and extending promoters using open-chromatin evidence. Do not
reconstruct these features from our own TSS windows.

The public GFF URL requested on 2026-09-08 was:

```
https://regulation.ensembl.org/api/annotation/v0.15/files/download/2025-12/homo_sapiens/GRCh38/Homo_sapiens.GRCh38.regulatory_features.gff.gz
```

The redirect headers explicitly resolve that request to **2025-05**, not
2025-12. Both dates are recorded: `requested_regulatory_release=2025-12`,
`regulatory_release=2025-05`. The durable public URL and actual compressed-file
SHA-256 identify the source; short-lived signed redirect URLs are not stored
as stable provenance. For future downloads, inspect the resolution headers
again, rather than assuming this alias continues to resolve identically.

| Export | SHA-256 |
|---|---|
| GFF.gz | `a5f5ef58ee7b3dfbc3667692d3cc6515c66789ad2de2c3a15784c5436367bb32` |
| BigBed | `0bc0285ef0c2d5f0e2ca5f398bf8015a20747c9e54e80a97f6a168ef6f13f389` |

The GFF contains 643,528 features: 35,983 promoters, 246,403 enhancers,
7,541 open-chromatin regions, 90,891 CTCF binding sites and 262,710 EMARs.
All promoters have `extended_start` and `extended_end`; the GFF interval
itself is the core. Gene IDs can be comma-separated. Every feature type and
the raw attributes are preserved, including records with no gene link.

### Coordinate audit: resolved

The initial rtracklayer check appeared to disagree for `ENSR1_958`:

| Interpretation, BED half-open | Core | Extended |
|---|---|---|
| GFF, standard `start - 1`, unchanged end | `[10935,11436)` | `[9952,11436)` |
| BigBed read by rtracklayer, GRanges start converted back to BED | `[10934,11436)` | `[9951,11436)` |
| BigBed read independently with UCSC bigBedToBed | `[10935,11436)` | `[9952,11436)` |

The GFF raw row is `10936..11436`, with `extended_start=9953` and
`extended_end=11436`. The independent reader establishes that the discrepancy
was in the earlier R import/conversion path, **not between Ensembl's exports**.
All 643,528 feature IDs, chromosomes, strands, types, core and extended bounds
agree with the GFF conversion. No source intervals need to be shifted again.
The exact payload/reader hashes and audit result are committed in
[`regulatory_GRCh38_2025-05_coordinate_audit.json`](../resources/ensembl/regulatory_GRCh38_2025-05_coordinate_audit.json).

Reproduce the whole-file streaming comparison using the official UCSC reader
for the local platform, without writing an intermediate BED file:

```bash
python3 scripts/build_regulatory_annotation.py audit \
  --gff SOURCE/Homo_sapiens.GRCh38.regulatory_features.gff.gz \
  --bigbed SOURCE/Homo_sapiens.GRCh38.regulatory_features.bb \
  --bigbed-to-bed TOOLS/bigBedToBed --assembly GRCh38 --output RUN/coordinate-audit
```

This is not a choice of regional coordinate conventions: Ensembl BED/BigBed
already use zero-based, half-open intervals. GENtle's genomic-region records
use `start_0based`, `end_0based_exclusive` and an explicit coordinate-convention
tag. Its genome-track projection adds one to the BED start, not the end, when
interfacing with a one-based genome anchor. Keep these boundaries explicit:
BED `[9952,11436)` displays as inclusive `9953-11436`; a one-base interval
`[100,101)` displays as `101-101`. No GENtle code changes were needed here.

An independent BigBed check, where Bioconductor rtracklayer is available:

```r
library(rtracklayer)
g <- import(BigBedFile("Homo_sapiens.GRCh38.regulatory_features.bb"),
            which=GRanges("1", IRanges(1,12000)),
            colnames=c("name", "thick", "featureType"))
d <- as.data.frame(g)
d <- d[d$name == "ENSR1_958", ]
data.frame(core_start=d$thick.start-1, core_end=d$thick.end,
           extended_start=d$start-1, extended_end=d$end)
```

Selecting columns avoids the export's unrelated AutoSQL `itemRgb` issue:
it declares a number while carrying comma-separated RGB values.

## Separate packages and keys

`scripts/build_regulatory_annotation.py` needs only Python's standard library
and the DuckDB CLI. Use its `import`, `annotate` and `select` subcommand help
for all options. Outputs are immutable directories: existing destinations are
refused, intermediates stay private, and only completed packages are promoted.
Source data are never removed. DuckDB defaults to two threads and a 1 GB
memory limit. Large generated packages remain outside Git.

The imported dimension is shared across chromosome membership packages,
rather than copied for every cofactor. It is sorted by chromosome and start
in ZSTD Parquet; annotation output is one package per chromosome. Every
package contains a file-size/SHA-256 inventory, builder hash and `schema.sql`.
Open DuckDB **from that package directory** and read `schema.sql` to create
the corresponding views. Views use exact file paths, not globbed input sets.

| Table | Grain and meaning |
|---|---|
| `regulatory_feature` | One physical Ensembl feature. ID hashes provider, assembly, resolved release, payload SHA-256 and native ID. Stores original interval, nullable core/extended bounds, strand, type and original attributes. |
| `regulatory_feature_gene` | Native gene ownership or a separately labelled TSS-overlap-derived link. `link_source` distinguishes native, core, extended and feature-extent links; annotation releases remain separate. |
| `regulatory_feature_tss` | Physical feature + selected bound definition + pinned physical TSS. Several genes and opposite-strand TSSs can share a feature. An enhancer overlapping a TSS is not thereby a proven target-gene link. |
| `tp73_anchor_regulatory_feature` | Physical anchor + feature + bound definition; positive overlap bases and a full-containment flag. Does not multiply the analysis cohort by genes. |
| `tp73_anchor_regulatory_membership` | Exactly one row per input physical anchor, including anchors with no feature. Boolean subset flags; unavailable legacy TSS-promoter annotation is NULL, not false. |

Missing extended bounds remain NULL in the dimension. Annotation refuses to
declare extended-promoter absence when promoter extension coordinates are
missing. An absent chromosome (including a `chr1`/`1` mismatch) is an error,
not an all-negative annotation. When supplied, TSS/ownership/promoter files
must match the explicitly pinned genome and GTF annotation release.

The existing `promoter`, `promoter_gene`, Q20--Q24 and gene-relation
precedence are unchanged. Regulatory feature IDs have no single owning TSS.

## Local commands

Download the pinned GFF into a dedicated source directory, retaining the
redirect headers. Import using its verified digest:

```bash
python3 scripts/build_regulatory_annotation.py import \
  --gff SOURCE/Homo_sapiens.GRCh38.regulatory_features.gff.gz \
  --assembly GRCh38 --requested-release 2025-12 --resolved-release 2025-05 \
  --source-uri https://regulation.ensembl.org/api/annotation/v0.15/files/download/2025-12/homo_sapiens/GRCh38/Homo_sapiens.GRCh38.regulatory_features.gff.gz \
  --expected-sha256 a5f5ef58ee7b3dfbc3667692d3cc6515c66789ad2de2c3a15784c5436367bb32 \
  --coordinate-audit resources/ensembl/regulatory_GRCh38_2025-05_coordinate_audit.json \
  --coordinate-audit-note 'All 643528 features agree with independently read BigBed' \
  --output RUN/regulatory_features

python3 scripts/build_regulatory_annotation.py annotate \
  --features RUN/regulatory_features --anchors TP73_EVIDENCE.parquet \
  --assembly GRCh38 --chrom 1 --output RUN/chrom-1
```

For gene links and the legacy-promoter comparator also supply `--tss`,
`--transcript-tss`, `--genome-id`, `--annotation-release`, `--promoters` and
`--promoter-definition-id`. These take the existing context-package tables,
not a newly inferred gene annotation.

The `select` command exports evidence rows, preserving their other columns:

```bash
python3 scripts/build_regulatory_annotation.py select \
  --anchors TP73_EVIDENCE.parquet --membership RUN/chrom-1 \
  --subset promoter_extended --output RUN/chrom-1-extended-evidence
```

Pass its `selected_anchors.parquet` to the existing TP73 distance-count
kernel, alongside the **unchanged** cofactor-position files. It recomputes
class frequencies and per-series/isoform block components on that cohort.
The annotation's n:m bridges are never directly joined into the model rows.

For H3K4me3, keep the original evidence, change and zero-complete maxima
inputs and add these options to the existing evaluator invocation:

```bash
--regulatory-membership RUN/chrom-1/tp73_anchor_regulatory_membership.parquet \
--regulatory-subset promoter_extended
```

Membership input is repeatable for multiple chromosomes and must cover the
exact physical-anchor universe of the signal input. Duplicate/missing/extra
keys, NULL flags and an empty selected subset fail before fitting. Each
result table carries `regulatory_subset`; run-config schema 7 records input
and selected anchor counts, exact sidecar paths and MD5 file identities.
The membership package additionally supplies SHA-256 provenance. Default
unfiltered runs retain schema 6 and their original statistical behaviour.

Every selected cohort is refitted, preserving TP73-score adjustment,
series-specific TA/DN contrasts, block-clustered uncertainty and per-family
BH adjustment. The existing all-zero H3K4me3 exclusion remains specific to
the intensity model; it is not an occupancy-selection rule. Intermediate
cofactor scores remain in frequency denominators. An entirely empty subset
fails explicitly; nonempty but underpowered contrasts keep existing
not-estimable statuses rather than being reported as no effect.

## Validation and pilot

```bash
python3 tests/test_regulatory_annotation.py
bash tests/test_tp73_distance_cofactor_counts.sh
bash tests/test_h3k4me3_cofactor_change.sh
bash tests/test_script_help.sh
```

Tests cover coordinates, nested bounds, native versus derived gene ownership,
bidirectional/shared promoters, unlinked enhancers, abutment, chromosome
mismatch, corrupted packages, immutable output and unique-anchor denominators.
The TP73 fixture retains an adjacent cofactor **outside** the selected promoter
and recomputes the subset TA odds ratio (3 instead of the original 6).
The H3K4me3 sidecar refit exactly matches independently filtered input data,
including effect, isoform-contrast and frequency-summary tables.

The real local pilot used the older
`dry_runs/h3k4me3_cofactor_change_chr1_20260809/tp73_anchor_evidence.parquet`,
which contains **310,782 chr1 anchors with score >= 0**. It is not the current
genome-wide, low-floor local-peak production cohort. The annotation
found 2,415 core-promoter anchors, 3,559 extended-promoter anchors, 20,068
enhancer anchors, 263 open-chromatin anchors and 282,712 anchors with no
regulatory overlap. These overlapping categories must not be added together.

Both inputs are restricted to the same explicit chromosome before range
joins, allowing DuckDB's inequality join rather than a chromosome hash bucket.
The optimized pilot took about 1.1 seconds; its full membership bridge was
identical to the slower equality-plus-range join. Pilot packages are in
`dry_runs/ensembl_regulation_20260908/`; they are validation artifacts, not new
enrichment or H3K4me3 results.

## Restart-safe chromosome production

`manage_regulatory_annotation.py prepare` pins exact chromosome evidence and
TSS/ownership/promoter files from completed catalogs, the clean source commit,
scientific source hashes and the passed coordinate audit. `setup` downloads
the small public GFF on a compute node, checks its digest and builds the
shared feature package once. No data are copied from the laptop to Haumea.

Each `run-task` stages a single chromosome's inputs and the shared annotation
on `/scratch`, verifies the staged bytes, annotates, then copies the small
result to a unique durable attempt directory before atomic promotion. Requeue
restages scratch and validates/reuses completed chromosomes; it never replaces
them. Haumea owns scratch cleanup. `finalize` requires every planned chromosome
and produces the combined, zero-complete membership Parquet and query schema.

The submission helper uses `requeue`, 2 CPUs, 8 GB and a 20-minute limit per
job, with 10 concurrent chromosome tasks. This is annotation, not a genome
scan; it has no FASTA dependency. SIGUSR1 reports worker phase and elapsed time.

```bash
SOURCE=/data/sm718/GitHub/jaspar-mapping
RUN=/data/sm718/jaspar_mapping_runs/ensembl_grch38_tp73_regulatory_20260909_v1
RUNTIME=/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_context_thresholds_v1/runtime
python3 "$SOURCE/scripts/manage_regulatory_annotation.py" --help
bash "$SOURCE/scripts/submit_regulatory_annotation_slurm.sh" \
  --source "$SOURCE" --run-root "$RUN" \
  --evidence-package /data/sm718/jaspar_mapping_runs/jaspar2026_grch38_tp73_cutandrun_evidence_v1/final/genome_evidence \
  --annotation-catalog /data/sm718/jaspar_mapping_runs/jaspar2026_grch38_tp73_annotation_v4_schema9/final \
  --duckdb "$RUNTIME/duckdb/bin/duckdb" --partition requeue --max-concurrent 10
```

Use an immutable source checkout for a submitted run. The optional `--dry-run`
creates a pinned plan and prints commands only; use another new run directory
for a subsequent fresh submission. Submitted IDs are recorded immediately in
`submissions.tsv` so a partial scheduler submission is diagnosable.

### Completed Haumea annotation, 2026-09-09

The production run used immutable source commit `940e54d`, fetched from the
public repository, and the audited GFF downloaded directly on Haumea:

```
/data/sm718/jaspar_mapping_runs/ensembl_grch38_tp73_regulatory_20260909_v1
```

Setup `5782607`, chromosome array `5782608` (22 tasks), and finalizer
`5782609` all completed with exit code `0:0`. Slurm recorded setup at
08:37:29 and finalization at 08:38:59 (cluster time), including dependency
and scheduling waits. Setup took 16 seconds; each chromosome task and the
finalizer took about one second. Maximum recorded job RSS was 1,191,336 KiB,
below the 8 GB allocation. Scratch staging was confirmed in the task log.

The combined `final/tp73_anchor_regulatory_membership.parquet` contains
**3,596,429 unique physical anchors across chromosomes 1-22**, occupying
15,481,758 bytes. Its SHA-256 is
`e8a8db032478a1af0a47348893bcbbae14d33657a7dd7c670485a54b8c43521c`.
The finalizer validated the chromosome packages, unique anchor keys and
core-to-extended nesting before publication. A bounded aggregate query on
the final Parquet produced:

| Anchor membership | Anchors |
|---|---:|
| Ensembl promoter core | 23,698 |
| Ensembl extended promoter | 34,035 |
| Ensembl enhancer | 190,540 |
| Ensembl open chromatin | 2,651 |
| No overlap with any imported regulatory feature | 3,333,743 |
| Legacy TSS upstream-2000/downstream-500 promoter | 436,537 |

These are overlapping memberships, not additive classes or enrichment
estimates. The no-overlap flag also considers CTCF sites and EMARs, whose
individual counts are not displayed above. This production cohort differs
from the older local pilot. Sex chromosomes and mitochondria were not part
of this autosomal run. No new motif scan or cofactor model fit was performed.

## Promoter-subset cofactor refits

Both statistical submission helpers accept:

```bash
--regulatory-package /data/sm718/jaspar_mapping_runs/ensembl_grch38_tp73_regulatory_20260909_v1/final \
--regulatory-subset promoter_extended
```

Use a **new immutable run directory for each subset and adjustment variant**.
The predeclared comparison is Ensembl `promoter_extended` (primary among
promoter-focused analyses), nested `promoter_core`, and parallel
`tss_window_promoter`. The existing unrestricted runs remain the whole-genome
comparators. Keep the positive cofactor score at 0 for compatibility with the
completed score-zero analyses; do not retune it within each promoter subset.
TP73 distance inference retains its existing strict `< -1` negative reference
and six exclusive bands. H3K4me3 retains both `< -1` and `< 0` references,
the six bands plus `all_150`, and the fixed `flank_150_1000` signal window.
Run H3K4me3 both with and without `--adjust-gfp-baseline` as separate packages.

The managers pin the regulatory manifest and payload SHA-256, audit, subset,
GTF release and legacy promoter definition. Preflight verifies exact physical
anchor coverage, Boolean flags, promoter-core nesting and a nonempty cohort.
TP73 workers consume newly selected chromosome evidence but the unchanged
low-floor cofactor positions. H3K4me3 workers stage the membership file with
the other fixed inputs; its evaluator validates the full cohort before
selecting anchors and refitting. No precomputed whole-genome odds ratio or
H3K4me3 coefficient is filtered after fitting.

Subset identity is checked on anchor splits, chromosome/motif checkpoints and
final packages. H3K4me3 validates evaluator schema 7 for subset results and
schema 6 for unrestricted results, including the selected/input anchor counts
from preflight. Its finalizer retains subset labels and recomputes all-motif BH
families within each immutable run. TP73 final schema 6 adds `regulatory_subset`
to result/frequency/contrast tables without changing the estimator. Completed
whole-cohort results cannot be reused as subset checkpoints.

Initial scheduling uses the `requeue` partition, chromosome-local scratch for
TP73 and batch-local scratch for H3K4me3. Start the extended-promoter runs at
modest concurrency, inspect a completed motif and memory use, then submit the
two comparators. With eight TP73 tasks and six H3K4me3 batches per adjustment
variant for each of three subsets, the combined concurrency ceiling is 60.
Use four sequential H3K4me3 motif checkpoints per batch, 64 GB per H3K4me3 job,
32 GB per TP73 job, and two-hour task limits initially. These are scheduling
choices, not part of the scientific cohort definition.

The manager tests exercise full prepare/preflight/run/reuse/finalize cycles,
reject a whole-cohort checkpoint presented as a promoter-subset result, and
reject changed membership bytes. In the TP73 fixture the positive frequency
changes from 4/8 to 3/7 after anchor selection; nearby cofactor positions remain
unchanged. The H3K4me3 fixture refits 198 of the original 264 anchors and checks
that the finalized tables and portable provenance retain that selection.

After completion, report TA and DN beside one another with frequency, support,
uncertainty and matched reference definitions. Differences between subset
estimates require a separate interaction test, not a comparison of significance
labels. Reference-tissue promoter membership does not demonstrate activity in
the experimental cell lines, and neither change model establishes causality.

No new cluster jobs or source synchronization are implied by the local pilot.

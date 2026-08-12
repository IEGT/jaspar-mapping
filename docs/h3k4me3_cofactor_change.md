# GFP-Referenced H3K4me3 Change at TP73 Motif Anchors

The completed schema-7 chromosome-1 production result is reported in
[`h3k4me3_chr1_production_20260810.md`](h3k4me3_chr1_production_20260810.md).

## Question and estimand

This analysis asks whether the sequence-defined context of a TP73 motif anchor
is associated with an additional H3K4me3 response after expression of TAp73 or
DNp73, relative to the matching GFP condition. H3K4me3 is an
activation-associated chromatin mark; it is not interpreted as transcription
factor occupancy.

For anchor `i`, experimental series `r`, and condition `c`, the windowed mark is

```text
H(i,r,c) = log2((H3K4me3_area(i,r,c) + alpha) /
                (input_area(i,r,c) + alpha))
```

and the two condition changes are

```text
delta_TA = H(TA) - H(GFP)
delta_DN = H(DN) - H(GFP).
```

The primary per-series estimand for a cofactor motif is

```text
E(delta | cofactor score >= calibrated threshold) -
E(delta | cofactor absent or score < negative reference).
```

This is a difference-in-differences association: static anchor signal and the
static main contribution of the cofactor cancel from the condition difference.
It is not a causal cofactor effect because sequence features can still be
associated with differential chromatin responsiveness.

## Included series

The source-controlled manifest is
[`config/h3k4me3_cutandrun_tracks_v1.tsv`](../config/h3k4me3_cutandrun_tracks_v1.tsv).
The primary analysis includes exactly:

- `saos2` (`SaOS-2`); and
- `skmel29_2` (`SK-Mel-29`).

`skmel29_1` is retained in the manifest but has
`analysis_included=false` and the investigator-specified exclusion reason that
the dataset is unresolved. The builder neither requires nor opens excluded
files. Filename discovery and wildcard inclusion are deliberately not used.
SaOS-2 and SK-Mel-29 are separate biological systems, not replicates. Results
are estimated separately and directional agreement is descriptive support;
there is no claim of three-cell-line replication or a pooled replication
P-value.

Every included series must provide the complete `GFP`, `TA`, and `DN` factorial
for four channels:

- `h3k4me3`: files historically named `pos_*`;
- `input`: the condition-matched input;
- `tp73`: anti-p73 CUT&RUN; and
- `negative_control`: files historically named `neg_*`.

The manifest currently marks the processed BigWig scale as pending provenance.
Exact normalization and clipping metadata must be resolved before publication.
The run sidecar records source size, modification time, and SHA-256.

## Window selection

H3K4me3 can flank a nucleosome-depleted regulatory centre, so the 150 bp motif
context is not automatically an appropriate signal window. Window selection is
blinded to the outcome:

1. use only `GFP` H3K4me3 and input tracks;
2. draw an aggregate profile from -2 kb to +2 kb around a deterministic,
   genome-ordered sample of TP73 motif-alignment centres;
3. inspect both included series; and
4. predeclare one shared primary window before examining TA/DN or cofactor
   effects.

Candidate durable windows are the motif span, central +/-150, +/-500, +/-1000
bp, and the two flanks from 150 to 1000 bp. Narrower windows remain sensitivity
analyses. A common window is preferred; a serious disagreement between the two
GFP profiles is reported rather than hidden by series-specific optimization.

The blinded 100,000-anchor chromosome-1 calibration selected
`flank_150_1000` as the primary window. Both accepted series showed a shallow,
broad central depletion and nearly identical relative profiles. The evidence
and decision are recorded in
[`h3k4me3_gfp_window_calibration_chr1_20260809.md`](h3k4me3_gfp_window_calibration_chr1_20260809.md).
The subsequent nine-motif real-data pilot is reported separately in
[`h3k4me3_chr1_cofactor_pilot_20260809.md`](h3k4me3_chr1_cofactor_pilot_20260809.md).

The profile's per-base pseudocount is only a display aid. The analysis uses an
integrated-signal pseudocount, default `alpha=1`, and records it explicitly.

## Signal and TP73 evidence

[`build_h3k4me3_anchor_signal.py`](../scripts/build_h3k4me3_anchor_signal.py)
uses BigWig as the source contract. BigWig-to-bedGraph conversion occurs only
inside a temporary directory and is deleted when the process exits; the
persistent data layer is Parquet, not BED or TSV.

The two durable outputs are:

- TP73 anchor evidence: one row per orientation-collapsed motif-alignment span,
  carrying strict TP73 and negative-control immersion/depth for every included
  series and condition; and
- H3K4me3 anchor signal: one row per anchor, series, condition, and named window,
  carrying integrated area, full-window mean, maximum, positive-signal bases,
  and positive-bp fraction for H3K4me3 and input.

Strict TP73 support retains the established rule: one merged positive-coverage
component must start before and end after the complete motif-alignment span.
H3K4me3 is summarized over a window and is not subjected to strict immersion.

The chromosome-1 production build passed its then-current schema-7
`tp73_context_anchor` Parquet with `--anchor-source`. Those rows had already
passed the declared TP73 local-peak selection. Tied orientation records at one
alignment span were collapsed to one physical anchor using their maximum score,
without recomputing source selection.

The whole-genome signal build instead passes each completed chromosome's
TP73/control evidence with `--tp73-evidence-input`. This signal-only mode keeps
the exact physical anchor set and reads only H3K4me3/input tracks; it does not
repeat the strict-immersion calculation. `--anchor-source` and
`--anchor-plus`/`--anchor-minus` remain available for older full and score-floor
pilots, and all source modes are mutually exclusive.

## TP73 binding state

Condition-specific confirmed TP73 support is

```text
strict TP73 immersion AND NOT strict negative-control immersion.
```

For each isoform comparison, anchors are also described as `gained`, `shared`,
`lost`, or `unsupported` relative to GFP. TP73 confirmation is post-treatment
and can mediate a cofactor-associated H3K4me3 response. It is therefore omitted
from the primary total-association model. The secondary model

```text
delta ~ spline(TP73 motif score) + cofactor_present * TP73_confirmed
```

answers the requested effect-modification question descriptively. It is emitted
only when all four cofactor-positive/negative by TP73-confirmed/unconfirmed
cells satisfy the declared minimum support.

## Cofactor and uncertainty contract

- Positive anchors use the motif-specific calibrated threshold.
- The primary negative class is absent or score `< -1`.
- The historical comparison uses absent or score `< 0` only when it does not
  overlap the positive class.
- Each negative reference is evaluated independently. A reference below the
  cofactor scan's retained score floor is emitted as
  `negative_reference_below_source_floor`; an observable reference in the same
  run remains analyzable. Comparisons whose positive and negative score classes
  overlap are marked invalid rather than summarized as disjoint.
- Intermediate anchors are excluded from that comparison, not relabelled.
- Anchors with H3K4me3 area zero in GFP, TA, and DN are omitted from intensity
  and continuous-score fits because they carry no mark-intensity information.
  They remain in the secondary gained/lost occurrence summaries.
- Linear-model uncertainty is clustered by 5 Mb genomic block.
- Models are fit separately for each included series and for TA and DN.
- Benjamini-Hochberg correction is applied separately by series, isoform, and
  negative reference. TP73-interaction contrasts form separate families.

With schema-9 annotation, the primary adjustment includes a TP73-score spline,
chromosome, compact genomic context, unsigned genomic distances to the nearest
TSS and CDS, and explicit upstream/downstream/overlap/mixed-strand direction
classes. The evaluator also writes context-stratified intensity
effects (including strict intergenic) and a continuous cofactor-score
sensitivity. The deterministic nearest-feature summaries are adjustment
covariates; the normalized tied-nearest tables remain authoritative for
biological interpretation.

The evaluator additionally writes fixed four-way gene-relation results using
`promoter > downstream > gene_body > intergenic` precedence. Promoter and
downstream membership come from versioned many-to-many physical interval
bridges. The downstream interval mirrors the promoter geometry around the
transcript end and follows transcript strand: its primary definition spans
500 bp toward the transcript body and 2,000 bp beyond the end. An anchor that
overlaps both a promoter and a downstream region is classified as promoter,
while both independent overlap flags and bridge rows remain available.
`gene_body` means transcript overlap outside those higher-precedence regions;
it is not synonymous with CDS. Within that pooled stratum the model continues
to adjust for the finer CDS/exonic/intronic context.

Two relation-stratified outputs are emitted. The H3K4me3 table has 32 rows per
motif: two negative references by two isoforms by two series by four relation
classes. The matched TP73-occupancy table has eight rows per motif: two
negative references by four relation classes. It refits the established
discordant anti-TP73 versus matched-control logistic model inside each class,
with sample and chromosome fixed effects, a TP73-score spline, and 5 Mb
block-clustered uncertainty. Consequently, each relation-specific figure uses
the same anchor stratum on both axes. The occupancy fit retains anchors with
zero H3K4me3 because its outcome is independent of the mark; only the H3K4me3
fit applies the all-mark-zero exclusion. Underpowered classes remain explicit
rows with a non-`ok` status rather than borrowing the genome-wide estimate.

The evaluator writes intensity effects, TP73 interactions, binding-state
summaries, occurrence summaries, cross-series directional summaries, and the
complete run configuration. PATZ1, TFAP2C, E2F1, and SP1 are sentinel biological
motifs, not null controls. Density-matched motif labels or block-preserving
permutations are required for a null validation.

The reproducible chromosome-1 pilot panel and its positive thresholds are in
[`config/h3k4me3_chr1_pilot_cofactors_v1.tsv`](../config/h3k4me3_chr1_pilot_cofactors_v1.tsv).
It covers E2F1, SP1, REST, POU2F2, KLF14, TCF7, POU4F1, TFAP2C, and PATZ1 and
points back to the immutable chromosome-1 TP73-context threshold set.

## Local GFP preflight

The following uses the checked local chromosome-1 TP73 scan and the processed
CUT&RUN resource directory. It reads only the four GFP H3K4me3/input BigWigs and
does not touch `skmel29_1`.

```bash
RUN=dry_runs/h3k4me3_cofactor_change_chr1_20260809
mkdir -p "$RUN"

scripts/build_h3k4me3_anchor_signal.py \
  --anchor-plus dry_runs/tp73_chr1_direct_sparse_parquet_benchmark_20260717/threshold_zero/tables/jaspar2026/motif_hit/motif_id=MA0861.2/score_mode=log2_relative_risk/pseudocount=1/minimum_score=0/minimum_pwm_relative_score=none/maximum_pwm_relative_score=none/chrom=1/strand=plus/n_policy=skip/matched_sequence=omitted/part-from=0-to=end-000000.parquet \
  --anchor-minus dry_runs/tp73_chr1_direct_sparse_parquet_benchmark_20260717/threshold_zero/tables/jaspar2026/motif_hit/motif_id=MA0861.2/score_mode=log2_relative_risk/pseudocount=1/minimum_score=0/minimum_pwm_relative_score=none/maximum_pwm_relative_score=none/chrom=1/strand=minus/n_policy=skip/matched_sequence=omitted/part-from=0-to=end-000000.parquet \
  --track-manifest config/h3k4me3_cutandrun_tracks_v1.tsv \
  --track-root ../gentle_rs/data/resources/cutandrun_20250602_noDuplicates \
  --profile-only --profile-output "$RUN/gfp_metaprofile.tsv" \
  --chrom 1 --chrom-length 248956422 --minimum-anchor-score 0 \
  --profile-flank 2000 --profile-bin-size 50 --profile-sample-size 20000 \
  --tmpdir /private/tmp

scripts/plot_h3k4me3_metaprofile.R \
  --input "$RUN/gfp_metaprofile.tsv" \
  --output "$RUN/gfp_metaprofile.png" \
  --table-output "$RUN/gfp_metaprofile_ratio.tsv"
```

After predeclaring the window, omit `--profile-only` and add
`--tp73-output`, `--signal-output`, and any explicit `--window` definitions.
The nine-motif local context file can then exercise the complete evaluator
before a per-motif all-JASPAR Slurm run is designed.

Full builds report separate timed phases for strict TP73/control evidence and
H3K4me3/input window aggregation. Outputs are promoted from temporary scratch
only after the corresponding DuckDB validation succeeds.

## Chromosome-1 production run

The production driver resolves the chromosome-1 `tp73_context_anchor` Parquet
from the finalized annotation catalog rather than from a filename pattern. The
partition's chromosome is supplied explicitly as `--chrom 1`; it is not
assumed to be redundantly stored inside the Parquet payload. The selected
schema-7 input contains 305,528 strand-aware records representing 305,492
physical motif-alignment spans, all marked `anchor_selection_class =
'local_peak'` and all at score `>= -1`.

The same job builds the strict TP73/negative-control evidence, the complete
H3K4me3/input signal factorial, and the nine-cofactor maxima before fitting the
predeclared `flank_150_1000` analysis. Motif Parquets are selected through the
scan catalog and copied to local scratch. BigWig conversion is scratch-only.
`skmel29_1` remains excluded by the source-controlled track manifest.

```sh
SOURCE=/data/sm718/GitHub/jaspar-mapping
ANNOTATION=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_tp73_annotation_v2
SCAN=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_sparse_v3/package
TRACKS=$SOURCE/cutandrun_20250602_noDuplicates
RUNTIME=/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_context_thresholds_v1/runtime
RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_h3k4me3_production_v4

"$SOURCE/scripts/submit_h3k4me3_chr1_production_slurm.sh" \
  --run-root "$RUN" --annotation-run "$ANNOTATION" \
  --scan-package "$SCAN" --track-root "$TRACKS" \
  --runtime-prefix "$RUNTIME" --partition requeue \
  --cpus 4 --memory 32G --time 08:00:00
```

The worker handles `SIGUSR1` without terminating and reports the current phase,
elapsed time, durable bytes, and scratch bytes. Attempts are immutable. Only a
validated attempt is atomically promoted to `$RUN/final`, whose `manifest.json`
pins the annotation catalog, anchor file, scan inventory, track manifest,
cofactor thresholds, source commit, output checksums, and cardinality checks.
Scratch-to-`/data` publication copies first to a hidden file beside the durable
target and renames that file atomically, because a direct cross-filesystem
`rename(2)` is not valid.

## Whole-genome production

The restart-safe whole-genome measurement, schema-9 annotation dependency, and
all-JASPAR inference commands are specified in
[`h3k4me3_whole_genome_production.md`](h3k4me3_whole_genome_production.md).
Autosomes are primary; X/Y are retained as sensitivity data; mitochondria are
excluded. The evaluator accepts repeated exact `--change`, `--annotation`, and
`--cofactor-maxima` files, so finalized inventories are consumed without a
wildcard or copied aggregate.

One-motif checkpoints necessarily have a one-motif local BH family. Their
`q_value_bh` is diagnostic. The whole-genome finalizer retains it as
`q_value_bh_task` and computes `q_value_bh_all_motifs` across every planned
non-TP73 motif within each declared series/isoform/reference/result family.

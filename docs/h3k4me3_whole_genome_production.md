# Whole-Genome H3K4me3 Change Production

This run extends the chromosome-1 GFP-referenced H3K4me3 analysis to the 24
GRCh38 nuclear chromosomes. Autosomes 1-22 form the primary inference set. X
and Y are retained as a separate sensitivity set. Mitochondrial labels (`M`,
`MT`, and `25`, with or without a `chr` prefix) are excluded because
mitochondrial abundance can produce a CUT&RUN bystander artifact.

The workflow has three independently restartable products. Schema-9 annotation
and H3K4me3 signal measurement can run in parallel. All-motif inference starts
only after both have finalized.

## Immutable inputs

- completed TP73/control evidence package;
- completed all-autosome cofactor-maxima package and applied threshold registry;
- finalized JASPAR 2026 GRCh38 scan package;
- Ensembl release-113 GTF;
- source-controlled H3K4me3/input track manifest; and
- a tracked-clean source commit.

No wildcard discovers scientific payloads. Each manager resolves exact files
from a finalized inventory, pins their size and checksum, and refuses changed
inputs. Durable output is below `/data/sm718`; temporary copies and DuckDB spill
are below the job's `/scratch` directory.

The H3K4me3 signal plan hashes every included H3K4me3/input BigWig once and
stores path, size, nanosecond mtime, and SHA-256 in an immutable track-file
inventory. Chromosome workers validate the pinned metadata and reuse that
checksum in their provenance sidecars. They therefore perform chromosome-local
BigWig reads without redundantly hashing every complete track in every task.

## A. Schema-9 annotation

Schema 9 retains the schema-8 normalized gene/CDS layer and adds physical TES,
transcript-to-TES, versioned downstream-region, downstream-to-gene, and
anchor-to-downstream-region tables. Its compact TP73 anchor surface records the
exhaustive `promoter > downstream > gene_body > intergenic` projection while
retaining every independent overlap flag and many-to-many relationship. The
normalized tables remain authoritative when several starts, ends, CDS
segments, transcripts, or genes are tied.

```sh
SOURCE=/data/sm718/GitHub/jaspar-mapping
SCAN=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_sparse_v3/package
GTF=/data/sm718/resources/ensembl/113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
ANNOTATION_RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_tp73_annotation_v4_schema9

"$SOURCE/scripts/submit_motif_context_slurm.sh" \
  --run-root "$ANNOTATION_RUN" --scan-package "$SCAN" --gtf "$GTF" \
  --anchor-only --output-tier summary \
  --chrom-file "$SOURCE/config/grch38_primary_nuclear_chromosomes.txt" \
  --annotation-release ensembl_113 \
  --promoter-definition tss_upstream_2000_downstream_500_v1 \
  --promoter-upstream-bp 2000 --promoter-downstream-bp 500 \
  --downstream-definition tes_upstream_500_downstream_2000_v1 \
  --downstream-upstream-bp 500 --downstream-downstream-bp 2000 \
  --account cluster --partition requeue --max-concurrent 20 \
  --cpus 4 --memory 32G --memory-limit 24GB --max-temp-size 100GB \
  --time 02:00:00 --dry-run
```

Remove `--dry-run` only after inspecting the immutable task plan. Each
chromosome task stages the two TP73 orientation files once and is safe to
requeue. The final catalog is `$ANNOTATION_RUN/final/context.duckdb`.

## B. H3K4me3 signal and GFP change

The signal manager reuses the completed TP73 evidence anchor set rather than
recomputing TP73/control immersion. One task per nuclear chromosome extracts
only the included H3K4me3/input tracks, preserving explicit zero signal when a
BigWig has no entries for that chromosome.

```sh
EVIDENCE=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_tp73_cutandrun_evidence_v1/final/genome_evidence
TRACKS=$SOURCE/cutandrun_20250602_noDuplicates
RUNTIME=/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_context_thresholds_v1/runtime
H3_RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_h3k4me3_anchor_signal_v1

"$SOURCE/scripts/submit_h3k4me3_genome_signal_slurm.sh" \
  --run-root "$H3_RUN" --evidence-package "$EVIDENCE" \
  --track-root "$TRACKS" --runtime-prefix "$RUNTIME" \
  --source "$SOURCE" --partition requeue --max-concurrent 12 \
  --cpus 4 --memory 32G --time 04:00:00 --dry-run
```

The durable signal contains five predeclared windows. The compact change table
uses the chromosome-1-selected `flank_150_1000` window and stores, per physical
anchor and accepted experimental series:

```text
log2((H3K4me3 area + 1)/(input area + 1))
delta_TA = TA - GFP
delta_DN = DN - GFP
```

The final catalog exposes complete, autosome-only, and sex-chromosome views.

## C. All-JASPAR cofactor inference

This stage consumes the finalized H3 change, schema-9 annotation, autosomal
TP73 evidence, and one zero-complete 150 bp context-maximum file per non-TP73
motif. A preflight job proves that all three fixed anchor layers have identical
physical keys before any model starts.

```sh
H3_PACKAGE=$H3_RUN/final/genome_h3k4me3_signal
ANNOTATION=$ANNOTATION_RUN/final
CONTEXT=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_tp73_context_maxima_autosomes_v1/final/context_maxima
ANALYSIS_RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_h3k4me3_cofactor_analysis_v3_schema9

"$SOURCE/scripts/submit_h3k4me3_cofactor_analysis_slurm.sh" \
  --run-root "$ANALYSIS_RUN" --h3-package "$H3_PACKAGE" \
  --evidence-package "$EVIDENCE" --context-maxima-package "$CONTEXT" \
  --annotation-catalog "$ANNOTATION" --runtime-prefix "$RUNTIME" \
  --run-id jaspar2026_grch38_h3k4me3_cofactor_analysis_v3_schema9 \
  --source "$SOURCE" --rscript "$RUNTIME/r/bin/Rscript" \
  --partition requeue --max-concurrent 20 \
  --motifs-per-batch 8 --cpus 4 --memory 32G --time 04:00:00 \
  --dry-run
```

Each batch stages the fixed autosomal layers once and evaluates motifs
sequentially. Every completed motif is validated and atomically published, so
a timeout or requeue repeats only its unfinished motif. `SIGUSR1` prints the
phase, batch, active motif, elapsed time, and child PID without terminating the
job.

The primary model is fit separately for SaOS-2 and SK-Mel-29 series 2 and for
TA-GFP and DN-GFP. It adjusts for TP73 motif score, chromosome, compact genomic
context, unsigned genomic distances to the nearest TSS and CDS, and explicit
upstream/downstream/overlap/mixed-strand direction classes. The
positive cofactor class uses the applied convenient threshold; strict `< -1`
and `< 0` references remain separate and are marked censored when the source
scan floor makes them unobservable. All-H3-zero anchors are omitted from
intensity fits but retained in occurrence summaries.

`cofactor_present x confirmed_TP73` is a secondary descriptive interaction:
TP73 confirmation is post-treatment and is deliberately omitted from the
primary total-association model. A continuous cofactor-score sensitivity and
context-stratified effects, including `strict_intergenic`, are emitted beside
the binary contrast.

The schema-9 evaluator also emits
`gene_relation_stratified_intensity_effect`: exactly 32 rows per motif (two
negative references by two isoforms by two series by four relation classes).
Its classes are promoter, downstream, gene body outside the two higher-
precedence regions, and intergenic. Underpowered classes remain explicit rows
with a non-`ok` status.

Array tasks contain one motif at a time for inferential output, so their local
`q_value_bh` is diagnostic only. Finalization recomputes
`q_value_bh_all_motifs` over every planned non-TP73 JASPAR motif within each
declared series/isoform/reference/result family. Scientific interpretation must
use the final value.

After finalization, build the compact joined tables and all six evidence plots
(the overall and detailed-context plots plus one four-panel plot for each
gene-relation class) with:

```sh
ENRICHMENT=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_tp73_cofactor_enrichment_autosomes_v1/final/cofactor_enrichment
INTERPRETATION=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_h3k4me3_cofactor_interpretation_v2_schema9
SOURCE_COMMIT=$(git -C "$SOURCE" rev-parse HEAD)

"$SOURCE/scripts/run_h3k4me3_genome_cofactor_interpretation.sh" \
  --analysis-run "$ANALYSIS_RUN" --enrichment-package "$ENRICHMENT" \
  --h3-package "$H3_PACKAGE" --output-dir "$INTERPRETATION" \
  --source "$SOURCE" --source-commit "$SOURCE_COMMIT" \
  --duckdb "$RUNTIME/duckdb/bin/duckdb" \
  --rscript "$RUNTIME/r/bin/Rscript"
```

The wrapper is restart-safe at its two publication boundaries. The compact
summary is an atomic directory produced by the Python extractor; all figures
are rendered into a job-specific staging directory and promoted together.
Existing complete products are validated and reused.

## Extensibility

The analysis grain remains one physical TP73 anchor by experimental series.
Future pre-treatment covariates such as repeat/ALU overlap, GC, mappability,
enhancer state, conserved syntenic identity, and detailed promoter membership
can join on the anchor key without changing the H3 outcome table. They should
be added as named feature blocks with pinned source versions. CUT&RUN-derived
TP73 confirmation remains a post-treatment evidence block, never silently
mixed into the primary sequence/annotation adjustment.

Use these status commands without traversing payload directories:

```sh
python3 "$SOURCE/scripts/manage_h3k4me3_genome_signal.py" status \
  --run-root "$H3_RUN"
python3 "$SOURCE/scripts/manage_h3k4me3_cofactor_analysis.py" status \
  --run-root "$ANALYSIS_RUN"
```

Contract tests are `tests/test_h3k4me3_genome_signal_manager.sh` and
`tests/test_h3k4me3_cofactor_analysis_manager.sh`.

# JASPAR 2026 negative-threshold sensitivity run

## Question

The first chromosome-1 TP73-context calibration only tested non-negative
integer thresholds. That was operationally consistent with the available
score floor, but it left the 438 motifs whose selected threshold was exactly
zero left-censored: a negative threshold could not win because it was never
offered to the model.

This sensitivity run asks whether that candidate-grid boundary induced a
meaningful filtering bias. It selects exactly those 438 accessions from
`thresholds/jaspar2026_grch38_chr1_tp73_context_v1.tsv`, rescans them around the
same 310,782 physical TP73 anchors, and evaluates every inclusive integer
threshold from -20 through 0.

The result is diagnostic. It does not modify or supersede the v1 threshold
list. A later decision may use the evidence to define a separately versioned
registry.

## Fixed contract

- Genome: GRCh38 Ensembl release 113 primary assembly, chromosome 1.
- Motifs: JASPAR 2026 CORE non-redundant accessions whose v1 recommendation is
  exactly 0; TP73 is not treated as its own cofactor.
- Score: `log2_relative_risk`, pseudocount 1, uniform A/C/G/T background.
- Anchors and outcomes: the immutable TP73/CUT&RUN evidence Parquet from the v1
  run.
- Context: signed interval-edge distance at most 150 bp.
- Recovered score floor: -20, inclusive.
- Candidate thresholds: every integer in `[-20, 0]`, inclusive.
- Feature: whether the maximum recoverable score around an anchor reaches the
  candidate threshold.
- Model, folds, and minimum class support: unchanged from v1.

Scores below -20 are not needed to test this grid. An anchor with no alignment
at or above -20 has `context_score = NULL` and is absent at every candidate.
Constant, under-supported, and nonconvergent candidates remain explicit metric
rows with an `evaluation_status`; they are not silently discarded.

## Direct Parquet scoring

`pssm_scan_parquet --context-maxima` reads a BED-coordinate anchor index and
writes one row per anchor directly to Zstandard-compressed Parquet. It does not
materialize chromosome-wide BED/TSV hits. The physical fields are:

```text
schema_version, chrom, anchor_start, anchor_end, motif_id, context_score,
source_score_floor, context_flank_bp, capture_prefilter_center_bp,
observed_max_anchor_span_bp, observed_max_context_span_bp,
context_distance_metric
```

The scanner checks every complete alignment on both strands whose motif
interval has signed edge distance at most 150 bp from the TP73 interval. The
same shared scoring function is used by sparse and dense scans. A regression
test compares the direct result with the established sparse-hit plus DuckDB
geometry on a synthetic chromosome.

On the local full-anchor pilot (`MA0007.4`), direct scoring took about 32 seconds
and produced a 3 MB Parquet with 310,782 rows. This is a planning observation,
not a cluster performance guarantee.

## Restart-safe execution

The manager writes an immutable 438-row task plan, a checksum-pinned run
configuration, and a 310,782-row anchor index. One requeueable Slurm array
element handles one motif:

1. validate the clean scanner build and planned inputs;
2. check durable free space before doing work;
3. stage chromosome 1 from the indexed FASTA to node-local `/scratch`;
4. write direct context maxima at source floor -20;
5. fit only the 21 thresholds from -20 through 0;
6. validate row counts, score floor, geometry, and candidate completeness;
7. atomically promote the result below `/data/sm718`.

An existing valid completion marker is reused after requeue. SIGUSR1 is
forwarded to the scanner during scoring; during other phases the wrapper emits
its own one-line phase report.

Example submission after building a clean Arrow scanner:

```sh
SOURCE=/data/sm718/GitHub/jaspar-mapping
INPUT=/data/sm718/codex/jaspar2026_chr1_dense_5motifs_5684241/input/public
V1_RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_context_thresholds_v1
RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_context_negative_sensitivity_v1

"$SOURCE/scripts/submit_negative_threshold_sensitivity_slurm.sh" \
  --run-root "$RUN" \
  --threshold-list "$SOURCE/thresholds/jaspar2026_grch38_chr1_tp73_context_v1.tsv" \
  --jaspar "$INPUT/JASPAR2026_CORE_non-redundant_pfms_jaspar.txt" \
  --genome "$INPUT/Homo_sapiens.GRCh38.dna.primary_assembly.fasta" \
  --anchor-evidence "$V1_RUN/input/tp73_chr1_anchor_evidence.parquet" \
  --runtime-prefix "$V1_RUN/runtime" \
  --source "$SOURCE" --account cluster --partition requeue \
  --max-concurrent 20 --memory 16G --time 00:30:00
```

Finalization creates consolidated threshold metrics and one comparison row per
motif. It reports the zero-threshold AUC gain, the strongest evaluable negative
candidate, their difference, support fractions, and a categorical sensitivity
conclusion. `negative_threshold_higher_auc` means only that the negative
candidate performed better in this exploratory chromosome-1 analysis; it is
not an automatic recommendation to change the production storage floor.

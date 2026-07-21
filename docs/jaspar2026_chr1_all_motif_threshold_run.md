# JASPAR 2026 chromosome-1 TP73-context threshold run

## Purpose and boundary

This run calibrates one convenient, inclusive score threshold for every JASPAR
2026 CORE non-redundant motif other than the TP73 target motif. The question is:

> Does the strongest alignment of this motif within signed interval-edge
> distance 150 bp of a TP73 alignment improve held-out prediction of strict
> anti-p73 CUT&RUN immersion over the TP73-score baseline?

The 2,632 neighboring motifs are independent Slurm tasks. `MA0861.2` is not
treated as its own cofactor: the anchor span would otherwise be the trivial
local maximum. TP73 direct support and TP73 tandem architecture keep their
separate threshold policies.

These are exploratory chromosome-1 transformations, not universal biochemical
binding cutoffs. Continuous scores and every production hit down to the -1
storage floor remain authoritative.

## Fixed scientific contract

- genome: GRCh38 Ensembl release 113 primary assembly, chromosome 1;
- motif set: JASPAR 2026 CORE non-redundant;
- score: `log2_relative_risk`, pseudocount 1, uniform A/C/G/T background;
- target anchors: physical `MA0861.2` alignment spans with maximum orientation
  score at least 0;
- neighboring candidates: production sparse records with score at least -1;
- context: signed interval-edge distance at most 150 bp;
- motif feature: maximum neighboring score per TP73 anchor;
- candidate thresholds: every observed non-negative integer, inclusive;
- fitting support: retained and absent classes must each contain at least 1% of
  physical anchors;
- outcome: strict immersion in anti-p73 only versus matched control only, for
  six TA/DN sample pairs;
- validation: five equal-width contiguous chromosome-1 folds;
- selection: maximum held-out macro ROC-AUC gain over the TP73-score baseline,
  with the lower threshold winning a tie.

No fit failure is converted to a threshold. Constant, under-supported, or
nonconvergent candidates receive explicit unevaluable metric rows. A motif with
no eligible positive gain retains a null recommendation and a declared status.

## Reconstructing the shared labels

The run does not depend on a file copied from a workstation. The setup job
reconstructs the 310,782-row anchor table from:

1. the exact plus/minus `MA0861.2` chromosome-1 paths in the production
   `scan_file_inventory`;
2. the six existing `tp73_*.clipped.bedGraph` tracks; and
3. chromosome-1 exports of the six matched `neg_*.bigWig` controls.

Plus and minus records at an identical alignment span are one physical anchor;
their maximum score is retained. Adjacent positive bedGraph segments are merged.
A merged component supports `[anchor_start,anchor_end)` only when its start is
strictly lower and its end strictly greater. `build_tp73_anchor_evidence.py`
processes one track at a time, bounding memory, and records source paths and
file statistics beside the Parquet output.

## Task and storage layout

The submission helper creates an immutable plan from the scan catalog rather
than discovering files with a package-wide glob:

```text
RUN/
  plan/
    calibration_tasks.tsv
    target_anchor_files.tsv
    motifs.txt
    run_config.json
    runtime-explicit.txt
    slurm_submission.tsv
  input/
    control_bedgraph/
    tp73_chr1_anchor_evidence.parquet
  tasks/
    task-000000-MA....../
      cofactor_maxima.parquet
      threshold_metrics.tsv
      baseline_metrics.tsv
      fold_manifest.tsv
      evaluator_run_config.tsv
      complete.json
  final/threshold_calibration/
    threshold_metrics.parquet
    tables/jaspar2026/motif_score_threshold/part-000000.parquet
    manifest.json
```

Each task reads exactly its two inventory paths. DuckDB spill and incomplete
working output live on node-local `/scratch`; a validated task directory is
atomically promoted on `/data`. Requeued tasks reuse a compatible completed
directory. Source scan files are opened read-only and are never rewritten.

## Haumea submission

Use a dedicated run below `/data/sm718`. The helper creates a recorded
Micromamba runtime when one is absent, submits the shared-label setup, submits
the dependent 2,632-element array with at most 20 live tasks, and submits a
finalizer dependent on the complete array:

```sh
SOURCE=/data/sm718/GitHub/jaspar-mapping
SCAN=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_sparse_v3/package
RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_context_thresholds_v1

"$SOURCE/scripts/submit_motif_threshold_calibration_slurm.sh" \
  --run-root "$RUN" \
  --scan-package "$SCAN" \
  --jaspar /data/sm718/jaspar_mapping_runs/jaspar2026_grch38_sparse_v3/input/JASPAR2026_CORE_non-redundant_pfms_jaspar.txt \
  --cutandrun-dir "$SOURCE/cutandrun_20250602_noDuplicates" \
  --source "$SOURCE" \
  --account cluster --partition requeue \
  --max-concurrent 20 --cpus 2 --memory 16G --time 02:00:00
```

The array's `--requeue` policy protects node interruptions; the scheduled USR1
report five minutes before the time limit provides a final phase indication.
The July 22 Haumea maintenance remains an operational reason to inspect job
state rather than assume uninterrupted execution.

## Progress and completion

Request one status line from a running array element without interrupting its
DuckDB/R child:

```sh
scancel --signal=USR1 --batch ARRAY_JOB_ID_TASK_INDEX
```

Inspect durable progress without scanning the large result tree:

```sh
python3 "$SOURCE/scripts/manage_motif_threshold_calibration.py" status \
  --run-root "$RUN"
```

The manager reads the immutable task plan and each exact expected completion
marker. Finalization refuses missing tasks and verifies that the registry has
exactly 2,632 rows and no accidental `pending` rows. `pending` is reserved for
missing metric input; motifs that genuinely have no usable threshold instead
receive `insufficient_class_support`, `no_finite_metric`, or
`no_positive_gain`.

The local KLF14 reference sweep took 174 seconds before compact-output and
class-support shortcuts. That makes several hours of wall time at 20-way
parallelism a reasonable planning estimate, but motif density and filesystem
load can widen it. The task logs and exact completion markers, not that estimate,
determine progress.

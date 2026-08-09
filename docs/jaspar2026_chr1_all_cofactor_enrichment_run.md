# JASPAR 2026 TP73 Cofactor Enrichment Run

## Purpose

This run asks, for every non-TP73 motif in JASPAR 2026, whether a
threshold-qualified motif occurrence within the interval-defined 150 bp TP73
context is enriched or depleted among strict anti-p73 CUT&RUN immersion events
relative to the matched control track. It also retains descriptive effects for
eight TP73 score strata and six track-normalized CUT&RUN depth tiers.

The result is exploratory. The operating points and enrichment measurements
both come from chromosome 1, and the anchor package lacks GC, mappability,
repeat/Alu, accessibility, and GFP covariates. The primary result is therefore
an all-motif overview and hypothesis generator, not independent validation or
evidence of causal binding.

Read
[`tp73_cofactor_cutandrun_enrichment.md`](tp73_cofactor_cutandrun_enrichment.md)
for the exact positive, negative, intermediate, immersion, depth, and matched
model definitions.

## Reused Inputs

The enrichment run does not rescan chromosome 1 and does not rebuild motif
contexts. It reuses the completed threshold-calibration package:

```text
/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_context_thresholds_v1
```

That package owns:

- one shared TP73/CUT&RUN anchor Parquet;
- one exact schema-v4 `cofactor_maxima.parquet` for each of 2,632 non-TP73
  motifs; and
- the finalized motif-score threshold registry.

The immutable enrichment plan records the exact source task directory, size,
and SHA-256 for every maxima file. A worker reads the path named in its plan
row; it does not use `*`, recursive discovery, or a directory listing.

For 2,615 motifs, the positive operating point is the historical threshold
recommendation. Seventeen registry rows have no recommendation. They remain
`NULL` in `recommended_threshold`; to keep the requested all-motif overview
complete, their tasks use threshold zero with
`positive_threshold_role = 'descriptive_fallback'`. This fallback is explicit
and is not presented as a newly inferred recommendation.

## Slurm Layout

One array task evaluates one motif. It copies only the shared 2.5 MB anchor
file and that motif's maxima file to node-local `$SLURM_TMPDIR` or `/scratch`,
runs the R evaluator in memory, validates all six compact tables, and atomically
promotes a task directory below `/data/sm718`. Source maxima are never copied
into the new durable run.

The submission is split into scheduler-compatible chunks. Chunks depend on the
preceding chunk, so at most 20 tasks are live across the complete submission.
Tasks request one CPU, 16 GiB, 30 minutes, and the `requeue` partition. A
dependent finalizer runs only after every task succeeds.

Interrupted attempts have unique names and do not count as complete. A
requeued task reuses an existing final directory only after its identity and
file inventory validate. No script deletes durable source, task, or final data.

## Submit

From the clean source checkout on Haumea:

```sh
SOURCE=/data/sm718/GitHub/jaspar-mapping
SOURCE_RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_context_thresholds_v1
RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_cofactor_enrichment_v1

"$SOURCE/scripts/submit_tp73_cofactor_enrichment_slurm.sh" \
  --run-root "$RUN" \
  --source-threshold-run "$SOURCE_RUN" \
  --source "$SOURCE" \
  --partition requeue \
  --max-concurrent 20 \
  --memory 16G \
  --time 00:30:00
```

The source run's pinned DuckDB and R runtimes are reused by default. The plan,
source hashes, submission IDs, and execution commit are recorded below
`$RUN/plan/`.

## Progress And Restart

Exact completion counts do not require scanning generated data:

```sh
python3 "$SOURCE/scripts/manage_tp73_cofactor_enrichment.py" status \
  --run-root "$RUN"
```

Slurm and recent logs provide live state:

```sh
squeue -u "$USER"
tail -n 30 "$RUN"/logs/enrichment-JOB_TASK.err
```

Slurm is configured to send `SIGUSR1` two minutes before the time limit. A
manual signal to a selected live array element also produces exactly one line
containing its phase, global task index, motif ID, elapsed time, and child PID:

```sh
scancel --signal=USR1 --batch ARRAY_JOB_ID_ARRAY_TASK_ID
```

If submission itself is interrupted before `slurm_submission.tsv` is written,
inspect `squeue` before resubmitting. Once the submission record exists, do not
submit a second array. Failed or preempted elements can be requeued normally;
completed elements are idempotent.

## Final Query Surface

The finalizer verifies every exact output checksum, confirms that all workers
used the same depth-tier manifest, and applies Benjamini-Hochberg once with the
complete 2,632-motif family size. Per-task one-row q-values are discarded.

The published directory is:

```text
$RUN/final/cofactor_enrichment/
```

It contains Zstandard-compressed Parquet tables, a compact TSV overview, a
manifest, and `tp73_cofactor_enrichment.duckdb`. The main tables are:

- `cofactor_overview`: one row per motif with class prevalence, primary matched
  odds ratio and confidence interval, all-JASPAR q-value, direction, and the
  strict-immersion sample-macro specificity ratio;
- `cofactor_primary_occupancy`: the complete primary model result;
- `cofactor_macro_summary`: motif by negative reference, TP73 score stratum,
  and CUT&RUN depth tier;
- `cofactor_descriptive`: the corresponding per-sample results;
- `cofactor_class_count`: positive, negative, and intermediate anchor counts;
- `cofactor_depth_tier`: the shared normalized depth thresholds; and
- `cofactor_task`: the exact threshold and source provenance for every motif.

Useful queries include:

```sql
-- Strongest estimable anti-p73 enrichments and depletions.
SELECT motif_id, factor_name, positive_threshold, positive_anchor_fraction,
       adjusted_odds_ratio, confidence_interval_95_lower,
       confidence_interval_95_upper, q_value_bh_all_jaspar,
       association_direction
FROM cofactor_overview
WHERE evaluation_status = 'ok'
ORDER BY q_value_bh_all_jaspar,
         abs(ln(adjusted_odds_ratio)) DESC;

-- Effect progression over TP73-score and CUT&RUN-depth strata.
SELECT motif_id, factor_name, tp73_score_stratum, depth_tier,
       mean_log2_anti_control_specificity_ratio_jeffreys,
       samples_anti_p73_enriched, samples_anti_p73_depleted
FROM cofactor_macro_summary
WHERE negative_reference_threshold = -1
ORDER BY motif_id, tp73_score_stratum, depth_tier_order;

-- Motifs for which threshold zero was only an explicit coverage fallback.
SELECT motif_id, factor_name, recommended_threshold, positive_threshold,
       source_calibration_status
FROM cofactor_overview
WHERE positive_threshold_role = 'descriptive_fallback';
```

## Local Contract Test

The small synthetic test creates two source tasks, including one null threshold,
runs and restarts both workers, finalizes the package, and queries its DuckDB:

```sh
bash tests/test_tp73_cofactor_enrichment_manager.sh
```

# Chromosome 1 PATZ1/TP73 Dense Dry Run

This dry run exercises the direct Parquet and DuckDB path for two JASPAR 2026
motifs:

| Factor | Motif ID | Model length |
| --- | --- | ---: |
| PATZ1 | `MA1961.2` | 11 bp |
| TP73 | `MA0861.2` | 16 bp |

For a five-base introduction where every expected score can be checked by
hand, begin with
[`../test_files/synthetic_dense/README.md`](../test_files/synthetic_dense/README.md).

Every run uses pseudocount 1, both reference orientations, skip-N handling,
and both `log2_relative_risk` and `log_odds`. This gives eight dense score
configurations. A stored element represents one possible motif-model alignment
start. It is not an assertion that a TF binds there and its `start`/`end`
interval is not an asserted physical footprint of the TF complex.

## Bounded Default

The default interval is chr1:3,600,000-3,700,000, which contains the TP73
locus and promoter reference region while keeping the end-to-end test small:

```sh
make dry_run_chr1_patz1_tp73
```

The target builds `pssm_scan_parquet`, creates the FASTA `.fai` if needed,
runs all eight configurations, validates their Parquet inventory, writes a
manifest, and creates:

```text
dry_runs/chr1_patz1_tp73_from-3600000-to-3700000/
  manifest.json
  jaspar2026_chr1_patz1_tp73.duckdb
  tables/jaspar2026/motif_metadata/...
  tables/jaspar2026/motif_score_dense/
    motif_id=.../score_mode=.../pseudocount=1/chrom=1/strand=.../...
```

Generated dry-run packages are ignored by Git. The runner never cleans an
output directory. It reuses complete expected Parquet pairs and rejects stale,
extra, or overlapping parts when their inventory disagrees with the manifest;
use a new output directory for a different run.

## Inspection

The overview is block-level and does not expand score vectors:

```sh
make inspect_chr1_patz1_tp73
```

The read-only inspector also supports focused operations:

```sh
scripts/inspect_chr1_dense_dry_run.sh overview
scripts/inspect_chr1_dense_dry_run.sh files
scripts/inspect_chr1_dense_dry_run.sh \
  region MA0861.2 log_odds + 3651800 3651820
scripts/inspect_chr1_dense_dry_run.sh \
  summary MA1961.2 log2_relative_risk - 3600000 3700000
scripts/inspect_chr1_dense_dry_run.sh \
  histogram MA0861.2 log_odds + 0.5 3600000 3700000
scripts/inspect_chr1_dense_dry_run.sh \
  bins MA0861.2 log2_relative_risk + 3600000 3700000
scripts/inspect_chr1_dense_dry_run.sh shell
```

Use `--package DIR` before the command for a non-default package. The `shell`
command opens DuckDB read-only from the package root, which is required because
the persisted views deliberately use package-relative Parquet paths.

The schema in [`../sql/chr1_dense_dry_run_schema.sql`](../sql/chr1_dense_dry_run_schema.sql)
provides:

- `run_manifest`: source commit and dirty flag, scanner/JASPAR/index hashes,
  chromosome/range, and run policy.
- `motif_metadata`: motif identity and model length.
- `motif_score_dense_block`: physical Parquet blocks without list expansion.
- `dense_run_inventory`: cheap file/block/window validation.
- `dense_scores_region(...)`: row expansion for a bounded start range.
- `dense_score_summary(...)`: valid/skipped counts and score quantiles.
- `dense_score_histogram(...)`: arbitrary positive fixed-width bins.
- `dense_score_calibration_bins(...)`: the shared biologically motivated bins.

Ready-to-run SQL examples are in
[`../sql/chr1_dense_dry_run_examples.sql`](../sql/chr1_dense_dry_run_examples.sql).

## Dynamic BED Export

`scripts/export_dense_bed.py` regenerates sparse, thresholded BED rows from a
dense package without rescanning the genome or creating an intermediate file.
First inspect the available motif/configuration partitions:

```sh
scripts/export_dense_bed.py \
  --package dry_runs/chr1_patz1_tp73_from-3600000-to-3700000 \
  --list-configs
```

Count a proposed TP73 selection before writing it:

```sh
scripts/export_dense_bed.py \
  --package dry_runs/chr1_patz1_tp73_from-3600000-to-3700000 \
  --motif MA0861.2 --score-mode log_odds --pseudocount 1 \
  --chrom 1 --strand both --from 3651800 --to 3652600 \
  --min-score 0 --count-only
```

Then stream the same selection directly to compressed BED:

```sh
scripts/export_dense_bed.py \
  --package dry_runs/chr1_patz1_tp73_from-3600000-to-3700000 \
  --motif MA0861.2 --score-mode log_odds --pseudocount 1 \
  --chrom 1 --strand both --from 3651800 --to 3652600 \
  --min-score 0 --columns provenance \
  --output TP73_MA0861.2_log_odds_chr1_3651800-3652600.bed.gz
```

`--chrom` is repeatable and accepts comma-separated values. `--min-score` and
`--max-score` are inclusive; the historical lower threshold of zero is the
default, while `--all-scores` explicitly removes it. `--no-header` supports
strict BED consumers, and `--score-decimals 3` reproduces the old display
precision. Existing output files are never replaced unless `--force` is given.

The range is a BED-style, 0-based, half-open selection of motif alignment
starts. Consequently, the reported sequence interval can extend beyond
`--to`; that interval records the sequence used for scoring, not the physical
footprint of a bound TF complex. NULL scores from skipped sequence are never
exported.

The default output is the compatible first six columns `Chromosome`, `From`,
`To`, `Name`, `Score`, and `Strand`. `--columns provenance` appends `MotifID`,
`ScoreMode`, `Pseudocount`, and `CoordinateMode`. Dense storage does not contain
the matched sequence or the motif-specific finite score range needed to recover
`PWMRelativeScore`, so the exporter deliberately does not fabricate those two
columns.

## Full Chromosome 1

The full run is deliberately gated because it contains approximately two
billion float scores across two motifs, two modes, and two strands: about 8 GB
of float payload before Parquet compression and metadata. It requires at least
12 GiB free in the output filesystem:

```sh
make dry_run_chr1_patz1_tp73 \
  DRY_RUN_FULL_CHR1=1 \
  DRY_RUN_OUTPUT=dry_runs/chr1_patz1_tp73_full
```

Send `SIGUSR1` to the active scanner process for an on-demand progress line.
The bounded dry run should pass before launching the full version.

## Interpretation Boundary

This package is sequence layer 1 only. It contains neither CUT&RUN signal nor
promoter/expression annotations. It can verify score calculation, orientation,
partitioning, query latency, and score distributions. It cannot establish
occupancy, PATZ1-p73 interaction at a locus, promoter regulation, or a useful
genome-wide storage threshold. Those require the later CUT&RUN calibration and
annotation layers described in
[`promoter_architecture_ml_schema.md`](promoter_architecture_ml_schema.md).

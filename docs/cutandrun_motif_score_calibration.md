# CUT&RUN Coverage Immersion And Motif-Score Calibration

## Question

For TP73 motif-model alignment spans on chromosome 1, estimate how the sequence
score changes:

1. the chance that positive CUT&RUN coverage completely immerses the scored DNA
   span; and
2. the maximum CUT&RUN depth within that span.

The completed real-data chromosome-1 threshold analysis is recorded in
[`tp73_chr1_cutandrun_threshold_20260713.md`](tp73_chr1_cutandrun_threshold_20260713.md).
The follow-up using the highest local PATZ1 score for each TP73 alignment is
recorded in
[`tp73_patz1_chr1_cutandrun_20260714.md`](tp73_patz1_chr1_cutandrun_20260714.md).
The corresponding E2F1 and three-model TFAP2C comparison is recorded in
[`tp73_e2f1_tfap2c_chr1_cutandrun_20260714.md`](tp73_e2f1_tfap2c_chr1_cutandrun_20260714.md).
The negative-predictor analysis using the strongest local POU2F2 score is
recorded in
[`tp73_pou2f2_chr1_cutandrun_20260717.md`](tp73_pou2f2_chr1_cutandrun_20260717.md).
The follow-up stratifying TP73 candidates by nearby TP73 pair architecture and
combining PATZ1, TFAP2C, POU2F2, and E2F1 is reported in
[`tp73_pair_stratified_chr1_cutandrun_20260718.md`](tp73_pair_stratified_chr1_cutandrun_20260718.md).

The scored span is the PSSM alignment interval `[start,end)`. It remains a
computational interval rather than an assertion about the complete physical
footprint of the TP73 complex. Here it is the requested minimum interval that
positive coverage must exceed on both sides.

## Coverage-union rule

All coordinates are BED 0-based half-open. First merge overlapping and directly
adjacent positive-coverage intervals. A merged component `[c_start,c_end)`
strictly immerses a motif-model span `[m_start,m_end)` only when:

```text
c_start < m_start AND c_end > m_end
```

Equality on either boundary is not sufficient. A gap inside the motif is not
sufficient. Directly adjacent intervals such as `[8,12)` and `[12,16)` are
treated as the continuous union `[8,16)`. This deliberately accepts the stated
error in which two reads jointly support a motif although neither read spans it
alone.

The local audit table records:

- `ordinary_max_depth`: maximum input depth within the motif (overlapping
  fragment count or bedGraph column-4 signal);
- `effective_max_depth`: ordinary maximum depth when a merged component
  strictly immerses the motif, otherwise 0; and
- the supporting component coordinates and source interval names.

Keeping ordinary and effective depth side by side exposes boundary mistakes.
The synthetic fixture includes two spans with ordinary depth 1 but effective
depth 0, plus one span immersed only after two adjacent intervals are joined.

## Accepted coverage inputs

Read/fragment BED, bedGraph, and bigWig are all conceptually sufficient because
read identity is not part of the label:

- BED rows are merged into the union of covered bases. Their overlap count can
  provide integer depth.
- bedGraph or bigWig segments with signal greater than zero define the covered
  union, while their signal supplies maximum depth.
- BAM may be converted to BED/bedGraph or streamed into the same union.

Both analyzers accept BED3+ and bedGraph, plain or gzip-compressed. Use
`--coverage-format fragments` when each BED row is one fragment and contributes
unit depth. Use `--coverage-format bedgraph` when column 4 is a positive finite
depth or signal. The default `auto` selects bedGraph semantics only for
`.bedGraph[.gz]` or `.bdg[.gz]`; all other names retain fragment semantics.
This conservative rule avoids mistaking a numeric BED name for depth.

BigWig is conceptually sufficient but is not read directly by these scripts;
convert it to bedGraph first. Prefer the supplied `.clipped.bedGraph` files over
the larger legacy `.clipped.clean.bed` copies because the original bedGraph
keeps depth unambiguously in column 4.

## Threshold metrics

At each motif-score threshold `t`, retain motif spans with `score >= t` and
report:

- **Coverage-component recall**: fraction of merged positive-coverage
  components that immerse at least one retained motif span. This replaces the
  earlier read-centric fraction once read identity is discarded.
- **Motif precision**: fraction of retained motif spans immersed in positive
  coverage. This is precision (positive predictive value), not sensitivity.
- **Motif recall / true-positive rate**: fraction of all immersed candidate
  spans retained at `t`.
- **Motif false-positive rate**: fraction of all non-immersed candidate spans
  retained at `t`.
- **Mean effective depth**: effective maximum depth averaged over all retained
  spans, including zero-depth spans.
- **Precision enrichment**: motif precision divided by the unthresholded
  prevalence of immersed spans.

Motif recall versus motif precision gives the precision-recall curve and
average precision. Motif recall versus false-positive rate gives the ROC curve
and ROC AUC. Because chromosome-wide candidate spans are extremely imbalanced,
average precision and precision enrichment are more informative for the
storage decision than ROC AUC alone.

CUT&RUN absence is an operational negative label, not proof that TP73 cannot
bind. Accessibility, antibody efficiency, cell state, recovery, and sequencing
depth can all create false negatives.

## Synthetic local test

The human-readable inputs are:

- `test_files/synthetic_cutandrun/tp73_motifs.bed`
- `test_files/synthetic_cutandrun/tp73_fragments.bed`
- `test_files/synthetic_cutandrun/tp73_coverage.bedGraph`

Run:

```sh
make check_cutandrun_containment
make synthetic_cutandrun_example
```

The retained example is written under
`dry_runs/synthetic_cutandrun_coverage_union/` with:

- `motif_evidence.tsv`: ordinary/effective depth per motif span;
- `coverage_component_evidence.tsv`: merged coverage components and their best
  immersed motif score;
- `threshold_curve.tsv`: exact score-threshold metrics; and
- `summary.json`: ROC AUC, average precision, and candidate thresholds.

The synthetic labels are intentionally imperfect. ROC AUC is 0.7333, average
precision is 0.8762, maximum Youden J occurs at score 6, and maximum F1 occurs
at score -2. These values validate the implementation and have no biological
meaning.

## True-data landing directory

Place true coverage files below the dedicated run root:

```text
/data/sm718/codex/jaspar2026_chr1_dense_5motifs_5684241/
  input/cutandrun/tp73/
    samples.tsv
    skmel29.TA.R1.GRCh38.coverage.bed.gz
    skmel29.DN.R1.GRCh38.coverage.bigWig
    skmel29.GFP.R1.GRCh38.coverage.bedGraph.gz
    ...
```

Keep original source files and checksums; do not place them in the Git
checkout. `samples.tsv` should record at least:

```text
sample_id filename format cell_line isoform antibody replicate duplicate_policy reference_build
```

Use GRCh38 coordinates matching the Ensembl release-113 FASTA used by the
dense run. Analyze biological replicates separately before any declared pooled
analysis. Existing files need not be renamed if `samples.tsv` maps them.

## Full chromosome implementation

For threshold calibration, `cutandrun_score_calibration` is the production
full-chromosome path. It streams paired plus/minus dense Parquet blocks and any
number of sorted, non-overlapping positive-depth bedGraph tracks. It advances
through each score and coverage stream once, writes only compact histograms and
curves, and never materializes chromosome-wide per-position evidence.

For an explicitly row-level downstream model, add `--feature-parquet FILE` to
a joint anchor/context run. The option writes one Zstandard-compressed Parquet
row for each anchor at or above `--minimum-anchor-score`; it is disabled by
default. Each row contains the exact anchor and best local context score,
0-based half-open spans, strand codes (`1` for plus and `-1` for minus), signed
center distance, orientation agreement, and one strict-support/maximum-depth
column pair per coverage track. The context is still the highest
orientation-collapsed score within `--context-flank`, with the existing
earliest-start/plus-orientation tie rule. Use the exact table when several
cofactor maxima must be joined by anchor coordinate; the compact joint
histogram cannot recover that co-occurrence after aggregation.

For TP73 motif length `L`, one merged component `[c_start,c_end)` can immerse
only alignment starts in:

```text
[c_start + 1, c_end - L)
```

The streaming calibrator tests this condition while maintaining the maximum
active bedGraph depth across each motif span. Combining the resulting compact
summaries is sufficient for fixed-bin thresholds, ROC/PR metrics, coverage
component recall, and depth summaries.

`scripts/analyze_dense_cutandrun_coverage.py` remains useful when a bounded
region needs a row-level `immersed_motif_evidence` audit. Its detailed
per-position output and base expansion are not the production path for a full,
deeply covered chromosome.

### Pseudocount comparison

For the first smoothing comparison, analyze only `log2_relative_risk` with
per-nucleotide pseudocounts 0 and 1. The latter is add-one smoothing and adds a
total count of 4 to each motif column. In the unsmoothed run, an alignment that
selects any zero-count nucleotide has score `-Inf`: favorable positions cannot
compensate for that mismatch. Keep these alignments as one lowest-score tie;
only genomic `N` under `--skip-N` is `NULL`. Both runs must therefore report the
same total candidate and immersed-alignment counts.

The calibration script accepts a repeatable `--score-mode`; use
`--score-mode log2_relative_risk` for this comparison. Its `-Inf` threshold row
is the complete ROC endpoint, with recall and false-positive rate both equal to
1. Unsmoothed per-position logistic odds are deferred because mixed positive
and negative infinities require a separate score definition.

Generate the two dense partitions with identical arguments apart from the
pseudocount; the hive path keeps them separate even in one output package:

```sh
./pssm_scan_parquet --dense-scores --outdir RUN/dense \
  --genome Homo_sapiens.GRCh38.dna.primary_assembly.fasta \
  --pssm JASPAR2026_CORE_non-redundant_pfms_jaspar.txt \
  --motif MA0861.2 --chr 1 --strand both \
  --motif-set-id jaspar2026_core_nonredundant \
  --genome-id homo_sapiens_grch38_ensembl113_primary \
  --score-mode log2_relative_risk --pseudocount 0

./pssm_scan_parquet --dense-scores --outdir RUN/dense \
  --genome Homo_sapiens.GRCh38.dna.primary_assembly.fasta \
  --pssm JASPAR2026_CORE_non-redundant_pfms_jaspar.txt \
  --motif MA0861.2 --chr 1 --strand both \
  --motif-set-id jaspar2026_core_nonredundant \
  --genome-id homo_sapiens_grch38_ensembl113_primary \
  --score-mode log2_relative_risk --pseudocount 1
```

Run `cutandrun_score_calibration` once per score definition and pass all
coverage tracks for that comparison as repeated `--coverage ID=FILE`
arguments. For a bounded row-level audit, run
`scripts/analyze_dense_cutandrun_coverage.py` once per pseudocount with the
same `--package`, `--coverage-bed`, and
`--score-mode log2_relative_risk`.

### Memory safety

The streaming C++ calibrator keeps coverage runs, components, compact
histograms, and one Parquet record batch in memory. It does not need a DuckDB
spill directory. The detailed Python/DuckDB auditor defaults to one DuckDB
thread, a `6GB` memory limit, disabled insertion-order preservation, and at
most `40GB` of temporary spill data. Those limits remain deliberate for
bounded audits on a 16 GB workstation.

On macOS, inspect system pressure and the DuckDB worker while a run is active:

```sh
memory_pressure -Q
pgrep -fl 'duckdb.*calibration'
ps -p PID -o pid=,rss=,%mem=,etime=,state=
```

Treat a sharply falling free percentage, swap growth, or an unresponsive system
as the stop signal. A large but bounded RSS alone is not sufficient if
`memory_pressure` still reports ample reclaimable memory. Check free disk space
beforehand because `--max-temp-size` is a ceiling rather than a reservation.

The primary unstranded analysis collapses the two orientation records to one
genomic span using their maximum score. `log2_relative_risk` and `log_odds`
are evaluated independently with explicit pseudocount provenance.

Canonical production outputs are compact score histograms, threshold curves,
and one calibration summary per sample and score mode. Detailed positive-span
Parquet remains an optional bounded audit artifact. The final genome-wide
storage threshold must not be selected from Youden J alone. Report Youden and
F1 candidates, then declare an effect-size and storage rule against matched
controls, as in the completed TP73 chromosome-1 analysis.

### Optional second-motif score

The streaming calibrator can also attach the single highest orientation-
collapsed context-motif score to each anchor without materializing a genomic
join. Supply paired context Parquet files plus the context motif length and
center-distance flank:

```sh
./cutandrun_score_calibration \
  --plus-parquet TP73_PLUS.parquet --minus-parquet TP73_MINUS.parquet \
  --context-plus-parquet PATZ1_PLUS.parquet \
  --context-minus-parquet PATZ1_MINUS.parquet \
  --context-motif-id MA1961.2 --context-motif-length 11 \
  --context-flank 150 --minimum-anchor-score -10 \
  --coverage SAMPLE=TRACK.bedGraph --output-dir RUN/joint_risk_p1
```

`joint_score_histogram.tsv` records the two score bins, support/depth evidence,
orientation agreement, and aggregate center distances. `joint_run_config.tsv`
pins both Parquet inputs and the deterministic tie rule. The ordinary TP73
histogram and curve are still emitted and should agree with a TP73-only pass.

After matched anti-p73 and IgG passes, use
`summarize_tp73_patz1_cutandrun_threshold.R`. Despite its historical filename,
`--context-label` makes its tables and plots cofactor-generic. It compares a
context-motif gate and assisted lower-TP73 rescue rule with a TP73-only
expectation at exactly the same storage count. The exact count match
interpolates within one tied 0.2-score bin; it does not invent a finer
deterministic TP73 threshold. Summary tables retain the minimum, median, and
maximum across the six TA/DN comparisons and nominate an exploratory beneficial
or detrimental rule only when all six support and depth effects agree on its
direction. Use `compare_tp73_cofactor_summaries.R` to compare the resulting
selected-policy tables at approximately equal retention.

# CUT&RUN Coverage Immersion And Motif-Score Calibration

## Question

For TP73 motif-model alignment spans on chromosome 1, estimate how the sequence
score changes:

1. the chance that positive CUT&RUN coverage completely immerses the scored DNA
   span; and
2. the maximum CUT&RUN depth within that span.

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

- `ordinary_max_depth`: maximum input-interval depth within the motif;
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

The synthetic analyzer accepts BED3+ (plain or gzip-compressed). The full
DuckDB run will consume BED/bedGraph directly and use a bigWig reader when the
true source is bigWig; an intermediate motif BED or TSV is not required.

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

The full run consumes dense Parquet blocks and CUT&RUN coverage directly. For
TP73 motif length `L`, one merged component `[c_start,c_end)` can immerse only
alignment starts in:

```text
[c_start + 1, c_end - L)
```

DuckDB can generate those starts, look scores up directly in dense block
vectors, and aggregate positive evidence without expanding the complete dense
chromosome into a stored row table. The all-candidate score histogram comes
directly from per-block list histograms; requested score modes are processed
sequentially. Combining those compact summaries is sufficient for exact or
user-selected thresholds, ROC/PR metrics, and depth summaries.

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
  --score-mode log2_relative_risk --pseudocount 0

./pssm_scan_parquet --dense-scores --outdir RUN/dense \
  --genome Homo_sapiens.GRCh38.dna.primary_assembly.fasta \
  --pssm JASPAR2026_CORE_non-redundant_pfms_jaspar.txt \
  --motif MA0861.2 --chr 1 --strand both \
  --score-mode log2_relative_risk --pseudocount 1
```

Run `scripts/analyze_dense_cutandrun_coverage.py` once per pseudocount with the
same `--package`, `--coverage-bed`, and `--score-mode log2_relative_risk`.

### Memory safety

`scripts/analyze_dense_cutandrun_coverage.py` defaults to one DuckDB thread, a
`6GB` DuckDB memory limit, disabled insertion-order preservation, and at most
`40GB` of temporary spill data below its output directory. These are deliberate
defaults for a 16 GB workstation. A DuckDB `Out of Memory Error` at the stated
limit is a controlled query failure; it prevents the query from consuming all
system memory. Do not remove the limit merely to make such a query continue.

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
genomic span using their maximum score. A secondary table retains
orientation-specific results. `log2_relative_risk` and `log_odds` are evaluated
independently with pseudocount 1 and explicit provenance.

Canonical outputs are Parquet tables for positive span evidence, threshold
curves, and one compact calibration summary per sample and score mode. The
final genome-wide storage threshold should not be selected from Youden J alone.
Report Youden and F1 candidates, then choose the lowest score that preserves an
agreed coverage-component recall while maintaining useful motif-precision
enrichment across biological replicates.

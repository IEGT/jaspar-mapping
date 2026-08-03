# Scanner maintenance contract

This document defines the operational boundary that must remain sound before
motif context, gene annotation, CUT&RUN calibration, or machine learning is
built from a genome scan. It covers `pssm_scan_parquet`, the immutable scan
plan, Slurm execution, finalization, integrity verification, and the DuckDB
inventory catalog. It does not define downstream biological features.

## Stable scientific contract

- Coordinates are BED 0-based half-open.
- Every retained row is one motif-model alignment start and orientation.
- Genome, motif set, score mode, pseudocount semantics, threshold policy,
  ambiguous-base policy, and source commit are explicit provenance.
- A durable task is one chromosome and one non-overlapping motif batch.
- A task is visible below `package/task_data/` only after both orientations
  have closed, passed row/configuration checks, and received SHA-256 digests.
- Existing promoted tasks and payload files are never overwritten.

The completed GRCh38 v3 package remains valid under this contract. Maintenance
does not rescan, rewrite, consolidate, or otherwise modify its hit payloads.

## Why two FASTA indexes appear

The index has no biological or scoring role. A `.fai` records byte offsets and
line geometry for FASTA records.

1. The source index lets `stage_fasta_region.py` seek directly to one
   chromosome in the immutable multi-record FASTA on `/data`. Without it, every
   chromosome worker would have to stream the complete assembly merely to find
   its input sequence.
2. The staged chromosome is written as a new one-record FASTA with regular
   60-base lines. Its small local `.fai` describes that new file so the existing
   indexed scanner can read it without a separate code path.

Staging verifies the sequence name, length, and the per-chromosome SHA-256 that
was recorded when the immutable plan was prepared. Thus the local representation
may have different line wrapping while the sequence identity remains exact.

## Chromosome-local Slurm execution

`run_genome_scan_slurm_chromosome.sh` maps one Slurm array element to one
planned chromosome. It creates a unique node-local job directory, stages and
verifies that chromosome once, and then runs all still-incomplete motif batches
against the local copy. On requeue, scratch is staged again; already promoted
batches are recognized from their exact task records and reused.

The default writes scanner output to unique durable staging below the run root
on `/data`. `--scratch-task-output` is optional: it writes and validates each
batch locally, streams every verified file once into a new durable attempt, and
then atomically promotes that attempt within `/data`. Local copies are not
removed by these scripts; they remain until Haumea applies its normal post-job
scratch cleanup. This option therefore requires enough scratch for all output
created during the allocation.

`--batch-workers` permits several independent batch processes to share the
staged FASTA. Start with one. Benchmark two and four on representative large
and small chromosomes before increasing production concurrency, and increase
`--cpus-per-task` and memory correspondingly. The submission helper does the
CPU adjustment automatically.

```sh
scripts/submit_genome_scan_slurm.sh \
  --run-root "$RUN" \
  --scanner "$SOURCE/pssm_scan_parquet" \
  --max-concurrent 20 --batch-workers 1 \
  --memory 16G --time 0-08:00:00 --finalize --dry-run
```

Remove `--dry-run` only after inspecting the rendered submission. Durable data
must remain below `/data/sm718`; `/scratch` is disposable input/work storage.

## Progress reporting

`pssm_scan` and `pssm_scan_parquet` retain the signal-safe scanner status
contract:

```sh
scancel --signal=USR1 --batch JOB_ID_ARRAY_INDEX
```

The scanner reports phase, global motif and chromosome indices, orientation,
position, emitted hits, and elapsed time. It also reports requests received
around indexed FASTA loading and Parquet closing. A chromosome supervisor sends
one request to one active batch, preventing a single signal from producing a
burst of competing lines.

For score-floor sensitivity analyses, `--context-maxima OUTPUT.parquet` with
`--regions ANCHORS.tsv` writes one maximum per anchor directly to Parquet. The
mode requires one motif, BED coordinates, an explicit recoverable score floor,
and an Arrow build. Its `--context-flank` uses signed interval-edge distance,
so overlapping motif intervals have negative distance and abutting intervals
have distance zero. SIGUSR1 reports anchor progress and the count of recoverable
maxima without interrupting the scan.

The finalizer and checksum auditor use their own one-line manager report with
phase, files and bytes complete/total, throughput, elapsed time, ETA, task, and
path. `status` also reports the most recently persisted checksum-audit state.

## Fast finalization and full verification

Task execution already parses every newly written Parquet file, validates its
rows and configuration, and computes SHA-256 before promotion. Normal
`finalize` therefore checks the exact expected task paths, provenance, byte
sizes, and recorded modification times. It does not reread all payload bytes.
Metadata and the DuckDB catalog are published first; `manifest.json` is still
published last.

A full payload reread is explicit:

```sh
scripts/manage_genome_scan.py verify \
  --run-root "$RUN" --checksums \
  --max-read-mib-per-second 125
```

Verification uses only inventory paths, writes append-only checkpoint segments,
and resumes completed files. `--max-files N` intentionally stops after at most
N new files and leaves a valid partial state. The final
`verification/checksum_audit.json` is written only after every current file has
the recorded digest. A same-size mutation is detected even if its timestamp is
restored.

## DuckDB catalog and exact hit queries

Small metadata tables are copied into the DuckDB catalog. Consequently metadata
queries work from any current directory and the catalog remains valid when the
package is moved as a unit. Catalog construction does not open hit payloads.

The former package-wide `motif_hit` wildcard was removed because creating that
view caused DuckDB to inspect every Parquet footer. The catalog now provides
the table macro `motif_hit_files(file_paths)`. `query_genome_scan.py` first
selects exact paths from `scan_file_inventory`, resolves them relative to the
package, and supplies only that bounded list to the macro.

```sh
scripts/query_genome_scan.py summary --package "$RUN/package"

scripts/query_genome_scan.py hits --package "$RUN/package" \
  --motif MA0861.2 --chrom 1 --minimum-score 0 --limit 20
```

An alternative catalog can be constructed without opening or changing any hit
file:

```sh
scripts/manage_genome_scan.py build-catalog \
  --run-root "$RUN" \
  --output "$RUN/package/jaspar_genome_scan.rebuilt.duckdb"
```

For an older manifest, pass that alternate database explicitly to the query
tool with `--database`.

## Build provenance

`pssm_scan --version-json` reports the program version, source commit, dirty
state, compiler, C++ language level, effective compile flags, LTO state,
assertion mode, Parquet support, and Arrow/Parquet versions. Every new task
records this object plus the executable SHA-256. By default the coordinator
rejects a dirty executable or a commit that differs from the immutable plan.
The explicit `--allow-scanner-provenance-mismatch` escape hatch records the
mismatch and is intended for diagnostics, not production.

Build the production scanner only after committing or checking out the planned
source revision:

```sh
make LTO=1 CXXOPTFLAGS="-march=x86-64-v3" pssm_scan_parquet
./pssm_scan_parquet --version-json
```

Scanner targets deliberately compile all of their translation units on every
explicit build. A change in optimization flags therefore cannot reuse objects
from a development build, and embedded Git provenance cannot remain stale.

## Physical-layout benchmark

No existing human payload is repacked. Before choosing a layout for mouse or
rat, compare the current two strand-partitioned files with a temporary combined
file for representative motif/chromosome pairs:

```sh
scripts/benchmark_sparse_layout.py \
  --package "$RUN/package" --motif MA0861.2 --chrom 1 \
  --output-dir "$SCRATCH/layout-MA0861.2-chr1"
```

The tool reports compressed bytes and median simple-query latency. A future
layout change requires a material improvement across several motif densities;
one favorable tiny example is insufficient.

## Maintenance acceptance gate

- `make check` and `make check-duckdb` pass.
- Scratch staging reproduces the planned sequence hash and builds a valid local
  index.
- Parallel batch children retain independent atomic promotion and retry.
- SIGUSR1 leaves scanner and checksum-audit processes alive and emits one line.
- Catalog rebuilding succeeds while hit payloads are unavailable.
- Exact hit queries work from outside the package directory.
- Partial checksum verification resumes and detects same-size corruption.
- The generated catalog inventory and hit totals equal the immutable manifest.

OpenMP inside one scanner, compressed BED output, arbitrary sparse region lists,
context features, gene annotation, CUT&RUN joins, and ML remain outside this
maintenance gate.

The downstream TP73 pair-layer contract in `tests/test_motif_context.sh` also
asserts unique configuration-aware `anchor_hit_id` values, exact additive
partitioning of partner loci by orientation, exhaustive pair classes, and
singleton partner-score nullability. Anchor and partner score floors are
recorded in both `context_run_config` and `tp73_pair_feature`.

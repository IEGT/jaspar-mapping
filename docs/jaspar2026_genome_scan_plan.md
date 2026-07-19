# JASPAR 2026 Whole-Genome Scan Preflight

## Scope and frozen decisions

This is an execution plan, not an executed run. It scans canonical human
GRCh38 sequence regions `1`-`22`, `X`, `Y`, and `MT`, both orientations, with
JASPAR 2026 CORE non-redundant motifs and direct sparse Parquet output.

- TP73 `MA0861.2`: retain `score >= -5`.
- Every other motif: retain `score >= -1`.
- Score: `log2_relative_risk`, pseudocount `1`, uniform A/C/G/T background.
- Coordinates: BED 0-based half-open; windows crossing `N` are skipped.
- Matched sequence is omitted. It remains recoverable from the pinned genome.
- TP73 is excluded from every default batch, so it cannot also be scanned at
  `-1` by accident.

The thresholds are a storage policy, not declarations of binding. TP73 is kept
deeper because its CUT&RUN calibration supports downstream sensitivity analyses;
the `-1` floor for motifs without matched CUT&RUN remains provisional.

## Identity and completeness contract

`genome_id` and `motif_set_id` are explicit partition keys. `genome` records
taxon ID, assembly name/accession, Ensembl release, FASTA SHA-256, `.fai`
SHA-256, and a sequence-set SHA-256 derived from per-sequence checksums.
`sequence_region` records order, length, checksum, and whether the region was
included. This prevents human, mouse, and rat coordinates from sharing a key.
Before each task, the manager re-hashes the small JASPAR and FASTA-index files
and checks the planned FASTA byte size and nanosecond modification time. The
full FASTA SHA-256 is computed once during planning rather than rereading the
multi-gigabyte genome once per array task.

The final package also contains `scan_run`, `scan_threshold_policy`,
`scan_task`, and `scan_file_inventory`. The file inventory is authoritative:
for each motif/chromosome/orientation it records expected windows, skipped-N
windows, sentinel-score windows, threshold/PWM rejections, emitted hits, bytes,
SHA-256, state, Slurm IDs/restart count, scanner checksum, and source commit.
Zero-hit files are retained and inventoried.

## Haumea preparation

Use one dedicated tree under `/data/sm718`; no generated file belongs in `$HOME`
or the source checkout. Pull/build from the public repository on Haumea rather
than transferring a laptop build or data product.

```sh
readonly SOURCE=/data/sm718/GitHub/jaspar-mapping
readonly RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_sparse_v1
readonly GENOME=$SOURCE/Homo_sapiens.GRCh38.dna.primary_assembly.fasta
readonly FAI=$GENOME.fai
readonly JASPAR=$SOURCE/JASPAR2026_CORE_non-redundant_pfms_jaspar.txt

cd "$SOURCE"
git pull --ff-only https://github.com/IEGT/jaspar-mapping.git main
SOURCE_COMMIT=$(git rev-parse HEAD)
make pssm_scan_parquet
df -h /data/sm718
```

Prepare the immutable plan. The 500 GiB reserve is deliberately conservative;
the manager checks it before planning and before every task. The existing TP73
chr1 result (7.731 MiB at `-1`) scales to roughly 0.1 GiB genome-wide for that
motif, but 2,633 motifs have very different information content, so a chr1 pilot
over all motifs remains the better disk forecast before submission.

```sh
chrom_args=()
for chrom in {1..22} X Y MT; do chrom_args+=(--chrom "$chrom"); done

scripts/manage_genome_scan.py prepare \
  --run-root "$RUN" --run-id jaspar2026_grch38_sparse_v1 \
  --source-commit "$SOURCE_COMMIT" \
  --genome "$GENOME" --fasta-index "$FAI" --jaspar "$JASPAR" \
  --genome-id homo_sapiens_grch38_ensembl113_primary \
  --motif-set-id jaspar2026_core_nonredundant \
  --taxon-id 9606 --species 'Homo sapiens' \
  --assembly-name GRCh38 --assembly-accession GCA_000001405.29 \
  --ensembl-release 113 \
  --fasta-url 'https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz' \
  --jaspar-url 'https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_non-redundant_pfms_jaspar.txt' \
  --special-motif MA0861.2 --special-minimum-score -5 \
  --default-minimum-score -1 --motif-batch-size 64 \
  --score-mode log2_relative_risk --pseudocount 1 \
  --minimum-free-gib 500 "${chrom_args[@]}"
```

With the current 2,633-motif file this produces 43 disjoint motif batches and
1,075 chromosome/batch tasks. Tasks are ordered batch-first and then by
chromosome, so the first concurrent array wave does not make every process read
chromosome 1. Confirm `TASK_COUNT` from the command rather than hardcoding it if
the source file changes.

## Slurm execution

Run at most 20 concurrent one-core tasks. Request 16 GiB per task even though
the scanner's expected resident set is much smaller: this leaves ample room for
Arrow and a DuckDB validation database capped at 12 GiB in memory. The
validation database is transient; motif hits remain direct Parquet output, and
the final durable DuckDB file is a small catalog of views and inventories rather
than a materialized copy of all hits.

All durable staging and promoted results stay below `/data/sm718`. When a node
offers writable `/scratch`, the wrapper uses a unique job/restart directory
there only for `TMPDIR` and possible DuckDB spill. It never relies on scratch
for a completed task, and no upfront scratch distribution is required. Requeue
creates a new durable attempt directory. A task is promoted only after every
motif/orientation file in its batch passes row, coordinate, score, configuration,
checksum, and window-accounting validation.

```sh
mkdir -p "$RUN/logs"
TASK_COUNT=$(python3 - "$RUN/plan/scan_plan.json" <<'PY'
import json, sys
print(len(json.load(open(sys.argv[1], encoding='utf-8'))['tasks']))
PY
)

JOB_ID=$(sbatch --parsable \
  --account=cluster --partition=requeue --requeue \
  --array="0-$((TASK_COUNT - 1))%20" \
  --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G --time=1-00:00:00 \
  --chdir="$SOURCE" \
  --output="$RUN/logs/slurm-%A_%a.out" \
  --error="$RUN/logs/slurm-%A_%a.err" \
  scripts/run_genome_scan_slurm_task.sh \
  --run-root "$RUN" --scanner "$SOURCE/pssm_scan_parquet" \
  --duckdb-memory-limit 12GB)
printf 'Submitted %s\n' "$JOB_ID"
```

No task removes or replaces a prior attempt. Failed attempts remain below
`$RUN/staging/TASK_ID/`; successful task packages are atomically renamed below
`$RUN/package/task_data/task_id=TASK_ID/`. A retry recognizes and reuses only a
promoted task whose plan, completion record, file sizes, and SHA-256 checksums
all agree.

## Status and SIGUSR1

The coordinator provides aggregate status without opening Parquet:

```sh
scripts/manage_genome_scan.py status --run-root "$RUN"
```

Request one live progress line from an array element with Slurm's signal
delivery (the wrapper `exec`s the coordinator, which forwards `SIGUSR1` to the
scanner):

```sh
scancel --signal=USR1 --batch "${JOB_ID}_7"
```

The process remains alive. The line reports phase, motif and chromosome
index/count within that scanner invocation, motif/chromosome/orientation,
emitted-hit count, elapsed time, current alignment start, and percent complete.
The task ID and global task index remain in the Slurm log and immutable task
plan. Phases also identify FASTA loading and Parquet closing, where an ordinary
scan-position percentage cannot be reported.
The production coordinator deliberately omits `-v`, so it does not flood
Slurm logs with continuously redrawn progress bars.

## Finalization

Only finalize after `tasks_complete == tasks_total`:

```sh
scripts/manage_genome_scan.py finalize --run-root "$RUN"
cd "$RUN/package"
duckdb jaspar_genome_scan.duckdb \
  -c 'SELECT * FROM scan_inventory_summary;'
duckdb jaspar_genome_scan.duckdb \
  -c 'SELECT * FROM scan_task_completeness WHERE NOT complete;'
```

`manifest.json` is written last and is the package completion marker. If a
finalizer is interrupted before that marker, a retry accepts existing metadata
only after semantic equality checks and accepts an existing DuckDB index only
after its run/file/hit inventory agrees; it never deletes a partial attempt.
The raw query contract is
[`../sql/genome_scan_schema.sql`](../sql/genome_scan_schema.sql).

## Mouse, rat, and synteny

Run each species as its own package with a distinct `genome_id`, taxon ID,
assembly accession, Ensembl release, FASTA/sequence-set checksums, and
`sequence_region` rows. The same `motif_set_id` may be reused when the exact
JASPAR file checksum agrees. Raw `motif_hit` does not otherwise change.

Synteny is a derived layer, not a replacement coordinate system. It should keep
the source and target `genome_id`, intervals on both assemblies, mapping
orientation, chain/alignment source and checksum, mapping multiplicity,
reciprocal status, aligned-base fraction, and quality. An
`orthologous_motif_hit` bridge can then record whether the lifted interval
contains the same motif, the source and target `motif_set_id`, its target
score/orientation, score delta, and whether the motif sequence itself is
conserved. One-to-many and unmapped cases must remain rows rather than being
silently discarded.

Thresholds should initially be treated as species-specific policies even when
the numeric human floor is reused. Conservation is then an additional feature:
it must not be used to overwrite the species' own sequence score or CUT&RUN
calibration. Held-out-species and held-out-synteny-block evaluation are natural
tests once matched mouse or rat occupancy/expression data are available.

### Initial provider: Ensembl Compara

Use Ensembl Compara as the first and canonical mapping source. Pin the Compara
release independently in every imported row. First check the archived release
113 pairwise inventories because that matches the existing Ensembl-113 genome
identities. If a required pairwise product is absent or materially improved in
a later release, import that explicitly versioned release only after confirming
the source and target assembly accessions. Do not silently mix data from the
current REST server into the package merely because assemblies have the same
display names.

Import two related but distinct products:

1. Compara `SYNTENY` regions provide coarse conserved-order blocks for
   stratification and neighborhood context. They are not precise enough to
   declare a short motif orthologous.
2. Human-mouse and human-rat pairwise `LASTZ_NET` genomic alignments provide
   the coordinate mapping used to project individual motif intervals. Prefer
   the release-pinned Ensembl MAF/Compara dump for the bulk import. The REST
   `alignment/region` endpoint is useful for small fixtures and independent
   spot checks, but its 10 Mb request limit makes it a poor whole-genome
   transport.

Store the Compara release/division, method-link type, method-link species-set
ID, Ensembl synteny or genomic-alignment ID when available, source URL/file and
SHA-256, and both assembly-specific genome IDs. Convert Ensembl's 1-based
inclusive coordinates to BED 0-based half-open coordinates at the importer
boundary. Keep pairwise CIGAR information so projections crossing a gap can be
reported as split or gapped rather than forced into an apparently contiguous
orthologous hit. The corresponding optional views are defined in
[`../sql/cross_species_schema.sql`](../sql/cross_species_schema.sql).

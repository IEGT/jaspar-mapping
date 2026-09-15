# Regulatory TFBS Tags And Bounded Exports

## Scope

`scripts/export_regulatory_tfbs.py` provides the general-purpose filtering
route requested for GENtle. It selects **every retained motif/orientation
record** with positive interval overlap with an Ensembl regulatory feature
or a versioned TSS window. It does not require TP73, a TP73 local peak,
CUT&RUN support, or selection as the strongest cofactor in a distance band.
It does not rescan DNA, rewrite the original atlas or refit association models.

This is distinct from Glen's 119 MB promoter-cofactor collaboration package.
That package carries TP73-conditioned statistics and strongest-per-band
positions. The new route carries individually retained hits with geometric
tags, without inventing TP73 association statistics for them.

Implementation is locally tested and the eight-motif chromosome-1 pilot below
has completed on Haumea. The subsequent all-motif, genome-wide **intersection**
export has been submitted, as recorded below; its final size is not yet known.
The earlier union pilot is a separate package and must not be substituted for
this new delivery.

## Genome-Wide Intersection For Glen

The subsequent requested export changes **both** the geometry and the scope
of the eight-motif pilot. It is not an enlargement of that pilot's union:

- All 2,633 scanned JASPAR matrices, regardless of source species; no TP73
  anchor requirement, cofactor ranking, strongest-hit selection or CUT&RUN
  selection.
- Every canonical chromosome in the source scan: 1-22, X, Y and MT.
- All physical TSSs from the pinned Ensembl 113 GTF, including alternative
  transcripts and shared gene starts. No representative-transcript filter.
- Strand-aware **700 bp upstream / 300 bp downstream** TSS windows.
- A retained motif interval must positively overlap the **intersection** of
  one such window and an audited Ensembl regulatory extent. A motif bridging
  two merely abutting/disjoint annotation intervals does not qualify.
- Retain every stored score: TP73 MA0861.2 down to -5, other matrices down to
  -1 in the permissive v3 atlas. No informative/density cutoff is imposed.

As in the existing promoter definition, offsets include the TSS base itself.
For zero-based TSS base `t`, the half-open windows are `[t-700,t+301)` on `+`
and `[t-300,t+701)` on `-`: 1,001 bases before clipping to chromosome bounds.
The minus-strand TSS is GTF transcript end minus one, not the exclusive BED
end. Motif intervals are retained whole, not clipped; membership uses positive
overlap rather than full containment. Regulatory extents use annotated
extended promoters where available and the native intervals of other feature
types, retaining separate core/extension tags. This remains a reference
annotation-based selection, not proof of regulatory activity in a cell line.

`manage_regulatory_tfbs.py prepare` builds a separate, whole-genome annotation
package from the original GTF and the **complete** audited regulatory feature
package. It does not depend on the autosome-only TP73 membership run. The
annotation contains physical TSS/window dimensions and a separate
`transcript_tss` bridge with gene IDs/names. Shared promoters therefore never
multiply motif-hit records. All source SHA-256 digests, chromosome lengths,
annotation releases and the resolved window definition are pinned.

The regulatory source has features on X and Y, but none on MT. The plan
records MT as `known_empty_intersection` for every matrix, rather than silently
omitting that chromosome or pretending mitochondrial sequence was not scanned.
An annotated chromosome without GTF TSS coverage fails preparation. GTF
contigs outside the canonical source scan are not claimed as exported.

### Production Commands

Use a fresh dedicated directory under `/data/sm718` and an immutable source
checkout fetched through Git. On the login node, create the run directory and
pin the clean source once:

```bash
RUNS=/data/sm718/jaspar_mapping_runs
RUN="$RUNS/glen_genome_regulatory_tss700_300_v1"
mkdir -p "$RUN"
python3 "$SOURCE/scripts/manage_regulatory_tfbs.py" pin-source --output "$RUN/source_provenance.json"
```

Preparation parses the GTF and runs on a compute node. Compute nodes do not
need Git; they verify the pinned source hashes instead. From a `requeue` allocation:

```bash
RUNS=/data/sm718/jaspar_mapping_runs
RUN="$RUNS/glen_genome_regulatory_tss700_300_v1"
python3 "$SOURCE/scripts/manage_regulatory_tfbs.py" prepare \
  --run-root "$RUN" --package "$RUNS/jaspar2026_grch38_sparse_v3/package" \
  --features "$RUNS/ensembl_grch38_tp73_regulatory_20260909_v1/features" \
  --gtf /data/sm718/resources/ensembl/113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz \
  --upstream 700 --downstream 300 --batch-size 128 --duckdb "$DUCKDB" \
  --scratch-root /scratch/sm718 --source-provenance "$RUN/source_provenance.json"
```

After preparation and a real-data resource pilot, submit the frozen plan:

```bash
python3 "$SOURCE/scripts/manage_regulatory_tfbs.py" submit \
  --run-root "$RUN" --concurrent 20 --slurm-memory 24G \
  --memory-limit 16GB --slurm-time 04:00:00
```

With 25 chromosomes and 2,633 motifs this creates 525 chromosome/motif-batch
tasks, each with at most 128 matrices, plus an `afterok` finalizer. The source
checkout must be clean. Job IDs are recorded immediately in
`submissions.jsonl`; a repeated `submit` refuses to duplicate jobs. Preparation
and workers reuse verified completed state after preemption, with new private
attempt/scratch directories for incomplete work. No source or failed attempt
is deleted.

Each motif's two exact input files are copied to node-local scratch and
checksum-verified. The original atlas is unchanged. Inputs over ten million
orientation records are queried in **5 Mb start-owned tiles**, so a hit
crossing a tile boundary occurs exactly once. Only the filtered tile results
are combined into one sorted Parquet per chromosome/motif. No intermediate
BED/TSV is created. Completed motifs are independently reusable on requeue.
Haumea owns post-job scratch cleanup.

There is **no row or delivery-byte cap** in this production run:
`--max-rows 0 --max-output-bytes 0 --source-floor`. The 50 GiB free-space
reserve, DuckDB memory ceiling and bounded input working sets remain safety
constraints. Workers stop rather than publish truncated data when a safety
constraint fails. Full transfer to a laptop waits for measured output sizes
and sufficient destination capacity.

### Delivery Layout

The finalizer verifies every expected task and motif, then publishes `final/`:

- `hits/chrom=<chrom>/<motif>.parquet`: every qualifying orientation record.
- `file_inventory.json` and `regulatory_tfbs.duckdb`: exact relative file
  paths, counts, bytes, SHA-256 digests, original scan inventory, genome and
  motif metadata. Known-empty chromosomes have explicit zero rows, no payload.
- `annotation/`: shared physical TSS/windows, transcript/gene ownership and
  original regulatory dimensions/provenance.
- `schema.sql` and `manifest.json`: coordinate, selection and coverage
  contracts. Package kind is `genome_regulatory_tfbs_subset`.

Final hit files are hardlinks to completed exports within the same `/data`
filesystem, avoiding a second payload-sized copy. They become ordinary files
when transferred. The package index does not bind a genome-wide wildcard:
select exact paths from `file_inventory`, then pass them to its `motif_hits`
table macro. For example, from the final package directory:

```sql
-- Open regulatory_tfbs.duckdb read-only for metadata or file-based queries.
SELECT path, rows, bytes FROM file_inventory
WHERE chrom='1' AND motif_id='MA1961.2';
SELECT * FROM motif_hits(['hits/chrom=1/MA1961.2.parquet'])
WHERE start < 1790000 AND "end" > 1770000 AND score >= 0;
```

The full atlas of unrestricted genomic matches is still a different provider.
GENtle must respect this package's declared regulatory/TSS intersection and
original motif floors, not interpret absence outside that subset as absence of
a genomic sequence match. The original pilot and its checksum below remain
unchanged.

### Haumea Execution: 2026-09-15

Production source is immutable commit
`0cfb8f14ed29b9d8dcc4633f5771591fe415453f`, fetched from GitHub. The dedicated
run directory is:

```text
/data/sm718/jaspar_mapping_runs/glen_genome_regulatory_tss700_300_20260915_v1
```

Setup/pilot job **5795221** completed `0:0` in **2:03**, with Slurm maximum RSS
2,829,636 KiB under its 24 GB allocation. It produced **320,404 physical TSS
windows** across all 25 scanned chromosomes. X has 9,420 windows and 16,138
regulatory features; Y has 1,568 windows and 676 regulatory features. MT has
37 TSS windows but zero regulatory features, so its intersection is explicitly
empty. The earlier setup attempt 5795219 stopped at source-version capture
because compute nodes lack Git; no motif exports were produced by that attempt.
Its logs and staging remain available. Source pinning was moved to the login
node and the full synthetic pipeline tested with Git absent from `PATH`.

Four real-data pilots verified the final intersection geometry and source-floor
configuration before production submission:

| Chromosome | Matrix | Retained records | Parquet bytes | Source floor |
|---|---|---:|---:|---:|
| 1 | PATZ1 MA1961.2 | 154,173 | 1,073,515 | -1 |
| 2 | MA0013.1 (dense stress test) | 1,084,334 | 4,760,685 | -1 |
| Y | POU2F2 MA0507.3 | 26 | 5,791 | -1 |
| 1 | TP73 MA0861.2 | 34,185 | 360,872 | -5 |

The dense input contained 153,954,062 records. It was processed in 49 disjoint
5 Mb start-owned chunks before combining the selected records. These pilots
are not a statistical subsample from which to infer the final archive size.

The full export was submitted as array **5795222**, tasks `0-524`, with a
20-task concurrency limit, two CPUs and 24 GB per allocation, 16 GB DuckDB
memory and a four-hour job limit. Finalizer **5795223** has an `afterok`
dependency on the complete array. This record documents **submission**, not
completed genome-wide output. Authoritative completion will be the published
`final/manifest.json`, with every motif/chromosome accounted for.

## Tags And Coordinates

Each hit keeps its full BED 0-based half-open interval, original fractional
score and orientation. Positive overlap means `hit.start < region.end AND
hit.end > region.start`. Mere abutment is excluded, and neither query-region
nor annotation boundaries clip the stored hit. Overlapping features and
shared promoters do not duplicate hits. Both retained orientations at one
physical span remain separate records.

The `regulation_tags` integer is an OR of these independent bits; the export
also contains the corresponding readable Boolean columns:

| Bit | Boolean column | Definition |
|---:|---|---|
| 1 | `overlaps_promoter_core` | Audited Ensembl core promoter |
| 2 | `overlaps_promoter_extended` | Audited Ensembl extended promoter |
| 4 | `overlaps_enhancer` | Ensembl enhancer |
| 8 | `overlaps_open_chromatin` | Ensembl open chromatin region |
| 16 | `overlaps_ctcf` | Ensembl CTCF binding-site annotation |
| 32 | `overlaps_emar` | Ensembl EMAR annotation |
| 64 | `overlaps_other_regulatory` | Other feature type retained in the source |
| 128 | `overlaps_tss_window` | Existing versioned, strand-aware GTF promoter window |

`overlaps_regulatory` is any of the first seven bits. `regulation_candidate`
is any bit. These are **annotation-based selection tags, not measured
activity, protein binding, effect size or an ordered biological ranking**.
Enhancers are not assigned a target gene by proximity.
An absent optional promoter extension contributes no extension tag; its core
is still considered. Missing extensions are never invented from core bounds.

The standalone default `--scope regulatory_or_tss` is a **union**, not an intersection.
`--scope promoter_or_tss` retains the union of the promoter core, promoter
extension and TSS window. Any individual tag is also a selectable scope,
as is `--scope regulatory` without the TSS-only sites. All tags remain in
the output, even when selection uses only one of them.
`--scope regulatory_and_tss` additionally requires a positive shared interval
and exports `overlaps_regulatory_tss_intersection`. New individual exports use
schema version 2; the older pilot's schema-1 packages remain unchanged.

The current completed regulatory run supplies whole-chromosome TSS windows
under `tss_upstream_2000_downstream_500_v1` (Ensembl GTF 113). Their already
computed genomic bounds encode strand orientation; this tool does not
reconstruct them or restrict them to promoters with TP73. Ensembl regulatory
release 2025-05 is independently pinned by its passed GFF/BigBed audit.
The annotation plan, completion marker, genome identity, coordinate audit and
selected chromosome's dimension hashes must agree before any hit is queried.
The TP73 membership table is never read.

## Source Floors And Coverage

Use the completed permissive source:

```text
/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_sparse_v3/package
```

It contains 2,633 matrices from all source taxa, not 2,633 human TF genes.
For a biological human panel, choose exact accessions using the normalized
JASPAR catalog and retain source-species provenance. `--all-motifs` explicitly
selects the entire catalog, including non-vertebrate patterns.

The default inclusive query threshold is `-1`; score zero is an optional
selection for historical comparison. The source's per-motif retention floor
is kept separately as `source_minimum_score`. A source floor above the
requested query threshold, or additional PWM-relative filtering, is rejected.
Thus the density-limited atlas cannot silently replace the permissive source
and lose abundant factors. TP73 alone can be requested at `-5`; other matrices
retained only down to `-1` cannot support that request.

One standalone invocation selects one exact chromosome and an explicit motif panel
(repeat `--motif`) or `--all-motifs`. Supply BED `--start`/`--end`, or explicitly
request `--whole-chromosome`. The **legacy TP73** annotation plan covers autosomes
1-22. X/Y/MT and chromosome aliases are rejected rather than becoming false
negative annotations when using `--regulatory-run`. The independent
`--annotation-package` production path above covers all source chromosomes.
Empty valid subsets publish a zero-row Parquet plus
their scope manifest; absent motif/chromosome coverage is an error.

## Small Haumea Export

Make the tested source available through Git before running this command.
Use a Slurm `requeue` allocation, for example 2 CPUs and 16 GB for a bounded
panel. This Python/DuckDB CLI tool supports Linux and macOS (POSIX locks).
The exporter defaults to two DuckDB threads, 1 GB memory and no disk
spill. It reads only exact inventory paths and resolves legacy catalog metadata
views from the source package directory, independently of the caller's directory.
Whole-chromosome/all-motif counts
are real payload work and must not run on the login node.

From the allocation, with the installed DuckDB on `PATH`:

```sh
SOURCE=/data/sm718/GitHub/jaspar-mapping
RUNS=/data/sm718/jaspar_mapping_runs
python3 "$SOURCE/scripts/export_regulatory_tfbs.py" \
  --package "$RUNS/jaspar2026_grch38_sparse_v3/package" \
  --regulatory-run "$RUNS/ensembl_grch38_tp73_regulatory_20260909_v1" \
  --motif MA0653.1 --motif MA1961.2 --motif MA0079.5 \
  --chrom 1 --start 1770000 --end 1790000 \
  --scope regulatory_or_tss --minimum-score -1 \
  --max-rows 100000 --max-output-bytes 100000000 \
  --output "$RUNS/glen_general_regulatory_pilot_v1/chr1_irf9_patz1_sp1"
```

This is a region query, not a claim of an experimental cofactor effect in
that region. Outputs are direct ZSTD Parquet, never intermediate BED/TSV.
`manifest.json` records the source inventory and its recorded checksums,
verified annotation hashes, score threshold, genomic/motif selection, tag
definitions and measured bytes. Large scan payloads are checked by exact
inventory and byte size, not rehashed in full; this distinction is explicit.
The saved SQL is an audit record. `schema.sql` exposes the completed Parquet
as `regulatory_tfbs` when opened from the export directory.

The original scan remains on `/data`; no FASTA is needed. Add
`--scratch-directory "$JOB_SCRATCH"` inside an allocation to copy only the
selected hit files and annotation dimensions into a unique local directory.
Their bytes and SHA-256 digests are verified before querying; the manifest
records `scratch_staged_and_verified`. The small scan catalog remains read-only
on `/data`. The tool does not remove scratch copies: Haumea's post-job cleanup
owns their lifetime. It does not submit cluster jobs itself.

## Sizing And Further Filtering

Add `--count-only` and choose a new output directory to write a small
`counts.parquet`: per-motif/tag-combination orientation-record and distinct
physical-locus counts. This still reads the selected hit payloads, but does
not write the full hit list. Counts cannot predict compressed bytes exactly;
measure representative actual exports before choosing a delivery panel.

For an exported hit package, from its directory:

```sh
duckdb -no-init :memory: -c "
  SELECT motif_id, count(*) AS orientation_records,
         count(DISTINCT (chrom,start,\"end\")) AS physical_loci
  FROM read_parquet('motif_hits.parquet')
  WHERE score >= 0 AND (overlaps_promoter_extended OR overlaps_tss_window)
  GROUP BY motif_id ORDER BY motif_id;"
```

Default publication ceilings are one million rows and one billion Parquet
bytes, with a further 1 GiB free-space reserve. Limits fail explicitly; they
never publish a silently truncated selection. The writer is checked every
250 ms for byte-budget/free-space violations, with a final size check before
publication. Temporary overshoot is possible between checks and while closing
Parquet; these are application safeguards, not filesystem quotas.
Existing outputs are refused unless `--resume` verifies an identical completed
export: selection, inputs, inventory, builder hash, publication caps and output
checksums must agree. A requeued job can therefore reuse completed motif exports
and retry interrupted ones in new attempt directories. A kernel advisory lock
prevents concurrent writers and releases automatically when a process exits or
is killed; its lock file remains harmlessly in place. Failed attempts are
preserved for inspection; these tools delete nothing.

A delivered package is complete only for its declared motif/chromosome/span,
annotation selection and score threshold. It must not be attached as a full
genome-scan provider or opened as the different TP73 collaborator-package
schema. GENtle can query its tagged Parquet; native attachment of this new
package kind is a separate, small adapter task. Gene/TSS ownership remains in
the pinned dimensions rather than multiplying motif-hit rows.

## Completed Chromosome-1 Pilot: 2026-09-15

Slurm job **5795208** completed with exit status `0:0` in **28 seconds** on
node145 (`cluster` account, `requeue` partition, 2 CPUs, 16 GB requested).
Slurm recorded maximum RSS of 5,383,544 KiB. The immutable source checkout was
commit `a31ce6a170c850f6dc0792501450db1e7b6e6987`, fetched from GitHub on Haumea;
no local source or input files were copied to the cluster.

The dedicated durable run directory is:

```text
/data/sm718/jaspar_mapping_runs/general_regulatory_tfbs_chr1_pilot_20260915_v1
```

It contains `source_commit.txt`, `job_id.txt`, the Slurm-captured
`submitted_batch.sh`, `slurm-5795208.log`, and the eight completed packages at
`motifs/<accession>/`. Each motif was exported over **all chromosome 1** with
`--scope regulatory_or_tss --minimum-score -1`, not just CD44/TGFB1 regions.
The source score mode is `log2_relative_risk`, additive pseudocount 1 per base,
uniform A/C/G/T background, and `n_policy=skip`.

The sequential motif commands used `--resume`, verified scratch staging,
10,000,000-row and 1,000,000,000-byte publication caps per motif, a 10 GiB
free-space reserve, two DuckDB threads, 8 GB DuckDB memory with spill disabled,
and a 3,300-second per-call timeout. All eight manifests confirm complete
exports and checksum-verified scratch staging. The input inventory comprised
16 exact files, 13,140,420 orientation records and 108,666,855 Parquet bytes.

| Motif | Accession | Exported orientation records | Physical loci | Parquet bytes |
|---|---|---:|---:|---:|
| IRF9 | MA0653.1 | 9,281 | 9,281 | 117,539 |
| PATZ1 | MA1961.2 | 573,301 | 573,301 | 4,011,684 |
| SP1 | MA0079.5 | 352,823 | 352,823 | 1,898,070 |
| E2F1 | MA0024.3 | 286,931 | 172,024 | 2,378,334 |
| REST | MA0138.3 | 83,504 | 83,467 | 973,757 |
| POU2F2 | MA0507.3 | 81,645 | 81,483 | 892,547 |
| TFAP2C | MA0524.3 | 1,933,785 | 1,715,658 | 15,775,367 |
| TP53 | MA0106.3 | 57,021 | 36,595 | 586,980 |
| **Total** | | **3,378,291** | **3,024,632** | **26,634,278** |

Physical loci are distinct `(start,end)` spans **within each motif**; the total
does not collapse coinciding spans across different matrices. The original
orientations and scores remain present. Local inspection of all output records
found zero duplicate `(motif_id,start,end,strand)` keys and zero invalid
coordinates, chromosomes, strands, scores below -1, or empty annotation tags.
All manifest-listed file sizes and SHA-256 digests passed after transfer.

Two useful filters measured directly on these outputs:

- Promoter-core/extension **or** TSS-window membership: 2,417,266 records.
- Score >= 0 across the full regulatory/TSS union: 2,612,628 records.

These are separate filters, not a combined count. Both can be applied locally
without rescanning DNA or transferring more atlas data.

The shareable ZIP is **25,982,743 bytes** (25.98 MB), with SHA-256
`fd7d3108c28c1c182b2b864d7f843500c72f2b6496a1b7f7880843e84a7bcead`.
It is present in the remote run directory and locally at
`dry_runs/general_regulatory_tfbs_chr1_pilot_20260915_v1/Glen_regulatory_TFBS_chr1_pilot_20260915.zip`.
The local extracted packages are under the adjacent `package/motifs/` directory.
Each package's `schema.sql` is portable when opened from its own directory;
`query.sql` is the production audit and still references Haumea/scratch inputs.

This bundle is suitable for early Parquet inspection and a GENtle adapter for
the new `regulatory_tfbs_subset` kind. It is not the existing TP73-conditioned
collaborator schema or a complete whole-genome TFBS provider. The small size
of this selected panel must not be extrapolated to all 2,633 matrices without
broader measurements.

## Tests

```sh
python3 tests/test_regulatory_tfbs.py
bash tests/test_script_help.sh
```

The Python test constructs explicitly synthetic hit intervals (not an actual
PSSM scan), two chromosomes, two artificial matrices and strand-aware TSS
windows. It exercises real Parquet and the full exporter: abutment, overlap,
shared promoters, both orientations, TSS-only and non-promoter hits, score
floors, unknown coverage, caps, writer termination at the byte budget, empty
selections, input immutability, legacy relative catalog views, verified scratch
staging and restart-safe reuse.
No genome, JASPAR or experimental data download is required.

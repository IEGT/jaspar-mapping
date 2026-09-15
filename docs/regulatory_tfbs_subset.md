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
has completed on Haumea. No genome-wide export has been submitted by this
change, and the full regulatory subset's disk size remains unmeasured.
Broader production should use bounded exports and counting on Haumea, then a
deliberately sized collaborator panel, not a laptop download of the whole
permissive atlas.

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

The default `--scope regulatory_or_tss` is a **union**, not an intersection.
`--scope promoter_or_tss` retains the union of the promoter core, promoter
extension and TSS window. Any individual tag is also a selectable scope,
as is `--scope regulatory` without the TSS-only sites. All tags remain in
the output, even when selection uses only one of them.

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

One invocation selects one exact chromosome and an explicit motif panel
(repeat `--motif`) or `--all-motifs`. Supply BED `--start`/`--end`, or explicitly
request `--whole-chromosome`. The present annotation plan covers autosomes
1-22. X/Y/MT and chromosome aliases are rejected rather than becoming false
negative annotations. Empty valid subsets publish a zero-row Parquet plus
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

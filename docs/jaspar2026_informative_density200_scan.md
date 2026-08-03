# JASPAR 2026 Informative, Density-Capped Genome Scan

## Purpose

This run supersedes the provisional fixed-floor policy used by the completed
GRCh38 v3 scan. It retains a motif-specific score floor that is deliberately
permissive with respect to the chromosome-1 TP73-context analysis, but prevents
low-information motifs from dominating storage and downstream context counts.
It is an execution plan; no Slurm job is submitted by the preparation commands.

The scoring contract remains `log2_relative_risk`, pseudocount 1, uniform
A/C/G/T background, both orientations, skipped N-containing alignments, and
BED 0-based half-open coordinates.

## Threshold rule

For motif `m`, define:

```text
candidate(m) = min(informative_threshold(m), -1)
final(m)     = max(candidate(m), density_threshold_chr1(m))
```

An absent informative recommendation uses `-1`. TP73 `MA0861.2` is supplied as
an explicit `-5` override because its direct CUT&RUN calibration is separate
from the cofactor registry.

The density threshold is the lowest integer score at or above the candidate
for which chromosome 1 contains no more than
`floor(valid_alignment_starts / 200)` retained physical loci. At each genomic
alignment start the plus- and minus-orientation scores are collapsed to their
maximum before counting. Thus the 200 bp value is reciprocal locus density,
not the arithmetic mean of nearest-neighbor distances. It excludes starts
whose motif span crosses N and is motif-length-specific.

Production still writes orientation-specific rows. A physical locus can
therefore yield two rows when both orientations pass the same final threshold;
the density ceiling applies to distinct genomic alignment spans, not stored
row count. The chromosome-1 ceiling is transferred as a policy to the other
chromosomes. It does not assert that every chromosome independently has exactly
the same density.

## Stage 1: density calibration

Use one durable tree below `/data/sm718`. The worker copies chromosome 1 and the
small JASPAR matrix source once to its allocation's `/scratch`, then runs
checkpointed motif batches against those local copies. A requeued worker stages
the inputs again but verifies and reuses every already promoted batch. Failed
durable attempts are retained, and the cluster owns cleanup of the node-local
files.

```sh
readonly SOURCE=/data/sm718/GitHub/jaspar-mapping
readonly DENSITY_RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_locus_density200_v1
readonly GENOME=$SOURCE/Homo_sapiens.GRCh38.dna.primary_assembly.fasta
readonly FAI=$GENOME.fai
readonly JASPAR=$SOURCE/JASPAR2026_CORE_non-redundant_pfms_jaspar.txt
readonly INFORMATIVE=$SOURCE/thresholds/jaspar2026_grch38_chr1_tp73_context_v1.tsv
readonly NEGATIVE=/data/sm718/jaspar_mapping_runs/jaspar2026_chr1_tp73_context_negative_sensitivity_v1/final/negative_threshold_sensitivity/negative_threshold_sensitivity.tsv
readonly OVERRIDES=$SOURCE/thresholds/jaspar2026_grch38_scan_overrides_v1.tsv

cd "$SOURCE"
SOURCE_COMMIT=$(git rev-parse HEAD)
make LTO=1 CXXOPTFLAGS="-march=x86-64-v3" pssm_scan pssm_scan_parquet

scripts/manage_motif_density_calibration.py prepare \
  --run-root "$DENSITY_RUN" \
  --run-id jaspar2026_chr1_locus_density200_v1 \
  --source-commit "$SOURCE_COMMIT" \
  --genome "$GENOME" --fasta-index "$FAI" --jaspar "$JASPAR" \
  --informative-thresholds "$INFORMATIVE" \
  --negative-sensitivity "$NEGATIVE" \
  --override-thresholds "$OVERRIDES" \
  --threshold-set-id jaspar2026_grch38_chr1_informative_density200_v1 \
  --genome-id homo_sapiens_grch38_ensembl113_primary \
  --motif-set-id jaspar2026_core_nonredundant \
  --chrom 1 --motif-batch-size 64 \
  --default-minimum-score -1 --minimum-spacing-bp 200 \
  --score-mode log2_relative_risk --pseudocount 1 \
  --minimum-free-gib 500

scripts/submit_motif_density_calibration_slurm.sh \
  --run-root "$DENSITY_RUN" --scanner "$SOURCE/pssm_scan" \
  --source "$SOURCE" --account cluster --partition requeue \
  --batch-workers 4 --memory 32G --time 0-08:00:00 --dry-run
```

Inspect the two rendered `sbatch` commands before removing `--dry-run`. Check
progress without opening distribution files:

```sh
scripts/manage_motif_density_calibration.py status --run-root "$DENSITY_RUN"
scancel --signal=USR1 --batch DENSITY_JOB_ID
```

If a non-requeued failure or time limit stops the worker, inspect its log and
resubmit the same prepared run with `--resume`. The replacement worker verifies
and skips every completed batch before doing new work.

The finalizer writes:

```text
$DENSITY_RUN/final/distribution_manifest.tsv
$DENSITY_RUN/final/motif_scan_thresholds.tsv
$DENSITY_RUN/final/motif_scan_thresholds.json
$DENSITY_RUN/final/manifest.json
```

The TSV contains the informative, candidate, density, and final thresholds,
counts at both candidate and final scores, valid/skipped alignment starts,
estimated spacing, and the source distribution checksum for every motif. The
JSON pins its formula, input checksums, score configuration, and threshold-set
identity. It also reports aggregate chromosome-1 loci at the candidate and
final thresholds plus the conservative two-orientation row ceiling. Inspect
those values and current durable free space before preparing Stage 2:

```sh
python3 - "$DENSITY_RUN/final/motif_scan_thresholds.json" <<'PY'
import json, sys
value = json.load(open(sys.argv[1], encoding="utf-8"))
for key in (
    "density_limited_motifs",
    "total_loci_at_final_threshold_density_chrom",
    "maximum_orientation_rows_at_final_threshold_density_chrom",
):
    print(f"{key}\t{value[key]}")
PY
df -h /data/sm718
```

The row ceiling is not a byte forecast and chromosome 1 is not a guarantee for
other chromosomes. It is an early scale check. During Stage 2, the existing
free-space reserve is checked before every durable batch.

## Stage 2: whole-genome scan

Only prepare this plan after Stage 1 has finalized. The manager verifies exact
motif-set coverage and registry metadata, then groups motifs by identical final
score. A batch therefore always has one scalar scanner threshold; per-motif
thresholds are never approximated or mixed inside an invocation.

```sh
readonly RUN=/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_informative_density200_v1
readonly THRESHOLDS=$DENSITY_RUN/final/motif_scan_thresholds.tsv
readonly THRESHOLD_METADATA=$DENSITY_RUN/final/motif_scan_thresholds.json

chrom_args=()
for chrom in {1..22} X Y MT; do chrom_args+=(--chrom "$chrom"); done

scripts/manage_genome_scan.py prepare \
  --run-root "$RUN" \
  --run-id jaspar2026_grch38_informative_density200_v1 \
  --source-commit "$SOURCE_COMMIT" \
  --genome "$GENOME" --fasta-index "$FAI" --jaspar "$JASPAR" \
  --genome-id homo_sapiens_grch38_ensembl113_primary \
  --motif-set-id jaspar2026_core_nonredundant \
  --taxon-id 9606 --species 'Homo sapiens' \
  --assembly-name GRCh38 --assembly-accession GCA_000001405.29 \
  --ensembl-release 113 \
  --fasta-url 'https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz' \
  --jaspar-url 'https://jaspar.elixir.no/download/data/2026/CORE/JASPAR2026_CORE_non-redundant_pfms_jaspar.txt' \
  --motif-thresholds "$THRESHOLDS" \
  --motif-threshold-metadata "$THRESHOLD_METADATA" \
  --motif-batch-size 64 \
  --score-mode log2_relative_risk --pseudocount 1 \
  --minimum-free-gib 500 "${chrom_args[@]}"

scripts/submit_genome_scan_slurm.sh \
  --run-root "$RUN" --scanner "$SOURCE/pssm_scan_parquet" \
  --source "$SOURCE" --account cluster --partition requeue \
  --max-concurrent 20 --batch-workers 2 \
  --memory 32G --time 0-08:00:00 --finalize --dry-run
```

One array element owns one chromosome. It copies that chromosome and the
JASPAR source from `/data` to a unique `/scratch` directory once and processes
every still-incomplete score-homogeneous motif batch there. Durable packages are promoted
atomically below `/data/sm718`; requeue never repeats a completed batch. The
scanner reads the staged chromosome once per batch process, not once from the
shared filesystem per motif.

Finalized packages expose `scan_motif_threshold` in DuckDB. It keeps the full
per-motif decision beside `scan_threshold_policy`, while legacy packages get a
compatible threshold view reconstructed from their file inventory.

## Interpretation limits

The informative component optimizes an exploratory TP73-context prediction on
chromosome 1 and is not an assay-independent binding cutoff. The density rule
is an engineering/storage constraint. Keep both components, their source data,
and the continuous motif score visible in downstream analyses. A later
held-out-chromosome or independent-assay threshold set should receive a new
identifier rather than mutating this version.

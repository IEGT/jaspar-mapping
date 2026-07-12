# JASPAR 2026 Chromosome 1 Motif Panel

`scripts/run_chr1_2026_motif_panel.sh` performs the full chromosome 1 dense
scan used to calibrate later sparse genome-wide storage. It runs five JASPAR
2026 motifs in two score modes and both orientations, always with pseudocount
1 and skip-N handling:

| Factor | Motif | Length |
| --- | --- | ---: |
| TP73 | `MA0861.2` | 16 |
| E2F1 | `MA0024.3` | 12 |
| SP1 | `MA0079.5` | 9 |
| PATZ1 | `MA1961.2` | 11 |
| POU2F2 | `MA0507.3` | 13 |

The runner downloads the inputs directly from the versioned public JASPAR and
Ensembl URLs printed by `--print-config`. Given `--source-commit`, it also
downloads that exact commit from the public GitHub archive, so compute nodes do
not need Git. It creates a pinned Conda environment inside the run root, builds
the direct-Parquet scanner, and launches exactly 20 independent Slurm steps. No
local workstation files are inputs to the analysis.

On Haumea, prepare a dedicated directory and submit with the requeue resources:

```sh
RUN_ROOT=/data/sm718/codex/jaspar2026_chr1_dense_5motifs
REPOSITORY=/data/sm718/GitHub/jaspar-mapping
SOURCE_COMMIT=$(/usr/bin/git -C "$REPOSITORY" rev-parse HEAD)
mkdir -p "$RUN_ROOT/logs"
sbatch \
  --account=cluster \
  --partition=requeue \
  --ntasks=20 \
  --cpus-per-task=1 \
  --mem=64G \
  --time=1-00:00:00 \
  --requeue \
  --chdir="$RUN_ROOT" \
  --output="$RUN_ROOT/logs/slurm-%j.out" \
  --error="$RUN_ROOT/logs/slurm-%j.err" \
  "$REPOSITORY/scripts/run_chr1_2026_motif_panel.sh" \
  --run-root "$RUN_ROOT" \
  --source-commit "$SOURCE_COMMIT"
```

The 64 GiB request includes 3 GiB for each scanner plus 4 GiB for the batch
coordinator and headroom. This is based on observed complete-chromosome peak
RSS values of approximately 1.8--2.1 GiB per scanner. The runner passes the
3 GiB reservation to every `srun` step explicitly; otherwise Slurm may assign the
whole job-memory request to each step and serialize the nominally parallel
scans.

Each task first writes to a job-specific staging directory. A Parquet file is
promoted to the final package only after its block extent and total number of
motif alignment starts match chromosome 1. A requeued job reuses only final
files that pass the same validation and never overwrites or removes an existing
file. The resulting package contains the Parquet partitions, motif metadata,
`manifest.json`, and a rebuildable DuckDB query index suitable for
`scripts/export_dense_bed.py`.

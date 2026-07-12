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
Ensembl URLs printed by `--print-config`. It creates a pinned Conda environment
inside the run root, archives the current public Git commit into that root,
builds the direct-Parquet scanner, and launches exactly 20 independent Slurm
steps. No local workstation files are inputs to the analysis.

On Haumea, prepare a dedicated directory and submit with the requeue resources:

```sh
RUN_ROOT=/data/sm718/codex/jaspar2026_chr1_dense_5motifs
mkdir -p "$RUN_ROOT/logs"
sbatch \
  --account=cluster \
  --partition=requeue \
  --ntasks=20 \
  --cpus-per-task=1 \
  --mem=20G \
  --time=1-00:00:00 \
  --requeue \
  --chdir="$RUN_ROOT" \
  --output="$RUN_ROOT/logs/slurm-%j.out" \
  --error="$RUN_ROOT/logs/slurm-%j.err" \
  scripts/run_chr1_2026_motif_panel.sh --run-root "$RUN_ROOT"
```

Each task first writes to a job-specific staging directory. A Parquet file is
promoted to the final package only after its block extent and total number of
motif alignment starts match chromosome 1. A requeued job reuses only final
files that pass the same validation and never overwrites or removes an existing
file. The resulting package contains the Parquet partitions, motif metadata,
`manifest.json`, and a rebuildable DuckDB query index suitable for
`scripts/export_dense_bed.py`.

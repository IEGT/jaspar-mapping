# AGENTS.md

## Promoter Architecture And ML Query Surface

For work on promoter architecture, expression-linked ML features, OpenClaw
integration, or agent-facing query surfaces, read
[`docs/promoter_architecture_ml_schema.md`](docs/promoter_architecture_ml_schema.md)
before designing or changing code.

For motif-score calibration against raw CUT&RUN fragments, also read
[`docs/cutandrun_motif_score_calibration.md`](docs/cutandrun_motif_score_calibration.md).

For TP73-local motif spacing, tandem TP73 sites, or transcript/intron context,
also read [`docs/tp73_motif_context.md`](docs/tp73_motif_context.md).

For whole-genome JASPAR 2026 scans, Slurm batching/requeue behavior, production
inventories, or cross-species identity, also read
[`docs/jaspar2026_genome_scan_plan.md`](docs/jaspar2026_genome_scan_plan.md).

Treat [`sql/schema.sql`](sql/schema.sql) and [`sql/queries.sql`](sql/queries.sql)
as the draft query contract:

- Use BED 0-based half-open coordinates for exported motif and promoter tables.
- Keep large generated data out of Git; package it as versioned Parquet plus a
  rebuildable DuckDB query index.
- Prefer the materialized `promoter_card` for low-latency agent lookups, bind
  its `score_mode` and `pseudocount` explicitly, and use named parameterized
  queries instead of ad hoc joins from agents.
- Keep ML features in long/tidy tables; pivot to wide matrices only at the ML
  export boundary.
- Keep sequence-derived promoter architecture separate from CUT&RUN-derived
  evidence unless an analysis explicitly declares an "architecture + binding"
  feature block.

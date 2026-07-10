# AGENTS.md

## Promoter Architecture And ML Query Surface

For work on promoter architecture, expression-linked ML features, OpenClaw
integration, or agent-facing query surfaces, read
[`docs/promoter_architecture_ml_schema.md`](docs/promoter_architecture_ml_schema.md)
before designing or changing code.

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

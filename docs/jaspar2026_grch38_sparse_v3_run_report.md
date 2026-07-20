# JASPAR 2026 GRCh38 sparse v3 production report

This report records the observed production package that motivated the scanner
maintenance work. It is an operational benchmark, not a reinterpretation of
the biological results.

## Identity

- Source commit: `42eb83df8855980e77002b9055558c307922608a`
- Durable package:
  `/data/sm718/jaspar_mapping_runs/jaspar2026_grch38_sparse_v3/package`
- Motif collection: JASPAR 2026 CORE non-redundant
- Motifs: 2,633
- Sequence regions: 25
- Thresholds: TP73 `MA0861.2` at `-5`; all other motifs at `-1`
- Score mode: `log2_relative_risk`
- Pseudocount: 1 per nucleotide count
- Strands: both
- Coordinates: BED 0-based half-open

## Completed package

- Planned/completed tasks: 1,075 / 1,075
- Failed tasks: 0
- Parquet files: 131,650
- Retained hits: 340,672,877,861
- Parquet bytes: 1,858,587,985,628
- Approximate binary size: 1.69 TiB

The task count is 25 chromosomes multiplied by 43 non-overlapping motif
batches. Each motif/chromosome/orientation has one inventory entry and one
Parquet file.

## Original finalizer observation

The successful finalizer Slurm job `5435555` used:

- elapsed time: 6:31:42;
- maximum resident set: 5,619,356 KiB (about 5.36 GiB);
- maximum disk read: 1,773,279.03 MiB;
- maximum disk write: 206.15 MiB.

Most time was spent recomputing checksums over essentially the entire payload.
Afterward, binding the package-wide DuckDB `motif_hit` wildcard with
`union_by_name=true` took approximately 45 minutes and raised memory from about
1.1 GiB to about 5.6 GiB despite little biological payload being queried.

The measured behavior led to two maintenance decisions:

1. ordinary finalization reuses task-time checksums and performs exact
   path/size/stat/provenance checks, while a full reread is an explicit,
   resumable audit; and
2. the DuckDB catalog contains metadata and an exact-file table macro rather
   than binding every hit file at catalog creation.

No maintenance operation should rewrite this production payload. A rebuilt
metadata-only catalog can be placed beside it under a new filename and compared
against the existing manifest before adoption.

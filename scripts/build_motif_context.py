#!/usr/bin/env python3

"""Build long-form TP73 motif-pair and transcript-context Parquet tables.

The pair table preserves every configured motif occurrence within the broad
capture radius, including overlapping TP73 alignments and opposite-strand
alignments of the same sequence span.  Summary features use a narrower,
provisional context radius.  A tandem TP73 partner is intentionally stricter:
it must be a distinct, non-overlapping TP73 span whose edge-to-edge gap is at
most ``--tandem-flank``.  Distances describe centers of scored sequence spans,
not inferred physical footprints of bound proteins.

When a GTF is supplied, coordinates are converted directly from GTF 1-based
inclusive to BED-style 0-based half-open coordinates.  The output is a movable
DuckDB-plus-Parquet package; no BED or TSV intermediate is produced.
"""

from __future__ import annotations

import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


class ContextBuildError(RuntimeError):
    pass


def nonnegative_integer(value: str) -> int:
    try:
        result = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected a non-negative integer") from error
    if result < 0:
        raise argparse.ArgumentTypeError("expected a non-negative integer")
    return result


def positive_integer(value: str) -> int:
    result = nonnegative_integer(value)
    if result == 0:
        raise argparse.ArgumentTypeError("expected an integer greater than zero")
    return result


def nonnegative_float(value: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected a non-negative number") from error
    if result < 0:
        raise argparse.ArgumentTypeError("expected a non-negative number")
    return result


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a de novo, long-form Parquet context around TP73 motif "
            "occurrences without an intermediate BED/TSV table."
        ),
        epilog=(
            "Distances are between BED interval centers. Genomic distance increases "
            "with the reference coordinate; anchor-oriented distance is sign-flipped "
            "for minus-strand TP73 anchors."
        ),
    )
    parser.add_argument(
        "--motif-hits",
        required=True,
        help="motif-hit Parquet file, directory, or glob",
    )
    parser.add_argument("--output", required=True, type=Path,
                        help="new output package directory")
    parser.add_argument("--gtf", type=Path,
                        help="optional Ensembl-compatible GTF or GTF.gz for TSS/introns")
    parser.add_argument("--anchor-motif", default="MA0861.2",
                        help="anchor motif accession (default: MA0861.2)")
    parser.add_argument(
        "--score-mode",
        default="log2_relative_risk",
        choices=("log2_relative_risk", "log_odds"),
    )
    parser.add_argument("--pseudocount", type=nonnegative_float, default=1.0)
    parser.add_argument(
        "--chrom",
        action="append",
        default=[],
        metavar="CHROM[,CHROM...]",
        help="restrict chromosomes; repeat or comma-separate (default: all)",
    )
    parser.add_argument(
        "--capture-flank",
        type=nonnegative_integer,
        default=500,
        help="maximum center-distance retained in the pair table (default: 500)",
    )
    parser.add_argument(
        "--context-flank",
        type=nonnegative_integer,
        default=150,
        help="provisional center-distance used for context features (default: 150)",
    )
    parser.add_argument(
        "--tandem-flank",
        type=nonnegative_integer,
        default=20,
        help="maximum edge-to-edge gap for a non-overlapping TP73 tandem partner (default: 20)",
    )
    parser.add_argument("--threads", type=positive_integer, default=1)
    parser.add_argument("--memory-limit", default="8GB")
    parser.add_argument("--max-temp-size", default="40GB")
    parser.add_argument("--duckdb", default="duckdb")
    parser.add_argument("--force", action="store_true",
                        help="replace an existing output package atomically")
    return parser


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_number(value: float) -> str:
    return format(value, ".17g")


def parse_chromosomes(values: list[str]) -> list[str]:
    result: list[str] = []
    for value in values:
        for chromosome in value.split(","):
            chromosome = chromosome.strip()
            if not chromosome:
                raise ContextBuildError("--chrom contains an empty chromosome")
            if chromosome not in result:
                result.append(chromosome)
    return result


def resolve_parquet_glob(value: str) -> str:
    expanded = os.path.expanduser(value)
    candidate = Path(expanded)
    if candidate.is_dir():
        pattern = str(candidate.resolve() / "**" / "*.parquet")
    else:
        pattern = expanded if os.path.isabs(expanded) else os.path.abspath(expanded)
    if not glob.glob(pattern, recursive=True):
        raise ContextBuildError(f"no motif-hit Parquet files match: {value}")
    return pattern


def report_disk_preflight(output_parent: Path, motif_glob: str) -> None:
    source_files = [Path(filename) for filename in glob.glob(motif_glob, recursive=True)]
    source_bytes = sum(filename.stat().st_size for filename in source_files)
    free_bytes = shutil.disk_usage(output_parent).free
    gib = 1024 ** 3
    print(
        f"I: Motif context input is {source_bytes / gib:.2f} GiB; "
        f"output filesystem has {free_bytes / gib:.2f} GiB free.",
        file=sys.stderr,
    )
    if free_bytes < source_bytes:
        print(
            "W: Free space is smaller than the motif-hit input; the pair table "
            "can be larger than its input.",
            file=sys.stderr,
        )


def chromosome_filter(chromosomes: list[str], column: str = "chrom") -> str:
    if not chromosomes:
        return ""
    values = ", ".join(sql_string(chromosome) for chromosome in chromosomes)
    return f" AND CAST({column} AS VARCHAR) IN ({values})"


def gtf_sql(gtf: Path, chromosomes: list[str]) -> str:
    chrom_clause = chromosome_filter(chromosomes, "chrom")
    return f"""
-- Convert GTF coordinates once at ingestion. GTF end is already the correct
-- exclusive BED end after subtracting one only from the GTF start.
CREATE TEMP VIEW gtf_raw AS
SELECT
    column0::VARCHAR AS chrom,
    column2::VARCHAR AS feature,
    column3::BIGINT - 1 AS start,
    column4::BIGINT AS "end",
    column6::VARCHAR AS strand,
    trim(regexp_extract(column8, 'gene_id "?([^";]+)"?', 1)) AS gene_id,
    trim(regexp_extract(column8, 'transcript_id "?([^";]+)"?', 1)) AS transcript_id,
    trim(regexp_extract(column8, 'gene_name "?([^";]+)"?', 1)) AS gene_name,
    trim(regexp_extract(column8, '(?:transcript_biotype|transcript_type) "?([^";]+)"?', 1))
        AS biotype
FROM read_csv(
    {sql_string(gtf)},
    delim='\\t', header=false, comment='#', quote='', auto_detect=false,
    strict_mode=false,
    columns={{
        'column0':'VARCHAR', 'column1':'VARCHAR', 'column2':'VARCHAR',
        'column3':'BIGINT', 'column4':'BIGINT', 'column5':'VARCHAR',
        'column6':'VARCHAR', 'column7':'VARCHAR', 'column8':'VARCHAR'
    }}
)
WHERE column2 IN ('transcript', 'exon')
  AND regexp_extract(column8, 'transcript_id "?([^";]+)"?', 1) <> ''
  {chrom_clause};

CREATE TEMP TABLE transcript_annotation AS
SELECT
    gene_id,
    transcript_id,
    CAST(chrom AS VARCHAR) AS chrom,
    strand,
    MIN(start)::BIGINT AS transcript_start,
    MAX("end")::BIGINT AS transcript_end,
    CASE WHEN strand = '+' THEN MIN(start) ELSE MAX("end") END::BIGINT AS tss,
    MAX(gene_name) AS gene_name,
    MAX(biotype) AS biotype
FROM gtf_raw
GROUP BY gene_id, transcript_id, chrom, strand;

CREATE TEMP TABLE intron_annotation AS
WITH ordered_exon AS (
    SELECT
        gene_id,
        transcript_id,
        CAST(chrom AS VARCHAR) AS chrom,
        strand,
        start AS exon_start,
        "end" AS exon_end,
        LAG("end") OVER (
            PARTITION BY gene_id, transcript_id
            ORDER BY start, "end"
        ) AS previous_exon_end
    FROM gtf_raw
    WHERE feature = 'exon'
), gaps AS (
    SELECT gene_id, transcript_id, chrom, strand,
           previous_exon_end::BIGINT AS start,
           exon_start::BIGINT AS "end"
    FROM ordered_exon
    WHERE previous_exon_end < exon_start
)
SELECT
    gene_id,
    transcript_id,
    chrom,
    strand,
    start,
    "end",
    CASE WHEN strand = '+' THEN
        ROW_NUMBER() OVER (PARTITION BY gene_id, transcript_id ORDER BY start)
    ELSE
        ROW_NUMBER() OVER (PARTITION BY gene_id, transcript_id ORDER BY start DESC)
    END::INTEGER AS intron_number
FROM gaps;

COPY transcript_annotation TO 'tables/jaspar2026/transcript.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY intron_annotation TO 'tables/jaspar2026/intron.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);

CREATE TEMP TABLE anchor_transcript_context AS
WITH overlap AS (
    SELECT
        a.anchor_hit_id,
        a.chrom,
        a.start AS anchor_start,
        a."end" AS anchor_end,
        t.gene_id,
        t.gene_name,
        t.transcript_id,
        t.strand AS transcript_strand,
        t.tss,
        CASE WHEN t.strand = '+' THEN a.center_bp - t.tss
             ELSE t.tss - a.center_bp END AS signed_tss_distance_bp,
        COALESCE(BOOL_OR(a.start >= i.start AND a."end" <= i."end"), false)
            AS fully_within_intron,
        COALESCE(BOOL_OR(a.start < i."end" AND a."end" > i.start), false)
            AS overlaps_intron
    FROM anchor_hit a
    JOIN transcript_annotation t
      ON a.chrom = t.chrom
     AND a.start >= t.transcript_start
     AND a."end" <= t.transcript_end
    LEFT JOIN intron_annotation i
      ON i.gene_id = t.gene_id
     AND i.transcript_id = t.transcript_id
     AND a.start < i."end"
     AND a."end" > i.start
    GROUP BY a.anchor_hit_id, a.chrom, a.start, a."end", t.gene_id,
             t.gene_name, t.transcript_id, t.strand, t.tss, a.center_bp
)
SELECT *,
       CASE WHEN fully_within_intron THEN 'intron'
            WHEN overlaps_intron THEN 'intron_boundary'
            ELSE 'transcribed_non_intron' END AS transcript_region
FROM overlap;

COPY anchor_transcript_context TO 'tables/jaspar2026/motif_transcript_context.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);

CREATE TEMP TABLE unique_tss AS
SELECT * EXCLUDE (tss_rank)
FROM (
    SELECT *, ROW_NUMBER() OVER (
        PARTITION BY chrom, tss
        ORDER BY gene_id, transcript_id
    ) AS tss_rank
    FROM transcript_annotation
)
WHERE tss_rank = 1;

CREATE TEMP TABLE nearest_left AS
SELECT a.anchor_hit_id, a.center_bp, t.gene_id, t.gene_name,
       t.transcript_id, t.strand, t.tss
FROM anchor_hit a
ASOF LEFT JOIN unique_tss t
  ON a.chrom = t.chrom AND a.center_bp >= t.tss;

CREATE TEMP TABLE nearest_right AS
SELECT a.anchor_hit_id, a.center_bp, t.gene_id, t.gene_name,
       t.transcript_id, t.strand, t.tss
FROM anchor_hit a
ASOF LEFT JOIN unique_tss t
  ON a.chrom = t.chrom AND a.center_bp <= t.tss;

CREATE TEMP TABLE nearest_tss AS
SELECT * EXCLUDE (nearest_rank)
FROM (
    SELECT *, ROW_NUMBER() OVER (
        PARTITION BY anchor_hit_id
        ORDER BY ABS(center_bp - tss), gene_id, transcript_id
    ) AS nearest_rank
    FROM (
        SELECT * FROM nearest_left WHERE transcript_id IS NOT NULL
        UNION ALL
        SELECT * FROM nearest_right WHERE transcript_id IS NOT NULL
    )
)
WHERE nearest_rank = 1;

CREATE TEMP TABLE transcript_context_summary AS
SELECT
    anchor_hit_id,
    true AS in_any_transcript,
    BOOL_OR(fully_within_intron) AS in_any_intron,
    BOOL_OR(overlaps_intron AND NOT fully_within_intron) AS overlaps_any_intron_boundary
FROM anchor_transcript_context
GROUP BY anchor_hit_id;

CREATE TEMP TABLE anchor_gene_context AS
SELECT
    a.anchor_hit_id,
    true AS gene_annotation_available,
    n.gene_id AS nearest_gene_id,
    n.gene_name AS nearest_gene_name,
    n.transcript_id AS nearest_transcript_id,
    CASE WHEN n.strand = '+' THEN a.center_bp - n.tss
         WHEN n.strand = '-' THEN n.tss - a.center_bp
         ELSE NULL END AS nearest_tss_distance_bp,
    COALESCE(s.in_any_transcript, false) AS in_any_transcript,
    COALESCE(s.in_any_intron, false) AS in_any_intron,
    COALESCE(s.overlaps_any_intron_boundary, false) AS overlaps_any_intron_boundary,
    CASE WHEN COALESCE(s.in_any_intron, false) THEN 'intron'
         WHEN COALESCE(s.overlaps_any_intron_boundary, false) THEN 'intron_boundary'
         WHEN COALESCE(s.in_any_transcript, false) THEN 'transcribed_non_intron'
         ELSE 'intergenic' END AS primary_transcript_region
FROM anchor_hit a
LEFT JOIN nearest_tss n USING (anchor_hit_id)
LEFT JOIN transcript_context_summary s USING (anchor_hit_id);
"""


def no_gtf_sql() -> str:
    return """
CREATE TEMP TABLE anchor_gene_context AS
SELECT
    anchor_hit_id,
    false AS gene_annotation_available,
    NULL::VARCHAR AS nearest_gene_id,
    NULL::VARCHAR AS nearest_gene_name,
    NULL::VARCHAR AS nearest_transcript_id,
    NULL::DOUBLE AS nearest_tss_distance_bp,
    NULL::BOOLEAN AS in_any_transcript,
    NULL::BOOLEAN AS in_any_intron,
    NULL::BOOLEAN AS overlaps_any_intron_boundary,
    'not_assessed'::VARCHAR AS primary_transcript_region
FROM anchor_hit;
"""


def build_sql(arguments: argparse.Namespace, motif_glob: str,
              chromosomes: list[str], gtf: Path | None) -> str:
    chrom_clause = chromosome_filter(chromosomes)
    annotation_sql = gtf_sql(gtf, chromosomes) if gtf is not None else no_gtf_sql()
    gtf_source = sql_string(gtf) if gtf is not None else "NULL"
    return f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;
SET temp_directory='duckdb_tmp';
SET max_temp_directory_size={sql_string(arguments.max_temp_size)};
PRAGMA enable_progress_bar;

CREATE TEMP TABLE context_run_config AS SELECT
    1::INTEGER AS schema_version,
    {sql_string(arguments.anchor_motif)}::VARCHAR AS anchor_motif_id,
    {sql_string(arguments.score_mode)}::VARCHAR AS score_mode,
    {sql_number(arguments.pseudocount)}::DOUBLE AS pseudocount,
    {arguments.capture_flank}::INTEGER AS capture_flank_bp,
    {arguments.context_flank}::INTEGER AS context_flank_bp,
    {arguments.tandem_flank}::INTEGER AS tandem_flank_bp,
    'motif_center'::VARCHAR AS distance_metric,
    'nonoverlapping_edge_gap'::VARCHAR AS tandem_distance_metric,
    'anchor_minus_strand_flips_sign'::VARCHAR AS oriented_distance_rule,
    'same_motif_id_nonoverlapping_distinct_span'::VARCHAR AS tandem_identity_rule,
    {sql_string(motif_glob)}::VARCHAR AS motif_hit_source,
    {gtf_source}::VARCHAR AS gtf_source;

COPY context_run_config TO 'tables/jaspar2026/context_run_config.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);

-- A dense scan may be assembled from overlapping/retried Parquet parts. Keep
-- one deterministic row per scored orientation and model span before joining.
CREATE TEMP TABLE configured_hit AS
SELECT * EXCLUDE (duplicate_rank)
FROM (
    SELECT
        md5(concat_ws('|', CAST(chrom AS VARCHAR), start::VARCHAR, "end"::VARCHAR,
                      motif_id, strand, score_mode, printf('%.17g', pseudocount)))
            AS hit_id,
        CAST(chrom AS VARCHAR) AS chrom,
        start::BIGINT AS start,
        "end"::BIGINT AS "end",
        (start::DOUBLE + "end"::DOUBLE) / 2.0 AS center_bp,
        motif_id::VARCHAR AS motif_id,
        motif_name::VARCHAR AS motif_name,
        strand::VARCHAR AS strand,
        score::DOUBLE AS score,
        score_mode::VARCHAR AS score_mode,
        pseudocount::DOUBLE AS pseudocount,
        pwm_relative_score::DOUBLE AS pwm_relative_score,
        ROW_NUMBER() OVER (
            PARTITION BY CAST(chrom AS VARCHAR), start, "end", motif_id, strand,
                         score_mode, pseudocount
            ORDER BY score DESC NULLS LAST, pwm_relative_score DESC NULLS LAST
        ) AS duplicate_rank
    FROM read_parquet({sql_string(motif_glob)}, hive_partitioning=1)
    WHERE score_mode = {sql_string(arguments.score_mode)}
      AND pseudocount = {sql_number(arguments.pseudocount)}
      {chrom_clause}
)
WHERE duplicate_rank = 1;

CREATE TEMP TABLE anchor_hit AS
SELECT hit_id AS anchor_hit_id, * EXCLUDE (hit_id)
FROM configured_hit
WHERE motif_id = {sql_string(arguments.anchor_motif)};

SELECT CASE WHEN (SELECT COUNT(*) FROM anchor_hit) = 0
    THEN error('no anchor motif hits matched the requested configuration') END;

-- Keep the broad de novo pair layer lossless with respect to local alignments.
-- same_anchor_motif_span lets summaries suppress the opposite orientation of
-- an identical TP73 span without deleting that observation from the pair data.
CREATE TEMP TABLE motif_context_pair AS
SELECT
    a.anchor_hit_id,
    a.chrom,
    a.start AS anchor_start,
    a."end" AS anchor_end,
    a.center_bp AS anchor_center_bp,
    a.strand AS anchor_strand,
    a.score AS anchor_score,
    a.pwm_relative_score AS anchor_pwm_relative_score,
    n.hit_id AS neighbor_hit_id,
    n.start AS neighbor_start,
    n."end" AS neighbor_end,
    n.center_bp AS neighbor_center_bp,
    n.motif_id AS neighbor_motif_id,
    n.motif_name AS neighbor_motif_name,
    n.strand AS neighbor_strand,
    n.score AS neighbor_score,
    n.pwm_relative_score AS neighbor_pwm_relative_score,
    n.center_bp - a.center_bp AS genomic_center_distance_bp,
    CASE WHEN a.strand = '-' THEN a.center_bp - n.center_bp
         ELSE n.center_bp - a.center_bp END AS anchor_oriented_center_distance_bp,
    ABS(n.center_bp - a.center_bp) AS absolute_center_distance_bp,
    CASE WHEN a.strand = n.strand THEN 'same' ELSE 'opposite' END
        AS relative_orientation,
    CASE WHEN (CASE WHEN a.strand = '-' THEN a.center_bp - n.center_bp
                    ELSE n.center_bp - a.center_bp END) < 0 THEN 'upstream'
         WHEN (CASE WHEN a.strand = '-' THEN a.center_bp - n.center_bp
                    ELSE n.center_bp - a.center_bp END) > 0 THEN 'downstream'
         ELSE 'coincident_center' END AS anchor_oriented_side,
    n.start = a.start AND n."end" = a."end" AS same_alignment_span,
    -- Tandem is a conservative architectural label: distinct TP73 spans may
    -- touch, but shifted overlapping alignments are retained only as pairs.
    n.motif_id = a.motif_id
        AND n.start = a.start AND n."end" = a."end"
        AS same_anchor_motif_span,
    GREATEST(
        0,
        LEAST(a."end", n."end") - GREATEST(a.start, n.start)
    )::BIGINT AS interval_overlap_bp,
    CASE WHEN n.start >= a."end" THEN n.start - a."end"
         WHEN a.start >= n."end" THEN a.start - n."end"
         ELSE 0 END::BIGINT AS inter_motif_gap_bp,
    ABS(n.center_bp - a.center_bp) <= {arguments.context_flank}
        AS within_context_flank,
    n.motif_id = a.motif_id
        AND NOT (n.start = a.start AND n."end" = a."end")
        AND LEAST(a."end", n."end") <= GREATEST(a.start, n.start)
        AND (CASE WHEN n.start >= a."end" THEN n.start - a."end"
                  WHEN a.start >= n."end" THEN a.start - n."end"
                  ELSE 0 END) <= {arguments.tandem_flank}
        AS is_tandem_tp73,
    a.score_mode,
    a.pseudocount,
    {arguments.capture_flank}::INTEGER AS capture_flank_bp,
    {arguments.context_flank}::INTEGER AS context_flank_bp,
    {arguments.tandem_flank}::INTEGER AS tandem_flank_bp
FROM anchor_hit a
JOIN configured_hit n
  ON a.chrom = n.chrom
 AND n.center_bp BETWEEN a.center_bp - {arguments.capture_flank}
                     AND a.center_bp + {arguments.capture_flank}
WHERE n.hit_id <> a.anchor_hit_id;

SELECT CASE WHEN (SELECT COUNT(*) FROM motif_context_pair) = 0
    THEN error('anchor hits have no neighboring motif hits inside capture flank') END;

COPY motif_context_pair TO 'tables/jaspar2026/motif_context_pair'
    (FORMAT PARQUET, PARTITION_BY (chrom), COMPRESSION ZSTD);

CREATE TEMP TABLE pair_summary AS
SELECT
    anchor_hit_id,
    COUNT(*) FILTER (
        WHERE within_context_flank AND NOT same_anchor_motif_span
    )::BIGINT AS n_context_neighbor_hits,
    COUNT(DISTINCT concat_ws('|', neighbor_motif_id, neighbor_start::VARCHAR,
                             neighbor_end::VARCHAR))
        FILTER (
            WHERE within_context_flank AND NOT same_anchor_motif_span
        )::BIGINT AS n_context_neighbor_loci,
    COUNT(DISTINCT neighbor_motif_id)
        FILTER (
            WHERE within_context_flank AND NOT same_anchor_motif_span
        )::BIGINT AS n_context_motifs,
    COUNT(DISTINCT neighbor_hit_id)
        FILTER (WHERE is_tandem_tp73)::BIGINT AS n_tandem_tp73_partners
FROM motif_context_pair
GROUP BY anchor_hit_id;

CREATE TEMP TABLE nearest_tandem AS
SELECT * EXCLUDE (tandem_rank)
FROM (
    SELECT
        anchor_hit_id,
        neighbor_hit_id,
        anchor_oriented_center_distance_bp,
        genomic_center_distance_bp,
        absolute_center_distance_bp,
        relative_orientation,
        neighbor_score,
        ROW_NUMBER() OVER (
            PARTITION BY anchor_hit_id
            ORDER BY absolute_center_distance_bp, neighbor_start, neighbor_end,
                     neighbor_strand, neighbor_hit_id
        ) AS tandem_rank
    FROM motif_context_pair
    WHERE is_tandem_tp73
)
WHERE tandem_rank = 1;

{annotation_sql}

CREATE TEMP TABLE tp73_context_anchor AS
SELECT
    a.anchor_hit_id,
    a.chrom,
    a.start,
    a."end",
    a.center_bp,
    a.motif_id,
    a.motif_name,
    a.strand,
    a.score,
    a.score_mode,
    a.pseudocount,
    a.pwm_relative_score,
    COALESCE(s.n_context_neighbor_hits, 0)::BIGINT AS n_context_neighbor_hits,
    COALESCE(s.n_context_neighbor_loci, 0)::BIGINT AS n_context_neighbor_loci,
    COALESCE(s.n_context_motifs, 0)::BIGINT AS n_context_motifs,
    COALESCE(s.n_tandem_tp73_partners, 0) > 0 AS has_tandem_tp73,
    COALESCE(s.n_tandem_tp73_partners, 0)::BIGINT AS n_tandem_tp73_partners,
    t.neighbor_hit_id AS nearest_tandem_hit_id,
    t.anchor_oriented_center_distance_bp AS nearest_tandem_oriented_distance_bp,
    t.genomic_center_distance_bp AS nearest_tandem_genomic_distance_bp,
    t.absolute_center_distance_bp AS nearest_tandem_absolute_distance_bp,
    t.relative_orientation AS nearest_tandem_relative_orientation,
    t.neighbor_score AS nearest_tandem_score,
    g.gene_annotation_available,
    g.nearest_gene_id,
    g.nearest_gene_name,
    g.nearest_transcript_id,
    g.nearest_tss_distance_bp,
    g.in_any_transcript,
    g.in_any_intron,
    g.overlaps_any_intron_boundary,
    g.primary_transcript_region,
    {arguments.capture_flank}::INTEGER AS capture_flank_bp,
    {arguments.context_flank}::INTEGER AS context_flank_bp,
    {arguments.tandem_flank}::INTEGER AS tandem_flank_bp
FROM anchor_hit a
LEFT JOIN pair_summary s USING (anchor_hit_id)
LEFT JOIN nearest_tandem t USING (anchor_hit_id)
JOIN anchor_gene_context g USING (anchor_hit_id);

COPY tp73_context_anchor TO 'tables/jaspar2026/tp73_context_anchor'
    (FORMAT PARQUET, PARTITION_BY (chrom), COMPRESSION ZSTD);

DROP TABLE context_run_config;
DROP TABLE motif_context_pair;
DROP TABLE tp73_context_anchor;
CREATE VIEW context_run_config AS
SELECT * FROM read_parquet('tables/jaspar2026/context_run_config.parquet');
CREATE VIEW motif_context_pair AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/motif_context_pair/*/*.parquet', hive_partitioning=1
);
CREATE VIEW tp73_context_anchor AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/tp73_context_anchor/*/*.parquet', hive_partitioning=1
);
""" + (
        """
DROP TABLE anchor_transcript_context;
CREATE VIEW transcript AS
SELECT * FROM read_parquet('tables/jaspar2026/transcript.parquet');
CREATE VIEW intron AS
SELECT * FROM read_parquet('tables/jaspar2026/intron.parquet');
CREATE VIEW motif_transcript_context AS
SELECT * FROM read_parquet('tables/jaspar2026/motif_transcript_context.parquet');
""" if gtf is not None else ""
    )


def run_duckdb(executable: str, staging: Path, sql: str) -> None:
    process = subprocess.run(
        [executable, "context.duckdb"],
        input=sql,
        text=True,
        cwd=staging,
        check=False,
    )
    if process.returncode != 0:
        raise ContextBuildError(
            f"DuckDB context build failed with exit code {process.returncode}"
        )


def build_package(arguments: argparse.Namespace) -> None:
    if shutil.which(arguments.duckdb) is None:
        raise ContextBuildError(f"DuckDB executable not found: {arguments.duckdb}")
    if arguments.context_flank > arguments.capture_flank:
        raise ContextBuildError("--context-flank cannot exceed --capture-flank")
    if arguments.tandem_flank > arguments.capture_flank:
        raise ContextBuildError("--tandem-flank cannot exceed --capture-flank")

    motif_glob = resolve_parquet_glob(arguments.motif_hits)
    chromosomes = parse_chromosomes(arguments.chrom)
    gtf = arguments.gtf.expanduser().resolve() if arguments.gtf else None
    if gtf is not None and not gtf.is_file():
        raise ContextBuildError(f"GTF not found: {gtf}")

    output = arguments.output.expanduser().resolve()
    if output.exists() and not arguments.force:
        raise ContextBuildError(f"output already exists (use --force): {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    report_disk_preflight(output.parent, motif_glob)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        (staging / "tables" / "jaspar2026").mkdir(parents=True)
        sql = build_sql(arguments, motif_glob, chromosomes, gtf)
        run_duckdb(arguments.duckdb, staging, sql)
        shutil.rmtree(staging / "duckdb_tmp", ignore_errors=True)

        if output.exists():
            replacement = output.with_name(f".{output.name}.replaced")
            if replacement.exists():
                shutil.rmtree(replacement)
            os.replace(output, replacement)
            try:
                os.replace(staging, output)
            except Exception:
                os.replace(replacement, output)
                raise
            shutil.rmtree(replacement)
        else:
            os.replace(staging, output)
        print(f"I: Wrote motif context package: {output}", file=sys.stderr)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def main() -> int:
    parser = argument_parser()
    arguments = parser.parse_args()
    try:
        build_package(arguments)
        return 0
    except ContextBuildError as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

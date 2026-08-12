#!/usr/bin/env python3

"""Build long-form TP73 motif-pair and transcript-context Parquet tables.

The pair table preserves every configured motif occurrence whose signed
interval distance is within the capture radius. Directly abutting BED spans
have distance zero, overlaps have negative distance, and unoccupied bases
between spans give positive distance. Center offsets remain available for
direction and orientation, but do not decide context membership. A tandem
TP73 partner is intentionally stricter: it must be a distinct, non-overlapping
TP73 span whose edge-to-edge gap is at most ``--tandem-flank``.

The raw pair table preserves orientation-specific scoring records.  A separate
TP73 pair-feature table collapses records for the same neighboring alignment
span into one partner locus.  A locus represented on both strands is labelled
orientation-ambiguous instead of being counted as two physical partners.

When a GTF is supplied, coordinates are converted directly from GTF 1-based
inclusive to BED-style 0-based half-open coordinates.  The output is a movable
DuckDB-plus-Parquet package; no BED or TSV intermediate is produced.
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import math
import os
import re
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


def finite_float(value: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected a finite number") from error
    if not math.isfinite(result):
        raise argparse.ArgumentTypeError("expected a finite number")
    return result


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build a de novo, long-form Parquet context around TP73 motif "
            "occurrences without an intermediate BED/TSV table."
        ),
        epilog=(
            "Context distance is signed BED interval separation: zero means "
            "abutting, overlap is negative, and unoccupied bases are positive. "
            "Center offsets retain genomic and TP73-oriented direction."
        ),
    )
    parser.add_argument(
        "--motif-hits",
        required=True,
        help="motif-hit Parquet file, directory, or glob",
    )
    parser.add_argument(
        "--motif-hit-source-label",
        help=(
            "durable provenance label recorded instead of the physical input "
            "glob, useful when files are copied to node-local scratch"
        ),
    )
    parser.add_argument(
        "--input-uniqueness",
        choices=("deduplicate", "validated_scan_inventory"),
        default="deduplicate",
        help=(
            "deduplicate arbitrary input rows, or trust exact files selected "
            "from a finalized scan inventory (default: deduplicate)"
        ),
    )
    parser.add_argument("--output", required=True, type=Path,
                        help="new output package directory")
    parser.add_argument("--gtf", type=Path,
                        help="optional Ensembl-compatible GTF or GTF.gz for TSS/introns")
    parser.add_argument(
        "--annotation-release",
        default="unspecified",
        help=(
            "gene-annotation release recorded with TSS/promoter dimensions "
            "(default: unspecified; production builds should set this)"
        ),
    )
    parser.add_argument(
        "--promoter-definition-id",
        default="tss_upstream_2000_downstream_500_v1",
        help=(
            "versioned promoter-window identity (default: "
            "tss_upstream_2000_downstream_500_v1)"
        ),
    )
    parser.add_argument(
        "--promoter-upstream-bp",
        type=nonnegative_integer,
        default=2000,
        help="promoter bases upstream of the TSS, inclusive (default: 2000)",
    )
    parser.add_argument(
        "--promoter-downstream-bp",
        type=nonnegative_integer,
        default=500,
        help="promoter bases downstream of the TSS, inclusive (default: 500)",
    )
    parser.add_argument(
        "--downstream-definition-id",
        default="tes_upstream_500_downstream_2000_v1",
        help=(
            "versioned transcript-end downstream-window identity (default: "
            "tes_upstream_500_downstream_2000_v1)"
        ),
    )
    parser.add_argument(
        "--downstream-upstream-bp",
        type=nonnegative_integer,
        default=500,
        help=(
            "bases upstream of the TES, toward the transcript body, inclusive "
            "(default: 500)"
        ),
    )
    parser.add_argument(
        "--downstream-downstream-bp",
        type=nonnegative_integer,
        default=2000,
        help="bases downstream of the TES, inclusive (default: 2000)",
    )
    parser.add_argument(
        "--gtf-sha256",
        help=(
            "expected SHA-256 of the GTF bytes; verified before building and "
            "recorded in package provenance"
        ),
    )
    parser.add_argument(
        "--gtf-size-bytes",
        type=nonnegative_integer,
        help="expected GTF byte size; verified before building",
    )
    parser.add_argument("--anchor-motif", default="MA0861.2",
                        help="anchor motif accession (default: MA0861.2)")
    parser.add_argument(
        "--motif-set-id",
        required=True,
        help="explicit motif-set identity carried into every derived table",
    )
    parser.add_argument(
        "--source-commit",
        default="unknown",
        help="builder source commit recorded in package provenance",
    )
    parser.add_argument(
        "--genome-id",
        required=True,
        help="explicit reference-genome identity carried into every derived table",
    )
    parser.add_argument(
        "--score-mode",
        default="log2_relative_risk",
        choices=("log2_relative_risk", "log_odds"),
    )
    parser.add_argument("--pseudocount", type=nonnegative_float, default=1.0)
    parser.add_argument(
        "--background-model-id", default="uniform_acgt_v1",
        help="scoring background identity (default: uniform_acgt_v1)",
    )
    parser.add_argument(
        "--pseudocount-scheme", default="additive_per_base",
        help="pseudocount semantics (default: additive_per_base)",
    )
    parser.add_argument(
        "--anchor-minimum-score", required=True, type=finite_float,
        help="inclusive score floor for TP73 anchors",
    )
    parser.add_argument(
        "--tandem-minimum-score", "--partner-minimum-score",
        dest="partner_minimum_score", required=True, type=finite_float,
        help=(
            "inclusive score floor that both TP73 members must meet to be "
            "called a tandem (legacy alias: --partner-minimum-score)"
        ),
    )
    parser.add_argument(
        "--default-neighbor-minimum-score", type=finite_float, default=0.0,
        help=(
            "inclusive neighbor-locus threshold used when motif-hit input "
            "does not expose its minimum_score field (default: 0)"
        ),
    )
    parser.add_argument(
        "--anchor-selection-mode",
        choices=("threshold", "local_peak"),
        default="threshold",
        help=(
            "retain every anchor above the floor, or only physical-locus "
            "local maxima (default: threshold)"
        ),
    )
    parser.add_argument(
        "--anchor-local-peak-flank", type=positive_integer, default=150,
        help=(
            "maximum start-coordinate distance used to test whether a "
            "sub-threshold TP73 locus is a local maximum (default: 150)"
        ),
    )
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
        default=150,
        help="maximum signed interval distance retained (default: 150)",
    )
    parser.add_argument(
        "--context-flank",
        type=nonnegative_integer,
        default=150,
        help="signed interval distance used for context features (default: 150)",
    )
    parser.add_argument(
        "--tandem-flank",
        type=nonnegative_integer,
        default=20,
        help="maximum edge-to-edge gap for a non-overlapping TP73 tandem partner (default: 20)",
    )
    parser.add_argument(
        "--cofactor-pair-flank",
        type=nonnegative_integer,
        default=150,
        help=(
            "maximum signed interval distance between distinct same-motif "
            "cofactor loci (default: 150)"
        ),
    )
    parser.add_argument(
        "--output-tier",
        choices=("selected", "summary", "band"),
        default="selected",
        help=(
            "selected retains raw TP73/cofactor relationships; summary "
            "persists compact all-motif and shared TP73 features; band "
            "persists only non-TP73 per-band features (default: selected)"
        ),
    )
    parser.add_argument("--threads", type=positive_integer, default=1)
    parser.add_argument("--memory-limit", default="8GB")
    parser.add_argument("--max-temp-size", default="40GB")
    parser.add_argument(
        "--temp-directory",
        type=Path,
        help=(
            "external DuckDB spill directory, such as node-local /scratch; "
            "it must be empty and is not removed by this program"
        ),
    )
    parser.add_argument("--duckdb", default="duckdb")
    parser.add_argument("--force", action="store_true",
                        help="replace an existing output package atomically")
    return parser


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_number(value: float) -> str:
    return format(value, ".17g")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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


def gtf_sql(arguments: argparse.Namespace, gtf: Path,
            chromosomes: list[str]) -> str:
    chrom_clause = chromosome_filter(chromosomes, "chrom")
    genome_id = sql_string(arguments.genome_id)
    annotation_release = sql_string(arguments.annotation_release)
    promoter_definition_id = sql_string(arguments.promoter_definition_id)
    promoter_upstream_bp = arguments.promoter_upstream_bp
    promoter_downstream_bp = arguments.promoter_downstream_bp
    downstream_definition_id = sql_string(arguments.downstream_definition_id)
    downstream_upstream_bp = arguments.downstream_upstream_bp
    downstream_downstream_bp = arguments.downstream_downstream_bp
    return f"""
-- Convert GTF coordinates once at ingestion. GTF end is already the correct
-- exclusive BED end after subtracting one only from the GTF start.
CREATE TEMP VIEW gtf_raw AS
SELECT
    column0::VARCHAR AS chrom,
    lower(column2)::VARCHAR AS feature,
    column3::BIGINT - 1 AS start,
    column4::BIGINT AS "end",
    column6::VARCHAR AS strand,
    column7::VARCHAR AS phase,
    trim(regexp_extract(column8, 'gene_id "?([^";]+)"?', 1)) AS gene_id,
    trim(regexp_extract(column8, 'transcript_id "?([^";]+)"?', 1)) AS transcript_id,
    trim(regexp_extract(column8, 'gene_name "?([^";]+)"?', 1)) AS gene_name,
    trim(regexp_extract(column8, '(?:transcript_biotype|transcript_type) "?([^";]+)"?', 1))
        AS biotype,
    try_cast(regexp_extract(column8, 'exon_number "?([^";]+)"?', 1) AS INTEGER)
        AS exon_number
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
WHERE lower(column2) IN (
        'transcript', 'exon', 'cds', 'start_codon', 'stop_codon',
        'five_prime_utr', 'three_prime_utr', 'utr'
      )
  AND regexp_extract(column8, 'transcript_id "?([^";]+)"?', 1) <> ''
  {chrom_clause};

CREATE TEMP TABLE transcript_annotation AS
SELECT
    {genome_id}::VARCHAR AS genome_id,
    {annotation_release}::VARCHAR AS annotation_release,
    gene_id,
    transcript_id,
    CAST(chrom AS VARCHAR) AS chrom,
    strand,
    MIN(start)::BIGINT AS transcript_start,
    MAX("end")::BIGINT AS transcript_end,
    -- TSS is the genomic base, not the exclusive transcript boundary.
    CASE WHEN strand = '+' THEN MIN(start) ELSE MAX("end") - 1 END::BIGINT AS tss,
    -- TES is likewise the final transcribed genomic base, not the boundary
    -- immediately beyond the transcript.
    CASE WHEN strand = '+' THEN MAX("end") - 1 ELSE MIN(start) END::BIGINT AS tes,
    MAX(gene_name) AS gene_name,
    MAX(biotype) AS biotype
FROM gtf_raw
GROUP BY gene_id, transcript_id, chrom, strand;

-- Gene bounds and the compatibility representative TSS are derived from the
-- complete transcript set. Alternative physical starts remain authoritative
-- in transcription_start_site/transcript_tss below.
CREATE TEMP TABLE gene_annotation AS
SELECT
    genome_id,
    annotation_release,
    gene_id,
    MAX(gene_name) AS gene_name,
    chrom,
    strand,
    MIN(transcript_start)::BIGINT AS gene_start,
    MAX(transcript_end)::BIGINT AS gene_end,
    CASE WHEN strand = '+' THEN MIN(tss) ELSE MAX(tss) END::BIGINT AS tss,
    MAX(biotype) AS biotype,
    COUNT(DISTINCT transcript_id)::BIGINT AS n_transcripts
FROM transcript_annotation
GROUP BY genome_id, annotation_release, gene_id, chrom, strand;

CREATE TEMP TABLE exon_annotation AS
SELECT DISTINCT
    {genome_id}::VARCHAR AS genome_id,
    {annotation_release}::VARCHAR AS annotation_release,
    gene_id,
    gene_name,
    transcript_id,
    CAST(chrom AS VARCHAR) AS chrom,
    strand,
    start,
    "end",
    exon_number
FROM gtf_raw
WHERE feature = 'exon';

-- CDS identity is physical. Transcript ownership and reading phase live in a
-- bridge so shared coding segments are not multiplied in interval queries.
CREATE TEMP TABLE coding_sequence_segment AS
SELECT DISTINCT
    md5(concat_ws('|', {genome_id}, {annotation_release}, chrom,
                  start::VARCHAR, "end"::VARCHAR, strand)) AS cds_segment_id,
    {genome_id}::VARCHAR AS genome_id,
    {annotation_release}::VARCHAR AS annotation_release,
    CAST(chrom AS VARCHAR) AS chrom,
    start,
    "end",
    strand
FROM gtf_raw
WHERE feature = 'cds';

CREATE TEMP TABLE transcript_cds AS
SELECT DISTINCT
    s.genome_id,
    s.annotation_release,
    s.cds_segment_id,
    r.gene_id,
    r.gene_name,
    r.transcript_id,
    r.exon_number,
    try_cast(NULLIF(r.phase, '.') AS INTEGER) AS phase
FROM gtf_raw r
JOIN coding_sequence_segment s
  ON s.chrom = r.chrom
 AND s.start = r.start
 AND s."end" = r."end"
 AND s.strand = r.strand
WHERE r.feature = 'cds';

CREATE TEMP TABLE transcript_coding_span AS
SELECT
    s.genome_id,
    s.annotation_release,
    c.gene_id,
    MAX(c.gene_name) AS gene_name,
    c.transcript_id,
    s.chrom,
    s.strand,
    MIN(s.start)::BIGINT AS cds_start,
    MAX(s."end")::BIGINT AS cds_end,
    CASE WHEN s.strand = '+' THEN MIN(s.start) ELSE MAX(s."end") - 1 END::BIGINT
        AS cds_five_prime_base,
    CASE WHEN s.strand = '+' THEN MAX(s."end") - 1 ELSE MIN(s.start) END::BIGINT
        AS cds_three_prime_base,
    COUNT(DISTINCT s.cds_segment_id)::BIGINT AS n_cds_segments
FROM transcript_cds c
JOIN coding_sequence_segment s
  USING (genome_id, annotation_release, cds_segment_id)
GROUP BY s.genome_id, s.annotation_release, c.gene_id, c.transcript_id,
         s.chrom, s.strand;

CREATE TEMP TABLE coding_landmark AS
SELECT DISTINCT
    md5(concat_ws('|', {genome_id}, {annotation_release}, feature, chrom,
                  start::VARCHAR, "end"::VARCHAR, strand, gene_id,
                  transcript_id)) AS coding_landmark_id,
    {genome_id}::VARCHAR AS genome_id,
    {annotation_release}::VARCHAR AS annotation_release,
    feature AS landmark_type,
    gene_id,
    gene_name,
    transcript_id,
    CAST(chrom AS VARCHAR) AS chrom,
    start,
    "end",
    strand
FROM gtf_raw
WHERE feature IN ('start_codon', 'stop_codon');

-- TSS identity is physical. Genes and transcripts attach through
-- transcript_tss, so shared starts are represented once without losing their
-- biological associations. A TSS is a one-base BED interval on both strands.
CREATE TEMP TABLE transcription_start_site AS
SELECT DISTINCT
    md5(concat_ws('|', genome_id, annotation_release, chrom, tss::VARCHAR,
                  strand)) AS tss_id,
    genome_id,
    annotation_release,
    chrom,
    tss::BIGINT AS start,
    (tss + 1)::BIGINT AS "end",
    strand
FROM transcript_annotation;

CREATE TEMP TABLE transcript_tss AS
SELECT
    t.genome_id,
    t.annotation_release,
    t.gene_id,
    t.gene_name,
    t.transcript_id,
    s.tss_id
FROM transcript_annotation t
JOIN transcription_start_site s
  ON s.genome_id = t.genome_id
 AND s.annotation_release = t.annotation_release
 AND s.chrom = t.chrom
 AND s.start = t.tss
 AND s.strand = t.strand;

-- TES identity is physical and normalized exactly like the TSS dimension.
-- Transcript and gene ownership stays in transcript_tes because alternative
-- isoforms and distinct genes can share the same terminal base.
CREATE TEMP TABLE transcription_end_site AS
SELECT DISTINCT
    md5(concat_ws('|', genome_id, annotation_release, chrom, tes::VARCHAR,
                  strand)) AS tes_id,
    genome_id,
    annotation_release,
    chrom,
    tes::BIGINT AS start,
    (tes + 1)::BIGINT AS "end",
    strand
FROM transcript_annotation;

CREATE TEMP TABLE transcript_tes AS
SELECT
    t.genome_id,
    t.annotation_release,
    t.gene_id,
    t.gene_name,
    t.transcript_id,
    s.tes_id
FROM transcript_annotation t
JOIN transcription_end_site s
  ON s.genome_id = t.genome_id
 AND s.annotation_release = t.annotation_release
 AND s.chrom = t.chrom
 AND s.start = t.tes
 AND s.strand = t.strand;

-- A promoter is a resolved genomic interval for one physical TSS under one
-- versioned definition. Offsets include the TSS base and are oriented in the
-- direction of transcription.
CREATE TEMP TABLE promoter_annotation AS
SELECT
    md5(concat_ws('|', tss_id, {promoter_definition_id})) AS promoter_id,
    tss_id,
    genome_id,
    annotation_release,
    {promoter_definition_id}::VARCHAR AS promoter_definition_id,
    chrom,
    strand,
    CASE WHEN strand = '+'
         THEN GREATEST(0, start - {promoter_upstream_bp})
         ELSE GREATEST(0, start - {promoter_downstream_bp})
    END::BIGINT AS promoter_start,
    CASE WHEN strand = '+'
         THEN start + {promoter_downstream_bp} + 1
         ELSE start + {promoter_upstream_bp} + 1
    END::BIGINT AS promoter_end,
    start AS tss_start,
    "end" AS tss_end,
    {promoter_upstream_bp}::INTEGER AS upstream_bp,
    {promoter_downstream_bp}::INTEGER AS downstream_bp
FROM transcription_start_site;

CREATE TEMP TABLE promoter_gene AS
SELECT
    p.genome_id,
    p.annotation_release,
    p.promoter_definition_id,
    p.promoter_id,
    tt.gene_id,
    MAX(tt.gene_name) AS gene_name,
    COUNT(DISTINCT tt.transcript_id)::BIGINT AS n_transcripts
FROM promoter_annotation p
JOIN transcript_tss tt USING (genome_id, annotation_release, tss_id)
GROUP BY p.genome_id, p.annotation_release, p.promoter_definition_id,
         p.promoter_id, tt.gene_id;

-- This region mirrors the promoter window around a physical TES: 500 bases
-- toward the transcript body and 2,000 bases beyond it by default. Offsets are
-- transcription-oriented and include the TES base itself.
CREATE TEMP TABLE downstream_region_annotation AS
SELECT
    md5(concat_ws('|', tes_id, {downstream_definition_id}))
        AS downstream_region_id,
    tes_id,
    genome_id,
    annotation_release,
    {downstream_definition_id}::VARCHAR AS downstream_definition_id,
    chrom,
    strand,
    CASE WHEN strand = '+'
         THEN GREATEST(0, start - {downstream_upstream_bp})
         ELSE GREATEST(0, start - {downstream_downstream_bp})
    END::BIGINT AS downstream_start,
    CASE WHEN strand = '+'
         THEN start + {downstream_downstream_bp} + 1
         ELSE start + {downstream_upstream_bp} + 1
    END::BIGINT AS downstream_end,
    start AS tes_start,
    "end" AS tes_end,
    {downstream_upstream_bp}::INTEGER AS upstream_bp,
    {downstream_downstream_bp}::INTEGER AS downstream_bp
FROM transcription_end_site;

CREATE TEMP TABLE downstream_region_gene AS
SELECT
    d.genome_id,
    d.annotation_release,
    d.downstream_definition_id,
    d.downstream_region_id,
    tt.gene_id,
    MAX(tt.gene_name) AS gene_name,
    COUNT(DISTINCT tt.transcript_id)::BIGINT AS n_transcripts
FROM downstream_region_annotation d
JOIN transcript_tes tt USING (genome_id, annotation_release, tes_id)
GROUP BY d.genome_id, d.annotation_release, d.downstream_definition_id,
         d.downstream_region_id, tt.gene_id;

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
    FROM exon_annotation
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

COPY transcription_start_site
TO 'tables/jaspar2026/transcription_start_site.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY transcript_tss TO 'tables/jaspar2026/transcript_tss.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY transcription_end_site
TO 'tables/jaspar2026/transcription_end_site.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY transcript_tes TO 'tables/jaspar2026/transcript_tes.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY gene_annotation TO 'tables/jaspar2026/gene.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY exon_annotation TO 'tables/jaspar2026/exon.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY coding_sequence_segment
TO 'tables/jaspar2026/coding_sequence_segment.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY transcript_cds TO 'tables/jaspar2026/transcript_cds.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY transcript_coding_span
TO 'tables/jaspar2026/transcript_coding_span.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY coding_landmark TO 'tables/jaspar2026/coding_landmark.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY promoter_annotation TO 'tables/jaspar2026/promoter.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY promoter_gene TO 'tables/jaspar2026/promoter_gene.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY downstream_region_annotation
TO 'tables/jaspar2026/downstream_region.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY downstream_region_gene
TO 'tables/jaspar2026/downstream_region_gene.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
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
        tt.tss_id,
        t.strand AS transcript_strand,
        t.tss,
        CASE WHEN t.strand = '+' THEN a.center_bp - t.tss
             ELSE t.tss - a.center_bp END AS signed_tss_distance_bp,
        a.start >= t.transcript_start AND a."end" <= t.transcript_end
            AS fully_within_transcript,
        true AS overlaps_transcript,
        COALESCE(BOOL_OR(a.start >= e.start AND a."end" <= e."end"), false)
            AS fully_within_exon,
        COALESCE(BOOL_OR(a.start < e."end" AND a."end" > e.start), false)
            AS overlaps_exon,
        COALESCE(BOOL_OR(a.start >= cds.start AND a."end" <= cds."end"), false)
            AS fully_within_cds,
        COALESCE(BOOL_OR(a.start < cds."end" AND a."end" > cds.start), false)
            AS overlaps_cds,
        COALESCE(BOOL_OR(a.start >= i.start AND a."end" <= i."end"), false)
            AS fully_within_intron,
        COALESCE(BOOL_OR(a.start < i."end" AND a."end" > i.start), false)
            AS overlaps_intron
    FROM anchor_hit a
    JOIN transcript_annotation t
      ON a.chrom = t.chrom
     AND a.start < t.transcript_end
     AND a."end" > t.transcript_start
    JOIN transcript_tss tt
      ON tt.genome_id = t.genome_id
     AND tt.annotation_release = t.annotation_release
     AND tt.gene_id = t.gene_id
     AND tt.transcript_id = t.transcript_id
    LEFT JOIN exon_annotation e
      ON e.gene_id = t.gene_id
     AND e.transcript_id = t.transcript_id
     AND a.start < e."end"
     AND a."end" > e.start
    LEFT JOIN (
        SELECT tc.gene_id, tc.transcript_id, s.start, s."end"
        FROM transcript_cds tc
        JOIN coding_sequence_segment s
          USING (genome_id, annotation_release, cds_segment_id)
    ) cds
      ON cds.gene_id = t.gene_id
     AND cds.transcript_id = t.transcript_id
     AND a.start < cds."end"
     AND a."end" > cds.start
    LEFT JOIN intron_annotation i
      ON i.gene_id = t.gene_id
     AND i.transcript_id = t.transcript_id
     AND a.start < i."end"
     AND a."end" > i.start
    GROUP BY a.anchor_hit_id, a.chrom, a.start, a."end", t.gene_id,
             t.gene_name, t.transcript_id, tt.tss_id, t.strand, t.tss,
             t.transcript_start, t.transcript_end, a.center_bp
)
SELECT *,
       CASE WHEN fully_within_cds THEN 'cds'
            WHEN overlaps_cds THEN 'cds_boundary'
            WHEN fully_within_exon THEN 'exonic_non_cds'
            WHEN overlaps_exon THEN 'exon_boundary'
            WHEN fully_within_intron THEN 'intron'
            WHEN overlaps_intron THEN 'intron_boundary'
            WHEN fully_within_transcript THEN 'transcribed_non_exonic'
            ELSE 'transcript_boundary' END AS transcript_region
FROM overlap;

COPY anchor_transcript_context TO 'tables/jaspar2026/motif_transcript_context.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);

CREATE TEMP TABLE unique_tss_coordinate AS
SELECT DISTINCT chrom, start AS tss
FROM transcription_start_site;

CREATE TEMP TABLE nearest_left AS
SELECT a.anchor_hit_id, a.center_bp, t.tss
FROM anchor_hit a
ASOF LEFT JOIN unique_tss_coordinate t
  ON a.chrom = t.chrom AND a.center_bp >= t.tss;

CREATE TEMP TABLE nearest_right AS
SELECT a.anchor_hit_id, a.center_bp, t.tss
FROM anchor_hit a
ASOF LEFT JOIN unique_tss_coordinate t
  ON a.chrom = t.chrom AND a.center_bp <= t.tss;

CREATE TEMP TABLE nearest_tss_coordinate AS
SELECT * EXCLUDE (nearest_rank)
FROM (
    SELECT *, DENSE_RANK() OVER (
        PARTITION BY anchor_hit_id
        ORDER BY ABS(center_bp - tss)
    ) AS nearest_rank
    FROM (
        SELECT * FROM nearest_left WHERE tss IS NOT NULL
        UNION
        SELECT * FROM nearest_right WHERE tss IS NOT NULL
    )
)
WHERE nearest_rank = 1;

-- Authoritative nearest relation. Every physical TSS tied for the minimum
-- distance is retained, including opposite-strand starts at one coordinate.
CREATE TEMP TABLE tp73_anchor_nearest_tss AS
WITH geometry AS (
    SELECT
        a.anchor_hit_id,
        a.genome_id,
        a.motif_set_id,
        a.chrom,
        a.start AS anchor_start,
        a."end" AS anchor_end,
        a.center_bp AS anchor_center_bp,
        s.tss_id,
        s.annotation_release,
        s.start AS tss_start,
        s."end" AS tss_end,
        s.strand AS tss_strand,
        a.center_bp - s.start AS genomic_center_offset_bp,
        CASE WHEN s.strand = '+' THEN a.center_bp - s.start
             ELSE s.start - a.center_bp END
            AS transcription_oriented_center_offset_bp,
        CASE
            WHEN a."end" <= s.start THEN s.start - a."end"
            WHEN s."end" <= a.start THEN a.start - s."end"
            ELSE -(LEAST(a."end", s."end") - GREATEST(a.start, s.start))
        END::BIGINT AS tss_interval_distance_bp,
        a.start <= s.start AND a."end" > s.start AS anchor_spans_tss
    FROM anchor_hit a
    JOIN nearest_tss_coordinate n USING (anchor_hit_id)
    JOIN transcription_start_site s
      ON s.chrom = a.chrom AND s.start = n.tss
)
SELECT
    *,
    CASE
        WHEN anchor_spans_tss THEN 'spans_tss'
        WHEN transcription_oriented_center_offset_bp < 0 THEN 'upstream'
        WHEN transcription_oriented_center_offset_bp > 0 THEN 'downstream'
        ELSE 'coincident_center'
    END AS anchor_tss_relation,
    COUNT(*) OVER (PARTITION BY anchor_hit_id)::BIGINT
        AS nearest_tss_tie_count,
    COUNT(DISTINCT tss_strand) OVER (PARTITION BY anchor_hit_id) > 1
        AS nearest_tss_has_mixed_strands
FROM geometry;

COPY tp73_anchor_nearest_tss
TO 'tables/jaspar2026/tp73_anchor_nearest_tss.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);

-- Find the nearest physical CDS segment without collapsing its transcript or
-- gene ownership. Overlapping segments take precedence over abutting ones;
-- otherwise the established signed interval-distance convention applies.
CREATE TEMP TABLE unique_cds_end AS
SELECT DISTINCT chrom, "end" FROM coding_sequence_segment;

CREATE TEMP TABLE unique_cds_start AS
SELECT DISTINCT chrom, start FROM coding_sequence_segment;

CREATE TEMP TABLE nearest_left_cds_coordinate AS
SELECT a.anchor_hit_id, c."end" AS cds_end
FROM anchor_hit a
ASOF LEFT JOIN unique_cds_end c
  ON a.chrom = c.chrom AND a.start >= c."end";

CREATE TEMP TABLE nearest_right_cds_coordinate AS
SELECT a.anchor_hit_id, c.start AS cds_start
FROM anchor_hit a
ASOF LEFT JOIN unique_cds_start c
  ON a.chrom = c.chrom AND a."end" <= c.start;

CREATE TEMP TABLE tp73_anchor_nearest_cds AS
WITH candidate AS (
    SELECT a.anchor_hit_id, s.*
    FROM anchor_hit a
    JOIN coding_sequence_segment s
      ON a.chrom = s.chrom
     AND a.start < s."end"
     AND a."end" > s.start
    UNION
    SELECT a.anchor_hit_id, s.*
    FROM anchor_hit a
    JOIN nearest_left_cds_coordinate n USING (anchor_hit_id)
    JOIN coding_sequence_segment s
      ON s.chrom = a.chrom AND s."end" = n.cds_end
    WHERE n.cds_end IS NOT NULL
    UNION
    SELECT a.anchor_hit_id, s.*
    FROM anchor_hit a
    JOIN nearest_right_cds_coordinate n USING (anchor_hit_id)
    JOIN coding_sequence_segment s
      ON s.chrom = a.chrom AND s.start = n.cds_start
    WHERE n.cds_start IS NOT NULL
), geometry AS (
    SELECT
        a.anchor_hit_id,
        a.genome_id,
        a.motif_set_id,
        a.chrom,
        a.start AS anchor_start,
        a."end" AS anchor_end,
        a.center_bp AS anchor_center_bp,
        c.cds_segment_id,
        c.annotation_release,
        c.start AS cds_start,
        c."end" AS cds_end,
        c.strand AS cds_strand,
        a.center_bp - ((c.start + c."end") / 2.0)
            AS genomic_center_offset_bp,
        CASE WHEN c.strand = '+'
             THEN a.center_bp - ((c.start + c."end") / 2.0)
             ELSE ((c.start + c."end") / 2.0) - a.center_bp END
            AS coding_oriented_center_offset_bp,
        CASE
            WHEN a."end" <= c.start THEN c.start - a."end"
            WHEN c."end" <= a.start THEN a.start - c."end"
            ELSE -(LEAST(a."end", c."end") - GREATEST(a.start, c.start))
        END::BIGINT AS cds_interval_distance_bp,
        GREATEST(
            0,
            CASE
                WHEN a."end" <= c.start THEN c.start - a."end"
                WHEN c."end" <= a.start THEN a.start - c."end"
                ELSE 0
            END
        )::BIGINT AS cds_genomic_distance_bp,
        GREATEST(
            0, LEAST(a."end", c."end") - GREATEST(a.start, c.start)
        )::BIGINT AS cds_overlap_bp
    FROM candidate c
    JOIN anchor_hit a USING (anchor_hit_id)
), nearest AS (
    SELECT *,
        DENSE_RANK() OVER (
            PARTITION BY anchor_hit_id
            ORDER BY CASE WHEN cds_interval_distance_bp < 0 THEN 0
                          ELSE cds_interval_distance_bp + 1 END
        ) AS nearest_rank
    FROM geometry
)
SELECT
    * EXCLUDE (nearest_rank),
    CASE
        WHEN cds_overlap_bp > 0 THEN 'overlaps_cds'
        WHEN coding_oriented_center_offset_bp < 0 THEN 'upstream'
        WHEN coding_oriented_center_offset_bp > 0 THEN 'downstream'
        ELSE 'coincident_center'
    END AS anchor_cds_relation,
    COUNT(*) OVER (PARTITION BY anchor_hit_id)::BIGINT
        AS nearest_cds_tie_count,
    COUNT(DISTINCT cds_strand) OVER (PARTITION BY anchor_hit_id) > 1
        AS nearest_cds_has_mixed_strands
FROM nearest
WHERE nearest_rank = 1;

COPY tp73_anchor_nearest_cds
TO 'tables/jaspar2026/tp73_anchor_nearest_cds.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);

-- Canonical many-to-many membership. Gene association remains normalized via
-- promoter_gene so shared physical promoters are not duplicated here.
CREATE TEMP TABLE tp73_anchor_promoter AS
WITH geometry AS (
    SELECT
        a.anchor_hit_id,
        a.genome_id,
        a.motif_set_id,
        a.chrom,
        a.start AS anchor_start,
        a."end" AS anchor_end,
        a.center_bp AS anchor_center_bp,
        p.promoter_id,
        p.tss_id,
        p.annotation_release,
        p.promoter_definition_id,
        p.promoter_start,
        p.promoter_end,
        p.tss_start,
        p.tss_end,
        p.strand AS tss_strand,
        LEAST(a."end", p.promoter_end)
            - GREATEST(a.start, p.promoter_start) AS promoter_overlap_bp,
        a.start >= p.promoter_start AND a."end" <= p.promoter_end
            AS anchor_fully_within_promoter,
        CASE
            WHEN a."end" <= p.tss_start THEN p.tss_start - a."end"
            WHEN p.tss_end <= a.start THEN a.start - p.tss_end
            ELSE -(LEAST(a."end", p.tss_end) - GREATEST(a.start, p.tss_start))
        END::BIGINT AS tss_interval_distance_bp,
        CASE WHEN p.strand = '+' THEN a.center_bp - p.tss_start
             ELSE p.tss_start - a.center_bp END
            AS transcription_oriented_center_offset_bp,
        a.start <= p.tss_start AND a."end" > p.tss_start AS anchor_spans_tss
    FROM anchor_hit a
    JOIN promoter_annotation p
      ON p.chrom = a.chrom
     AND a.start < p.promoter_end
     AND a."end" > p.promoter_start
)
SELECT
    *,
    CASE
        WHEN anchor_spans_tss THEN 'spans_tss'
        WHEN transcription_oriented_center_offset_bp < 0 THEN 'upstream'
        WHEN transcription_oriented_center_offset_bp > 0 THEN 'downstream'
        ELSE 'coincident_center'
    END AS anchor_tss_relation
FROM geometry;

COPY tp73_anchor_promoter
TO 'tables/jaspar2026/tp73_anchor_promoter.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);

-- Canonical many-to-many membership in a versioned transcript-end downstream
-- region. Gene ownership remains normalized through downstream_region_gene.
CREATE TEMP TABLE tp73_anchor_downstream_region AS
WITH geometry AS (
    SELECT
        a.anchor_hit_id,
        a.genome_id,
        a.motif_set_id,
        a.chrom,
        a.start AS anchor_start,
        a."end" AS anchor_end,
        a.center_bp AS anchor_center_bp,
        d.downstream_region_id,
        d.tes_id,
        d.annotation_release,
        d.downstream_definition_id,
        d.downstream_start,
        d.downstream_end,
        d.tes_start,
        d.tes_end,
        d.strand AS tes_strand,
        LEAST(a."end", d.downstream_end)
            - GREATEST(a.start, d.downstream_start) AS downstream_overlap_bp,
        a.start >= d.downstream_start AND a."end" <= d.downstream_end
            AS anchor_fully_within_downstream_region,
        CASE
            WHEN a."end" <= d.tes_start THEN d.tes_start - a."end"
            WHEN d.tes_end <= a.start THEN a.start - d.tes_end
            ELSE -(LEAST(a."end", d.tes_end) - GREATEST(a.start, d.tes_start))
        END::BIGINT AS tes_interval_distance_bp,
        CASE WHEN d.strand = '+' THEN a.center_bp - d.tes_start
             ELSE d.tes_start - a.center_bp END
            AS transcription_oriented_center_offset_bp,
        a.start <= d.tes_start AND a."end" > d.tes_start AS anchor_spans_tes
    FROM anchor_hit a
    JOIN downstream_region_annotation d
      ON d.chrom = a.chrom
     AND a.start < d.downstream_end
     AND a."end" > d.downstream_start
)
SELECT
    *,
    CASE
        WHEN anchor_spans_tes THEN 'spans_tes'
        WHEN transcription_oriented_center_offset_bp < 0 THEN 'upstream'
        WHEN transcription_oriented_center_offset_bp > 0 THEN 'downstream'
        ELSE 'coincident_center'
    END AS anchor_tes_relation
FROM geometry;

COPY tp73_anchor_downstream_region
TO 'tables/jaspar2026/tp73_anchor_downstream_region.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);

CREATE TEMP TABLE nearest_tss_summary AS
SELECT * EXCLUDE (representative_rank, n_tss_strands)
FROM (
    SELECT
        n.anchor_hit_id,
        n.tss_id AS nearest_tss_id,
        n.tss_start AS nearest_tss_start,
        CASE WHEN COUNT(DISTINCT n.tss_strand) OVER anchor_partition = 1
             THEN n.tss_strand ELSE NULL END AS nearest_tss_strand,
        n.nearest_tss_tie_count,
        n.nearest_tss_has_mixed_strands,
        ABS(n.genomic_center_offset_bp) AS nearest_tss_genomic_distance_bp,
        n.transcription_oriented_center_offset_bp AS nearest_tss_distance_bp,
        n.tss_interval_distance_bp AS nearest_tss_interval_distance_bp,
        n.anchor_tss_relation AS nearest_tss_relation,
        tt.gene_id,
        tt.gene_name,
        tt.transcript_id,
        COUNT(DISTINCT n.tss_strand) OVER anchor_partition AS n_tss_strands,
        ROW_NUMBER() OVER (
            anchor_partition
            ORDER BY n.tss_start, n.tss_strand, n.tss_id,
                     tt.gene_id, tt.transcript_id
        ) AS representative_rank
    FROM tp73_anchor_nearest_tss n
    JOIN transcript_tss tt USING (genome_id, tss_id)
    WINDOW anchor_partition AS (PARTITION BY n.anchor_hit_id)
)
WHERE representative_rank = 1;

CREATE TEMP TABLE nearest_cds_summary AS
SELECT * EXCLUDE (representative_rank)
FROM (
    SELECT
        n.anchor_hit_id,
        n.cds_segment_id AS nearest_cds_segment_id,
        n.cds_start AS nearest_cds_start,
        n.cds_end AS nearest_cds_end,
        CASE WHEN n.nearest_cds_has_mixed_strands
             THEN NULL ELSE n.cds_strand END AS nearest_cds_strand,
        n.nearest_cds_tie_count,
        n.nearest_cds_has_mixed_strands,
        n.cds_genomic_distance_bp AS nearest_cds_genomic_distance_bp,
        n.cds_interval_distance_bp AS nearest_cds_interval_distance_bp,
        n.cds_overlap_bp AS nearest_cds_overlap_bp,
        n.coding_oriented_center_offset_bp AS nearest_cds_distance_bp,
        n.anchor_cds_relation AS nearest_cds_relation,
        tc.gene_id AS nearest_cds_gene_id,
        tc.gene_name AS nearest_cds_gene_name,
        tc.transcript_id AS nearest_cds_transcript_id,
        ROW_NUMBER() OVER (
            PARTITION BY n.anchor_hit_id
            ORDER BY n.cds_start, n.cds_end, n.cds_strand,
                     n.cds_segment_id, tc.gene_id, tc.transcript_id
        ) AS representative_rank
    FROM tp73_anchor_nearest_cds n
    LEFT JOIN transcript_cds tc
      ON tc.genome_id = n.genome_id
     AND tc.annotation_release = n.annotation_release
     AND tc.cds_segment_id = n.cds_segment_id
)
WHERE representative_rank = 1;

CREATE TEMP TABLE transcript_context_summary AS
SELECT
    anchor_hit_id,
    true AS in_any_transcript,
    BOOL_OR(fully_within_transcript) AS fully_within_any_transcript,
    BOOL_OR(fully_within_exon) AS in_any_exon,
    BOOL_OR(overlaps_exon AND NOT fully_within_exon)
        AS overlaps_any_exon_boundary,
    BOOL_OR(fully_within_cds) AS in_any_cds,
    BOOL_OR(overlaps_cds AND NOT fully_within_cds)
        AS overlaps_any_cds_boundary,
    BOOL_OR(fully_within_intron) AS in_any_intron,
    BOOL_OR(overlaps_intron AND NOT fully_within_intron) AS overlaps_any_intron_boundary
FROM anchor_transcript_context
GROUP BY anchor_hit_id;

CREATE TEMP TABLE promoter_context_summary AS
SELECT anchor_hit_id,
       true AS overlaps_any_promoter,
       COUNT(*)::BIGINT AS n_overlapping_promoters
FROM tp73_anchor_promoter
GROUP BY anchor_hit_id;

CREATE TEMP TABLE downstream_context_summary AS
SELECT anchor_hit_id,
       true AS overlaps_any_downstream_region,
       COUNT(*)::BIGINT AS n_overlapping_downstream_regions
FROM tp73_anchor_downstream_region
GROUP BY anchor_hit_id;

CREATE TEMP TABLE anchor_gene_context AS
SELECT
    a.anchor_hit_id,
    true AS gene_annotation_available,
    n.gene_id AS nearest_gene_id,
    n.gene_name AS nearest_gene_name,
    n.transcript_id AS nearest_transcript_id,
    n.nearest_tss_id,
    n.nearest_tss_start,
    n.nearest_tss_strand,
    n.nearest_tss_tie_count,
    n.nearest_tss_has_mixed_strands,
    n.nearest_tss_genomic_distance_bp,
    n.nearest_tss_distance_bp,
    n.nearest_tss_interval_distance_bp,
    n.nearest_tss_relation,
    c.nearest_cds_segment_id,
    c.nearest_cds_start,
    c.nearest_cds_end,
    c.nearest_cds_strand,
    c.nearest_cds_tie_count,
    c.nearest_cds_has_mixed_strands,
    c.nearest_cds_genomic_distance_bp,
    c.nearest_cds_interval_distance_bp,
    c.nearest_cds_overlap_bp,
    c.nearest_cds_distance_bp,
    c.nearest_cds_relation,
    c.nearest_cds_gene_id,
    c.nearest_cds_gene_name,
    c.nearest_cds_transcript_id,
    COALESCE(s.in_any_transcript, false) AS in_any_transcript,
    COALESCE(s.fully_within_any_transcript, false)
        AS fully_within_any_transcript,
    COALESCE(s.in_any_exon, false) AS in_any_exon,
    COALESCE(s.overlaps_any_exon_boundary, false)
        AS overlaps_any_exon_boundary,
    COALESCE(s.in_any_cds, false) AS in_any_cds,
    COALESCE(s.overlaps_any_cds_boundary, false)
        AS overlaps_any_cds_boundary,
    COALESCE(s.in_any_intron, false) AS in_any_intron,
    COALESCE(s.overlaps_any_intron_boundary, false) AS overlaps_any_intron_boundary,
    COALESCE(p.overlaps_any_promoter, false) AS overlaps_any_promoter,
    COALESCE(p.n_overlapping_promoters, 0)::BIGINT AS n_overlapping_promoters,
    COALESCE(d.overlaps_any_downstream_region, false)
        AS overlaps_any_downstream_region,
    COALESCE(d.n_overlapping_downstream_regions, 0)::BIGINT
        AS n_overlapping_downstream_regions,
    NOT COALESCE(s.in_any_transcript, false)
        AND NOT COALESCE(p.overlaps_any_promoter, false)
        AND NOT COALESCE(d.overlaps_any_downstream_region, false)
        AS strict_intergenic,
    CASE WHEN COALESCE(s.in_any_intron, false) THEN 'intron'
         WHEN COALESCE(s.overlaps_any_intron_boundary, false) THEN 'intron_boundary'
         WHEN COALESCE(s.in_any_transcript, false) THEN 'transcribed_non_intron'
         ELSE 'intergenic' END AS primary_transcript_region,
    CASE
        WHEN COALESCE(p.overlaps_any_promoter, false)
             AND COALESCE(s.in_any_cds, false) THEN 'promoter_and_cds'
        WHEN COALESCE(p.overlaps_any_promoter, false)
             AND COALESCE(s.in_any_transcript, false) THEN 'promoter_and_transcribed'
        WHEN COALESCE(p.overlaps_any_promoter, false) THEN 'promoter_only'
        WHEN COALESCE(d.overlaps_any_downstream_region, false)
             AND COALESCE(s.in_any_cds, false) THEN 'downstream_and_cds'
        WHEN COALESCE(d.overlaps_any_downstream_region, false)
             AND COALESCE(s.in_any_transcript, false)
            THEN 'downstream_and_transcribed'
        WHEN COALESCE(d.overlaps_any_downstream_region, false)
            THEN 'downstream_only'
        WHEN COALESCE(s.in_any_cds, false) THEN 'cds'
        WHEN COALESCE(s.overlaps_any_cds_boundary, false) THEN 'cds_boundary'
        WHEN COALESCE(s.in_any_exon, false) THEN 'exonic_non_cds'
        WHEN COALESCE(s.overlaps_any_exon_boundary, false) THEN 'exon_boundary'
        WHEN COALESCE(s.in_any_intron, false) THEN 'intron'
        WHEN COALESCE(s.overlaps_any_intron_boundary, false) THEN 'intron_boundary'
        WHEN COALESCE(s.in_any_transcript, false) THEN 'transcript_boundary'
        ELSE 'strict_intergenic'
    END AS primary_genomic_context,
    CASE
        WHEN COALESCE(p.overlaps_any_promoter, false) THEN 'promoter'
        WHEN COALESCE(d.overlaps_any_downstream_region, false) THEN 'downstream'
        WHEN COALESCE(s.in_any_transcript, false) THEN 'gene_body'
        ELSE 'intergenic'
    END AS gene_relation_class
FROM anchor_hit a
LEFT JOIN nearest_tss_summary n USING (anchor_hit_id)
LEFT JOIN nearest_cds_summary c USING (anchor_hit_id)
LEFT JOIN transcript_context_summary s USING (anchor_hit_id)
LEFT JOIN promoter_context_summary p USING (anchor_hit_id)
LEFT JOIN downstream_context_summary d USING (anchor_hit_id);
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
    NULL::VARCHAR AS nearest_tss_id,
    NULL::BIGINT AS nearest_tss_start,
    NULL::VARCHAR AS nearest_tss_strand,
    NULL::BIGINT AS nearest_tss_tie_count,
    NULL::BOOLEAN AS nearest_tss_has_mixed_strands,
    NULL::DOUBLE AS nearest_tss_genomic_distance_bp,
    NULL::DOUBLE AS nearest_tss_distance_bp,
    NULL::BIGINT AS nearest_tss_interval_distance_bp,
    NULL::VARCHAR AS nearest_tss_relation,
    NULL::VARCHAR AS nearest_cds_segment_id,
    NULL::BIGINT AS nearest_cds_start,
    NULL::BIGINT AS nearest_cds_end,
    NULL::VARCHAR AS nearest_cds_strand,
    NULL::BIGINT AS nearest_cds_tie_count,
    NULL::BOOLEAN AS nearest_cds_has_mixed_strands,
    NULL::BIGINT AS nearest_cds_genomic_distance_bp,
    NULL::BIGINT AS nearest_cds_interval_distance_bp,
    NULL::BIGINT AS nearest_cds_overlap_bp,
    NULL::DOUBLE AS nearest_cds_distance_bp,
    NULL::VARCHAR AS nearest_cds_relation,
    NULL::VARCHAR AS nearest_cds_gene_id,
    NULL::VARCHAR AS nearest_cds_gene_name,
    NULL::VARCHAR AS nearest_cds_transcript_id,
    NULL::BOOLEAN AS in_any_transcript,
    NULL::BOOLEAN AS fully_within_any_transcript,
    NULL::BOOLEAN AS in_any_exon,
    NULL::BOOLEAN AS overlaps_any_exon_boundary,
    NULL::BOOLEAN AS in_any_cds,
    NULL::BOOLEAN AS overlaps_any_cds_boundary,
    NULL::BOOLEAN AS in_any_intron,
    NULL::BOOLEAN AS overlaps_any_intron_boundary,
    NULL::BOOLEAN AS overlaps_any_promoter,
    NULL::BIGINT AS n_overlapping_promoters,
    NULL::BOOLEAN AS overlaps_any_downstream_region,
    NULL::BIGINT AS n_overlapping_downstream_regions,
    NULL::BOOLEAN AS strict_intergenic,
    'not_assessed'::VARCHAR AS primary_transcript_region,
    'not_assessed'::VARCHAR AS primary_genomic_context,
    'not_assessed'::VARCHAR AS gene_relation_class
FROM anchor_hit;
"""


def parquet_columns(executable: str, motif_glob: str) -> set[str]:
    query = (
        "DESCRIBE SELECT * FROM read_parquet("
        f"{sql_string(motif_glob)}, hive_partitioning=1);"
    )
    process = subprocess.run(
        [executable, "-csv", "-noheader", ":memory:", "-c", query],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        detail = process.stderr.strip() or "DuckDB DESCRIBE failed"
        raise ContextBuildError(f"could not inspect motif-hit schema: {detail}")
    return {row[0] for row in csv.reader(process.stdout.splitlines()) if row}


def build_sql(arguments: argparse.Namespace, motif_glob: str,
              chromosomes: list[str], gtf: Path | None,
              temp_directory: Path,
              source_columns: set[str]) -> str:
    chrom_clause = chromosome_filter(chromosomes)
    annotation_sql = (
        gtf_sql(arguments, gtf, chromosomes)
        if gtf is not None else no_gtf_sql()
    )
    gtf_source = sql_string(gtf) if gtf is not None else "NULL"
    gtf_sha256 = (
        sql_string(arguments.gtf_sha256)
        if arguments.gtf_sha256 is not None else "NULL"
    )
    gtf_size_bytes = (
        str(arguments.gtf_size_bytes)
        if arguments.gtf_size_bytes is not None else "NULL"
    )
    motif_name = (
        "COALESCE(NULLIF(motif_name::VARCHAR, ''), motif_id::VARCHAR)"
        if "motif_name" in source_columns else "motif_id::VARCHAR"
    )
    qualifying_threshold = (
        "COALESCE(TRY_CAST(minimum_score AS DOUBLE), "
        f"{sql_number(arguments.default_neighbor_minimum_score)})"
        if "minimum_score" in source_columns
        else f"{sql_number(arguments.default_neighbor_minimum_score)}::DOUBLE"
    )
    genome_id = (
        "genome_id::VARCHAR" if "genome_id" in source_columns
        else f"{sql_string(arguments.genome_id)}::VARCHAR"
    )
    motif_set_id = (
        "motif_set_id::VARCHAR" if "motif_set_id" in source_columns
        else f"{sql_string(arguments.motif_set_id)}::VARCHAR"
    )
    identity_filters = ""
    if "genome_id" in source_columns:
        identity_filters += (
            f"\n      AND genome_id = {sql_string(arguments.genome_id)}"
        )
    if "motif_set_id" in source_columns:
        identity_filters += (
            f"\n      AND motif_set_id = {sql_string(arguments.motif_set_id)}"
        )
    if "background_model_id" in source_columns:
        identity_filters += (
            "\n      AND background_model_id = "
            f"{sql_string(arguments.background_model_id)}"
        )
    if "pseudocount_scheme" in source_columns:
        identity_filters += (
            "\n      AND pseudocount_scheme = "
            f"{sql_string(arguments.pseudocount_scheme)}"
        )
    if arguments.input_uniqueness == "deduplicate":
        duplicate_rank_sql = """
        ROW_NUMBER() OVER (
            PARTITION BY CAST(chrom AS VARCHAR), start, "end", motif_id, strand,
                         score_mode, pseudocount
            ORDER BY score DESC NULLS LAST, pwm_relative_score DESC NULLS LAST
        )"""
    else:
        duplicate_rank_sql = "1::BIGINT"
    if arguments.output_tier == "selected":
        motif_context_copy_sql = """
COPY motif_context_pair TO 'tables/jaspar2026/motif_context_pair'
    (FORMAT PARQUET, PARTITION_BY (genome_id, chrom), COMPRESSION ZSTD);
"""
        cofactor_locus_copy_sql = """
COPY cofactor_motif_locus
TO 'tables/jaspar2026/cofactor_motif_locus.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        cofactor_pair_copy_sql = """
COPY cofactor_motif_pair
TO 'tables/jaspar2026/cofactor_motif_pair.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        cofactor_locus_feature_copy_sql = """
COPY cofactor_locus_pair_feature
TO 'tables/jaspar2026/cofactor_locus_pair_feature.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        cofactor_pair_context_copy_sql = """
COPY tp73_cofactor_pair_context
TO 'tables/jaspar2026/tp73_cofactor_pair_context.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
    else:
        motif_context_copy_sql = """
COPY (SELECT * FROM motif_context_pair WHERE false)
TO 'tables/jaspar2026/motif_context_pair/empty.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        cofactor_locus_copy_sql = """
COPY (SELECT * FROM cofactor_motif_locus WHERE false)
TO 'tables/jaspar2026/cofactor_motif_locus.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        cofactor_pair_copy_sql = """
COPY (SELECT * FROM cofactor_motif_pair WHERE false)
TO 'tables/jaspar2026/cofactor_motif_pair.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        cofactor_locus_feature_copy_sql = """
COPY (SELECT * FROM cofactor_locus_pair_feature WHERE false)
TO 'tables/jaspar2026/cofactor_locus_pair_feature.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        cofactor_pair_context_copy_sql = """
COPY (SELECT * FROM tp73_cofactor_pair_context WHERE false)
TO 'tables/jaspar2026/tp73_cofactor_pair_context.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
    # All-JASPAR production runs motifs in batches. Keep the one batch-specific
    # feature table and emit typed empty files for shared TP73 surfaces so every
    # package remains queryable without multiplying those rows across batches.
    if arguments.output_tier == "band":
        anchor_locus_copy_sql = """
COPY (SELECT * FROM tp73_anchor_locus WHERE false)
TO 'tables/jaspar2026/tp73_anchor_locus/empty.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        motif_summary_copy_sql = """
COPY (SELECT * FROM tp73_motif_context_summary WHERE false)
TO 'tables/jaspar2026/tp73_motif_context_summary/empty.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        band_qualifying_filter_sql = (
            "\n      AND neighbor_motif_id <> "
            + sql_string(arguments.anchor_motif)
        )
        band_feature_copy_sql = """
COPY anchor_motif_band_feature
TO 'tables/jaspar2026/anchor_motif_band_feature/data.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        cofactor_pair_summary_copy_sql = """
COPY (SELECT * FROM tp73_cofactor_pair_summary WHERE false)
TO 'tables/jaspar2026/tp73_cofactor_pair_summary.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        pair_feature_copy_sql = """
COPY (SELECT * FROM tp73_pair_feature WHERE false)
TO 'tables/jaspar2026/tp73_pair_feature/empty.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        context_anchor_copy_sql = """
COPY (SELECT * FROM tp73_context_anchor WHERE false)
TO 'tables/jaspar2026/tp73_context_anchor/empty.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
    else:
        anchor_locus_copy_sql = """
COPY tp73_anchor_locus TO 'tables/jaspar2026/tp73_anchor_locus'
    (FORMAT PARQUET, PARTITION_BY (genome_id, chrom), COMPRESSION ZSTD);
"""
        motif_summary_copy_sql = """
COPY tp73_motif_context_summary
TO 'tables/jaspar2026/tp73_motif_context_summary'
    (FORMAT PARQUET,
     PARTITION_BY (genome_id, chrom, neighbor_motif_id), COMPRESSION ZSTD);
"""
        band_qualifying_filter_sql = ""
        band_feature_copy_sql = """
COPY anchor_motif_band_feature
TO 'tables/jaspar2026/anchor_motif_band_feature'
    (FORMAT PARQUET,
     PARTITION_BY (genome_id, chrom, neighbor_motif_id), COMPRESSION ZSTD);
"""
        cofactor_pair_summary_copy_sql = """
COPY tp73_cofactor_pair_summary
TO 'tables/jaspar2026/tp73_cofactor_pair_summary.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        pair_feature_copy_sql = """
COPY tp73_pair_feature TO 'tables/jaspar2026/tp73_pair_feature'
    (FORMAT PARQUET, PARTITION_BY (genome_id, chrom), COMPRESSION ZSTD);
"""
        context_anchor_copy_sql = """
COPY tp73_context_anchor TO 'tables/jaspar2026/tp73_context_anchor'
    (FORMAT PARQUET, PARTITION_BY (genome_id, chrom), COMPRESSION ZSTD);
"""
    motif_hit_source = arguments.motif_hit_source_label or motif_glob
    return f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;
SET temp_directory={sql_string(temp_directory)};
SET max_temp_directory_size={sql_string(arguments.max_temp_size)};
PRAGMA enable_progress_bar;

-- Arbitrary input may contain overlapping/retried Parquet parts. Deduplicate
-- those rows by default; finalized inventory selections explicitly bypass the
-- global window sort because their orientation/model-span identity is unique.
CREATE TEMP TABLE configured_hit AS
SELECT * EXCLUDE (duplicate_rank)
FROM (
    SELECT
        md5(concat_ws('|', CAST(chrom AS VARCHAR), start::VARCHAR, "end"::VARCHAR,
                      {genome_id}, {motif_set_id}, motif_id, strand, score_mode,
                      printf('%.17g', pseudocount::DOUBLE)))
            AS hit_id,
        {genome_id} AS genome_id,
        {motif_set_id} AS motif_set_id,
        CAST(chrom AS VARCHAR) AS chrom,
        start::BIGINT AS start,
        "end"::BIGINT AS "end",
        ("end"::BIGINT - start::BIGINT)::BIGINT AS span_bp,
        (start::DOUBLE + "end"::DOUBLE) / 2.0 AS center_bp,
        motif_id::VARCHAR AS motif_id,
        {motif_name} AS motif_name,
        CASE lower(strand::VARCHAR)
            WHEN 'plus' THEN '+'
            WHEN 'minus' THEN '-'
            ELSE strand::VARCHAR
        END AS strand,
        score::DOUBLE AS score,
        score_mode::VARCHAR AS score_mode,
        pseudocount::DOUBLE AS pseudocount,
        pwm_relative_score::DOUBLE AS pwm_relative_score,
        {qualifying_threshold} AS qualifying_threshold,
        {duplicate_rank_sql} AS duplicate_rank
    FROM read_parquet({sql_string(motif_glob)}, hive_partitioning=1)
    WHERE score_mode = {sql_string(arguments.score_mode)}
      AND pseudocount = {sql_number(arguments.pseudocount)}
      {identity_filters}
      {chrom_clause}
)
WHERE duplicate_rank = 1;

SELECT CASE WHEN EXISTS (
    SELECT 1 FROM configured_hit WHERE strand NOT IN ('+', '-')
) THEN error('motif-hit strand must be +/-, plus/minus') END;

SELECT CASE WHEN EXISTS (
    SELECT motif_id
    FROM configured_hit
    GROUP BY motif_id
    HAVING COUNT(DISTINCT qualifying_threshold) <> 1
) THEN error('each motif must expose one neighbor qualifying threshold') END;

-- Collapse strand alternatives before local-maximum selection. In local-peak
-- mode, an orientation can represent a retained anchor only when it is the
-- best score at its physical span and that span has no strictly stronger TP73
-- locus in the configured flank.
CREATE TEMP TABLE anchor_locus_score AS
WITH collapsed AS (
    SELECT
        md5(concat_ws('|', genome_id, motif_set_id, chrom, start::VARCHAR,
                      "end"::VARCHAR, motif_id, score_mode,
                      printf('%.17g', pseudocount))) AS anchor_locus_id,
        genome_id,
        motif_set_id,
        chrom,
        start,
        "end",
        center_bp,
        motif_id,
        MAX(score) AS locus_best_score
    FROM configured_hit
    WHERE motif_id = {sql_string(arguments.anchor_motif)}
    GROUP BY genome_id, motif_set_id, chrom, start, "end", center_bp, motif_id,
             score_mode, pseudocount
), neighboring AS (
    SELECT
        *,
        MAX(locus_best_score) OVER (
            PARTITION BY genome_id, motif_set_id, chrom
            ORDER BY start
            RANGE BETWEEN {arguments.anchor_local_peak_flank} PRECEDING
                      AND 1 PRECEDING
        ) AS best_left_locus_score,
        MAX(locus_best_score) OVER (
            PARTITION BY genome_id, motif_set_id, chrom
            ORDER BY start
            RANGE BETWEEN 1 FOLLOWING
                      AND {arguments.anchor_local_peak_flank} FOLLOWING
        ) AS best_right_locus_score
    FROM collapsed
)
SELECT
    * EXCLUDE (best_left_locus_score, best_right_locus_score),
    GREATEST(best_left_locus_score, best_right_locus_score)
        AS best_other_locus_score,
    locus_best_score - GREATEST(best_left_locus_score, best_right_locus_score)
        AS locus_score_prominence,
    GREATEST(best_left_locus_score, best_right_locus_score) IS NULL
        OR locus_best_score >=
           GREATEST(best_left_locus_score, best_right_locus_score)
        AS is_local_peak
FROM neighboring;

CREATE TEMP TABLE anchor_hit AS
SELECT
    h.hit_id AS anchor_hit_id,
    l.anchor_locus_id,
    h.* EXCLUDE (hit_id),
    l.locus_best_score AS anchor_locus_best_score,
    l.best_other_locus_score AS best_other_anchor_locus_score,
    l.locus_score_prominence AS anchor_locus_score_prominence,
    l.is_local_peak AS anchor_locus_is_local_peak,
    CASE
        WHEN {sql_string(arguments.anchor_selection_mode)} = 'threshold'
            THEN 'threshold'
        ELSE 'local_peak'
    END AS anchor_selection_class
FROM configured_hit h
JOIN anchor_locus_score l
  ON h.genome_id = l.genome_id
 AND h.motif_set_id = l.motif_set_id
 AND h.chrom = l.chrom
 AND h.start = l.start
 AND h."end" = l."end"
 AND h.motif_id = l.motif_id
WHERE h.motif_id = {sql_string(arguments.anchor_motif)}
  AND h.score >= {sql_number(arguments.anchor_minimum_score)}
  AND (
      {sql_string(arguments.anchor_selection_mode)} = 'threshold'
      OR (
          h.score = l.locus_best_score
          AND l.is_local_peak
      )
  );

SELECT CASE WHEN (SELECT COUNT(*) FROM anchor_hit) = 0
    THEN error('no anchor motif hits matched the requested configuration') END;

DROP TABLE anchor_locus_score;

-- The center prefilter is only an acceleration device. Use the deliberately
-- conservative full-span bound agreed for schema v4, then apply exact interval
-- geometry below. One-bin expansion remains lossless because the bin width
-- equals the largest center distance admitted by this prefilter.
CREATE TEMP TABLE capture_parameter AS
WITH widths AS (
    SELECT
        (SELECT MAX(span_bp) FROM anchor_hit)::BIGINT AS max_anchor_span_bp,
        (SELECT MAX(span_bp) FROM configured_hit)::BIGINT AS max_neighbor_span_bp
)
SELECT
    max_anchor_span_bp,
    max_neighbor_span_bp,
    (
        {arguments.capture_flank}
        + max_anchor_span_bp
        + max_neighbor_span_bp
    )::BIGINT AS capture_prefilter_center_bp,
    (
        {arguments.cofactor_pair_flank}
        + (2 * max_neighbor_span_bp)
    )::BIGINT
        AS cofactor_pair_prefilter_center_bp
FROM widths;

CREATE TEMP TABLE context_run_config AS
SELECT
    9::INTEGER AS schema_version,
    {sql_string(arguments.source_commit)}::VARCHAR AS builder_source_commit,
    {sql_string(arguments.input_uniqueness)}::VARCHAR AS input_uniqueness,
    {sql_string(arguments.genome_id)}::VARCHAR AS genome_id,
    {sql_string(arguments.motif_set_id)}::VARCHAR AS motif_set_id,
    {sql_string(arguments.anchor_motif)}::VARCHAR AS anchor_motif_id,
    {sql_string(arguments.score_mode)}::VARCHAR AS score_mode,
    {sql_number(arguments.pseudocount)}::DOUBLE AS pseudocount,
    {sql_string(arguments.background_model_id)}::VARCHAR AS background_model_id,
    {sql_string(arguments.pseudocount_scheme)}::VARCHAR AS pseudocount_scheme,
    {sql_number(arguments.anchor_minimum_score)}::DOUBLE AS anchor_minimum_score,
    {sql_number(arguments.partner_minimum_score)}::DOUBLE AS partner_minimum_score,
    {sql_number(arguments.default_neighbor_minimum_score)}::DOUBLE
        AS default_neighbor_minimum_score,
    'source_minimum_score_else_default_inclusive'::VARCHAR
        AS neighbor_qualifying_threshold_rule,
    {sql_string(arguments.anchor_selection_mode)}::VARCHAR
        AS anchor_selection_mode,
    {arguments.anchor_local_peak_flank}::INTEGER AS anchor_local_peak_flank_bp,
    'physical_locus_best_score_no_stronger_start_within_flank'::VARCHAR
        AS anchor_local_peak_rule,
    'both_orientation_specific_scores_at_or_above_tandem_minimum'::VARCHAR
        AS tandem_score_rule,
    {arguments.capture_flank}::INTEGER AS capture_flank_bp,
    {arguments.context_flank}::INTEGER AS context_flank_bp,
    {arguments.tandem_flank}::INTEGER AS tandem_flank_bp,
    {arguments.cofactor_pair_flank}::INTEGER AS cofactor_pair_flank_bp,
    {sql_string(arguments.output_tier)}::VARCHAR AS output_tier,
    'interval'::VARCHAR AS capture_geometry,
    'signed_interval_edge_distance'::VARCHAR AS distance_metric,
    'motif_center'::VARCHAR AS center_distance_metric,
    'capture_plus_observed_anchor_and_neighbor_spans'::VARCHAR
        AS capture_prefilter_rule,
    p.capture_prefilter_center_bp,
    p.cofactor_pair_prefilter_center_bp,
    p.max_anchor_span_bp AS observed_max_anchor_span_bp,
    p.max_neighbor_span_bp AS observed_max_neighbor_span_bp,
    'overlap|adjacent_0_5|gap_6_20|gap_21_50|gap_51_100|gap_101_150'
        ::VARCHAR AS distance_band_rule,
    'nonoverlapping_edge_gap'::VARCHAR AS tandem_distance_metric,
    'anchor_minus_strand_flips_sign'::VARCHAR AS oriented_distance_rule,
    'same_motif_id_nonoverlapping_distinct_span'::VARCHAR AS tandem_identity_rule,
    'same_alignment_span_collapses_orientation_records'::VARCHAR
        AS partner_locus_identity_rule,
    'singleton_same_opposite_mixed_or_ambiguous'::VARCHAR AS pair_class_rule,
    'at_least_one_member_is_a_tp73_context_locus'::VARCHAR
        AS cofactor_pair_scope,
    'tp73_context_loci_plus_their_pair_partners'::VARCHAR
        AS cofactor_motif_locus_scope,
    'tp73_context_loci_only'::VARCHAR AS cofactor_locus_pair_feature_scope,
    'physical_anchor_locus_neighbor_motif_interval_distance_band'::VARCHAR
        AS anchor_motif_band_feature_grain,
    'highest_score_then_nearest_center_then_coordinates'::VARCHAR
        AS anchor_motif_band_winner_rule,
    'both_locus_best_scores_at_or_above_neighbor_qualifying_threshold'::VARCHAR
        AS cofactor_pair_score_rule,
    {sql_string(motif_hit_source)}::VARCHAR AS motif_hit_source,
    {sql_string(arguments.annotation_release)}::VARCHAR AS annotation_release,
    {sql_string(arguments.promoter_definition_id)}::VARCHAR
        AS promoter_definition_id,
    {arguments.promoter_upstream_bp}::INTEGER AS promoter_upstream_bp,
    {arguments.promoter_downstream_bp}::INTEGER AS promoter_downstream_bp,
    {sql_string(arguments.downstream_definition_id)}::VARCHAR
        AS downstream_definition_id,
    {arguments.downstream_upstream_bp}::INTEGER AS downstream_upstream_bp,
    {arguments.downstream_downstream_bp}::INTEGER AS downstream_downstream_bp,
    'one_base_bed_start_plus_end_minus_1'::VARCHAR AS tss_coordinate_rule,
    'one_base_bed_end_minus_1_plus_start'::VARCHAR AS tes_coordinate_rule,
    'all_physical_tss_at_minimum_center_distance'::VARCHAR AS nearest_tss_rule,
    'physical_cds_segments_overlap_first_then_minimum_interval_gap'::VARCHAR
        AS nearest_cds_rule,
    'gtf_cds_1_based_inclusive_to_bed_half_open'::VARCHAR
        AS cds_coordinate_rule,
    'anchor_overlaps_resolved_promoter_interval'::VARCHAR
        AS promoter_membership_rule,
    'anchor_overlaps_resolved_downstream_interval'::VARCHAR
        AS downstream_membership_rule,
    'promoter_then_downstream_then_gene_body_then_intergenic'::VARCHAR
        AS gene_relation_precedence,
    'no_transcript_promoter_or_downstream_region_overlap'::VARCHAR
        AS strict_intergenic_rule,
    {gtf_source}::VARCHAR AS gtf_source,
    {gtf_sha256}::VARCHAR AS gtf_sha256,
    {gtf_size_bytes}::BIGINT AS gtf_size_bytes
FROM capture_parameter p;

COPY context_run_config TO 'tables/jaspar2026/context_run_config.parquet'
    (FORMAT PARQUET, COMPRESSION ZSTD);

CREATE TEMP VIEW anchor_hit_binned AS
SELECT
    a.*,
    FLOOR(a.center_bp / GREATEST(p.capture_prefilter_center_bp, 1))::BIGINT
        AS capture_bin
FROM anchor_hit a
CROSS JOIN capture_parameter p;

CREATE TEMP VIEW candidate_neighbor_hit AS
SELECT
    n.*,
    FLOOR(n.center_bp / GREATEST(p.capture_prefilter_center_bp, 1))::BIGINT
        + offsets.bin_offset AS anchor_capture_bin
FROM configured_hit n
CROSS JOIN capture_parameter p
CROSS JOIN (VALUES (-1), (0), (1)) AS offsets(bin_offset);

-- Keep the raw pair layer orientation-specific. Exact signed interval distance
-- decides capture and feature membership; center offsets retain direction.
CREATE TEMP TABLE motif_context_pair AS
WITH candidate_geometry AS (
    SELECT
        a.*,
        n.hit_id AS neighbor_hit_id,
        n.start AS neighbor_start,
        n."end" AS neighbor_end,
        n.span_bp AS neighbor_span_bp,
        n.center_bp AS neighbor_center_bp,
        n.motif_id AS neighbor_motif_id,
        n.motif_name AS neighbor_motif_name,
        n.strand AS neighbor_strand,
        n.score AS neighbor_score,
        n.pwm_relative_score AS neighbor_pwm_relative_score,
        n.qualifying_threshold AS neighbor_qualifying_threshold,
        (
            GREATEST(a.start, n.start) - LEAST(a."end", n."end")
        )::BIGINT AS anchor_neighbor_interval_distance_bp
    FROM anchor_hit_binned a
    JOIN candidate_neighbor_hit n
      ON a.genome_id = n.genome_id
     AND a.motif_set_id = n.motif_set_id
     AND a.chrom = n.chrom
     AND a.capture_bin = n.anchor_capture_bin
    CROSS JOIN capture_parameter p
    WHERE n.hit_id <> a.anchor_hit_id
      AND n.center_bp BETWEEN a.center_bp - p.capture_prefilter_center_bp
                          AND a.center_bp + p.capture_prefilter_center_bp
)
SELECT
    a.anchor_hit_id,
    a.anchor_locus_id,
    a.genome_id,
    a.motif_set_id,
    a.chrom,
    a.start AS anchor_start,
    a."end" AS anchor_end,
    a.center_bp AS anchor_center_bp,
    a.strand AS anchor_strand,
    a.score AS anchor_score,
    a.pwm_relative_score AS anchor_pwm_relative_score,
    a.neighbor_hit_id,
    a.neighbor_start,
    a.neighbor_end,
    a.neighbor_center_bp,
    a.neighbor_motif_id,
    a.neighbor_motif_name,
    a.neighbor_strand,
    a.neighbor_score,
    a.neighbor_pwm_relative_score,
    a.neighbor_qualifying_threshold,
    a.neighbor_center_bp - a.center_bp AS genomic_center_distance_bp,
    CASE WHEN a.strand = '-' THEN a.center_bp - a.neighbor_center_bp
         ELSE a.neighbor_center_bp - a.center_bp END
        AS anchor_oriented_center_distance_bp,
    ABS(a.neighbor_center_bp - a.center_bp) AS absolute_center_distance_bp,
    CASE WHEN a.neighbor_center_bp < a.center_bp THEN 'left'
         WHEN a.neighbor_center_bp > a.center_bp THEN 'right'
         ELSE 'coincident_center' END AS genomic_side,
    CASE WHEN a.strand = a.neighbor_strand THEN 'same' ELSE 'opposite' END
        AS relative_orientation,
    CASE WHEN (CASE WHEN a.strand = '-' THEN a.center_bp - a.neighbor_center_bp
                    ELSE a.neighbor_center_bp - a.center_bp END) < 0
            THEN 'upstream'
         WHEN (CASE WHEN a.strand = '-' THEN a.center_bp - a.neighbor_center_bp
                    ELSE a.neighbor_center_bp - a.center_bp END) > 0
            THEN 'downstream'
         ELSE 'coincident_center' END AS anchor_oriented_side,
    a.neighbor_start = a.start AND a.neighbor_end = a."end"
        AS same_alignment_span,
    -- Tandem is a conservative architectural label: distinct TP73 spans may
    -- touch, but shifted overlapping alignments are retained only as pairs.
    a.neighbor_motif_id = a.motif_id
        AND a.neighbor_start = a.start AND a.neighbor_end = a."end"
        AS same_anchor_motif_span,
    a.anchor_neighbor_interval_distance_bp,
    GREATEST(0, -a.anchor_neighbor_interval_distance_bp)::BIGINT
        AS interval_overlap_bp,
    GREATEST(0, a.anchor_neighbor_interval_distance_bp)::BIGINT
        AS inter_motif_gap_bp,
    CASE
        WHEN a.neighbor_start = a.start AND a.neighbor_end = a."end"
            THEN 'identical_span'
        WHEN a.start <= a.neighbor_start AND a."end" >= a.neighbor_end
            THEN 'anchor_contains_neighbor'
        WHEN a.neighbor_start <= a.start AND a.neighbor_end >= a."end"
            THEN 'neighbor_contains_anchor'
        WHEN a.anchor_neighbor_interval_distance_bp < 0 THEN 'partial_overlap'
        WHEN a.anchor_neighbor_interval_distance_bp = 0 THEN 'abutting'
        ELSE 'disjoint'
    END AS interval_relation,
    CASE
        WHEN a.anchor_neighbor_interval_distance_bp < 0 THEN 'overlap'
        WHEN a.anchor_neighbor_interval_distance_bp <= 5 THEN 'adjacent_0_5'
        WHEN a.anchor_neighbor_interval_distance_bp <= 20 THEN 'gap_6_20'
        WHEN a.anchor_neighbor_interval_distance_bp <= 50 THEN 'gap_21_50'
        WHEN a.anchor_neighbor_interval_distance_bp <= 100 THEN 'gap_51_100'
        WHEN a.anchor_neighbor_interval_distance_bp <= 150 THEN 'gap_101_150'
        ELSE 'outside_150'
    END AS interval_distance_band,
    a.anchor_neighbor_interval_distance_bp <= 5 AS within_5,
    a.anchor_neighbor_interval_distance_bp <= 20 AS within_20,
    a.anchor_neighbor_interval_distance_bp <= 50 AS within_50,
    a.anchor_neighbor_interval_distance_bp <= 100 AS within_100,
    a.anchor_neighbor_interval_distance_bp <= 150 AS within_150,
    a.anchor_neighbor_interval_distance_bp <= {arguments.context_flank}
        AS within_context_flank,
    a.neighbor_motif_id = a.motif_id
        AND NOT (a.neighbor_start = a.start AND a.neighbor_end = a."end")
        AND a.score >= {sql_number(arguments.partner_minimum_score)}
        AND a.neighbor_score >= {sql_number(arguments.partner_minimum_score)}
        AND a.anchor_neighbor_interval_distance_bp >= 0
        AND a.anchor_neighbor_interval_distance_bp <= {arguments.tandem_flank}
        AS is_tandem_tp73,
    a.score_mode,
    a.pseudocount,
    {sql_string(arguments.background_model_id)}::VARCHAR AS background_model_id,
    {sql_string(arguments.pseudocount_scheme)}::VARCHAR AS pseudocount_scheme,
    {sql_number(arguments.anchor_minimum_score)}::DOUBLE AS anchor_minimum_score,
    {sql_number(arguments.partner_minimum_score)}::DOUBLE AS partner_minimum_score,
    {arguments.capture_flank}::INTEGER AS capture_flank_bp,
    {arguments.context_flank}::INTEGER AS context_flank_bp,
    {arguments.tandem_flank}::INTEGER AS tandem_flank_bp,
    {arguments.cofactor_pair_flank}::INTEGER AS cofactor_pair_flank_bp
FROM candidate_geometry a
WHERE a.anchor_neighbor_interval_distance_bp <= {arguments.capture_flank};

SELECT CASE WHEN (SELECT COUNT(*) FROM motif_context_pair) = 0
    THEN error('anchor hits have no neighboring motif hits inside capture flank') END;

{motif_context_copy_sql}

DROP VIEW candidate_neighbor_hit;
DROP VIEW anchor_hit_binned;

-- CUT&RUN and other physical-location labels live at locus grain. Keep every
-- scored TP73 orientation above, but also provide one row per alignment span
-- with explicit orientation ambiguity instead of duplicating an ML label.
CREATE TEMP TABLE tp73_anchor_locus AS
WITH collapsed AS (
    SELECT
        md5(concat_ws('|', genome_id, motif_set_id, chrom, start::VARCHAR,
                      "end"::VARCHAR, motif_id, score_mode,
                      printf('%.17g', pseudocount))) AS anchor_locus_id,
        genome_id,
        motif_set_id,
        chrom,
        start,
        "end",
        span_bp,
        center_bp,
        motif_id,
        MAX(motif_name) AS motif_name,
        BOOL_OR(strand = '+') AS has_plus_orientation,
        BOOL_OR(strand = '-') AS has_minus_orientation,
        COUNT(*)::BIGINT AS n_orientation_records,
        MAX(score) FILTER (WHERE strand = '+') AS plus_score,
        MAX(score) FILTER (WHERE strand = '-') AS minus_score,
        MAX(pwm_relative_score) FILTER (WHERE strand = '+')
            AS plus_pwm_relative_score,
        MAX(pwm_relative_score) FILTER (WHERE strand = '-')
            AS minus_pwm_relative_score,
        MAX(score) AS best_score,
        MAX(anchor_locus_best_score) AS anchor_locus_best_score,
        MAX(best_other_anchor_locus_score) AS best_other_anchor_locus_score,
        MAX(anchor_locus_score_prominence) AS anchor_locus_score_prominence,
        BOOL_OR(anchor_locus_is_local_peak) AS anchor_locus_is_local_peak,
        MAX(anchor_selection_class) AS anchor_selection_class,
        MAX(score_mode) AS score_mode,
        MAX(pseudocount) AS pseudocount
    FROM anchor_hit
    GROUP BY genome_id, motif_set_id, chrom, start, "end", span_bp,
             center_bp, motif_id, score_mode, pseudocount
)
SELECT
    *,
    CASE WHEN has_plus_orientation AND has_minus_orientation THEN 'ambiguous'
         WHEN has_plus_orientation THEN 'plus'
         ELSE 'minus' END AS orientation_state,
    CASE
        WHEN plus_score IS NULL THEN 'minus'
        WHEN minus_score IS NULL THEN 'plus'
        WHEN plus_score > minus_score THEN 'plus'
        WHEN minus_score > plus_score THEN 'minus'
        ELSE 'ambiguous'
    END AS best_orientation_state,
    CASE
        WHEN plus_score IS NULL THEN minus_pwm_relative_score
        WHEN minus_score IS NULL THEN plus_pwm_relative_score
        WHEN plus_score >= minus_score THEN plus_pwm_relative_score
        ELSE minus_pwm_relative_score
    END AS best_pwm_relative_score,
    {sql_number(arguments.anchor_minimum_score)}::DOUBLE AS anchor_minimum_score,
    {arguments.anchor_local_peak_flank}::INTEGER AS anchor_local_peak_flank_bp,
    {sql_string(arguments.background_model_id)}::VARCHAR AS background_model_id,
    {sql_string(arguments.pseudocount_scheme)}::VARCHAR AS pseudocount_scheme
FROM collapsed;

{anchor_locus_copy_sql}

-- Collapse neighbor strand alternatives at the same span before counting
-- context loci. The orientation state remains relative to this TP73 record.
CREATE TEMP TABLE context_neighbor_locus AS
WITH orientation_records AS (
    SELECT
        anchor_hit_id,
        anchor_locus_id,
        genome_id,
        motif_set_id,
        chrom,
        neighbor_start,
        neighbor_end,
        neighbor_center_bp,
        neighbor_motif_id,
        MAX(neighbor_motif_name) AS neighbor_motif_name,
        BOOL_OR(relative_orientation = 'same') AS has_same_orientation,
        BOOL_OR(relative_orientation = 'opposite') AS has_opposite_orientation,
        COUNT(*)::BIGINT AS n_orientation_records,
        MAX(neighbor_score) FILTER (WHERE neighbor_strand = '+')
            AS neighbor_plus_score,
        MAX(neighbor_score) FILTER (WHERE neighbor_strand = '-')
            AS neighbor_minus_score,
        MAX(neighbor_pwm_relative_score) FILTER (WHERE neighbor_strand = '+')
            AS neighbor_plus_pwm_relative_score,
        MAX(neighbor_pwm_relative_score) FILTER (WHERE neighbor_strand = '-')
            AS neighbor_minus_pwm_relative_score,
        MAX(neighbor_score) AS neighbor_score,
        MAX(genomic_center_distance_bp) AS genomic_center_distance_bp,
        MAX(anchor_oriented_center_distance_bp)
            AS anchor_oriented_center_distance_bp,
        MAX(absolute_center_distance_bp) AS absolute_center_distance_bp,
        MAX(genomic_side) AS genomic_side,
        MAX(anchor_oriented_side) AS anchor_oriented_side,
        MAX(anchor_neighbor_interval_distance_bp)
            AS anchor_neighbor_interval_distance_bp,
        MAX(interval_overlap_bp)::BIGINT AS interval_overlap_bp,
        MAX(inter_motif_gap_bp)::BIGINT AS inter_motif_gap_bp,
        MAX(interval_relation) AS interval_relation,
        MAX(interval_distance_band) AS interval_distance_band,
        MAX(neighbor_qualifying_threshold) AS neighbor_qualifying_threshold,
        MAX(score_mode) AS score_mode,
        MAX(pseudocount) AS pseudocount,
        MAX(background_model_id) AS background_model_id,
        MAX(pseudocount_scheme) AS pseudocount_scheme
    FROM motif_context_pair
    WHERE within_context_flank AND NOT same_anchor_motif_span
    GROUP BY anchor_hit_id, anchor_locus_id, genome_id, motif_set_id, chrom, neighbor_start,
             neighbor_end, neighbor_center_bp, neighbor_motif_id
)
SELECT
    md5(concat_ws('|', genome_id, motif_set_id, chrom,
                  neighbor_start::VARCHAR, neighbor_end::VARCHAR,
                  neighbor_motif_id)) AS neighbor_locus_id,
    *,
    CASE WHEN has_same_orientation AND has_opposite_orientation THEN 'ambiguous'
         WHEN has_same_orientation THEN 'same'
         ELSE 'opposite' END AS relative_orientation_state,
    CASE
        WHEN neighbor_plus_score IS NULL THEN 'minus'
        WHEN neighbor_minus_score IS NULL THEN 'plus'
        WHEN neighbor_plus_score > neighbor_minus_score THEN 'plus'
        WHEN neighbor_minus_score > neighbor_plus_score THEN 'minus'
        ELSE 'ambiguous'
    END AS strongest_neighbor_orientation_state,
    CASE
        WHEN neighbor_plus_score IS NULL
            THEN neighbor_minus_pwm_relative_score
        WHEN neighbor_minus_score IS NULL
            THEN neighbor_plus_pwm_relative_score
        WHEN neighbor_plus_score >= neighbor_minus_score
            THEN neighbor_plus_pwm_relative_score
        ELSE neighbor_minus_pwm_relative_score
    END AS neighbor_pwm_relative_score
FROM orientation_records;

CREATE TEMP TABLE tp73_motif_context_summary AS
SELECT
    anchor_hit_id,
    genome_id,
    motif_set_id,
    chrom,
    neighbor_motif_id,
    MAX(neighbor_motif_name) AS neighbor_motif_name,
    interval_distance_band,
    anchor_oriented_side,
    relative_orientation_state,
    COUNT(*)::BIGINT AS n_neighbor_loci,
    SUM(n_orientation_records)::BIGINT AS n_orientation_records,
    MAX(neighbor_score) AS best_neighbor_score,
    MAX(neighbor_pwm_relative_score) AS best_neighbor_pwm_relative_score,
    MIN(GREATEST(anchor_neighbor_interval_distance_bp, 0))::BIGINT
        AS nearest_nonoverlap_distance_bp,
    MAX(interval_overlap_bp)::BIGINT AS maximum_overlap_bp,
    MIN(absolute_center_distance_bp) AS nearest_absolute_center_distance_bp,
    MAX(score_mode) AS score_mode,
    MAX(pseudocount) AS pseudocount,
    MAX(background_model_id) AS background_model_id,
    MAX(pseudocount_scheme) AS pseudocount_scheme,
    {arguments.context_flank}::INTEGER AS context_flank_bp
FROM context_neighbor_locus
GROUP BY anchor_hit_id, genome_id, motif_set_id, chrom, neighbor_motif_id,
         interval_distance_band, anchor_oriented_side,
         relative_orientation_state;

{motif_summary_copy_sql}

DROP TABLE tp73_motif_context_summary;

-- Generic cofactor pairing starts at physical motif-locus grain. Identical
-- spans reported on both strands become one ambiguous locus and can never
-- create a phantom distance-zero self-pair.
CREATE TEMP TABLE cofactor_motif_locus_all AS
WITH collapsed AS (
    SELECT
        md5(concat_ws('|', genome_id, motif_set_id, chrom, start::VARCHAR,
                      "end"::VARCHAR, motif_id)) AS cofactor_locus_id,
        genome_id,
        motif_set_id,
        chrom,
        start,
        "end",
        span_bp,
        center_bp,
        motif_id,
        MAX(motif_name) AS motif_name,
        BOOL_OR(strand = '+') AS has_plus_orientation,
        BOOL_OR(strand = '-') AS has_minus_orientation,
        COUNT(*)::BIGINT AS n_orientation_records,
        MAX(score) FILTER (WHERE strand = '+') AS plus_score,
        MAX(score) FILTER (WHERE strand = '-') AS minus_score,
        MAX(pwm_relative_score) FILTER (WHERE strand = '+')
            AS plus_pwm_relative_score,
        MAX(pwm_relative_score) FILTER (WHERE strand = '-')
            AS minus_pwm_relative_score,
        MAX(score) AS best_score,
        MAX(qualifying_threshold) AS qualifying_threshold,
        MAX(score_mode) AS score_mode,
        MAX(pseudocount) AS pseudocount
    FROM configured_hit
    WHERE motif_id <> {sql_string(arguments.anchor_motif)}
    GROUP BY genome_id, motif_set_id, chrom, start, "end", span_bp,
             center_bp, motif_id
)
SELECT
    *,
    CASE WHEN has_plus_orientation AND has_minus_orientation THEN 'ambiguous'
         WHEN has_plus_orientation THEN 'plus'
         ELSE 'minus' END AS orientation_state,
    CASE
        WHEN plus_score IS NULL THEN minus_pwm_relative_score
        WHEN minus_score IS NULL THEN plus_pwm_relative_score
        WHEN plus_score >= minus_score THEN plus_pwm_relative_score
        ELSE minus_pwm_relative_score
    END AS best_pwm_relative_score,
    {sql_string(arguments.background_model_id)}::VARCHAR AS background_model_id,
    {sql_string(arguments.pseudocount_scheme)}::VARCHAR AS pseudocount_scheme
FROM collapsed;

DROP TABLE configured_hit;

-- Only cofactor loci captured around a TP73 anchor can seed a cofactor pair in
-- this package. A same-motif partner may lie outside the TP73 context flank;
-- retaining it preserves the one-member-in-context boundary semantics.
CREATE TEMP TABLE context_cofactor_locus_id AS
SELECT DISTINCT neighbor_locus_id AS cofactor_locus_id
FROM context_neighbor_locus
WHERE neighbor_motif_id <> {sql_string(arguments.anchor_motif)};

CREATE TEMP VIEW all_cofactor_locus_binned AS
SELECT
    l.*,
    FLOOR(l.center_bp / GREATEST(p.cofactor_pair_prefilter_center_bp, 1))
        ::BIGINT AS pair_bin
FROM cofactor_motif_locus_all l
CROSS JOIN capture_parameter p;

CREATE TEMP VIEW context_cofactor_locus_binned AS
SELECT l.*
FROM all_cofactor_locus_binned l
JOIN context_cofactor_locus_id c USING (cofactor_locus_id);

CREATE TEMP VIEW candidate_right_all_cofactor_locus AS
SELECT
    r.*,
    FLOOR(r.center_bp / GREATEST(p.cofactor_pair_prefilter_center_bp, 1))
        ::BIGINT + offsets.bin_offset AS left_pair_bin
FROM cofactor_motif_locus_all r
CROSS JOIN capture_parameter p
CROSS JOIN (VALUES (-1), (0), (1)) AS offsets(bin_offset);

CREATE TEMP VIEW candidate_right_context_cofactor_locus AS
SELECT r.*
FROM candidate_right_all_cofactor_locus r
JOIN context_cofactor_locus_id c
  ON r.cofactor_locus_id = c.cofactor_locus_id;

CREATE TEMP TABLE cofactor_motif_pair AS
WITH candidate_pair AS (
    -- The genomic-left member is itself in TP73 context.
    SELECT
        l.cofactor_locus_id AS left_locus_id,
        r.cofactor_locus_id AS right_locus_id
    FROM context_cofactor_locus_binned l
    JOIN candidate_right_all_cofactor_locus r
      ON l.genome_id = r.genome_id
     AND l.motif_set_id = r.motif_set_id
     AND l.chrom = r.chrom
     AND l.motif_id = r.motif_id
     AND l.pair_bin = r.left_pair_bin
    CROSS JOIN capture_parameter p
    WHERE (
          l.start < r.start
       OR (l.start = r.start AND l."end" < r."end")
       OR (l.start = r.start AND l."end" = r."end"
           AND l.cofactor_locus_id < r.cofactor_locus_id)
    )
      AND r.center_bp BETWEEN l.center_bp - p.cofactor_pair_prefilter_center_bp
                          AND l.center_bp + p.cofactor_pair_prefilter_center_bp

    UNION ALL

    -- The genomic-right member is in context and the left member is not. This
    -- second branch retains boundary-crossing pairs without duplicating pairs
    -- for which both members are context loci.
    SELECT
        l.cofactor_locus_id AS left_locus_id,
        r.cofactor_locus_id AS right_locus_id
    FROM all_cofactor_locus_binned l
    JOIN candidate_right_context_cofactor_locus r
      ON l.genome_id = r.genome_id
     AND l.motif_set_id = r.motif_set_id
     AND l.chrom = r.chrom
     AND l.motif_id = r.motif_id
     AND l.pair_bin = r.left_pair_bin
    LEFT JOIN context_cofactor_locus_id left_context
      ON l.cofactor_locus_id = left_context.cofactor_locus_id
    CROSS JOIN capture_parameter p
    WHERE left_context.cofactor_locus_id IS NULL
      AND (
          l.start < r.start
       OR (l.start = r.start AND l."end" < r."end")
       OR (l.start = r.start AND l."end" = r."end"
           AND l.cofactor_locus_id < r.cofactor_locus_id)
      )
      AND r.center_bp BETWEEN l.center_bp - p.cofactor_pair_prefilter_center_bp
                          AND l.center_bp + p.cofactor_pair_prefilter_center_bp
), candidate_geometry AS (
    SELECT
        l.*,
        r.cofactor_locus_id AS right_locus_id,
        r.start AS right_start,
        r."end" AS right_end,
        r.span_bp AS right_span_bp,
        r.center_bp AS right_center_bp,
        r.orientation_state AS right_orientation_state,
        r.best_score AS right_score,
        r.best_pwm_relative_score AS right_pwm_relative_score,
        r.qualifying_threshold AS right_qualifying_threshold,
        (
            GREATEST(l.start, r.start) - LEAST(l."end", r."end")
        )::BIGINT AS pair_member_interval_distance_bp
    FROM candidate_pair pair
    JOIN cofactor_motif_locus_all l
      ON pair.left_locus_id = l.cofactor_locus_id
    JOIN cofactor_motif_locus_all r
      ON pair.right_locus_id = r.cofactor_locus_id
)
SELECT
    md5(concat_ws('|', genome_id, motif_set_id, chrom, motif_id,
                  cofactor_locus_id, right_locus_id)) AS cofactor_pair_id,
    genome_id,
    motif_set_id,
    chrom,
    motif_id,
    motif_name,
    cofactor_locus_id AS left_locus_id,
    start AS left_start,
    "end" AS left_end,
    center_bp AS left_center_bp,
    orientation_state AS left_orientation_state,
    best_score AS left_score,
    best_pwm_relative_score AS left_pwm_relative_score,
    right_locus_id,
    right_start,
    right_end,
    right_center_bp,
    right_orientation_state,
    right_score,
    right_pwm_relative_score,
    qualifying_threshold,
    pair_member_interval_distance_bp,
    GREATEST(0, -pair_member_interval_distance_bp)::BIGINT
        AS pair_member_overlap_bp,
    GREATEST(0, pair_member_interval_distance_bp)::BIGINT
        AS pair_member_gap_bp,
    right_center_bp - center_bp AS pair_member_center_distance_bp,
    CASE
        WHEN start <= right_start AND "end" >= right_end
            THEN 'left_contains_right'
        WHEN right_start <= start AND right_end >= "end"
            THEN 'right_contains_left'
        WHEN pair_member_interval_distance_bp < 0 THEN 'partial_overlap'
        WHEN pair_member_interval_distance_bp = 0 THEN 'abutting'
        ELSE 'disjoint'
    END AS pair_member_interval_relation,
    CASE
        WHEN pair_member_interval_distance_bp < 0 THEN 'overlap'
        WHEN pair_member_interval_distance_bp <= 5 THEN 'adjacent_0_5'
        WHEN pair_member_interval_distance_bp <= 20 THEN 'gap_6_20'
        WHEN pair_member_interval_distance_bp <= 50 THEN 'gap_21_50'
        WHEN pair_member_interval_distance_bp <= 100 THEN 'gap_51_100'
        WHEN pair_member_interval_distance_bp <= 150 THEN 'gap_101_150'
        ELSE 'outside_150'
    END AS pair_member_distance_band,
    pair_member_interval_distance_bp <= 5 AS within_5,
    pair_member_interval_distance_bp <= 20 AS within_20,
    pair_member_interval_distance_bp <= 50 AS within_50,
    pair_member_interval_distance_bp <= 100 AS within_100,
    pair_member_interval_distance_bp <= 150 AS within_150,
    CASE
        WHEN orientation_state = 'ambiguous'
          OR right_orientation_state = 'ambiguous' THEN 'ambiguous'
        WHEN orientation_state = 'plus'
         AND right_orientation_state = 'plus' THEN 'codirectional_plus'
        WHEN orientation_state = 'minus'
         AND right_orientation_state = 'minus' THEN 'codirectional_minus'
        WHEN orientation_state = 'plus'
         AND right_orientation_state = 'minus' THEN 'convergent'
        ELSE 'divergent'
    END AS pair_arrangement,
    LEAST(best_score, right_score) AS pair_min_score,
    best_score + right_score AS pair_sum_score,
    LEAST(best_pwm_relative_score, right_pwm_relative_score)
        AS pair_min_pwm_relative_score,
    score_mode,
    pseudocount,
    background_model_id,
    pseudocount_scheme,
    {arguments.cofactor_pair_flank}::INTEGER AS cofactor_pair_flank_bp
FROM candidate_geometry
WHERE pair_member_interval_distance_bp <= {arguments.cofactor_pair_flank}
  AND best_score >= qualifying_threshold
  AND right_score >= right_qualifying_threshold;

CREATE TEMP TABLE relevant_cofactor_locus_id AS
SELECT cofactor_locus_id FROM context_cofactor_locus_id
UNION
SELECT left_locus_id FROM cofactor_motif_pair
UNION
SELECT right_locus_id FROM cofactor_motif_pair;

CREATE TEMP TABLE cofactor_motif_locus AS
SELECT
    l.*,
    c.cofactor_locus_id IS NOT NULL AS in_any_tp73_context
FROM cofactor_motif_locus_all l
JOIN relevant_cofactor_locus_id relevant USING (cofactor_locus_id)
LEFT JOIN context_cofactor_locus_id c USING (cofactor_locus_id);

{cofactor_locus_copy_sql}
{cofactor_pair_copy_sql}

DROP VIEW candidate_right_context_cofactor_locus;
DROP VIEW candidate_right_all_cofactor_locus;
DROP VIEW context_cofactor_locus_binned;
DROP VIEW all_cofactor_locus_binned;
DROP TABLE relevant_cofactor_locus_id;
DROP TABLE cofactor_motif_locus_all;

CREATE TEMP TABLE cofactor_locus_pair_feature AS
WITH directed AS (
    SELECT
        p.left_locus_id AS cofactor_locus_id,
        p.right_locus_id AS partner_locus_id,
        p.right_score AS partner_score,
        p.right_pwm_relative_score AS partner_pwm_relative_score,
        p.* EXCLUDE (left_locus_id, right_locus_id, right_score,
                     right_pwm_relative_score)
    FROM cofactor_motif_pair p
    JOIN context_cofactor_locus_id c
      ON p.left_locus_id = c.cofactor_locus_id
    UNION ALL
    SELECT
        p.right_locus_id AS cofactor_locus_id,
        p.left_locus_id AS partner_locus_id,
        p.left_score AS partner_score,
        p.left_pwm_relative_score AS partner_pwm_relative_score,
        p.* EXCLUDE (left_locus_id, right_locus_id, left_score,
                     left_pwm_relative_score)
    FROM cofactor_motif_pair p
    JOIN context_cofactor_locus_id c
      ON p.right_locus_id = c.cofactor_locus_id
), ranked AS (
    SELECT
        *,
        ROW_NUMBER() OVER (
            PARTITION BY cofactor_locus_id
            ORDER BY GREATEST(pair_member_interval_distance_bp, 0),
                     ABS(pair_member_center_distance_bp), cofactor_pair_id
        ) AS nearest_rank
    FROM directed
), summary AS (
    SELECT
        cofactor_locus_id,
        COUNT(*)::BIGINT AS n_same_motif_partner_loci,
        COUNT(*) FILTER (WHERE pair_arrangement = 'codirectional_plus')::BIGINT
            AS n_codirectional_plus_pairs,
        COUNT(*) FILTER (WHERE pair_arrangement = 'codirectional_minus')::BIGINT
            AS n_codirectional_minus_pairs,
        COUNT(*) FILTER (WHERE pair_arrangement = 'convergent')::BIGINT
            AS n_convergent_pairs,
        COUNT(*) FILTER (WHERE pair_arrangement = 'divergent')::BIGINT
            AS n_divergent_pairs,
        COUNT(*) FILTER (WHERE pair_arrangement = 'ambiguous')::BIGINT
            AS n_ambiguous_pairs,
        MAX(pair_min_score) AS best_pair_min_score,
        MAX(pair_sum_score) AS best_pair_sum_score
    FROM directed
    GROUP BY cofactor_locus_id
)
SELECT
    l.*,
    COALESCE(s.n_same_motif_partner_loci, 0)::BIGINT
        AS n_same_motif_partner_loci,
    COALESCE(s.n_codirectional_plus_pairs, 0)::BIGINT
        AS n_codirectional_plus_pairs,
    COALESCE(s.n_codirectional_minus_pairs, 0)::BIGINT
        AS n_codirectional_minus_pairs,
    COALESCE(s.n_convergent_pairs, 0)::BIGINT AS n_convergent_pairs,
    COALESCE(s.n_divergent_pairs, 0)::BIGINT AS n_divergent_pairs,
    COALESCE(s.n_ambiguous_pairs, 0)::BIGINT AS n_ambiguous_pairs,
    n.partner_locus_id AS nearest_partner_locus_id,
    n.pair_member_interval_distance_bp AS nearest_pair_member_distance_bp,
    n.pair_arrangement AS nearest_pair_arrangement,
    n.partner_score AS nearest_partner_score,
    n.partner_pwm_relative_score AS nearest_partner_pwm_relative_score,
    s.best_pair_min_score,
    s.best_pair_sum_score,
    COALESCE(s.n_same_motif_partner_loci, 0) > 0 AS has_same_motif_partner
FROM cofactor_motif_locus l
JOIN context_cofactor_locus_id context_locus USING (cofactor_locus_id)
LEFT JOIN summary s USING (cofactor_locus_id)
LEFT JOIN ranked n
  ON l.cofactor_locus_id = n.cofactor_locus_id AND n.nearest_rank = 1;

{cofactor_locus_feature_copy_sql}

-- The ML-facing context feature stays at physical TP73-locus grain. Collapse
-- TP73 orientation alternatives before selecting a neighboring locus, then
-- keep score, geometry, orientation, and pair state from one ranked winner.
CREATE TEMP TABLE physical_context_neighbor_locus AS
WITH collapsed AS (
    SELECT
        anchor_locus_id,
        genome_id,
        motif_set_id,
        chrom,
        neighbor_locus_id,
        neighbor_start,
        neighbor_end,
        neighbor_center_bp,
        neighbor_motif_id,
        MAX(neighbor_motif_name) AS neighbor_motif_name,
        MAX(neighbor_plus_score) AS neighbor_plus_score,
        MAX(neighbor_minus_score) AS neighbor_minus_score,
        MAX(neighbor_score) AS neighbor_score,
        MAX(neighbor_pwm_relative_score) AS neighbor_pwm_relative_score,
        MAX(n_orientation_records)::BIGINT AS n_orientation_records,
        MAX(genomic_center_distance_bp) AS genomic_center_distance_bp,
        MAX(absolute_center_distance_bp) AS absolute_center_distance_bp,
        MAX(genomic_side) AS genomic_side,
        MAX(anchor_neighbor_interval_distance_bp)::BIGINT
            AS anchor_neighbor_interval_distance_bp,
        MAX(interval_overlap_bp)::BIGINT AS interval_overlap_bp,
        MAX(inter_motif_gap_bp)::BIGINT AS inter_motif_gap_bp,
        MAX(interval_relation) AS interval_relation,
        MAX(interval_distance_band) AS interval_distance_band,
        MAX(neighbor_qualifying_threshold)::DOUBLE AS qualifying_threshold,
        MAX(score_mode) AS score_mode,
        MAX(pseudocount)::DOUBLE AS pseudocount,
        MAX(background_model_id) AS background_model_id,
        MAX(pseudocount_scheme) AS pseudocount_scheme
    FROM context_neighbor_locus
    GROUP BY anchor_locus_id, genome_id, motif_set_id, chrom,
             neighbor_locus_id, neighbor_start, neighbor_end,
             neighbor_center_bp, neighbor_motif_id
), oriented AS (
    SELECT
        c.*,
        a.start AS anchor_start,
        a."end" AS anchor_end,
        a.center_bp AS anchor_center_bp,
        a.motif_id AS anchor_motif_id,
        a.best_score AS anchor_best_score,
        a.best_pwm_relative_score AS anchor_best_pwm_relative_score,
        a.orientation_state AS anchor_orientation_state,
        a.best_orientation_state AS anchor_best_orientation_state,
        CASE
            WHEN c.neighbor_plus_score IS NULL THEN 'minus'
            WHEN c.neighbor_minus_score IS NULL THEN 'plus'
            WHEN c.neighbor_plus_score > c.neighbor_minus_score THEN 'plus'
            WHEN c.neighbor_minus_score > c.neighbor_plus_score THEN 'minus'
            ELSE 'ambiguous'
        END AS strongest_neighbor_orientation_state
    FROM collapsed c
    JOIN tp73_anchor_locus a USING (
        anchor_locus_id, genome_id, motif_set_id, chrom
    )
)
SELECT
    *,
    CASE
        WHEN anchor_best_orientation_state = 'plus'
            THEN genomic_center_distance_bp
        WHEN anchor_best_orientation_state = 'minus'
            THEN -genomic_center_distance_bp
        ELSE NULL::DOUBLE
    END AS anchor_oriented_center_distance_bp,
    CASE
        WHEN anchor_best_orientation_state = 'ambiguous' THEN NULL::VARCHAR
        WHEN anchor_best_orientation_state = 'plus'
         AND genomic_center_distance_bp < 0
            THEN 'upstream'
        WHEN anchor_best_orientation_state = 'plus'
         AND genomic_center_distance_bp > 0
            THEN 'downstream'
        WHEN anchor_best_orientation_state = 'minus'
         AND genomic_center_distance_bp > 0
            THEN 'upstream'
        WHEN anchor_best_orientation_state = 'minus'
         AND genomic_center_distance_bp < 0
            THEN 'downstream'
        ELSE 'coincident_center'
    END AS anchor_oriented_side,
    CASE
        WHEN anchor_best_orientation_state = 'ambiguous'
          OR strongest_neighbor_orientation_state = 'ambiguous'
            THEN 'ambiguous'
        WHEN anchor_best_orientation_state = strongest_neighbor_orientation_state
            THEN 'same'
        ELSE 'opposite'
    END AS strongest_relative_orientation_state
FROM oriented;

CREATE TEMP TABLE anchor_motif_band_feature AS
WITH qualifying AS (
    SELECT *
    FROM physical_context_neighbor_locus
    WHERE neighbor_score >= qualifying_threshold{band_qualifying_filter_sql}
), ranked AS (
    SELECT
        *,
        ROW_NUMBER() OVER (
            PARTITION BY anchor_locus_id, genome_id, motif_set_id,
                         neighbor_motif_id, interval_distance_band
            ORDER BY neighbor_score DESC,
                     absolute_center_distance_bp,
                     neighbor_start,
                     neighbor_end,
                     neighbor_locus_id
        ) AS best_rank
    FROM qualifying
), group_stat AS (
    SELECT
        anchor_locus_id,
        genome_id,
        motif_set_id,
        neighbor_motif_id,
        interval_distance_band,
        COUNT(*)::BIGINT AS n_neighbor_loci_above_threshold,
        COUNT(*) FILTER (
            WHERE neighbor_score = group_best_score
        )::BIGINT AS n_best_score_ties
    FROM (
        SELECT
            *,
            MAX(neighbor_score) OVER (
                PARTITION BY anchor_locus_id, genome_id, motif_set_id,
                             neighbor_motif_id, interval_distance_band
            ) AS group_best_score
        FROM qualifying
    ) counted
    GROUP BY anchor_locus_id, genome_id, motif_set_id,
             neighbor_motif_id, interval_distance_band
)
SELECT
    9::INTEGER AS schema_version,
    r.anchor_locus_id,
    r.genome_id,
    r.motif_set_id,
    r.chrom,
    r.anchor_start,
    r.anchor_end,
    r.anchor_center_bp,
    r.anchor_motif_id,
    r.anchor_best_score,
    r.anchor_best_pwm_relative_score,
    r.anchor_orientation_state,
    r.anchor_best_orientation_state,
    r.neighbor_motif_id,
    r.neighbor_motif_name,
    r.interval_distance_band,
    CASE r.interval_distance_band
        WHEN 'overlap' THEN 0
        WHEN 'adjacent_0_5' THEN 1
        WHEN 'gap_6_20' THEN 2
        WHEN 'gap_21_50' THEN 3
        WHEN 'gap_51_100' THEN 4
        WHEN 'gap_101_150' THEN 5
        ELSE 6
    END::INTEGER AS interval_distance_band_order,
    r.qualifying_threshold,
    TRUE AS threshold_inclusive,
    s.n_neighbor_loci_above_threshold,
    s.n_best_score_ties,
    r.neighbor_locus_id AS best_neighbor_locus_id,
    r.neighbor_start AS best_neighbor_start,
    r.neighbor_end AS best_neighbor_end,
    r.neighbor_score AS best_neighbor_score,
    r.neighbor_pwm_relative_score AS best_neighbor_pwm_relative_score,
    r.n_orientation_records AS best_neighbor_n_orientation_records,
    r.neighbor_plus_score AS best_neighbor_plus_score,
    r.neighbor_minus_score AS best_neighbor_minus_score,
    r.anchor_neighbor_interval_distance_bp AS best_interval_distance_bp,
    r.interval_overlap_bp AS best_interval_overlap_bp,
    r.inter_motif_gap_bp AS best_inter_motif_gap_bp,
    r.interval_relation AS best_interval_relation,
    r.genomic_center_distance_bp AS best_genomic_center_distance_bp,
    r.anchor_oriented_center_distance_bp
        AS best_anchor_oriented_center_distance_bp,
    r.genomic_side AS best_genomic_side,
    r.anchor_oriented_side AS best_anchor_oriented_side,
    r.strongest_neighbor_orientation_state
        AS best_neighbor_orientation_state,
    r.strongest_relative_orientation_state
        AS best_relative_orientation_state,
    p.cofactor_locus_id IS NOT NULL AS best_hit_pair_architecture_assessed,
    p.has_same_motif_partner AS best_hit_has_same_motif_partner,
    p.n_same_motif_partner_loci AS best_hit_n_same_motif_partner_loci,
    p.n_codirectional_plus_pairs AS best_hit_n_codirectional_plus_pairs,
    p.n_codirectional_minus_pairs AS best_hit_n_codirectional_minus_pairs,
    p.n_convergent_pairs AS best_hit_n_convergent_pairs,
    p.n_divergent_pairs AS best_hit_n_divergent_pairs,
    p.n_ambiguous_pairs AS best_hit_n_ambiguous_pairs,
    p.nearest_partner_locus_id AS best_hit_nearest_partner_locus_id,
    p.nearest_pair_member_distance_bp
        AS best_hit_nearest_pair_member_distance_bp,
    p.nearest_pair_arrangement AS best_hit_nearest_pair_arrangement,
    p.nearest_partner_score AS best_hit_nearest_partner_score,
    p.best_pair_min_score AS best_hit_best_pair_min_score,
    p.best_pair_sum_score AS best_hit_best_pair_sum_score,
    r.score_mode,
    r.pseudocount,
    r.background_model_id,
    r.pseudocount_scheme,
    {arguments.context_flank}::INTEGER AS context_flank_bp
FROM ranked r
JOIN group_stat s USING (
    anchor_locus_id, genome_id, motif_set_id,
    neighbor_motif_id, interval_distance_band
)
LEFT JOIN cofactor_locus_pair_feature p
  ON r.neighbor_locus_id = p.cofactor_locus_id
WHERE r.best_rank = 1;

{band_feature_copy_sql}

DROP TABLE anchor_motif_band_feature;
DROP TABLE physical_context_neighbor_locus;
DROP TABLE tp73_anchor_locus;
DROP TABLE cofactor_locus_pair_feature;
DROP TABLE cofactor_motif_locus;
DROP TABLE context_cofactor_locus_id;

-- Attribute each canonical cofactor pair once per TP73 anchor. Joining through
-- the collapsed context loci prevents strand alternatives and two in-context
-- members from multiplying the pair.
CREATE TEMP VIEW cofactor_pair_member AS
SELECT cofactor_pair_id, motif_id, left_locus_id AS cofactor_locus_id,
       'left'::VARCHAR AS pair_member
FROM cofactor_motif_pair
UNION ALL
SELECT cofactor_pair_id, motif_id, right_locus_id AS cofactor_locus_id,
       'right'::VARCHAR AS pair_member
FROM cofactor_motif_pair;

CREATE TEMP TABLE tp73_cofactor_pair_context AS
WITH member_context AS (
    SELECT
        c.anchor_hit_id,
        c.genome_id,
        c.motif_set_id,
        c.chrom,
        m.cofactor_pair_id,
        m.pair_member,
        c.neighbor_locus_id,
        c.anchor_neighbor_interval_distance_bp,
        c.interval_distance_band,
        c.absolute_center_distance_bp,
        c.anchor_oriented_side,
        c.relative_orientation_state
    FROM context_neighbor_locus c
    JOIN cofactor_pair_member m
      ON c.neighbor_motif_id = m.motif_id
     AND c.neighbor_locus_id = m.cofactor_locus_id
), ranked AS (
    SELECT
        *,
        ROW_NUMBER() OVER (
            PARTITION BY anchor_hit_id, cofactor_pair_id
            ORDER BY GREATEST(anchor_neighbor_interval_distance_bp, 0),
                     absolute_center_distance_bp, pair_member
        ) AS nearest_rank
    FROM member_context
), attribution AS (
    SELECT
        anchor_hit_id,
        cofactor_pair_id,
        MAX(genome_id) AS genome_id,
        MAX(motif_set_id) AS motif_set_id,
        MAX(chrom) AS chrom,
        COUNT(DISTINCT pair_member)::INTEGER AS n_pair_members_in_context,
        BOOL_OR(pair_member = 'left') AS left_member_in_context,
        BOOL_OR(pair_member = 'right') AS right_member_in_context,
        MAX(anchor_neighbor_interval_distance_bp)
            FILTER (WHERE pair_member = 'left') AS left_anchor_interval_distance_bp,
        MAX(anchor_neighbor_interval_distance_bp)
            FILTER (WHERE pair_member = 'right') AS right_anchor_interval_distance_bp
    FROM member_context
    GROUP BY anchor_hit_id, cofactor_pair_id
)
SELECT
    a.anchor_hit_id,
    a.genome_id,
    a.motif_set_id,
    a.chrom,
    p.cofactor_pair_id,
    p.motif_id AS cofactor_motif_id,
    p.motif_name AS cofactor_motif_name,
    p.left_locus_id,
    p.right_locus_id,
    p.pair_member_interval_distance_bp,
    p.pair_member_overlap_bp,
    p.pair_member_gap_bp,
    p.pair_member_center_distance_bp,
    p.pair_member_interval_relation,
    p.pair_member_distance_band,
    p.pair_arrangement,
    p.left_score,
    p.right_score,
    p.pair_min_score,
    p.pair_sum_score,
    p.pair_min_pwm_relative_score,
    p.score_mode,
    p.pseudocount,
    p.background_model_id,
    p.pseudocount_scheme,
    a.n_pair_members_in_context,
    a.left_member_in_context,
    a.right_member_in_context,
    a.n_pair_members_in_context = 2 AS pair_fully_within_context,
    a.left_anchor_interval_distance_bp,
    a.right_anchor_interval_distance_bp,
    n.neighbor_locus_id AS nearest_member_locus_id,
    n.pair_member AS nearest_pair_member,
    n.anchor_neighbor_interval_distance_bp
        AS nearest_member_anchor_neighbor_interval_distance_bp,
    n.interval_distance_band AS nearest_member_anchor_distance_band,
    n.anchor_oriented_side AS nearest_member_anchor_oriented_side,
    n.relative_orientation_state AS nearest_member_relative_orientation,
    {arguments.context_flank}::INTEGER AS context_flank_bp,
    {arguments.cofactor_pair_flank}::INTEGER AS cofactor_pair_flank_bp
FROM attribution a
JOIN cofactor_motif_pair p USING (cofactor_pair_id)
JOIN ranked n
  ON a.anchor_hit_id = n.anchor_hit_id
 AND a.cofactor_pair_id = n.cofactor_pair_id
 AND n.nearest_rank = 1;

{cofactor_pair_context_copy_sql}

CREATE TEMP TABLE tp73_cofactor_pair_summary AS
SELECT
    anchor_hit_id,
    genome_id,
    motif_set_id,
    chrom,
    cofactor_motif_id,
    MAX(cofactor_motif_name) AS cofactor_motif_name,
    pair_member_distance_band,
    pair_arrangement,
    nearest_member_anchor_distance_band,
    nearest_member_anchor_oriented_side,
    nearest_member_relative_orientation,
    COUNT(*)::BIGINT AS n_cofactor_pairs,
    COUNT(*) FILTER (WHERE n_pair_members_in_context = 1)::BIGINT
        AS n_pairs_with_one_member_in_context,
    COUNT(*) FILTER (WHERE n_pair_members_in_context = 2)::BIGINT
        AS n_pairs_with_both_members_in_context,
    MAX(pair_min_score) AS best_pair_min_score,
    MAX(pair_sum_score) AS best_pair_sum_score,
    MAX(score_mode) AS score_mode,
    MAX(pseudocount) AS pseudocount,
    MAX(background_model_id) AS background_model_id,
    MAX(pseudocount_scheme) AS pseudocount_scheme,
    MIN(GREATEST(nearest_member_anchor_neighbor_interval_distance_bp, 0))
        ::BIGINT AS nearest_pair_nonoverlap_distance_bp,
    MAX(GREATEST(-nearest_member_anchor_neighbor_interval_distance_bp, 0))
        ::BIGINT AS maximum_nearest_member_overlap_bp
FROM tp73_cofactor_pair_context
GROUP BY anchor_hit_id, genome_id, motif_set_id, chrom, cofactor_motif_id,
         pair_member_distance_band, pair_arrangement,
         nearest_member_anchor_distance_band,
         nearest_member_anchor_oriented_side,
         nearest_member_relative_orientation;

{cofactor_pair_summary_copy_sql}

DROP TABLE tp73_cofactor_pair_context;
DROP TABLE tp73_cofactor_pair_summary;
DROP VIEW cofactor_pair_member;
DROP TABLE cofactor_motif_pair;

-- A near-palindromic motif can produce two strand records for one neighboring
-- alignment span. Preserve both records above, but collapse them here so the
-- feature layer counts physical partner loci rather than orientation choices.
CREATE TEMP TABLE tandem_partner_locus AS
WITH orientation_records AS (
    SELECT
        anchor_hit_id,
        chrom,
        neighbor_start,
        neighbor_end,
        neighbor_center_bp,
        genomic_center_distance_bp,
        anchor_oriented_center_distance_bp,
        absolute_center_distance_bp,
        inter_motif_gap_bp,
        BOOL_OR(relative_orientation = 'same') AS has_same_orientation,
        BOOL_OR(relative_orientation = 'opposite') AS has_opposite_orientation,
        COUNT(*)::BIGINT AS n_orientation_records,
        MAX(neighbor_score) AS partner_score,
        MAX(neighbor_pwm_relative_score) AS partner_pwm_relative_score
    FROM motif_context_pair
    WHERE is_tandem_tp73
    GROUP BY anchor_hit_id, chrom, neighbor_start, neighbor_end,
             neighbor_center_bp, genomic_center_distance_bp,
             anchor_oriented_center_distance_bp,
             absolute_center_distance_bp, inter_motif_gap_bp
)
SELECT
    *,
    CASE WHEN has_same_orientation AND has_opposite_orientation
             THEN 'ambiguous'
         WHEN has_same_orientation THEN 'same'
         ELSE 'opposite' END AS partner_orientation
FROM orientation_records;

CREATE TEMP TABLE tandem_partner_summary AS
SELECT
    anchor_hit_id,
    COUNT(*)::BIGINT AS n_tandem_tp73_partner_loci,
    COUNT(*) FILTER (WHERE partner_orientation = 'same')::BIGINT
        AS n_same_orientation_partner_loci,
    COUNT(*) FILTER (WHERE partner_orientation = 'opposite')::BIGINT
        AS n_opposite_orientation_partner_loci,
    COUNT(*) FILTER (WHERE partner_orientation = 'ambiguous')::BIGINT
        AS n_ambiguous_orientation_partner_loci,
    MIN(inter_motif_gap_bp)::BIGINT AS nearest_tandem_inter_motif_gap_bp,
    MIN(absolute_center_distance_bp) AS nearest_tandem_absolute_center_distance_bp,
    MIN(inter_motif_gap_bp) FILTER (WHERE partner_orientation = 'same')::BIGINT
        AS nearest_same_orientation_gap_bp,
    MIN(inter_motif_gap_bp) FILTER (WHERE partner_orientation = 'opposite')::BIGINT
        AS nearest_opposite_orientation_gap_bp,
    MIN(inter_motif_gap_bp) FILTER (WHERE partner_orientation = 'ambiguous')::BIGINT
        AS nearest_ambiguous_orientation_gap_bp,
    MAX(partner_score) AS best_partner_score,
    MAX(partner_pwm_relative_score) AS best_partner_pwm_relative_score,
    MAX(partner_score) FILTER (WHERE partner_orientation = 'same')
        AS best_same_orientation_partner_score,
    MAX(partner_score) FILTER (WHERE partner_orientation = 'opposite')
        AS best_opposite_orientation_partner_score,
    MAX(partner_score) FILTER (WHERE partner_orientation = 'ambiguous')
        AS best_ambiguous_orientation_partner_score
FROM tandem_partner_locus
GROUP BY anchor_hit_id;

-- These are sequence-compatible architecture features, not observations of a
-- protein complex or its quaternary structure.
CREATE TEMP TABLE tp73_pair_feature AS
SELECT
    a.anchor_hit_id,
    a.genome_id,
    a.motif_set_id,
    a.chrom,
    a.start,
    a."end",
    a.center_bp,
    a.motif_id,
    a.motif_name,
    a.strand,
    a.score,
    a.anchor_locus_id,
    a.anchor_selection_class,
    a.anchor_locus_best_score,
    a.best_other_anchor_locus_score,
    a.anchor_locus_score_prominence,
    a.anchor_locus_is_local_peak,
    a.score_mode,
    a.pseudocount,
    {sql_string(arguments.background_model_id)}::VARCHAR AS background_model_id,
    {sql_string(arguments.pseudocount_scheme)}::VARCHAR AS pseudocount_scheme,
    {sql_number(arguments.anchor_minimum_score)}::DOUBLE AS anchor_minimum_score,
    {sql_number(arguments.partner_minimum_score)}::DOUBLE AS partner_minimum_score,
    {arguments.anchor_local_peak_flank}::INTEGER AS anchor_local_peak_flank_bp,
    a.pwm_relative_score,
    CASE
        WHEN COALESCE(s.n_tandem_tp73_partner_loci, 0) = 0
            THEN 'singleton'
        WHEN COALESCE(s.n_ambiguous_orientation_partner_loci, 0) > 0
            THEN 'tandem_orientation_ambiguous'
        WHEN COALESCE(s.n_same_orientation_partner_loci, 0) > 0
         AND COALESCE(s.n_opposite_orientation_partner_loci, 0) > 0
            THEN 'tandem_mixed_orientation'
        WHEN COALESCE(s.n_same_orientation_partner_loci, 0) > 0
            THEN 'tandem_same_orientation'
        ELSE 'tandem_opposite_orientation'
    END AS pair_class,
    COALESCE(s.n_tandem_tp73_partner_loci, 0)::BIGINT
        AS n_tandem_tp73_partner_loci,
    COALESCE(s.n_same_orientation_partner_loci, 0)::BIGINT
        AS n_same_orientation_partner_loci,
    COALESCE(s.n_opposite_orientation_partner_loci, 0)::BIGINT
        AS n_opposite_orientation_partner_loci,
    COALESCE(s.n_ambiguous_orientation_partner_loci, 0)::BIGINT
        AS n_ambiguous_orientation_partner_loci,
    COALESCE(s.n_tandem_tp73_partner_loci, 0) > 1
        AS has_multiple_tandem_partner_loci,
    s.nearest_tandem_inter_motif_gap_bp,
    s.nearest_tandem_absolute_center_distance_bp,
    s.nearest_same_orientation_gap_bp,
    s.nearest_opposite_orientation_gap_bp,
    s.nearest_ambiguous_orientation_gap_bp,
    s.best_partner_score,
    s.best_partner_pwm_relative_score,
    s.best_same_orientation_partner_score,
    s.best_opposite_orientation_partner_score,
    s.best_ambiguous_orientation_partner_score,
    CASE WHEN s.best_partner_score IS NULL THEN NULL
         ELSE LEAST(a.score, s.best_partner_score) END AS best_pair_min_score,
    a.score + s.best_partner_score AS best_pair_sum_score,
    CASE WHEN s.best_partner_pwm_relative_score IS NULL THEN NULL
         ELSE LEAST(a.pwm_relative_score, s.best_partner_pwm_relative_score)
         END AS best_pair_min_pwm_relative_score,
    {arguments.tandem_flank}::INTEGER AS tandem_flank_bp
FROM anchor_hit a
LEFT JOIN tandem_partner_summary s USING (anchor_hit_id);

{pair_feature_copy_sql}

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
    a.genome_id,
    a.motif_set_id,
    a.chrom,
    a.start,
    a."end",
    a.center_bp,
    a.motif_id,
    a.motif_name,
    a.strand,
    a.score,
    a.anchor_locus_id,
    a.anchor_selection_class,
    a.anchor_locus_best_score,
    a.best_other_anchor_locus_score,
    a.anchor_locus_score_prominence,
    a.anchor_locus_is_local_peak,
    a.score_mode,
    a.pseudocount,
    {sql_string(arguments.background_model_id)}::VARCHAR AS background_model_id,
    {sql_string(arguments.pseudocount_scheme)}::VARCHAR AS pseudocount_scheme,
    {sql_number(arguments.anchor_minimum_score)}::DOUBLE AS anchor_minimum_score,
    {sql_number(arguments.partner_minimum_score)}::DOUBLE AS partner_minimum_score,
    {arguments.anchor_local_peak_flank}::INTEGER AS anchor_local_peak_flank_bp,
    a.pwm_relative_score,
    COALESCE(s.n_context_neighbor_hits, 0)::BIGINT AS n_context_neighbor_hits,
    COALESCE(s.n_context_neighbor_loci, 0)::BIGINT AS n_context_neighbor_loci,
    COALESCE(s.n_context_motifs, 0)::BIGINT AS n_context_motifs,
    COALESCE(s.n_tandem_tp73_partners, 0) > 0 AS has_tandem_tp73,
    COALESCE(s.n_tandem_tp73_partners, 0)::BIGINT AS n_tandem_tp73_partners,
    f.pair_class,
    f.n_tandem_tp73_partner_loci,
    f.n_same_orientation_partner_loci,
    f.n_opposite_orientation_partner_loci,
    f.n_ambiguous_orientation_partner_loci,
    f.has_multiple_tandem_partner_loci,
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
    g.nearest_tss_id,
    g.nearest_tss_start,
    g.nearest_tss_strand,
    g.nearest_tss_tie_count,
    g.nearest_tss_has_mixed_strands,
    g.nearest_tss_genomic_distance_bp,
    g.nearest_tss_distance_bp,
    g.nearest_tss_interval_distance_bp,
    g.nearest_tss_relation,
    g.nearest_cds_segment_id,
    g.nearest_cds_start,
    g.nearest_cds_end,
    g.nearest_cds_strand,
    g.nearest_cds_tie_count,
    g.nearest_cds_has_mixed_strands,
    g.nearest_cds_genomic_distance_bp,
    g.nearest_cds_interval_distance_bp,
    g.nearest_cds_overlap_bp,
    g.nearest_cds_distance_bp,
    g.nearest_cds_relation,
    g.nearest_cds_gene_id,
    g.nearest_cds_gene_name,
    g.nearest_cds_transcript_id,
    g.in_any_transcript,
    g.fully_within_any_transcript,
    g.in_any_exon,
    g.overlaps_any_exon_boundary,
    g.in_any_cds,
    g.overlaps_any_cds_boundary,
    g.in_any_intron,
    g.overlaps_any_intron_boundary,
    g.overlaps_any_promoter,
    g.n_overlapping_promoters,
    g.overlaps_any_downstream_region,
    g.n_overlapping_downstream_regions,
    g.strict_intergenic,
    g.primary_transcript_region,
    g.primary_genomic_context,
    g.gene_relation_class,
    {arguments.capture_flank}::INTEGER AS capture_flank_bp,
    {arguments.context_flank}::INTEGER AS context_flank_bp,
    {arguments.tandem_flank}::INTEGER AS tandem_flank_bp
FROM anchor_hit a
LEFT JOIN pair_summary s USING (anchor_hit_id)
LEFT JOIN nearest_tandem t USING (anchor_hit_id)
JOIN tp73_pair_feature f USING (anchor_hit_id)
JOIN anchor_gene_context g USING (anchor_hit_id);

{context_anchor_copy_sql}

DROP TABLE context_run_config;
DROP TABLE motif_context_pair;
DROP TABLE tp73_pair_feature;
DROP TABLE tp73_context_anchor;
CREATE VIEW motif_context_run_config AS
SELECT * FROM read_parquet('tables/jaspar2026/context_run_config.parquet');
-- Retain the original short name for packages and queries written before the
-- canonical schema-contract name was introduced.
CREATE VIEW context_run_config AS
SELECT * FROM motif_context_run_config;
CREATE VIEW motif_context_pair AS
SELECT * REPLACE (
    CAST(genome_id AS VARCHAR) AS genome_id,
    CAST(chrom AS VARCHAR) AS chrom
)
FROM read_parquet(
    'tables/jaspar2026/motif_context_pair/**/*.parquet',
    hive_partitioning=1
);
CREATE VIEW tp73_anchor_locus AS
SELECT * REPLACE (
    CAST(genome_id AS VARCHAR) AS genome_id,
    CAST(chrom AS VARCHAR) AS chrom
)
FROM read_parquet(
    'tables/jaspar2026/tp73_anchor_locus/**/*.parquet',
    hive_partitioning=1
);
CREATE VIEW tp73_motif_context_summary AS
SELECT * REPLACE (
    CAST(genome_id AS VARCHAR) AS genome_id,
    CAST(chrom AS VARCHAR) AS chrom,
    CAST(neighbor_motif_id AS VARCHAR) AS neighbor_motif_id
)
FROM read_parquet(
    'tables/jaspar2026/tp73_motif_context_summary/**/*.parquet',
    hive_partitioning=1
);
CREATE VIEW anchor_motif_band_feature AS
SELECT * REPLACE (
    CAST(genome_id AS VARCHAR) AS genome_id,
    CAST(chrom AS VARCHAR) AS chrom,
    CAST(neighbor_motif_id AS VARCHAR) AS neighbor_motif_id
)
FROM read_parquet(
    'tables/jaspar2026/anchor_motif_band_feature/**/*.parquet',
    hive_partitioning=1
);
CREATE VIEW cofactor_motif_locus AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/cofactor_motif_locus.parquet'
);
CREATE VIEW cofactor_motif_pair AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/cofactor_motif_pair.parquet'
);
CREATE VIEW cofactor_locus_pair_feature AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/cofactor_locus_pair_feature.parquet'
);
CREATE VIEW tp73_cofactor_pair_context AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/tp73_cofactor_pair_context.parquet'
);
CREATE VIEW tp73_cofactor_pair_summary AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/tp73_cofactor_pair_summary.parquet'
);
CREATE VIEW tp73_pair_feature AS
SELECT * REPLACE (
    CAST(genome_id AS VARCHAR) AS genome_id,
    CAST(chrom AS VARCHAR) AS chrom
)
FROM read_parquet(
    'tables/jaspar2026/tp73_pair_feature/**/*.parquet',
    hive_partitioning=1
);
CREATE VIEW tp73_context_pair_feature AS
SELECT
    p.*,
    f.motif_id AS anchor_motif_id,
    f.pair_class AS anchor_pair_class,
    f.n_tandem_tp73_partner_loci,
    f.n_same_orientation_partner_loci,
    f.n_opposite_orientation_partner_loci,
    f.n_ambiguous_orientation_partner_loci,
    CASE
        WHEN p.is_tandem_tp73 THEN 'same_motif_tandem'
        WHEN p.neighbor_motif_id = f.motif_id AND p.interval_overlap_bp > 0
            THEN 'same_motif_overlapping_alignment'
        WHEN p.neighbor_motif_id = f.motif_id THEN 'same_motif_context'
        ELSE 'heterotypic_context'
    END AS pair_relation,
    FLOOR(p.anchor_oriented_center_distance_bp / 5.0) * 5.0
        AS oriented_distance_bin_start_bp
FROM motif_context_pair p
JOIN tp73_pair_feature f USING (anchor_hit_id);
CREATE VIEW legacy_tp73_context_100 AS
SELECT
    p.*,
    p.neighbor_start - p.anchor_start AS legacy_genomic_start_distance_bp,
    CASE WHEN p.anchor_strand = '-'
         THEN p.anchor_start - p.neighbor_start
         ELSE p.neighbor_start - p.anchor_start END
        AS legacy_anchor_oriented_start_distance_bp
FROM motif_context_pair p
WHERE ABS(p.neighbor_start - p.anchor_start) <= 100;
CREATE VIEW tp73_context_anchor AS
SELECT * REPLACE (
    CAST(genome_id AS VARCHAR) AS genome_id,
    CAST(chrom AS VARCHAR) AS chrom
)
FROM read_parquet(
    'tables/jaspar2026/tp73_context_anchor/**/*.parquet',
    hive_partitioning=1
);
""" + (
        """
DROP TABLE anchor_transcript_context;
DROP TABLE transcription_start_site;
DROP TABLE transcript_tss;
DROP TABLE transcription_end_site;
DROP TABLE transcript_tes;
DROP TABLE promoter_gene;
DROP TABLE downstream_region_gene;
DROP TABLE tp73_anchor_nearest_tss;
DROP TABLE tp73_anchor_nearest_cds;
DROP TABLE tp73_anchor_promoter;
DROP TABLE tp73_anchor_downstream_region;
CREATE VIEW transcription_start_site AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/transcription_start_site.parquet'
);
CREATE VIEW transcript_tss AS
SELECT * FROM read_parquet('tables/jaspar2026/transcript_tss.parquet');
CREATE VIEW transcription_end_site AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/transcription_end_site.parquet'
);
CREATE VIEW transcript_tes AS
SELECT * FROM read_parquet('tables/jaspar2026/transcript_tes.parquet');
CREATE VIEW gene AS
SELECT * FROM read_parquet('tables/jaspar2026/gene.parquet');
CREATE VIEW exon AS
SELECT * FROM read_parquet('tables/jaspar2026/exon.parquet');
CREATE VIEW coding_sequence_segment AS
SELECT * FROM read_parquet('tables/jaspar2026/coding_sequence_segment.parquet');
CREATE VIEW transcript_cds AS
SELECT * FROM read_parquet('tables/jaspar2026/transcript_cds.parquet');
CREATE VIEW transcript_coding_span AS
SELECT * FROM read_parquet('tables/jaspar2026/transcript_coding_span.parquet');
CREATE VIEW coding_landmark AS
SELECT * FROM read_parquet('tables/jaspar2026/coding_landmark.parquet');
CREATE VIEW promoter AS
SELECT * FROM read_parquet('tables/jaspar2026/promoter.parquet');
CREATE VIEW promoter_gene AS
SELECT * FROM read_parquet('tables/jaspar2026/promoter_gene.parquet');
CREATE VIEW downstream_region AS
SELECT * FROM read_parquet('tables/jaspar2026/downstream_region.parquet');
CREATE VIEW downstream_region_gene AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/downstream_region_gene.parquet'
);
CREATE VIEW tp73_anchor_nearest_tss AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/tp73_anchor_nearest_tss.parquet'
);
CREATE VIEW tp73_anchor_nearest_cds AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/tp73_anchor_nearest_cds.parquet'
);
CREATE VIEW tp73_anchor_promoter AS
SELECT * FROM read_parquet('tables/jaspar2026/tp73_anchor_promoter.parquet');
CREATE VIEW tp73_anchor_downstream_region AS
SELECT * FROM read_parquet(
    'tables/jaspar2026/tp73_anchor_downstream_region.parquet'
);
CREATE VIEW transcript AS
SELECT * FROM read_parquet('tables/jaspar2026/transcript.parquet');
CREATE VIEW intron AS
SELECT * FROM read_parquet('tables/jaspar2026/intron.parquet');
CREATE VIEW motif_transcript_context AS
SELECT * FROM read_parquet('tables/jaspar2026/motif_transcript_context.parquet');
CREATE VIEW motif_transcript_context_pair AS
SELECT
    p.*,
    t.gene_id,
    t.gene_name,
    t.transcript_id,
    t.transcript_strand,
    t.tss,
    t.signed_tss_distance_bp AS anchor_signed_tss_distance_bp,
    CASE WHEN t.transcript_strand = '+' THEN p.neighbor_center_bp - t.tss
         ELSE t.tss - p.neighbor_center_bp END
        AS neighbor_signed_tss_distance_bp,
    CASE WHEN t.transcript_strand = '+' THEN p.genomic_center_distance_bp
         ELSE -p.genomic_center_distance_bp END
        AS transcript_oriented_center_distance_bp,
    CASE
        WHEN (CASE WHEN t.transcript_strand = '+'
                        THEN p.genomic_center_distance_bp
                   ELSE -p.genomic_center_distance_bp END) < 0 THEN 'upstream'
        WHEN (CASE WHEN t.transcript_strand = '+'
                        THEN p.genomic_center_distance_bp
                   ELSE -p.genomic_center_distance_bp END) > 0 THEN 'downstream'
        ELSE 'coincident_center'
    END AS transcript_oriented_side,
    p.anchor_start <= t.tss AND p.anchor_end > t.tss AS anchor_spans_tss,
    p.neighbor_start <= t.tss AND p.neighbor_end > t.tss AS neighbor_spans_tss
FROM motif_context_pair p
JOIN motif_transcript_context t USING (anchor_hit_id);
""" if gtf is not None else ""
    )


def run_duckdb(executable: str, staging: Path, sql: str) -> None:
    process = subprocess.run(
        [executable, "-bail", "context.duckdb"],
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
    identifier_pattern = re.compile(r"^[A-Za-z0-9._-]+$")
    for option, value in (
        ("--motif-set-id", arguments.motif_set_id),
        ("--genome-id", arguments.genome_id),
        ("--background-model-id", arguments.background_model_id),
        ("--pseudocount-scheme", arguments.pseudocount_scheme),
        ("--source-commit", arguments.source_commit),
        ("--annotation-release", arguments.annotation_release),
        ("--promoter-definition-id", arguments.promoter_definition_id),
        ("--downstream-definition-id", arguments.downstream_definition_id),
    ):
        if not identifier_pattern.fullmatch(value):
            raise ContextBuildError(
                f"{option} must contain only letters, digits, '.', '_', or '-'"
            )

    motif_glob = resolve_parquet_glob(arguments.motif_hits)
    columns = parquet_columns(arguments.duckdb, motif_glob)
    required_columns = {
        "chrom", "start", "end", "motif_id", "strand", "score",
        "score_mode", "pseudocount", "pwm_relative_score",
    }
    missing_columns = sorted(required_columns - columns)
    if missing_columns:
        raise ContextBuildError(
            "motif-hit Parquet lacks required columns: "
            + ", ".join(missing_columns)
        )
    chromosomes = parse_chromosomes(arguments.chrom)
    gtf = arguments.gtf.expanduser().resolve() if arguments.gtf else None
    if gtf is not None and not gtf.is_file():
        raise ContextBuildError(f"GTF not found: {gtf}")
    if gtf is None:
        if arguments.gtf_sha256 is not None or arguments.gtf_size_bytes is not None:
            raise ContextBuildError(
                "--gtf-sha256/--gtf-size-bytes require --gtf"
            )
    else:
        actual_size = gtf.stat().st_size
        if (arguments.gtf_size_bytes is not None
                and arguments.gtf_size_bytes != actual_size):
            raise ContextBuildError(
                "GTF size differs from --gtf-size-bytes: "
                f"expected {arguments.gtf_size_bytes}, observed {actual_size}"
            )
        expected_sha256 = arguments.gtf_sha256
        if expected_sha256 is not None and not re.fullmatch(
            r"[0-9A-Fa-f]{64}", expected_sha256
        ):
            raise ContextBuildError("--gtf-sha256 must contain 64 hexadecimal digits")
        actual_sha256 = sha256_file(gtf)
        if (expected_sha256 is not None
                and expected_sha256.lower() != actual_sha256):
            raise ContextBuildError(
                "GTF content differs from --gtf-sha256: "
                f"expected {expected_sha256.lower()}, observed {actual_sha256}"
            )
        arguments.gtf_sha256 = actual_sha256
        arguments.gtf_size_bytes = actual_size

    output = arguments.output.expanduser().resolve()
    if output.exists() and not arguments.force:
        raise ContextBuildError(f"output already exists (use --force): {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    report_disk_preflight(output.parent, motif_glob)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        (staging / "tables" / "jaspar2026").mkdir(parents=True)
        if arguments.output_tier != "selected":
            (staging / "tables" / "jaspar2026" / "motif_context_pair").mkdir()
        if arguments.output_tier == "band":
            for table in (
                "tp73_anchor_locus",
                "tp73_motif_context_summary",
                "anchor_motif_band_feature",
                "tp73_pair_feature",
                "tp73_context_anchor",
            ):
                (staging / "tables" / "jaspar2026" / table).mkdir()

        external_temp_directory = arguments.temp_directory is not None
        if external_temp_directory:
            temp_directory = arguments.temp_directory.expanduser().resolve()
            if temp_directory.exists() and not temp_directory.is_dir():
                raise ContextBuildError(
                    f"--temp-directory is not a directory: {temp_directory}"
                )
            temp_directory.mkdir(parents=True, exist_ok=True)
            if any(temp_directory.iterdir()):
                raise ContextBuildError(
                    f"--temp-directory must be empty: {temp_directory}"
                )
            print(
                f"I: DuckDB spill files will use {temp_directory}",
                file=sys.stderr,
            )
        else:
            temp_directory = staging / "duckdb_tmp"

        sql = build_sql(
            arguments, motif_glob, chromosomes, gtf, temp_directory,
            source_columns=columns,
        )
        run_duckdb(arguments.duckdb, staging, sql)
        if not external_temp_directory:
            shutil.rmtree(temp_directory, ignore_errors=True)

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

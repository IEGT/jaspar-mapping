#!/usr/bin/env python3

"""Calibrate dense motif scores against merged CUT&RUN coverage.

Fragment BED rows contribute unit depth; bedGraph rows preserve the signal in
column 4 when maximum depth is calculated inside a strictly immersed motif.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from pathlib import Path


NEGATIVE_INFINITY_BIN = -9223372036854775807


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a memory-bounded chromosome-wide ROC/PR calibration from "
            "dense motif-score Parquet and positive CUT&RUN BED coverage."
        )
    )
    parser.add_argument("--package", required=True, type=Path)
    parser.add_argument("--coverage-bed", required=True, type=Path)
    parser.add_argument(
        "--coverage-format",
        choices=("auto", "fragments", "bedgraph"),
        default="auto",
        help=(
            "depth interpretation (default: auto; .bedGraph[.gz] uses column "
            "4, other BED files contribute one per row)"
        ),
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--motif-id", default="MA0861.2")
    parser.add_argument("--motif-length", type=int, default=16)
    parser.add_argument("--chrom", default="1")
    parser.add_argument("--pseudocount", type=float, default=1.0)
    parser.add_argument(
        "--score-mode",
        action="append",
        choices=("log2_relative_risk", "log_odds"),
        dest="score_modes",
        help="Score mode to analyze; repeat for both (default: both)",
    )
    parser.add_argument("--bin-width", type=float, default=0.1)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--memory-limit", default="6GB")
    parser.add_argument("--max-temp-size", default="40GB")
    parser.add_argument("--duckdb", default="duckdb")
    arguments = parser.parse_args()
    if arguments.score_modes is None:
        arguments.score_modes = ["log2_relative_risk", "log_odds"]
    arguments.score_modes = list(dict.fromkeys(arguments.score_modes))
    return arguments


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def resolve_coverage_format(requested: str, path: Path) -> str:
    if requested != "auto":
        return requested
    filename = path.name.lower()
    if filename.endswith(".gz"):
        filename = filename[:-3]
    return "bedgraph" if filename.endswith((".bedgraph", ".bdg")) else "fragments"


def coverage_depth_semantics(coverage_format: str) -> str:
    return (
        "column_4_signal"
        if coverage_format == "bedgraph"
        else "unit_per_fragment_interval"
    )


def coverage_input_sql(arguments: argparse.Namespace) -> str:
    source = sql_string(arguments.coverage_bed.resolve())
    normalized_chrom = (
        "CASE WHEN LOWER(LEFT(chrom, 3)) = 'chr' "
        "THEN SUBSTR(chrom, 4) ELSE chrom END"
    )
    if arguments.coverage_format == "bedgraph":
        return f"""
SELECT
    {normalized_chrom} AS chrom,
    start::BIGINT AS start,
    "end"::BIGINT AS "end",
    depth::DOUBLE AS depth,
    'run_' || ROW_NUMBER() OVER () AS source_name
FROM read_csv(
    {source}, delim='\\t', header=false, comment='#', null_padding=true,
    columns={{
        'chrom':'VARCHAR', 'start':'BIGINT', 'end':'BIGINT', 'depth':'DOUBLE'
    }}
)"""
    return f"""
SELECT
    {normalized_chrom} AS chrom,
    start::BIGINT AS start,
    "end"::BIGINT AS "end",
    1.0::DOUBLE AS depth,
    COALESCE(NULLIF(name, ''), 'interval_' || ROW_NUMBER() OVER ()) AS source_name
FROM read_csv(
    {source}, delim='\\t', header=false, comment='#', null_padding=true,
    columns={{
        'chrom':'VARCHAR', 'start':'BIGINT', 'end':'BIGINT',
        'name':'VARCHAR', 'score':'VARCHAR', 'strand':'VARCHAR'
    }}
)"""


def build_sql(arguments: argparse.Namespace, output_dir: Path) -> str:
    dense_glob = (
        arguments.package
        / "tables/jaspar2026/motif_score_dense/*/*/*/*/*/*.parquet"
    )
    database = output_dir / "calibration.duckdb"
    temp_dir = output_dir / "duckdb_tmp"
    threshold_parquet = output_dir / "threshold_curve.parquet"
    threshold_tsv = output_dir / "threshold_curve.tsv"
    immersed_parquet = output_dir / "immersed_motif_evidence.parquet"
    immersed_tsv = output_dir / "immersed_motif_evidence.tsv"
    component_tsv = output_dir / "coverage_component_evidence.tsv"
    summary_tsv = output_dir / "summary.tsv"
    score_mode_sql = ", ".join(sql_string(mode) for mode in arguments.score_modes)
    coverage_sql = coverage_input_sql(arguments)
    depth_semantics = coverage_depth_semantics(arguments.coverage_format)
    histogram_statements = []
    for index, score_mode in enumerate(arguments.score_modes):
        create_or_insert = (
            "CREATE TABLE all_score_histogram AS"
            if index == 0
            else "INSERT INTO all_score_histogram"
        )
        histogram_statements.append(
            f"""
{create_or_insert}
WITH block_histograms AS (
    SELECT b.score_mode,
           list_histogram(list_transform(
               list_zip(b.plus_scores, b.minus_scores),
               lambda pair: CASE
                   WHEN isinf(GREATEST(pair[1], pair[2]))
                        AND GREATEST(pair[1], pair[2]) < 0
                   THEN {NEGATIVE_INFINITY_BIN}::BIGINT
                   ELSE FLOOR(
                       GREATEST(pair[1], pair[2]) / {arguments.bin_width}
                   )::BIGINT
               END
           )) AS bin_counts
    FROM paired_block b
    WHERE b.score_mode = {sql_string(score_mode)}
)
SELECT score_mode,
       entry.key::BIGINT AS bin_index,
       SUM(entry.value)::BIGINT AS n_total
FROM block_histograms
CROSS JOIN UNNEST(map_entries(bin_counts)) AS u(entry)
GROUP BY score_mode, entry.key;
"""
        )
    histogram_sql = "\n".join(histogram_statements)

    return f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;
SET temp_directory={sql_string(temp_dir)};
SET max_temp_directory_size={sql_string(arguments.max_temp_size)};
PRAGMA enable_progress_bar;

CREATE TABLE run_config AS SELECT
    {sql_string(arguments.sample_id)}::VARCHAR AS sample_id,
    {sql_string(arguments.motif_id)}::VARCHAR AS motif_id,
    {arguments.motif_length}::INTEGER AS motif_length,
    {sql_string(arguments.chrom)}::VARCHAR AS chrom,
    {arguments.pseudocount}::DOUBLE AS pseudocount,
    {sql_string(",".join(arguments.score_modes))}::VARCHAR AS score_modes,
    {arguments.bin_width}::DOUBLE AS score_bin_width,
    {sql_string(arguments.coverage_bed.resolve())}::VARCHAR AS coverage_source,
    {sql_string(arguments.coverage_format)}::VARCHAR AS coverage_format,
    {sql_string(depth_semantics)}::VARCHAR AS coverage_depth_semantics,
    {sql_string(dense_glob.resolve())}::VARCHAR AS dense_source_glob,
    'overlap_or_adjacent_union'::VARCHAR AS coverage_merge_rule,
    'component_start < motif_start AND component_end > motif_end'::VARCHAR
        AS immersion_rule,
    'max_orientation_score'::VARCHAR AS orientation_handling;

CREATE TEMP VIEW coverage_input AS
{coverage_sql};

CREATE TABLE coverage_interval AS
SELECT
    ROW_NUMBER() OVER ()::BIGINT AS interval_id,
    chrom,
    start,
    "end",
    depth,
    source_name
FROM coverage_input
WHERE chrom = {sql_string(arguments.chrom)};

SELECT CASE
    WHEN (SELECT COUNT(*) FROM coverage_interval) = 0
    THEN error('coverage BED has no intervals on the selected chromosome')
    WHEN EXISTS (
        SELECT 1 FROM coverage_interval WHERE start < 0 OR "end" <= start
    ) THEN error('coverage BED contains an invalid interval')
    WHEN EXISTS (
        SELECT 1 FROM coverage_interval WHERE NOT isfinite(depth) OR depth <= 0
    ) THEN error('coverage depth must be finite and positive')
END;

CREATE TABLE coverage_component AS
WITH ordered AS (
    SELECT *,
        MAX("end") OVER (
            PARTITION BY chrom ORDER BY start, "end", interval_id
            ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING
        ) AS prior_max_end
    FROM coverage_interval
), marked AS (
    SELECT *, CASE WHEN prior_max_end IS NULL OR start > prior_max_end
                   THEN 1 ELSE 0 END AS starts_component
    FROM ordered
), grouped AS (
    SELECT *,
        SUM(starts_component) OVER (
            PARTITION BY chrom ORDER BY start, "end", interval_id
        ) AS component_group
    FROM marked
), aggregated AS (
    SELECT chrom, component_group,
           MIN(start)::BIGINT AS start,
           MAX("end")::BIGINT AS "end",
           COUNT(*)::BIGINT AS source_interval_count,
           STRING_AGG(source_name, ',' ORDER BY interval_id) AS source_intervals
    FROM grouped
    GROUP BY chrom, component_group
)
SELECT ROW_NUMBER() OVER (ORDER BY chrom, start, "end")::BIGINT AS component_id,
       chrom, start, "end", source_interval_count, source_intervals
FROM aggregated;

CREATE VIEW dense_block AS
SELECT
    score_mode,
    CAST(pseudocount AS DOUBLE) AS pseudocount,
    CAST(chrom AS VARCHAR) AS chrom,
    strand,
    block_start,
    scores
FROM read_parquet(
    {sql_string(dense_glob.resolve())},
    hive_partitioning=true,
    union_by_name=true
)
WHERE motif_id = {sql_string(arguments.motif_id)}
  AND CAST(pseudocount AS DOUBLE) = {arguments.pseudocount}
  AND CAST(chrom AS VARCHAR) = {sql_string(arguments.chrom)}
  AND score_mode IN ({score_mode_sql});

CREATE VIEW paired_block AS
SELECT p.score_mode, p.chrom, p.block_start,
       p.scores AS plus_scores, m.scores AS minus_scores
FROM dense_block p
JOIN dense_block m
  ON p.score_mode = m.score_mode
 AND p.chrom = m.chrom
 AND p.block_start = m.block_start
WHERE p.strand = 'plus' AND m.strand = 'minus';

SELECT CASE
    WHEN (SELECT COUNT(DISTINCT score_mode) FROM paired_block) <> {len(arguments.score_modes)}
    THEN error('dense input does not contain every requested score mode and orientation')
END;

-- Process modes sequentially, and collapse each stored block to a histogram
-- before creating rows. This keeps the configured DuckDB memory limit useful on
-- small machines instead of expanding chromosome-wide score lists at once.
{histogram_sql}

CREATE TABLE immersed_start AS
SELECT c.component_id, c.chrom,
       UNNEST(RANGE(c.start + 1, c."end" - {arguments.motif_length}))::BIGINT AS start
FROM coverage_component c
WHERE c."end" - c.start > {arguments.motif_length} + 1;

CREATE TABLE immersed_motif_raw AS
SELECT b.score_mode, s.component_id, s.chrom, s.start,
       s.start + {arguments.motif_length} AS "end",
       GREATEST(
           LIST_EXTRACT(b.plus_scores, s.start - b.block_start + 1),
           LIST_EXTRACT(b.minus_scores, s.start - b.block_start + 1)
       ) AS score
FROM immersed_start s
JOIN paired_block b
  ON s.chrom = b.chrom
 AND s.start >= b.block_start
 AND s.start < b.block_start + LEN(b.plus_scores);

CREATE TABLE immersed_motif AS
WITH base_depth AS (
    SELECT m.score_mode, m.component_id, m.chrom, m.start, m."end", m.score,
           base,
           COALESCE(SUM(i.depth), 0.0)::DOUBLE AS depth
    FROM immersed_motif_raw m
    CROSS JOIN UNNEST(RANGE(m.start, m."end")) AS b(base)
    LEFT JOIN coverage_interval i
      ON i.chrom = m.chrom AND i.start <= base AND i."end" > base
    WHERE m.score IS NOT NULL
    GROUP BY m.score_mode, m.component_id, m.chrom, m.start, m."end", m.score, base
)
SELECT score_mode, component_id, chrom, start, "end", score,
       CASE WHEN isinf(score) AND score < 0
            THEN {NEGATIVE_INFINITY_BIN}::BIGINT
            ELSE FLOOR(score / {arguments.bin_width})::BIGINT
       END AS bin_index,
       MAX(depth)::DOUBLE AS effective_max_depth
FROM base_depth
GROUP BY score_mode, component_id, chrom, start, "end", score;

CREATE TABLE immersed_score_histogram AS
SELECT score_mode, bin_index,
       COUNT(*)::BIGINT AS n_supported,
       SUM(effective_max_depth)::DOUBLE AS effective_depth_sum
FROM immersed_motif
GROUP BY score_mode, bin_index;

CREATE TABLE component_best_score AS
SELECT score_mode, component_id, MAX(score) AS best_score
FROM immersed_motif
GROUP BY score_mode, component_id;

CREATE TABLE threshold_curve AS
WITH histogram AS (
    SELECT a.score_mode, a.bin_index, a.n_total,
           COALESCE(s.n_supported, 0)::BIGINT AS n_supported,
           COALESCE(s.effective_depth_sum, 0.0)::DOUBLE AS effective_depth_sum
    FROM all_score_histogram a
    LEFT JOIN immersed_score_histogram s USING (score_mode, bin_index)
), cumulative AS (
    SELECT *,
        CASE WHEN bin_index = {NEGATIVE_INFINITY_BIN}
             THEN '-Infinity'::DOUBLE
             ELSE bin_index * {arguments.bin_width}
        END AS threshold,
        SUM(n_total) OVER (
            PARTITION BY score_mode ORDER BY bin_index DESC
        )::BIGINT AS selected_motifs,
        SUM(n_supported) OVER (
            PARTITION BY score_mode ORDER BY bin_index DESC
        )::BIGINT AS supported_selected_motifs,
        SUM(effective_depth_sum) OVER (
            PARTITION BY score_mode ORDER BY bin_index DESC
        )::DOUBLE AS selected_effective_depth
    FROM histogram
), totals AS (
    SELECT score_mode, SUM(n_total)::BIGINT AS total_motifs,
           SUM(n_supported)::BIGINT AS total_supported
    FROM histogram
    GROUP BY score_mode
), rates AS (
    SELECT c.score_mode,
           c.bin_index,
           c.threshold,
           c.selected_motifs,
           c.supported_selected_motifs,
           c.selected_motifs - c.supported_selected_motifs
               AS unsupported_selected_motifs,
           c.supported_selected_motifs::DOUBLE / c.selected_motifs
               AS motif_precision,
           c.supported_selected_motifs::DOUBLE / t.total_supported
               AS motif_recall,
           (c.selected_motifs - c.supported_selected_motifs)::DOUBLE
               / (t.total_motifs - t.total_supported)
               AS motif_false_positive_rate,
           (SELECT COUNT(*)::DOUBLE FROM component_best_score b
            WHERE b.score_mode = c.score_mode
              AND b.best_score >= c.threshold)
               / (SELECT COUNT(*) FROM coverage_component)
               AS coverage_component_recall,
           c.selected_effective_depth::DOUBLE / c.selected_motifs
               AS mean_effective_depth,
           t.total_supported::DOUBLE / t.total_motifs AS support_prevalence
    FROM cumulative c
    JOIN totals t USING (score_mode)
)
SELECT *,
       motif_precision / support_prevalence AS precision_enrichment,
       2.0 * motif_precision * motif_recall
           / NULLIF(motif_precision + motif_recall, 0.0) AS f1,
       motif_recall - motif_false_positive_rate AS youden_j
FROM rates
ORDER BY score_mode, threshold DESC;

CREATE TABLE calibration_summary AS
WITH histogram AS (
    SELECT a.score_mode, a.bin_index, a.n_total,
           COALESCE(s.n_supported, 0)::BIGINT AS n_supported,
           a.n_total - COALESCE(s.n_supported, 0) AS n_unsupported
    FROM all_score_histogram a
    LEFT JOIN immersed_score_histogram s USING (score_mode, bin_index)
), ranked AS (
    SELECT *, COALESCE(SUM(n_unsupported) OVER (
        PARTITION BY score_mode ORDER BY bin_index
        ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING
    ), 0) AS lower_unsupported
    FROM histogram
), auc AS (
    SELECT score_mode,
           SUM(n_supported * (lower_unsupported + 0.5 * n_unsupported))
             / (SUM(n_supported) * SUM(n_unsupported)) AS roc_auc_binned
    FROM ranked GROUP BY score_mode
), ap AS (
    SELECT c.score_mode,
           SUM(h.n_supported::DOUBLE / t.total_supported * c.motif_precision)
               AS average_precision
    FROM threshold_curve c
    JOIN histogram h
      ON h.score_mode = c.score_mode
     AND h.bin_index = c.bin_index
    JOIN (
        SELECT score_mode, SUM(n_supported)::DOUBLE AS total_supported
        FROM histogram GROUP BY score_mode
    ) t ON t.score_mode = c.score_mode
    GROUP BY c.score_mode
), best_youden AS (
    SELECT score_mode, threshold AS youden_threshold, youden_j
    FROM threshold_curve
    QUALIFY ROW_NUMBER() OVER (
        PARTITION BY score_mode ORDER BY youden_j DESC, threshold DESC
    ) = 1
), best_f1 AS (
    SELECT score_mode, threshold AS f1_threshold, f1
    FROM threshold_curve
    QUALIFY ROW_NUMBER() OVER (
        PARTITION BY score_mode ORDER BY f1 DESC, threshold DESC
    ) = 1
), totals AS (
    SELECT score_mode, SUM(n_total)::BIGINT AS n_valid_motifs,
           SUM(n_supported)::BIGINT AS n_supported_motifs,
           COUNT(*)::BIGINT AS n_score_bins
    FROM histogram GROUP BY score_mode
)
SELECT t.score_mode,
       {sql_string(arguments.coverage_format)}::VARCHAR AS coverage_format,
       {sql_string(depth_semantics)}::VARCHAR AS coverage_depth_semantics,
       {arguments.bin_width}::DOUBLE AS score_bin_width,
       t.n_valid_motifs, t.n_supported_motifs,
       t.n_supported_motifs::DOUBLE / t.n_valid_motifs AS support_prevalence,
       (SELECT COUNT(*) FROM coverage_interval) AS n_coverage_intervals,
       (SELECT COUNT(*) FROM coverage_component) AS n_coverage_components,
       t.n_score_bins, a.roc_auc_binned, ap.average_precision,
       y.youden_threshold, y.youden_j, f.f1_threshold, f.f1
FROM totals t
JOIN auc a USING (score_mode)
JOIN ap USING (score_mode)
JOIN best_youden y USING (score_mode)
JOIN best_f1 f USING (score_mode)
ORDER BY t.score_mode;

COPY threshold_curve TO {sql_string(threshold_parquet)}
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY threshold_curve TO {sql_string(threshold_tsv)}
    (FORMAT CSV, DELIMITER '\t', HEADER true);
COPY immersed_motif TO {sql_string(immersed_parquet)}
    (FORMAT PARQUET, COMPRESSION ZSTD);
COPY immersed_motif TO {sql_string(immersed_tsv)}
    (FORMAT CSV, DELIMITER '\t', HEADER true);
COPY coverage_component TO {sql_string(component_tsv)}
    (FORMAT CSV, DELIMITER '\t', HEADER true);
COPY calibration_summary TO {sql_string(summary_tsv)}
    (FORMAT CSV, DELIMITER '\t', HEADER true);
"""


def main() -> int:
    arguments = parse_arguments()
    arguments.coverage_format = resolve_coverage_format(
        arguments.coverage_format, arguments.coverage_bed
    )
    if arguments.bin_width <= 0:
        raise SystemExit("E: --bin-width must be greater than zero")
    if arguments.threads < 1:
        raise SystemExit("E: --threads must be at least one")
    if not arguments.package.is_dir():
        raise SystemExit(f"E: package directory not found: {arguments.package}")
    if not arguments.coverage_bed.is_file():
        raise SystemExit(f"E: coverage BED not found: {arguments.coverage_bed}")
    duckdb = shutil.which(arguments.duckdb)
    if duckdb is None:
        raise SystemExit(f"E: DuckDB CLI not found: {arguments.duckdb}")

    output_dir = arguments.output_dir.resolve()
    expected = [
        output_dir / "analysis.sql",
        output_dir / "calibration.duckdb",
        output_dir / "threshold_curve.parquet",
        output_dir / "threshold_curve.tsv",
        output_dir / "immersed_motif_evidence.parquet",
        output_dir / "immersed_motif_evidence.tsv",
        output_dir / "coverage_component_evidence.tsv",
        output_dir / "summary.tsv",
        output_dir / "summary.json",
    ]
    existing = [path for path in expected if path.exists()]
    if existing:
        raise SystemExit(
            "E: refusing to replace existing output(s): "
            + ", ".join(str(path) for path in existing)
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "duckdb_tmp").mkdir(exist_ok=True)

    analysis_sql = output_dir / "analysis.sql"
    analysis_sql.write_text(build_sql(arguments, output_dir), encoding="utf-8")
    database = output_dir / "calibration.duckdb"
    subprocess.run(
        [duckdb, str(database), "-bail", "-f", str(analysis_sql)],
        check=True,
    )
    summary_output = subprocess.check_output(
        [
            duckdb,
            "-json",
            "-readonly",
            str(database),
            "-c",
            "SELECT * FROM calibration_summary ORDER BY score_mode;",
        ],
        text=True,
    )
    summary = json.loads(summary_output)
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

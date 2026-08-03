#!/usr/bin/env python3

"""Select versioned motif-score thresholds from a calibration metric table."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


class ThresholdBuildError(RuntimeError):
    pass


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_number(value: float) -> str:
    return format(value, ".17g")


def nullable_metric_cast(column: str, sql_type: str) -> str:
    """Cast an optional evaluator field while preserving malformed-value errors."""
    return f"CAST(NULLIF(TRIM(CAST({column} AS VARCHAR)), '') AS {sql_type})"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_jaspar_names(path: Path) -> dict[str, str]:
    names: dict[str, str] = {}
    header = re.compile(r"^>(\S+)\s+(.+?)\s*$")
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            match = header.match(line)
            if match:
                names[match.group(1)] = match.group(2)
    if not names:
        raise ThresholdBuildError(f"no motif headers found in JASPAR file: {path}")
    return names


def requested_motifs(arguments: argparse.Namespace) -> list[str]:
    motifs = list(arguments.motif or [])
    if arguments.motif_list is not None:
        with arguments.motif_list.open(encoding="ascii") as handle:
            for line_number, line in enumerate(handle, 1):
                value = line.strip()
                if not value or value.startswith("#"):
                    continue
                if len(value.split()) != 1:
                    raise ThresholdBuildError(
                        f"invalid motif-list line {line_number}: {line.rstrip()}"
                    )
                motifs.append(value)
    motifs = list(dict.fromkeys(motifs))
    if not motifs:
        raise ThresholdBuildError("at least one --motif or --motif-list is required")
    return motifs


def metrics_sql(path: Path) -> str:
    if path.suffix.lower() == ".parquet":
        return f"read_parquet({sql_string(path)})"
    return (
        f"read_csv_auto({sql_string(path)}, delim='\\t', header=true, "
        "nullstr='NA')"
    )


def build_sql(
    arguments: argparse.Namespace,
    motif_names: list[tuple[str, str]],
    staging_output: Path,
    metrics_sha256: str,
    jaspar_sha256: str,
) -> str:
    minimum_distance_sql = (
        "NULL::BIGINT"
        if arguments.context_min_interval_distance is None
        else f"{arguments.context_min_interval_distance}::BIGINT"
    )
    motif_values = ",\n        ".join(
        f"({sql_string(motif_id)}, {sql_string(motif_name)})"
        for motif_id, motif_name in motif_names
    )
    metric = arguments.selection_metric
    augmented_metric = {
        "delta_macro_roc_auc": "augmented_macro_roc_auc",
        "delta_macro_average_precision": "augmented_macro_average_precision",
        "delta_macro_log_loss": "augmented_macro_log_loss",
    }[metric]
    baseline_metric = {
        "delta_macro_roc_auc": "baseline_macro_roc_auc",
        "delta_macro_average_precision": "baseline_macro_average_precision",
        "delta_macro_log_loss": "baseline_macro_log_loss",
    }[metric]
    order = "ASC" if metric == "delta_macro_log_loss" else "DESC"
    positive_test = "< 0" if metric == "delta_macro_log_loss" else "> 0"
    near_test = (
        f"m.metric_gain <= b.metric_gain * {sql_number(arguments.near_optimal_fraction)}"
        if metric == "delta_macro_log_loss"
        else f"m.metric_gain >= b.metric_gain * {sql_number(arguments.near_optimal_fraction)}"
    )
    return f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};

COPY (
WITH expected(motif_id, motif_name) AS (
    VALUES
        {motif_values}
),
metric_input AS (
    SELECT
        CAST(motif_id AS VARCHAR) AS motif_id,
        CAST(threshold AS DOUBLE) AS threshold,
        CAST(anchors_total AS BIGINT) AS anchors_total,
        CAST(anchors_retained AS BIGINT) AS anchors_retained,
        CAST(retained_fraction AS DOUBLE) AS retained_fraction,
        CAST(discordant_observations AS BIGINT) AS discordant_observations,
        {nullable_metric_cast(baseline_metric, "DOUBLE")} AS baseline_metric_value,
        {nullable_metric_cast(augmented_metric, "DOUBLE")} AS augmented_metric_value,
        {nullable_metric_cast(metric, "DOUBLE")} AS metric_gain,
        {nullable_metric_cast("median_adjusted_odds_ratio", "DOUBLE")}
            AS adjusted_odds_ratio,
        {nullable_metric_cast("samples_with_raw_odds_ratio_below_one", "BIGINT")}
            AS samples_with_raw_odds_ratio_below_one,
        {nullable_metric_cast("samples_total", "BIGINT")} AS samples_total,
        {nullable_metric_cast("sample_fold_cells", "BIGINT")} AS sample_fold_cells,
        {nullable_metric_cast("sample_fold_cells_with_roc_auc_gain", "BIGINT")}
            AS sample_fold_cells_with_roc_auc_gain,
        {nullable_metric_cast("sample_fold_cells_with_log_loss_gain", "BIGINT")}
            AS sample_fold_cells_with_log_loss_gain
    FROM {metrics_sql(arguments.metrics)}
),
candidate_stats AS (
    SELECT
        motif_id,
        MIN(threshold) AS candidate_threshold_min,
        MAX(threshold) AS candidate_threshold_max,
        COUNT(*)::BIGINT AS candidate_threshold_count,
        COUNT(*) FILTER (
            WHERE retained_fraction >= {sql_number(arguments.minimum_class_fraction)}
              AND retained_fraction <= {sql_number(1 - arguments.minimum_class_fraction)}
        )::BIGINT AS eligible_threshold_count
    FROM metric_input
    GROUP BY motif_id
),
ranked AS (
    SELECT
        *,
        ROW_NUMBER() OVER (
            PARTITION BY motif_id
            ORDER BY metric_gain {order}, threshold ASC
        ) AS metric_rank
    FROM metric_input
    WHERE isfinite(metric_gain)
      AND retained_fraction >= {sql_number(arguments.minimum_class_fraction)}
      AND retained_fraction <= {sql_number(1 - arguments.minimum_class_fraction)}
),
best AS (
    SELECT * FROM ranked WHERE metric_rank = 1
),
near_optimal AS (
    SELECT
        m.motif_id,
        MIN(m.threshold) AS useful_threshold_min,
        MAX(m.threshold) AS useful_threshold_max
    FROM metric_input m
    JOIN best b USING (motif_id)
    WHERE b.metric_gain {positive_test}
      AND {near_test}
      AND m.retained_fraction >= {sql_number(arguments.minimum_class_fraction)}
      AND m.retained_fraction <= {sql_number(1 - arguments.minimum_class_fraction)}
    GROUP BY m.motif_id
)
SELECT
    1::INTEGER AS schema_version,
    {sql_string(arguments.threshold_set_id)}::VARCHAR AS threshold_set_id,
    {sql_string(arguments.genome_id)}::VARCHAR AS genome_id,
    {sql_string(arguments.motif_set_id)}::VARCHAR AS motif_set_id,
    e.motif_id,
    e.motif_name,
    {sql_string(arguments.threshold_role)}::VARCHAR AS threshold_role,
    {sql_string(arguments.target_motif_id)}::VARCHAR AS target_motif_id,
    {sql_string(arguments.target_motif_name)}::VARCHAR AS target_motif_name,
    {sql_string(arguments.score_mode)}::VARCHAR AS score_mode,
    {sql_number(arguments.pseudocount)}::DOUBLE AS pseudocount,
    {sql_string(arguments.background_model_id)}::VARCHAR AS background_model_id,
    {sql_string(arguments.pseudocount_scheme)}::VARCHAR AS pseudocount_scheme,
    {sql_number(arguments.source_minimum_score)}::DOUBLE AS source_minimum_score,
    TRUE AS threshold_inclusive,
    {arguments.context_flank}::BIGINT AS context_flank_bp,
    {sql_string(arguments.context_distance_metric)}::VARCHAR
        AS context_distance_metric,
    {minimum_distance_sql} AS context_min_interval_distance_bp,
    {arguments.context_max_interval_distance}::BIGINT
        AS context_max_interval_distance_bp,
    {sql_string(arguments.context_relation_filter)}::VARCHAR
        AS context_relation_filter,
    {sql_string(arguments.calibration_stratum_id)}::VARCHAR
        AS calibration_stratum_id,
    {sql_string(arguments.calibration_stratum_json)}::JSON
        AS calibration_stratum,
    {sql_string(arguments.calibration_scope)}::VARCHAR AS calibration_scope,
    {sql_string(arguments.evidence_dataset_id)}::VARCHAR AS evidence_dataset_id,
    {sql_string(arguments.outcome_id)}::VARCHAR AS outcome_id,
    {sql_string(arguments.fold_definition)}::VARCHAR AS fold_definition,
    {arguments.folds}::INTEGER AS n_folds,
    {sql_string(arguments.candidate_grid)}::VARCHAR AS candidate_grid,
    c.candidate_threshold_min,
    c.candidate_threshold_max,
    c.candidate_threshold_count,
    c.eligible_threshold_count,
    {sql_number(arguments.minimum_class_fraction)}::DOUBLE
        AS minimum_class_fraction,
    {sql_string(metric)}::VARCHAR AS selection_metric,
    'optimize_metric_then_lower_threshold'::VARCHAR AS selection_rule,
    {sql_number(arguments.near_optimal_fraction)}::DOUBLE AS near_optimal_fraction,
    CASE WHEN b.metric_gain {positive_test} THEN b.threshold END
        AS recommended_threshold,
    CASE WHEN b.metric_gain {positive_test} THEN n.useful_threshold_min END
        AS useful_threshold_min,
    CASE WHEN b.metric_gain {positive_test} THEN n.useful_threshold_max END
        AS useful_threshold_max,
    b.baseline_metric_value,
    b.augmented_metric_value AS selected_metric_value,
    b.metric_gain AS selected_metric_gain,
    b.anchors_total AS n_anchors,
    b.discordant_observations AS n_observations,
    b.anchors_retained AS selected_retained_anchor_count,
    b.retained_fraction AS selected_retained_fraction,
    b.adjusted_odds_ratio AS selected_adjusted_odds_ratio,
    CASE
        WHEN b.adjusted_odds_ratio < 1 THEN 'inhibitory_association'
        WHEN b.adjusted_odds_ratio > 1 THEN 'facilitative_association'
        WHEN b.adjusted_odds_ratio = 1 THEN 'neutral_association'
    END AS association_direction,
    CASE
        WHEN b.adjusted_odds_ratio < 1
            THEN b.samples_with_raw_odds_ratio_below_one
        WHEN b.adjusted_odds_ratio > 1
            THEN b.samples_total - b.samples_with_raw_odds_ratio_below_one
    END AS samples_consistent_with_direction,
    b.samples_total,
    b.sample_fold_cells,
    b.sample_fold_cells_with_roc_auc_gain,
    b.sample_fold_cells_with_log_loss_gain,
    CASE
        WHEN c.motif_id IS NULL THEN 'pending'
        WHEN c.eligible_threshold_count = 0 THEN 'insufficient_class_support'
        WHEN b.motif_id IS NULL THEN 'no_finite_metric'
        WHEN b.metric_gain {positive_test} THEN 'exploratory_positive_gain'
        ELSE 'no_positive_gain'
    END AS calibration_status,
    {sql_string(arguments.calibration_run_id)}::VARCHAR AS calibration_run_id,
    {sql_string(arguments.source_metrics_uri)}::VARCHAR AS source_metrics_uri,
    {sql_string(metrics_sha256)}::VARCHAR AS source_metrics_sha256,
    {sql_string(arguments.jaspar_uri)}::VARCHAR AS jaspar_uri,
    {sql_string(jaspar_sha256)}::VARCHAR AS jaspar_sha256,
    {sql_string(arguments.source_commit)}::VARCHAR AS source_commit,
    {str(arguments.source_dirty).upper()}::BOOLEAN AS source_dirty,
    {sql_string(arguments.notes)}::VARCHAR AS notes
FROM expected e
LEFT JOIN candidate_stats c USING (motif_id)
LEFT JOIN best b USING (motif_id)
LEFT JOIN near_optimal n USING (motif_id)
ORDER BY e.motif_id
) TO {sql_string(staging_output)} (FORMAT PARQUET, COMPRESSION ZSTD);
"""


def run_duckdb(duckdb: str, sql: str) -> None:
    process = subprocess.run(
        [duckdb, "-light-mode", ":memory:", "-c", sql],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ThresholdBuildError(
            process.stderr.strip() or "DuckDB threshold build failed"
        )


def validate_output(duckdb: str, output: Path, expected_rows: int) -> None:
    query = f"""
WITH t AS (SELECT * FROM read_parquet({sql_string(output)}))
SELECT
    COUNT(*)::BIGINT AS rows,
    COUNT(DISTINCT (
        threshold_set_id, genome_id, motif_set_id, motif_id, threshold_role,
        target_motif_id, score_mode, pseudocount, background_model_id,
        pseudocount_scheme, calibration_stratum_id
    ))::BIGINT AS keys,
    COUNT(*) FILTER (
        WHERE recommended_threshold < source_minimum_score
           OR useful_threshold_min > recommended_threshold
           OR useful_threshold_max < recommended_threshold
           OR threshold_inclusive IS DISTINCT FROM TRUE
    )::BIGINT AS invalid_thresholds,
    COUNT(*) FILTER (
        WHERE (calibration_status = 'exploratory_positive_gain'
               AND recommended_threshold IS NULL)
           OR (calibration_status <> 'exploratory_positive_gain'
               AND recommended_threshold IS NOT NULL)
    )::BIGINT AS invalid_status_nullability
FROM t;
"""
    process = subprocess.run(
        [duckdb, "-light-mode", "-json", ":memory:", "-c", query],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ThresholdBuildError(
            process.stderr.strip() or "DuckDB threshold validation failed"
        )
    try:
        rows = json.loads(process.stdout)
    except json.JSONDecodeError as error:
        raise ThresholdBuildError("DuckDB returned invalid validation JSON") from error
    if len(rows) != 1:
        raise ThresholdBuildError("DuckDB returned an invalid validation row count")
    result = {key: int(value) for key, value in rows[0].items()}
    if result["rows"] != expected_rows or result["keys"] != expected_rows:
        raise ThresholdBuildError("threshold output is incomplete or has duplicate keys")
    if result["invalid_thresholds"] or result["invalid_status_nullability"]:
        raise ThresholdBuildError("threshold output violates nullability or range rules")


def build(arguments: argparse.Namespace) -> None:
    if shutil.which(arguments.duckdb) is None:
        raise ThresholdBuildError(f"DuckDB executable not found: {arguments.duckdb}")
    if arguments.source_metrics_uri is None:
        arguments.source_metrics_uri = str(arguments.metrics)
    if arguments.jaspar_uri is None:
        arguments.jaspar_uri = str(arguments.jaspar)
    arguments.metrics = arguments.metrics.expanduser().resolve()
    arguments.jaspar = arguments.jaspar.expanduser().resolve()
    arguments.output = arguments.output.expanduser().resolve()
    if not arguments.metrics.is_file():
        raise ThresholdBuildError(f"metrics file not found: {arguments.metrics}")
    if not arguments.jaspar.is_file():
        raise ThresholdBuildError(f"JASPAR file not found: {arguments.jaspar}")
    if arguments.output.exists():
        raise ThresholdBuildError(f"output already exists: {arguments.output}")

    motif_ids = requested_motifs(arguments)
    names = parse_jaspar_names(arguments.jaspar)
    missing = [motif_id for motif_id in motif_ids if motif_id not in names]
    if missing:
        raise ThresholdBuildError(
            "motifs absent from JASPAR: " + ", ".join(missing)
        )
    if arguments.target_motif_id not in names:
        raise ThresholdBuildError(
            f"target motif absent from JASPAR: {arguments.target_motif_id}"
        )
    arguments.target_motif_name = names[arguments.target_motif_id]
    motif_names = [(motif_id, names[motif_id]) for motif_id in motif_ids]

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, staging_name = tempfile.mkstemp(
        prefix=f".{arguments.output.name}.", suffix=".tmp",
        dir=arguments.output.parent,
    )
    os.close(descriptor)
    staging_output = Path(staging_name)
    staging_output.unlink()
    try:
        sql = build_sql(
            arguments, motif_names, staging_output,
            sha256_file(arguments.metrics), sha256_file(arguments.jaspar),
        )
        run_duckdb(arguments.duckdb, sql)
        validate_output(arguments.duckdb, staging_output, len(motif_names))
        os.replace(staging_output, arguments.output)
    finally:
        if staging_output.exists():
            staging_output.unlink()

    print(
        f"I: Wrote {len(motif_names)} motif threshold rows: {arguments.output}",
        file=sys.stderr,
    )


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Select one provenance-rich, inclusive motif-score threshold per "
            "requested motif from an evaluator metric table."
        )
    )
    parser.add_argument(
        "--metrics", required=True, type=Path,
        help="threshold-metric TSV or Parquet produced by the evaluator",
    )
    parser.add_argument(
        "--source-metrics-uri",
        help="portable provenance URI; defaults to the --metrics argument",
    )
    parser.add_argument(
        "--jaspar", required=True, type=Path,
        help="JASPAR-format motif file used to resolve stable motif names",
    )
    parser.add_argument(
        "--jaspar-uri",
        help="portable JASPAR provenance URI; defaults to the --jaspar argument",
    )
    parser.add_argument(
        "--motif", action="append",
        help="expected motif accession; repeat to preserve pending motifs",
    )
    parser.add_argument(
        "--motif-list", type=Path,
        help="file containing one expected motif accession per line",
    )
    parser.add_argument(
        "--output", required=True, type=Path,
        help="new Parquet registry file; existing files are refused",
    )
    parser.add_argument(
        "--threshold-set-id", required=True,
        help="immutable, versioned identifier for this threshold set",
    )
    parser.add_argument(
        "--calibration-run-id", required=True,
        help="identifier of the metric-producing calibration run",
    )
    parser.add_argument("--genome-id", required=True, help="stable genome identity")
    parser.add_argument(
        "--motif-set-id", required=True, help="stable motif-set identity"
    )
    parser.add_argument(
        "--target-motif-id", default="MA0861.2",
        help="motif whose experimental support is predicted (default: MA0861.2)",
    )
    parser.add_argument(
        "--threshold-role", default="tp73_context_binary_feature",
        help="declared downstream use of the threshold",
    )
    parser.add_argument(
        "--score-mode", default="log2_relative_risk",
        help="score semantics used by the source scan",
    )
    parser.add_argument(
        "--pseudocount", type=float, default=1.0,
        help="per-base pseudocount used by the source scan (default: 1)",
    )
    parser.add_argument(
        "--background-model-id", default="uniform_acgt_v1",
        help="scoring-background identity",
    )
    parser.add_argument(
        "--pseudocount-scheme", default="additive_per_base",
        help="pseudocount semantics",
    )
    parser.add_argument(
        "--source-minimum-score", type=float, default=-1.0,
        help="inclusive storage floor of the source motif scans (default: -1)",
    )
    parser.add_argument(
        "--context-flank", type=int, default=150,
        help="reported context radius in bp (default: 150)",
    )
    parser.add_argument(
        "--context-distance-metric", default="signed_interval_edge_distance",
        help="distance geometry used to select neighboring hits",
    )
    parser.add_argument(
        "--context-min-interval-distance", type=int,
        help="optional inclusive lower signed interval-distance bound",
    )
    parser.add_argument(
        "--context-max-interval-distance", type=int,
        help="inclusive upper bound; defaults to --context-flank",
    )
    parser.add_argument(
        "--context-relation-filter", default="any",
        help="interval relation subset represented by the metrics (default: any)",
    )
    parser.add_argument(
        "--calibration-stratum-id", default="all_tp73_anchors",
        help="stable ID for the biological subset represented by each row",
    )
    parser.add_argument(
        "--calibration-stratum-json",
        default='{"anchor_pair_class":"all","genomic_region":"all"}',
        help="JSON object defining the declared biological subset",
    )
    parser.add_argument(
        "--calibration-scope", required=True,
        help="genomic scope used to fit/select thresholds, e.g. grch38_chr1",
    )
    parser.add_argument(
        "--evidence-dataset-id", required=True,
        help="stable identity of the experimental evidence collection",
    )
    parser.add_argument(
        "--outcome-id",
        default="discordant_anti_p73_only_vs_matched_control_only_immersion",
        help="modeled outcome definition",
    )
    parser.add_argument(
        "--fold-definition", default="equal_width_contiguous_chromosome_spans",
        help="cross-validation partition semantics",
    )
    parser.add_argument("--folds", type=int, default=5, help="number of folds")
    parser.add_argument(
        "--candidate-grid", default="observed_integer_grid",
        help="human-readable candidate-grid provenance",
    )
    parser.add_argument(
        "--minimum-class-fraction", type=float, default=0.01,
        help="minimum retained and absent anchor fraction (default: 0.01)",
    )
    parser.add_argument(
        "--selection-metric",
        choices=(
            "delta_macro_roc_auc", "delta_macro_average_precision",
            "delta_macro_log_loss",
        ),
        default="delta_macro_roc_auc",
        help="held-out metric gain optimized by the selector",
    )
    parser.add_argument(
        "--near-optimal-fraction", type=float, default=0.9,
        help="fraction of best positive gain defining the useful range",
    )
    parser.add_argument(
        "--source-commit", required=True,
        help="full Git commit underlying the calibration implementation",
    )
    parser.add_argument(
        "--source-dirty", action="store_true",
        help="record that the calibration used uncommitted source changes",
    )
    parser.add_argument("--notes", default="", help="free-text limitations")
    parser.add_argument("--duckdb", default="duckdb", help="DuckDB CLI command")
    parser.add_argument("--threads", type=int, default=2, help="DuckDB threads")
    parser.add_argument(
        "--memory-limit", default="1GB", help="DuckDB memory limit"
    )
    return parser


def main() -> int:
    parser = argument_parser()
    arguments = parser.parse_args()
    if arguments.motif is None and arguments.motif_list is None:
        parser.error("at least one --motif or --motif-list is required")
    if arguments.motif_list is not None:
        arguments.motif_list = arguments.motif_list.expanduser().resolve()
        if not arguments.motif_list.is_file():
            parser.error(f"--motif-list not found: {arguments.motif_list}")
    if re.fullmatch(r"[0-9a-f]{40}", arguments.source_commit) is None:
        parser.error("--source-commit must be a full lowercase Git SHA")
    if not (0 < arguments.near_optimal_fraction <= 1):
        parser.error("--near-optimal-fraction must be in (0,1]")
    if not (0 < arguments.minimum_class_fraction < 0.5):
        parser.error("--minimum-class-fraction must be in (0,0.5)")
    if arguments.context_flank < 0 or arguments.folds < 2:
        parser.error("--context-flank must be non-negative and --folds at least 2")
    if arguments.context_max_interval_distance is None:
        arguments.context_max_interval_distance = arguments.context_flank
    if (arguments.context_min_interval_distance is not None and
            arguments.context_min_interval_distance >
            arguments.context_max_interval_distance):
        parser.error(
            "--context-min-interval-distance cannot exceed the maximum"
        )
    try:
        stratum = json.loads(arguments.calibration_stratum_json)
    except json.JSONDecodeError as error:
        parser.error(f"--calibration-stratum-json is invalid: {error}")
    if not isinstance(stratum, dict):
        parser.error("--calibration-stratum-json must describe a JSON object")
    arguments.calibration_stratum_json = json.dumps(
        stratum, sort_keys=True, separators=(",", ":")
    )
    for option in (arguments.pseudocount, arguments.source_minimum_score):
        if not (-float("inf") < option < float("inf")):
            parser.error("numeric scoring arguments must be finite")
    try:
        build(arguments)
        return 0
    except (ThresholdBuildError, OSError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

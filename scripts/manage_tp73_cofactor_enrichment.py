#!/usr/bin/env python3
"""Plan, inspect, and finalize all-JASPAR TP73 cofactor enrichment runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any


OUTPUT_FILES = (
    "class_counts.tsv",
    "depth_tier_manifest.tsv",
    "descriptive.tsv",
    "macro_summary.tsv",
    "primary_occupancy.tsv",
    "evaluator_run_config.tsv",
)


class EnrichmentError(RuntimeError):
    """Raised when an immutable run contract is not satisfied."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_path_list(paths: list[Path]) -> str:
    return "[" + ",".join(sql_string(path) for path in paths) + "]"


def run_json(command: list[str]) -> list[dict[str, Any]]:
    process = subprocess.run(
        command, text=True, capture_output=True, check=False
    )
    if process.returncode != 0:
        raise EnrichmentError(process.stderr.strip() or "command failed")
    try:
        value = json.loads(process.stdout)
    except json.JSONDecodeError as error:
        raise EnrichmentError(f"invalid JSON command output: {error}") from error
    if not isinstance(value, list):
        raise EnrichmentError("JSON command output is not a row array")
    return value


def immutable_write(path: Path, content: str) -> None:
    encoded = content.encode("utf-8")
    if path.exists():
        if path.read_bytes() != encoded:
            raise EnrichmentError(f"immutable file differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        temporary.write_bytes(encoded)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def git_identity(source: Path) -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        text=True, capture_output=True, check=False,
    )
    if commit.returncode != 0:
        raise EnrichmentError(commit.stderr.strip() or "cannot read Git commit")
    dirty = False
    for arguments in (
        ["diff", "--quiet", "--ignore-submodules", "--"],
        ["diff", "--cached", "--quiet", "--ignore-submodules", "--"],
    ):
        status = subprocess.run(
            ["git", "-C", str(source), *arguments],
            text=True, capture_output=True, check=False,
        )
        if status.returncode not in (0, 1):
            raise EnrichmentError(status.stderr.strip() or "cannot read Git status")
        dirty = dirty or status.returncode == 1
    return commit.stdout.strip(), dirty


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise EnrichmentError(f"JSON file not found: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise EnrichmentError(f"JSON value is not an object: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise EnrichmentError(f"TSV file not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if reader.fieldnames is None or not rows:
        raise EnrichmentError(f"TSV is empty or lacks a header: {path}")
    return rows


def source_task_directory(source_run: Path, row: dict[str, str]) -> Path:
    return source_run / "tasks" / (
        f"task-{int(row['task_index']):06d}-{row['motif_id']}"
    )


def task_directory(run_root: Path, row: dict[str, str]) -> Path:
    return run_root / "tasks" / (
        f"task-{int(row['task_index']):06d}-{row['motif_id']}"
    )


def planned_tasks(run_root: Path) -> list[dict[str, str]]:
    rows = read_tsv(run_root / "plan" / "enrichment_tasks.tsv")
    expected = list(range(len(rows)))
    observed = [int(row["task_index"]) for row in rows]
    if observed != expected or len({row["motif_id"] for row in rows}) != len(rows):
        raise EnrichmentError("task indices are not contiguous or motifs are duplicated")
    return rows


def validate_source_task(
    source_run: Path, row: dict[str, str], verify_checksum: bool = False
) -> tuple[Path, dict[str, Any], Path]:
    directory = source_task_directory(source_run, row)
    marker_path = directory / "complete.json"
    marker = load_json(marker_path)
    if marker.get("motif_id") != row["motif_id"]:
        raise EnrichmentError(f"source marker motif mismatch: {marker_path}")
    record = marker.get("files", {}).get("cofactor_maxima.parquet")
    if not isinstance(record, dict):
        raise EnrichmentError(f"source marker lacks cofactor maxima: {marker_path}")
    maxima = directory / "cofactor_maxima.parquet"
    if not maxima.is_file() or maxima.stat().st_size != int(record.get("bytes", -1)):
        raise EnrichmentError(f"source maxima size mismatch: {maxima}")
    expected_digest = str(record.get("sha256", ""))
    if not re.fullmatch(r"[0-9a-f]{64}", expected_digest):
        raise EnrichmentError(f"invalid source maxima checksum: {marker_path}")
    if verify_checksum and sha256(maxima) != expected_digest:
        raise EnrichmentError(f"source maxima checksum mismatch: {maxima}")
    return maxima, marker, marker_path


def validate_output_marker(
    directory: Path, row: dict[str, str], verify_checksums: bool = False
) -> dict[str, Any]:
    marker_path = directory / "complete.json"
    marker = load_json(marker_path)
    if (marker.get("schema_version") != 1
            or marker.get("task_index") != int(row["task_index"])
            or marker.get("motif_id") != row["motif_id"]
            or marker.get("source_maxima_sha256") !=
                row["source_maxima_sha256"]):
        raise EnrichmentError(f"task marker identity mismatch: {marker_path}")
    records = marker.get("files")
    if not isinstance(records, dict) or set(records) != set(OUTPUT_FILES):
        raise EnrichmentError(f"task marker output inventory differs: {marker_path}")
    for name in OUTPUT_FILES:
        record = records[name]
        path = directory / name
        if not isinstance(record, dict) or not path.is_file():
            raise EnrichmentError(f"task output is missing: {path}")
        if path.stat().st_size != int(record.get("bytes", -1)):
            raise EnrichmentError(f"task output size differs: {path}")
        digest = str(record.get("sha256", ""))
        if not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise EnrichmentError(f"task output checksum is invalid: {path}")
        if verify_checksums and sha256(path) != digest:
            raise EnrichmentError(f"task output checksum differs: {path}")
    return marker


def prepare(arguments: argparse.Namespace) -> None:
    duckdb = shutil.which(arguments.duckdb)
    if duckdb is None:
        raise EnrichmentError(f"DuckDB executable not found: {arguments.duckdb}")
    run_root = arguments.run_root.expanduser().resolve()
    source_run = arguments.source_threshold_run.expanduser().resolve()
    source = arguments.source.expanduser().resolve()
    if not source_run.is_dir() or not source.is_dir():
        raise EnrichmentError("source threshold run or repository is missing")
    for name in ("plan", "logs", "tasks", "final"):
        (run_root / name).mkdir(parents=True, exist_ok=True)

    source_config_path = source_run / "plan" / "run_config.json"
    source_config = load_json(source_config_path)
    source_plan_path = source_run / "plan" / "calibration_tasks.tsv"
    source_rows = read_tsv(source_plan_path)
    anchor = Path(str(source_config.get("anchor_evidence", ""))).resolve()
    if not anchor.is_file():
        raise EnrichmentError(f"source TP73 anchor evidence is missing: {anchor}")
    registry = (source_run / "final" / "threshold_calibration" / "tables" /
                "jaspar2026" / "motif_score_threshold" / "part-000000.parquet")
    if not registry.is_file():
        raise EnrichmentError(f"final threshold registry is missing: {registry}")

    registry_rows = run_json([
        duckdb, "-light-mode", "-json", ":memory:", "-c",
        "SELECT motif_id::VARCHAR AS motif_id, motif_name::VARCHAR AS motif_name, "
        "recommended_threshold::DOUBLE AS recommended_threshold, "
        "calibration_status::VARCHAR AS calibration_status "
        f"FROM read_parquet({sql_string(registry)}) ORDER BY motif_id;",
    ])
    registry_by_motif = {str(row["motif_id"]): row for row in registry_rows}
    source_ids = [row["motif_id"] for row in source_rows]
    if len(source_ids) != len(set(source_ids)) or set(source_ids) != set(registry_by_motif):
        raise EnrichmentError("source task plan and final threshold registry differ")

    fields = [
        "task_index", "motif_id", "factor_name", "source_task_index",
        "source_task_relative_path", "source_marker_sha256",
        "source_maxima_bytes", "source_maxima_sha256",
        "recommended_threshold", "positive_threshold",
        "positive_threshold_source", "positive_threshold_role",
        "selection_semantics", "source_calibration_status",
    ]
    output = io.StringIO(newline="")
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    fallback_count = 0
    for task_index, source_row in enumerate(source_rows):
        motif_id = source_row["motif_id"]
        maxima, marker, marker_path = validate_source_task(source_run, source_row)
        registry_row = registry_by_motif[motif_id]
        recommendation = registry_row.get("recommended_threshold")
        if recommendation is None:
            fallback_count += 1
            recommended_text = "NA"
            positive_threshold = 0.0
            threshold_source = "fallback_zero_no_registry_recommendation"
            threshold_role = "descriptive_fallback"
            semantics = (
                "fallback_zero_for_all_motif_coverage;"
                "registry_recommendation_remains_null"
            )
        else:
            recommended_text = format(float(recommendation), ".17g")
            positive_threshold = float(recommendation)
            threshold_source = "historical_threshold_complement_recommendation"
            threshold_role = "historical_operating_point"
            semantics = (
                "fixed_historical_operating_point;exploratory_reuse_on_chr1"
            )
        relative = maxima.relative_to(source_run)
        writer.writerow({
            "task_index": task_index,
            "motif_id": motif_id,
            "factor_name": registry_row.get("motif_name") or source_row["motif_name"],
            "source_task_index": source_row["task_index"],
            "source_task_relative_path": str(relative.parent),
            "source_marker_sha256": sha256(marker_path),
            "source_maxima_bytes": marker["files"]["cofactor_maxima.parquet"]["bytes"],
            "source_maxima_sha256": marker["files"]["cofactor_maxima.parquet"]["sha256"],
            "recommended_threshold": recommended_text,
            "positive_threshold": format(positive_threshold, ".17g"),
            "positive_threshold_source": threshold_source,
            "positive_threshold_role": threshold_role,
            "selection_semantics": semantics,
            "source_calibration_status": registry_row.get("calibration_status") or "unknown",
        })
    task_path = run_root / "plan" / "enrichment_tasks.tsv"
    immutable_write(task_path, output.getvalue())

    commit, dirty = git_identity(source)
    config = {
        "schema_version": 1,
        "run_id": arguments.run_id,
        "analysis_scope": "all_non_target_jaspar2026_motifs_vs_tp73_chr1_cutrun",
        "task_count": len(source_rows),
        "fallback_zero_task_count": fallback_count,
        "source_threshold_run": str(source_run),
        "source_threshold_run_config": str(source_config_path),
        "source_threshold_run_config_sha256": sha256(source_config_path),
        "source_threshold_task_plan": str(source_plan_path),
        "source_threshold_task_plan_sha256": sha256(source_plan_path),
        "source_threshold_registry": str(registry),
        "source_threshold_registry_sha256": sha256(registry),
        "anchor_evidence": str(anchor),
        "anchor_evidence_sha256": sha256(anchor),
        "source": str(source),
        "source_commit": commit,
        "source_dirty": dirty,
        "source_dirty_scope": "tracked_and_staged_files_only",
        "primary_negative_reference": -1,
        "secondary_negative_reference": 0,
        "negative_reference_semantics": "strict context_score < N or absent",
        "positive_semantics": "inclusive context_score >= T",
        "intermediate_semantics": "excluded N <= context_score < T",
        "null_recommendation_policy": "explicit threshold-zero descriptive fallback",
        "context_flank_bp": 150,
        "context_geometry": "signed_interval_edge_distance",
        "tp73_score_breaks": [-5, -1, 0, 1, 2, 5, 10, 15, "Inf"],
        "depth_quantiles": [0.5, 0.75, 0.9, 0.95, 0.99],
        "block_size_bp": arguments.block_size,
        "spline_df": arguments.spline_df,
        "minimum_class_fraction": arguments.minimum_class_fraction,
        "multiple_testing_family_size": len(source_rows),
        "multiple_testing_method": "Benjamini-Hochberg across all planned motifs",
        "inference_status": (
            "exploratory chr1; operating points selected on chr1; "
            "no independent validation"
        ),
        "enrichment_task_plan_sha256": sha256(task_path),
    }
    immutable_write(
        run_root / "plan" / "run_config.json",
        json.dumps(config, indent=2, sort_keys=True) + "\n",
    )
    print(len(source_rows))


def status(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    rows = planned_tasks(run_root)
    complete = 0
    invalid = 0
    for row in rows:
        directory = task_directory(run_root, row)
        if not (directory / "complete.json").is_file():
            continue
        try:
            validate_output_marker(directory, row)
            complete += 1
        except (EnrichmentError, OSError, ValueError, json.JSONDecodeError):
            invalid += 1
    print("key\tvalue")
    print(f"planned\t{len(rows)}")
    print(f"complete\t{complete}")
    print(f"pending\t{len(rows) - complete - invalid}")
    print(f"invalid_complete_markers\t{invalid}")


def finalize(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(run_root / "plan" / "run_config.json")
    rows = planned_tasks(run_root)
    if int(config.get("task_count", -1)) != len(rows):
        raise EnrichmentError("run configuration task count differs from its plan")
    finalization_commit = arguments.finalization_source_commit
    if re.fullmatch(r"[0-9a-f]{40}", finalization_commit) is None:
        raise EnrichmentError("finalization source commit must be a full Git object ID")

    files: dict[str, list[Path]] = {name: [] for name in OUTPUT_FILES}
    depth_manifest_digest: str | None = None
    incomplete: list[str] = []
    for row in rows:
        directory = task_directory(run_root, row)
        if not (directory / "complete.json").is_file():
            incomplete.append(row["motif_id"])
            continue
        marker = validate_output_marker(directory, row, verify_checksums=True)
        current_depth_digest = marker["files"]["depth_tier_manifest.tsv"]["sha256"]
        if depth_manifest_digest is None:
            depth_manifest_digest = current_depth_digest
        elif current_depth_digest != depth_manifest_digest:
            raise EnrichmentError(
                f"depth-tier manifests differ at motif {row['motif_id']}"
            )
        for name in OUTPUT_FILES:
            files[name].append(directory / name)
    if incomplete:
        raise EnrichmentError(
            f"cannot finalize; {len(incomplete)} tasks are incomplete, "
            f"first={incomplete[:5]}"
        )

    duckdb = shutil.which(arguments.duckdb)
    if duckdb is None:
        raise EnrichmentError(f"DuckDB executable not found: {arguments.duckdb}")
    final = run_root / "final" / "cofactor_enrichment"
    if final.exists():
        if (final / "manifest.json").is_file():
            print(f"I: Reusing finalized cofactor enrichment: {final}", file=sys.stderr)
            return
        raise EnrichmentError(f"incomplete final output already exists: {final}")
    final.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".cofactor-enrichment.", dir=final.parent))
    try:
        database = staging / "tp73_cofactor_enrichment.duckdb"
        table_root = staging / "tables" / "jaspar2026"
        table_root.mkdir(parents=True)
        table_paths = {
            "cofactor_class_count": table_root / "cofactor_class_count" / "part-000000.parquet",
            "cofactor_depth_tier": table_root / "cofactor_depth_tier" / "part-000000.parquet",
            "cofactor_descriptive": table_root / "cofactor_descriptive" / "part-000000.parquet",
            "cofactor_macro_summary": table_root / "cofactor_macro_summary" / "part-000000.parquet",
            "cofactor_primary_occupancy": table_root / "cofactor_primary_occupancy" / "part-000000.parquet",
            "cofactor_overview": table_root / "cofactor_overview" / "part-000000.parquet",
            "cofactor_task": table_root / "cofactor_task" / "part-000000.parquet",
        }
        for path in table_paths.values():
            path.parent.mkdir(parents=True, exist_ok=True)
        task_count = len(rows)
        csv_options = "delim='\\t', header=true, nullstr=['NA',''], union_by_name=true"
        sql = f"""
SET preserve_insertion_order=false;
CREATE TABLE cofactor_task AS
SELECT * FROM read_csv_auto(
    {sql_string(run_root / 'plan' / 'enrichment_tasks.tsv')},
    delim='\\t', header=true, nullstr=['NA','']
);
CREATE TABLE cofactor_class_count AS
SELECT * REPLACE (
    try_cast(positive_threshold AS DOUBLE) AS positive_threshold,
    try_cast(negative_reference_threshold AS DOUBLE)
        AS negative_reference_threshold,
    try_cast(anchors_total AS BIGINT) AS anchors_total,
    try_cast(anchors_positive AS BIGINT) AS anchors_positive,
    try_cast(anchors_negative_reference AS BIGINT) AS anchors_negative_reference,
    try_cast(anchors_intermediate AS BIGINT) AS anchors_intermediate
)
FROM read_csv_auto({sql_path_list(files['class_counts.tsv'])}, {csv_options});
CREATE TABLE cofactor_depth_tier AS
SELECT * FROM read_csv_auto({sql_string(files['depth_tier_manifest.tsv'][0])},
    delim='\\t', header=true, nullstr=['NA','']);
CREATE TABLE cofactor_descriptive AS
SELECT * FROM read_csv_auto({sql_path_list(files['descriptive.tsv'])}, {csv_options});
CREATE TABLE cofactor_macro_summary AS
SELECT * REPLACE (
    try_cast(positive_threshold AS DOUBLE) AS positive_threshold,
    try_cast(negative_reference_threshold AS DOUBLE)
        AS negative_reference_threshold,
    try_cast(depth_tier_order AS INTEGER) AS depth_tier_order,
    try_cast(samples_total AS BIGINT) AS samples_total,
    try_cast(samples_estimable AS BIGINT) AS samples_estimable,
    try_cast(mean_log2_anti_control_specificity_ratio_jeffreys AS DOUBLE)
        AS mean_log2_anti_control_specificity_ratio_jeffreys,
    try_cast(median_log2_anti_control_specificity_ratio_jeffreys AS DOUBLE)
        AS median_log2_anti_control_specificity_ratio_jeffreys,
    try_cast(minimum_log2_anti_control_specificity_ratio_jeffreys AS DOUBLE)
        AS minimum_log2_anti_control_specificity_ratio_jeffreys,
    try_cast(maximum_log2_anti_control_specificity_ratio_jeffreys AS DOUBLE)
        AS maximum_log2_anti_control_specificity_ratio_jeffreys,
    try_cast(mean_anti_control_risk_difference_difference AS DOUBLE)
        AS mean_anti_control_risk_difference_difference,
    try_cast(samples_anti_p73_enriched AS BIGINT) AS samples_anti_p73_enriched,
    try_cast(samples_anti_p73_depleted AS BIGINT) AS samples_anti_p73_depleted
)
FROM read_csv_auto({sql_path_list(files['macro_summary.tsv'])}, {csv_options});
CREATE TEMP TABLE primary_raw AS
SELECT * REPLACE (
    try_cast(positive_threshold AS DOUBLE) AS positive_threshold,
    try_cast(negative_reference_threshold AS DOUBLE)
        AS negative_reference_threshold,
    try_cast(anchors_total AS BIGINT) AS anchors_total,
    try_cast(anchors_positive AS BIGINT) AS anchors_positive,
    try_cast(anchors_negative_reference AS BIGINT) AS anchors_negative_reference,
    try_cast(anchors_intermediate AS BIGINT) AS anchors_intermediate,
    try_cast(discordant_observations AS BIGINT) AS discordant_observations,
    try_cast(genomic_blocks AS BIGINT) AS genomic_blocks,
    try_cast(adjusted_log_odds AS DOUBLE) AS adjusted_log_odds,
    try_cast(adjusted_odds_ratio AS DOUBLE) AS adjusted_odds_ratio,
    try_cast(block_clustered_standard_error AS DOUBLE)
        AS block_clustered_standard_error,
    try_cast(confidence_interval_95_lower AS DOUBLE)
        AS confidence_interval_95_lower,
    try_cast(confidence_interval_95_upper AS DOUBLE)
        AS confidence_interval_95_upper,
    try_cast(p_value AS DOUBLE) AS p_value,
    try_cast(q_value_bh_panel AS DOUBLE) AS q_value_bh_panel,
    try_cast(requested_motifs_in_multiple_testing_scope AS BIGINT)
        AS requested_motifs_in_multiple_testing_scope,
    try_cast(estimable_motifs_in_multiple_testing_scope AS BIGINT)
        AS estimable_motifs_in_multiple_testing_scope
)
FROM read_csv_auto({sql_path_list(files['primary_occupancy.tsv'])}, {csv_options});
CREATE TABLE cofactor_primary_occupancy AS
WITH valid AS (
    SELECT motif_id, p_value,
           count(*) OVER (
               ORDER BY p_value
               RANGE BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
           ) AS bh_rank
    FROM primary_raw
    WHERE evaluation_status = 'ok' AND isfinite(p_value)
), adjusted AS (
    SELECT motif_id,
           least(1.0, min(p_value * {task_count} / bh_rank) OVER (
               ORDER BY p_value DESC
               ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
           )) AS q_value_bh_all_jaspar
    FROM valid
), counts AS (
    SELECT count(*)::BIGINT AS estimable_motifs FROM valid
)
SELECT p.* EXCLUDE (
           q_value_bh_panel,
           requested_motifs_in_multiple_testing_scope,
           estimable_motifs_in_multiple_testing_scope,
           multiple_testing_scope
       ),
       t.recommended_threshold,
       t.positive_threshold_role,
       t.source_calibration_status,
       a.q_value_bh_all_jaspar,
       {task_count}::BIGINT AS requested_motifs_in_multiple_testing_scope,
       c.estimable_motifs AS estimable_motifs_in_multiple_testing_scope,
       'all planned non-TP73 JASPAR 2026 motifs on chromosome 1; exploratory'
           AS multiple_testing_scope,
       CASE
         WHEN p.evaluation_status <> 'ok' THEN 'not_estimable'
         WHEN p.adjusted_odds_ratio > 1 THEN 'anti_p73_enriched'
         WHEN p.adjusted_odds_ratio < 1 THEN 'anti_p73_depleted'
         ELSE 'neutral'
       END AS association_direction
FROM primary_raw p
JOIN cofactor_task t USING (motif_id)
LEFT JOIN adjusted a USING (motif_id)
CROSS JOIN counts c
ORDER BY p.motif_id;

CREATE VIEW cofactor_overview AS
SELECT p.motif_id, p.factor_name, p.recommended_threshold,
       p.positive_threshold, p.positive_threshold_source,
       p.positive_threshold_role, p.source_calibration_status,
       c.anchors_total, c.anchors_positive, c.anchors_negative_reference,
       c.anchors_intermediate,
       c.anchors_positive / c.anchors_total::DOUBLE AS positive_anchor_fraction,
       c.anchors_negative_reference / c.anchors_total::DOUBLE
           AS negative_anchor_fraction,
       p.adjusted_odds_ratio, p.confidence_interval_95_lower,
       p.confidence_interval_95_upper, p.p_value,
       p.q_value_bh_all_jaspar, p.association_direction,
       p.evaluation_status, p.evaluation_note,
       m.mean_log2_anti_control_specificity_ratio_jeffreys
           AS mean_log2_specificity_ratio_strict_immersion,
       m.samples_anti_p73_enriched, m.samples_anti_p73_depleted
FROM cofactor_primary_occupancy p
JOIN cofactor_class_count c
  ON c.motif_id = p.motif_id
 AND c.negative_reference_threshold = -1
 AND c.tp73_score_stratum = 'all'
JOIN cofactor_macro_summary m
  ON m.motif_id = p.motif_id
 AND m.negative_reference_threshold = -1
 AND m.tp73_score_stratum = 'all'
 AND m.depth_tier = 'strict_immersion';

CREATE INDEX cofactor_primary_motif_idx
    ON cofactor_primary_occupancy(motif_id);
CREATE INDEX cofactor_class_motif_idx
    ON cofactor_class_count(motif_id, tp73_score_stratum);
CREATE INDEX cofactor_macro_motif_idx
    ON cofactor_macro_summary(motif_id, tp73_score_stratum, depth_tier);
CREATE INDEX cofactor_descriptive_motif_idx
    ON cofactor_descriptive(motif_id, tp73_score_stratum, sample_id, depth_tier);

COPY cofactor_task TO {sql_string(table_paths['cofactor_task'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_class_count TO {sql_string(table_paths['cofactor_class_count'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_depth_tier TO {sql_string(table_paths['cofactor_depth_tier'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_descriptive TO {sql_string(table_paths['cofactor_descriptive'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_macro_summary TO {sql_string(table_paths['cofactor_macro_summary'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_primary_occupancy TO {sql_string(table_paths['cofactor_primary_occupancy'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * FROM cofactor_overview ORDER BY motif_id)
  TO {sql_string(table_paths['cofactor_overview'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * FROM cofactor_overview ORDER BY motif_id)
  TO {sql_string(staging / 'cofactor_overview.tsv')}
  (FORMAT CSV, DELIMITER '\\t', HEADER, NULL 'NA');
COPY cofactor_primary_occupancy
  TO {sql_string(staging / 'primary_occupancy.tsv')}
  (FORMAT CSV, DELIMITER '\\t', HEADER, NULL 'NA');
CHECKPOINT;
"""
        sql_path = staging / "build.sql"
        sql_path.write_text(sql, encoding="utf-8")
        with sql_path.open("r", encoding="utf-8") as sql_input:
            process = subprocess.run(
                [duckdb, "-light-mode", "-batch", str(database)],
                stdin=sql_input, text=True, capture_output=True, check=False,
            )
        if process.returncode != 0:
            raise EnrichmentError(
                process.stderr.strip() or "DuckDB consolidation failed"
            )
        sql_path.unlink()

        counts = run_json([
            duckdb, "-light-mode", "-json", str(database), "-c",
            "SELECT "
            "(SELECT count(*) FROM cofactor_task) AS tasks, "
            "(SELECT count(*) FROM cofactor_primary_occupancy) AS primary_rows, "
            "(SELECT count(*) FROM cofactor_class_count) AS class_rows, "
            "(SELECT count(*) FROM cofactor_depth_tier) AS depth_rows, "
            "(SELECT count(*) FROM cofactor_descriptive) AS descriptive_rows, "
            "(SELECT count(*) FROM cofactor_macro_summary) AS macro_rows, "
            "(SELECT count(*) FROM cofactor_overview) AS overview_rows, "
            "(SELECT count(*) FROM cofactor_primary_occupancy "
            " WHERE evaluation_status='ok') AS estimable, "
            "(SELECT count(*) FROM cofactor_primary_occupancy "
            " WHERE positive_threshold_role='descriptive_fallback') AS fallbacks;",
        ])
        if len(counts) != 1:
            raise EnrichmentError("final row-count query did not return one row")
        count = counts[0]
        sample_count = len(read_tsv(files["depth_tier_manifest.tsv"][0])) // 12
        expected = {
            "tasks": task_count,
            "primary_rows": task_count,
            "class_rows": task_count * 18,
            "depth_rows": sample_count * 12,
            "descriptive_rows": task_count * sample_count * 108,
            "macro_rows": task_count * 108,
            "overview_rows": task_count,
            "fallbacks": int(config["fallback_zero_task_count"]),
        }
        for name, value in expected.items():
            if int(count[name]) != value:
                raise EnrichmentError(
                    f"final {name} count differs: {count[name]} != {value}"
                )
        manifest = {
            "schema_version": 1,
            "run_id": config["run_id"],
            "analysis_scope": config["analysis_scope"],
            "motifs": task_count,
            "estimable_primary_models": count["estimable"],
            "fallback_zero_motifs": count["fallbacks"],
            "multiple_testing_family_size": task_count,
            "multiple_testing_method": config["multiple_testing_method"],
            "source_threshold_run": config["source_threshold_run"],
            "source_threshold_registry_sha256":
                config["source_threshold_registry_sha256"],
            "anchor_evidence_sha256": config["anchor_evidence_sha256"],
            "metric_source_commit": config["source_commit"],
            "metric_source_dirty": config["source_dirty"],
            "finalization_source_commit": finalization_commit,
            "finalization_source_dirty": arguments.finalization_source_dirty,
            "source_task_payloads_copied": False,
            "row_counts": count,
            "database_sha256": sha256(database),
            "tables": {
                name: {
                    "path": str(path.relative_to(staging)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for name, path in table_paths.items()
            },
            "inference_status": config["inference_status"],
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, final)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    print(
        f"I: Finalized {len(rows)} cofactor enrichment results: {final}",
        file=sys.stderr,
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser(
        "prepare", help="write an immutable all-motif task plan"
    )
    prepare_parser.add_argument("--run-root", type=Path, required=True)
    prepare_parser.add_argument("--source-threshold-run", type=Path, required=True)
    prepare_parser.add_argument("--source", type=Path, required=True)
    prepare_parser.add_argument("--duckdb", default="duckdb")
    prepare_parser.add_argument("--block-size", type=int, default=5_000_000)
    prepare_parser.add_argument("--spline-df", type=int, default=4)
    prepare_parser.add_argument("--minimum-class-fraction", type=float, default=0.01)
    prepare_parser.add_argument(
        "--run-id", default="jaspar2026_chr1_tp73_cofactor_enrichment_v1"
    )
    prepare_parser.set_defaults(function=prepare)

    status_parser = subparsers.add_parser(
        "status", help="report exact planned/completed task counts"
    )
    status_parser.add_argument("--run-root", type=Path, required=True)
    status_parser.set_defaults(function=status)

    finalize_parser = subparsers.add_parser(
        "finalize", help="validate every task and build Parquet/DuckDB outputs"
    )
    finalize_parser.add_argument("--run-root", type=Path, required=True)
    finalize_parser.add_argument("--duckdb", default="duckdb")
    finalize_parser.add_argument("--finalization-source-commit", required=True)
    finalize_parser.add_argument(
        "--finalization-source-dirty", action="store_true"
    )
    finalize_parser.set_defaults(function=finalize)
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        arguments.function(arguments)
    except (EnrichmentError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

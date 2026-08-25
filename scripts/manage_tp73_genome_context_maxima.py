#!/usr/bin/env python3

"""Plan, inspect, and finalize autosomal TP73-context maxima production."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class GenomeContextError(RuntimeError):
    pass


AUTOSOMES = tuple(str(chrom) for chrom in range(1, 23))
SCIENTIFIC_SOURCE_FILES = (
    "scripts/build_sparse_context_maxima.py",
    "scripts/manage_tp73_genome_context_maxima.py",
    "scripts/run_tp73_genome_context_maxima_batch.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_string_list(values: list[str] | tuple[str, ...]) -> str:
    return "[" + ",".join(sql_string(value) for value in values) + "]"


def executable(value: str) -> str:
    resolved = shutil.which(value)
    if resolved is None and Path(value).is_file():
        resolved = str(Path(value).resolve())
    if resolved is None:
        raise GenomeContextError(f"executable is unavailable: {value}")
    return resolved


def run_json(command: list[str], *, cwd: Path | None = None) -> list[dict[str, Any]]:
    process = subprocess.run(
        command, cwd=cwd, text=True, capture_output=True, check=False
    )
    if process.returncode != 0:
        raise GenomeContextError(process.stderr.strip() or "command failed")
    try:
        value = json.loads(process.stdout or "[]")
    except json.JSONDecodeError as error:
        raise GenomeContextError("command returned invalid JSON") from error
    if not isinstance(value, list):
        raise GenomeContextError("JSON command result is not a row array")
    return value


def run_sql(duckdb: str, database: str | Path, sql: str,
            *, cwd: Path | None = None, readonly: bool = False) -> None:
    mode = ["-readonly"] if readonly and str(database) != ":memory:" else []
    process = subprocess.run(
        [duckdb, "-light-mode", *mode, "-batch", str(database), "-c", sql],
        cwd=cwd, text=True, capture_output=True, check=False,
    )
    if process.returncode != 0:
        raise GenomeContextError(process.stderr.strip() or "DuckDB command failed")


def immutable_write(path: Path, content: str) -> None:
    encoded = content.encode("utf-8")
    if path.exists():
        if path.read_bytes() != encoded:
            raise GenomeContextError(f"immutable file differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        temporary.write_bytes(encoded)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def install_immutable_file(staged: Path, destination: Path) -> None:
    if destination.exists():
        if sha256(destination) != sha256(staged):
            raise GenomeContextError(f"immutable file differs: {destination}")
        staged.unlink()
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    os.replace(staged, destination)


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise GenomeContextError(f"JSON file is missing: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise GenomeContextError(f"JSON value is not an object: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise GenomeContextError(f"TSV file is missing: {path}")
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        rows = list(reader)
    if reader.fieldnames is None or not rows:
        raise GenomeContextError(f"TSV is empty or lacks a header: {path}")
    return rows


def package_database(package: Path) -> tuple[Path, Path]:
    manifest_path = package / "manifest.json"
    manifest = load_json(manifest_path)
    database_value = manifest.get("database")
    if not database_value:
        raise GenomeContextError(f"package manifest lacks database: {manifest_path}")
    database = Path(str(database_value))
    if not database.is_absolute():
        database = package / database
    database = database.resolve()
    if not database.is_file():
        raise GenomeContextError(f"package database is missing: {database}")
    return manifest_path, database


def scan_threshold_source_sql(
    duckdb: str, database: Path, package: Path, run_id: str
) -> tuple[str, str]:
    rows = run_json([
        duckdb, "-light-mode", "-readonly", "-json", str(database), "-c", """
SELECT COUNT(DISTINCT column_name) FILTER (
         WHERE column_name IN (
           'motif_id', 'final_minimum_score', 'threshold_set_id',
           'density_limited'
         )
       )::BIGINT AS required_columns
FROM information_schema.columns
WHERE table_schema = 'main' AND table_name = 'scan_motif_threshold';
""",
    ], cwd=package)
    if len(rows) != 1:
        raise GenomeContextError("could not inspect scan threshold metadata")
    if int(rows[0]["required_columns"]) == 4:
        return """
SELECT motif_id::VARCHAR AS motif_id,
       final_minimum_score::DOUBLE AS final_minimum_score,
       threshold_set_id::VARCHAR AS threshold_set_id,
       density_limited::BOOLEAN AS density_limited
FROM scan_motif_threshold
""", "scan_motif_threshold"

    inventory_rows = run_json([
        duckdb, "-light-mode", "-readonly", "-json", str(database), "-c", """
SELECT COUNT(DISTINCT column_name) FILTER (
         WHERE column_name IN ('motif_id', 'minimum_score')
       )::BIGINT AS required_columns
FROM information_schema.columns
WHERE table_schema = 'main' AND table_name = 'scan_file_inventory';
""",
    ], cwd=package)
    if (len(inventory_rows) != 1
            or int(inventory_rows[0]["required_columns"]) != 2):
        raise GenomeContextError(
            "scan catalog has neither per-motif thresholds nor inventory floors"
        )
    legacy_set_id = f"legacy_inventory_{run_id}"
    return f"""
SELECT motif_id::VARCHAR AS motif_id,
       MIN(minimum_score)::DOUBLE AS final_minimum_score,
       {sql_string(legacy_set_id)}::VARCHAR AS threshold_set_id,
       NULL::BOOLEAN AS density_limited
FROM scan_file_inventory
GROUP BY motif_id
HAVING COUNT(DISTINCT minimum_score) = 1
""", "scan_file_inventory_unique_floor"


def git_identity(source: Path) -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        text=True, capture_output=True, check=False,
    )
    if commit.returncode != 0:
        raise GenomeContextError(commit.stderr.strip() or "cannot read Git commit")
    dirty = False
    for arguments in (
        ["diff", "--quiet", "--ignore-submodules", "--"],
        ["diff", "--cached", "--quiet", "--ignore-submodules", "--"],
    ):
        result = subprocess.run(
            ["git", "-C", str(source), *arguments],
            text=True, capture_output=True, check=False,
        )
        if result.returncode not in (0, 1):
            raise GenomeContextError(result.stderr.strip() or "cannot read Git state")
        dirty = dirty or result.returncode == 1
    return commit.stdout.strip(), dirty


def scientific_hashes(source: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for relative in SCIENTIFIC_SOURCE_FILES:
        path = source / relative
        if not path.is_file():
            raise GenomeContextError(f"scientific source is missing: {path}")
        result[relative] = sha256(path)
    return result


def verify_scientific_hashes(config: dict[str, Any]) -> None:
    source = Path(str(config["source"]))
    expected = config.get("scientific_source_file_sha256")
    if not isinstance(expected, dict) or not expected:
        raise GenomeContextError("run plan lacks scientific source hashes")
    for relative, digest in expected.items():
        path = source / relative
        if not path.is_file() or sha256(path) != digest:
            raise GenomeContextError(f"scientific source changed: {path}")


def task_directory(run_root: Path, row: dict[str, str]) -> Path:
    return run_root / "tasks" / (
        f"task-{int(row['task_index']):06d}-{row['motif_id']}"
    )


def planned_tasks(run_root: Path) -> list[dict[str, str]]:
    rows = read_tsv(run_root / "plan" / "calibration_tasks.tsv")
    indices = [int(row["task_index"]) for row in rows]
    if indices != list(range(len(rows))):
        raise GenomeContextError("task indices are not contiguous")
    if len({row["motif_id"] for row in rows}) != len(rows):
        raise GenomeContextError("task plan contains duplicate motifs")
    return rows


def validate_task_marker(directory: Path, row: dict[str, str],
                         verify_checksum: bool = False) -> dict[str, Any]:
    marker_path = directory / "complete.json"
    marker = load_json(marker_path)
    if (marker.get("schema_version") != 1
            or marker.get("task_index") != int(row["task_index"])
            or marker.get("motif_id") != row["motif_id"]
            or marker.get("feature_schema_version") != 3
            or marker.get("distance_band_schema_id") != (
                "exclusive_overlap_adjacent_0_5_gap_6_20_gap_21_50_"
                "gap_51_100_gap_101_150_v1"
            )
            or marker.get("analysis_partition") != "autosome"):
        raise GenomeContextError(f"task marker identity differs: {marker_path}")
    record = marker.get("files", {}).get("cofactor_maxima.parquet")
    if not isinstance(record, dict):
        raise GenomeContextError(f"task marker lacks maxima: {marker_path}")
    maxima = directory / "cofactor_maxima.parquet"
    if not maxima.is_file() or maxima.stat().st_size != int(record.get("bytes", -1)):
        raise GenomeContextError(f"task maxima size differs: {maxima}")
    digest = str(record.get("sha256", ""))
    if re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        raise GenomeContextError(f"task maxima digest is invalid: {marker_path}")
    if verify_checksum and sha256(maxima) != digest:
        raise GenomeContextError(f"task maxima checksum differs: {maxima}")
    validation = marker.get("validation")
    if (not isinstance(validation, dict)
            or int(validation.get("chromosomes", -1)) != 22
            or int(validation.get("duplicate_keys", -1)) != 0
            or int(validation.get("feature_schema_version", -1)) != 3
            or int(validation.get("distance_band_schema_values", -1)) != 1
            or int(validation.get("wrong_distance_band_schema_rows", -1)) != 0):
        raise GenomeContextError(f"task validation is incomplete: {marker_path}")
    return marker


def evidence_inputs(package: Path) -> tuple[Path, list[dict[str, Any]], dict[str, Any]]:
    manifest = load_json(package / "manifest.json")
    if (manifest.get("state") != "complete"
            or manifest.get("primary_inference_partition") != "autosome"):
        raise GenomeContextError("evidence package is not a complete autosome contract")
    aggregate = package / "tables" / "tp73_anchor_evidence_autosome.parquet"
    inventory_path = package / "chromosome_file_inventory.tsv"
    if not aggregate.is_file() or not inventory_path.is_file():
        raise GenomeContextError("evidence package lacks aggregate or chromosome inventory")
    rows = read_tsv(inventory_path)
    selected: list[dict[str, Any]] = []
    for row in rows:
        if row.get("analysis_partition") != "autosome":
            continue
        relative = Path(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise GenomeContextError(f"unsafe evidence path: {relative}")
        path = package / relative
        if not path.is_file() or path.stat().st_size != int(row["bytes"]):
            raise GenomeContextError(f"evidence chromosome file differs: {path}")
        selected.append({
            "sequence_order": int(row["sequence_order"]),
            "chrom": row["chrom"],
            "chrom_length": int(row["chrom_length"]),
            "path": str(path.resolve()),
            "bytes": int(row["bytes"]),
            "sha256": row["sha256"],
        })
    if sorted(row["chrom"] for row in selected) != sorted(AUTOSOMES):
        raise GenomeContextError("evidence package does not contain exactly autosomes 1-22")
    selected.sort(key=lambda row: row["sequence_order"])
    return aggregate.resolve(), selected, manifest


def build_applied_registry(
    duckdb: str,
    scan_database: Path,
    scan_package: Path,
    source_registry: Path,
    output: Path,
    threshold_set_id: str,
    target_motif: str,
    scan_threshold_sql: str,
) -> None:
    staged = output.with_name(f".{output.name}.tmp-{os.getpid()}")
    if staged.exists():
        staged.unlink()
    sql = f"""
COPY (
WITH registry AS (
  SELECT * FROM read_parquet({sql_string(source_registry)})
), scan_threshold AS (
  {scan_threshold_sql}
)
SELECT
  r.motif_id,
  r.motif_name,
  {sql_string(threshold_set_id)}::VARCHAR AS threshold_set_id,
  r.genome_id,
  r.motif_set_id,
  r.threshold_role,
  r.target_motif_id,
  r.score_mode,
  r.pseudocount,
  r.background_model_id,
  r.pseudocount_scheme,
  r.calibration_stratum_id,
  GREATEST(COALESCE(r.recommended_threshold, 0),
           s.final_minimum_score)::DOUBLE AS recommended_threshold,
  CASE
    WHEN r.recommended_threshold IS NULL AND s.final_minimum_score > 0
      THEN 'fallback_raised_to_scan_floor'
    WHEN r.recommended_threshold IS NULL
      THEN 'fallback_zero_no_source_recommendation'
    WHEN r.recommended_threshold < s.final_minimum_score
      THEN 'raised_to_scan_retention_floor'
    ELSE 'source_recommendation_reachable'
  END::VARCHAR AS calibration_status,
  r.context_min_interval_distance_bp,
  r.context_max_interval_distance_bp,
  r.context_relation_filter,
  true::BOOLEAN AS threshold_inclusive,
  'signed_interval_edge_distance'::VARCHAR AS context_distance_metric,
  s.final_minimum_score::DOUBLE AS source_minimum_score,
  r.threshold_set_id::VARCHAR AS source_threshold_set_id,
  r.recommended_threshold::DOUBLE AS source_recommended_threshold,
  r.calibration_status::VARCHAR AS source_calibration_status,
  s.threshold_set_id::VARCHAR AS scan_threshold_set_id,
  s.final_minimum_score::DOUBLE AS scan_retention_floor,
  s.density_limited::BOOLEAN AS scan_density_limited,
  (r.recommended_threshold IS NULL)::BOOLEAN AS used_fallback,
  (COALESCE(r.recommended_threshold, 0) < s.final_minimum_score)::BOOLEAN
      AS raised_to_scan_floor
FROM registry r
JOIN scan_threshold s USING (motif_id)
WHERE r.motif_id <> {sql_string(target_motif)}
ORDER BY r.motif_id
) TO {sql_string(staged)} (FORMAT PARQUET, COMPRESSION ZSTD);
"""
    run_sql(duckdb, scan_database, sql, cwd=scan_package, readonly=True)
    install_immutable_file(staged, output)


def prepare(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    source = arguments.source.expanduser().resolve()
    scan_package = arguments.scan_package.expanduser().resolve()
    evidence_package = arguments.evidence_package.expanduser().resolve()
    source_registry = arguments.threshold_registry.expanduser().resolve()
    runtime_prefix = arguments.runtime_prefix.expanduser().resolve()
    for path, label in (
        (source, "source"), (scan_package, "scan package"),
        (evidence_package, "evidence package"),
        (runtime_prefix, "runtime prefix"),
    ):
        if not path.is_dir():
            raise GenomeContextError(f"{label} is missing: {path}")
    if not source_registry.is_file():
        raise GenomeContextError(f"threshold registry is missing: {source_registry}")
    duckdb = executable(arguments.duckdb)
    scan_manifest, scan_database = package_database(scan_package)
    scan_manifest_value = load_json(scan_manifest)
    scan_threshold_sql, scan_threshold_source = scan_threshold_source_sql(
        duckdb, scan_database, scan_package,
        str(scan_manifest_value.get("run_id") or "unknown_scan"),
    )
    anchor, anchor_files, evidence_manifest = evidence_inputs(evidence_package)
    for name in ("plan", "logs", "tasks", "final"):
        (run_root / name).mkdir(parents=True, exist_ok=True)

    anchor_summary = run_json([
        duckdb, "-light-mode", "-json", ":memory:", "-c", f"""
SELECT count(*)::BIGINT AS anchors,
       count(DISTINCT chrom)::BIGINT AS chromosomes,
       count(*) - count(DISTINCT (chrom, anchor_start, anchor_end))::BIGINT
           AS duplicate_keys,
       min(anchor_score)::DOUBLE AS minimum_score
FROM read_parquet({sql_string(anchor)});
""",
    ])
    if len(anchor_summary) != 1:
        raise GenomeContextError("anchor validation returned no summary")
    anchor_values = anchor_summary[0]
    if (int(anchor_values["anchors"]) <= 0
            or int(anchor_values["chromosomes"]) != 22
            or int(anchor_values["duplicate_keys"]) != 0
            or float(anchor_values["minimum_score"]) < -1 - 1e-9):
        raise GenomeContextError(f"autosome anchor validation failed: {anchor_values}")

    applied_registry = run_root / "plan" / "applied_context_thresholds.parquet"
    build_applied_registry(
        duckdb, scan_database, scan_package, source_registry,
        applied_registry, arguments.applied_threshold_set_id,
        arguments.target_motif, scan_threshold_sql,
    )
    registry_summary = run_json([
        duckdb, "-light-mode", "-readonly", "-json", str(scan_database), "-c",
        f"""
WITH registry AS (
  SELECT * FROM read_parquet({sql_string(applied_registry)})
), source_registry AS (
  SELECT * FROM read_parquet({sql_string(source_registry)})
), scan_threshold AS (
  {scan_threshold_sql}
), expected AS (
  SELECT DISTINCT motif_id FROM scan_file_inventory
  WHERE motif_id <> {sql_string(arguments.target_motif)}
)
SELECT count(*)::BIGINT AS rows,
       count(DISTINCT motif_id)::BIGINT AS motifs,
       count(*) FILTER (WHERE motif_id = {sql_string(arguments.target_motif)})::BIGINT
           AS target_rows,
       count(*) FILTER (WHERE recommended_threshold IS NULL
                         OR recommended_threshold < source_minimum_score
                         OR context_max_interval_distance_bp > {arguments.flank}
                         OR threshold_inclusive IS DISTINCT FROM true)::BIGINT
           AS invalid_rows,
       count(*) FILTER (WHERE used_fallback)::BIGINT AS fallback_rows,
       count(*) FILTER (WHERE raised_to_scan_floor)::BIGINT AS raised_rows,
       count(*) FILTER (
         WHERE source_minimum_score > {arguments.maximum_source_score_floor}
       )::BIGINT AS source_floors_above_contract,
       count(DISTINCT threshold_set_id)::BIGINT AS threshold_sets,
       count(DISTINCT threshold_role)::BIGINT AS threshold_roles,
       count(DISTINCT calibration_stratum_id)::BIGINT AS strata
       ,(SELECT count(*) FROM expected)::BIGINT AS expected_scan_motifs
       ,(SELECT count(*) FROM scan_threshold
         WHERE motif_id = {sql_string(arguments.target_motif)})::BIGINT
           AS scan_target_rows
       ,(SELECT count(*) FROM expected e
         LEFT JOIN source_registry s USING (motif_id)
         WHERE s.motif_id IS NULL)::BIGINT AS missing_source_motifs
       ,(SELECT count(*) FROM source_registry s
         LEFT JOIN scan_threshold t USING (motif_id)
         WHERE s.motif_id <> {sql_string(arguments.target_motif)}
           AND t.motif_id IS NULL)::BIGINT AS extra_source_motifs
FROM registry;
""",
    ], cwd=scan_package)
    values = {key: int(value) for key, value in registry_summary[0].items()}
    if (values["rows"] <= 0 or values["rows"] != values["motifs"]
            or values["rows"] != values["expected_scan_motifs"]
            or values["target_rows"] != 0 or values["invalid_rows"] != 0
            or values["scan_target_rows"] > 1
            or values["missing_source_motifs"] != 0
            or values["extra_source_motifs"] != 0
            or values["source_floors_above_contract"] != 0
            or values["threshold_sets"] != 1
            or values["threshold_roles"] != 1 or values["strata"] != 1):
        raise GenomeContextError(f"applied threshold registry is invalid: {values}")

    identity = run_json([
        duckdb, "-light-mode", "-readonly", "-json", str(scan_database), "-c",
        f"""
SELECT min(threshold_role)::VARCHAR AS threshold_role,
       min(calibration_stratum_id)::VARCHAR AS calibration_stratum_id
FROM read_parquet({sql_string(applied_registry)});
""",
    ], cwd=scan_package)[0]

    task_temporary = run_root / "plan" / f".calibration_tasks.tmp-{os.getpid()}.tsv"
    run_sql(duckdb, scan_database, f"""
COPY (
  SELECT (row_number() OVER (ORDER BY r.motif_id) - 1)::BIGINT AS task_index,
         r.motif_id, r.motif_name, m.motif_length::BIGINT AS motif_length,
         r.source_minimum_score::DOUBLE AS scan_minimum_score,
         r.recommended_threshold::DOUBLE AS applied_context_threshold,
         r.source_recommended_threshold::DOUBLE AS source_recommended_threshold,
         r.calibration_status AS threshold_application_status
  FROM read_parquet({sql_string(applied_registry)}) r
  JOIN motif_metadata m USING (motif_id)
  ORDER BY r.motif_id
) TO {sql_string(task_temporary)}
  (FORMAT CSV, DELIMITER '\t', HEADER true, NULL 'NA');
""", cwd=scan_package, readonly=True)
    task_path = run_root / "plan" / "calibration_tasks.tsv"
    install_immutable_file(task_temporary, task_path)
    immutable_write(
        run_root / "plan" / "motif_tasks.tsv",
        task_path.read_text(encoding="utf-8"),
    )
    task_rows = planned_tasks(run_root)
    if len(task_rows) != values["motifs"]:
        raise GenomeContextError("task plan and applied registry differ")

    anchor_output = io.StringIO(newline="")
    anchor_writer = csv.DictWriter(
        anchor_output, fieldnames=tuple(anchor_files[0]), delimiter="\t",
        lineterminator="\n",
    )
    anchor_writer.writeheader()
    anchor_writer.writerows(anchor_files)
    anchor_plan = run_root / "plan" / "anchor_files.tsv"
    immutable_write(anchor_plan, anchor_output.getvalue())

    scan_temporary = run_root / "plan" / f".scan_files.tmp-{os.getpid()}.tsv"
    run_sql(duckdb, scan_database, f"""
COPY (
  SELECT i.motif_id, s.sequence_order::BIGINT AS sequence_order,
         CAST(i.chrom AS VARCHAR) AS chrom, i.strand,
         ('task_data/task_id=' || CAST(i.task_id AS VARCHAR) || '/' ||
          i.output_relative_path)::VARCHAR AS relative_path,
         i.bytes::BIGINT AS bytes, i.sha256::VARCHAR AS sha256,
         i.emitted_hits::BIGINT AS emitted_hits,
         i.minimum_score::DOUBLE AS minimum_score
  FROM scan_file_inventory i
  JOIN sequence_region s ON CAST(s.chrom AS VARCHAR) = CAST(i.chrom AS VARCHAR)
  JOIN read_parquet({sql_string(applied_registry)}) r USING (motif_id)
  WHERE CAST(i.chrom AS VARCHAR) IN {sql_string_list(AUTOSOMES)}
  ORDER BY i.motif_id, s.sequence_order, i.strand
) TO {sql_string(scan_temporary)}
  (FORMAT CSV, DELIMITER '\t', HEADER true);
""", cwd=scan_package, readonly=True)
    scan_plan = run_root / "plan" / "scan_files.tsv"
    install_immutable_file(scan_temporary, scan_plan)

    task_floor = {
        row["motif_id"]: float(row["scan_minimum_score"])
        for row in task_rows
    }
    observed: dict[tuple[str, str], set[str]] = {}
    scan_rows = 0
    with scan_plan.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            scan_rows += 1
            motif_id = row["motif_id"]
            chrom = row["chrom"]
            strand = row["strand"]
            relative = Path(row["relative_path"])
            if (motif_id not in task_floor or chrom not in AUTOSOMES
                    or strand not in {"+", "-"}
                    or relative.is_absolute() or ".." in relative.parts):
                raise GenomeContextError(f"invalid scan-plan row: {row}")
            if not math.isclose(
                float(row["minimum_score"]), task_floor[motif_id],
                rel_tol=1e-12, abs_tol=1e-12,
            ):
                raise GenomeContextError(f"scan floor differs for {motif_id}")
            key = (motif_id, chrom)
            strands = observed.setdefault(key, set())
            if strand in strands:
                raise GenomeContextError(f"duplicate scan file for {key} {strand}")
            strands.add(strand)
    expected_scan_rows = len(task_rows) * len(AUTOSOMES) * 2
    if (scan_rows != expected_scan_rows or len(observed) != len(task_rows) * 22
            or any(strands != {"+", "-"} for strands in observed.values())):
        raise GenomeContextError("scan plan is not motif x autosome x strand complete")

    commit, dirty = git_identity(source)
    batch_count = math.ceil(len(task_rows) / arguments.motifs_per_batch)
    config = {
        "schema_version": 1,
        "run_id": arguments.run_id,
        "analysis": "whole_autosome_tp73_context_maxima",
        "analysis_scope": "all_non_target_jaspar2026_motifs_vs_tp73_autosomes",
        "inference_status": (
            "autosomes 1-22; operating points calibrated on chromosome 1; "
            "whole-autosome association remains exploratory"
        ),
        "source": str(source),
        "source_commit": commit,
        "source_dirty": dirty,
        "scientific_source_file_sha256": scientific_hashes(source),
        "scan_package": str(scan_package),
        "scan_manifest_sha256": sha256(scan_manifest),
        "scan_database": str(scan_database),
        "scan_threshold_source": scan_threshold_source,
        "evidence_package": str(evidence_package),
        "evidence_manifest_sha256": sha256(evidence_package / "manifest.json"),
        "anchor_evidence": str(anchor),
        "anchor_evidence_sha256": sha256(anchor),
        "anchor_count": int(anchor_values["anchors"]),
        "anchor_minimum_score": float(anchor_values["minimum_score"]),
        "source_threshold_registry": str(source_registry),
        "source_threshold_registry_sha256": sha256(source_registry),
        "threshold_registry": str(applied_registry),
        "threshold_registry_sha256": sha256(applied_registry),
        "threshold_set_id": arguments.applied_threshold_set_id,
        "threshold_role": identity["threshold_role"],
        "calibration_stratum_id": identity["calibration_stratum_id"],
        "target_motif_id": arguments.target_motif,
        "analysis_partition": "autosome",
        "chromosomes": list(AUTOSOMES),
        "context_flank_bp": arguments.flank,
        "context_distance_metric": "signed_interval_edge_distance",
        "task_count": len(task_rows),
        "motifs_per_batch": arguments.motifs_per_batch,
        "batch_count": batch_count,
        "fallback_threshold_count": values["fallback_rows"],
        "raised_to_scan_floor_count": values["raised_rows"],
        "maximum_source_score_floor": arguments.maximum_source_score_floor,
        "runtime_prefix": str(runtime_prefix),
        "scratch_root": str(arguments.scratch_root),
        "threads": arguments.threads,
        "memory_limit": arguments.memory_limit,
        "max_temp_size": arguments.max_temp_size,
        "minimum_free_run_gb": arguments.minimum_free_run_gb,
        "minimum_free_scratch_gb": arguments.minimum_free_scratch_gb,
        "task_plan_sha256": sha256(task_path),
        "anchor_file_plan_sha256": sha256(anchor_plan),
        "scan_file_plan_sha256": sha256(scan_plan),
        "threshold_application_rule": (
            "max(coalesce(chr1_context_recommendation,0),scan_retention_floor)"
        ),
        "source_retention_contract": (
            "every_non_target_motif_scan_floor_at_or_below_"
            f"{arguments.maximum_source_score_floor:g}"
        ),
        "payload_discovery_policy": "exact_scan_and_anchor_inventory_paths_only",
    }
    immutable_write(
        run_root / "plan" / "run_config.json",
        json.dumps(config, indent=2, sort_keys=True) + "\n",
    )
    print(batch_count)


def status(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(run_root / "plan" / "run_config.json")
    rows = planned_tasks(run_root)
    complete = invalid = 0
    for row in rows:
        directory = task_directory(run_root, row)
        if not (directory / "complete.json").is_file():
            continue
        try:
            validate_task_marker(directory, row)
            complete += 1
        except (GenomeContextError, OSError, ValueError, json.JSONDecodeError):
            invalid += 1
    print("key\tvalue")
    print(f"planned_motifs\t{len(rows)}")
    print(f"planned_batches\t{config['batch_count']}")
    print(f"complete_motifs\t{complete}")
    print(f"pending_motifs\t{len(rows) - complete - invalid}")
    print(f"invalid_complete_markers\t{invalid}")


def output_inventory(root: Path) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path.name == "manifest.json":
            continue
        result.append({
            "relative_path": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    return result


def finalize(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(run_root / "plan" / "run_config.json")
    verify_scientific_hashes(config)
    rows = planned_tasks(run_root)
    if int(config.get("task_count", -1)) != len(rows):
        raise GenomeContextError("run configuration task count differs")
    records: list[dict[str, Any]] = []
    incomplete: list[str] = []
    for row in rows:
        directory = task_directory(run_root, row)
        if not (directory / "complete.json").is_file():
            incomplete.append(row["motif_id"])
            continue
        marker = validate_task_marker(
            directory, row, verify_checksum=arguments.verify_checksums
        )
        file_record = marker["files"]["cofactor_maxima.parquet"]
        records.append({
            "task_index": int(row["task_index"]),
            "motif_id": row["motif_id"],
            "motif_name": row["motif_name"],
            "motif_length": int(row["motif_length"]),
            "analysis_partition": "autosome",
            "chromosome_count": int(marker["validation"]["chromosomes"]),
            "anchor_count": int(marker["validation"]["anchors"]),
            "feature_schema_version": int(
                marker["validation"]["feature_schema_version"]
            ),
            "distance_band_schema_id": marker["distance_band_schema_id"],
            "scan_minimum_score": float(row["scan_minimum_score"]),
            "applied_context_threshold": float(row["applied_context_threshold"]),
            "relative_path": str(
                (directory / "cofactor_maxima.parquet").relative_to(run_root)
            ),
            "absolute_path": str((directory / "cofactor_maxima.parquet").resolve()),
            "bytes": int(file_record["bytes"]),
            "sha256": file_record["sha256"],
        })
    if incomplete:
        raise GenomeContextError(
            f"cannot finalize; {len(incomplete)} motifs incomplete, "
            f"first={incomplete[:5]}"
        )

    final = run_root / "final" / "context_maxima"
    if final.exists():
        manifest = load_json(final / "manifest.json")
        if (manifest.get("state") == "complete"
                and manifest.get("run_id") == config.get("run_id")):
            print(f"I: reusing finalized context maxima: {final}", file=sys.stderr)
            return
        raise GenomeContextError(f"final output has another identity: {final}")
    final.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".context-maxima.", dir=final.parent))
    try:
        inventory_tsv = staging / "context_maxima_file_inventory.tsv"
        with inventory_tsv.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=tuple(records[0]), delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(records)
        table_root = staging / "tables" / "jaspar2026"
        inventory_parquet = (
            table_root / "context_maxima_file_inventory" / "part-000000.parquet"
        )
        registry_parquet = (
            table_root / "motif_score_threshold" / "part-000000.parquet"
        )
        inventory_parquet.parent.mkdir(parents=True)
        registry_parquet.parent.mkdir(parents=True)
        shutil.copy2(Path(str(config["threshold_registry"])), registry_parquet)
        duckdb = executable(arguments.duckdb)
        run_sql(duckdb, ":memory:", f"""
COPY (
  SELECT * FROM read_csv_auto(
    {sql_string(inventory_tsv)}, delim='\t', header=true
  ) ORDER BY task_index
) TO {sql_string(inventory_parquet)} (FORMAT PARQUET, COMPRESSION ZSTD);
""")
        database = staging / "tp73_genome_context_maxima.duckdb"
        run_sql(duckdb, database, """
CREATE TABLE context_maxima_file_inventory AS
SELECT * FROM read_parquet(
  'tables/jaspar2026/context_maxima_file_inventory/part-000000.parquet'
);
CREATE TABLE motif_score_threshold AS
SELECT * FROM read_parquet(
  'tables/jaspar2026/motif_score_threshold/part-000000.parquet'
);
CREATE MACRO context_maxima_files(file_paths) AS TABLE
SELECT * FROM read_parquet(file_paths, hive_partitioning=false);
CREATE MACRO tp73_motif_threshold_count_files(file_paths) AS TABLE
SELECT * FROM read_parquet(file_paths, hive_partitioning=false);
""", cwd=staging)
        manifest = {
            "schema_version": 1,
            "state": "complete",
            "run_id": config["run_id"],
            "analysis": config["analysis"],
            "analysis_scope": config["analysis_scope"],
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "source_commit": config["source_commit"],
            "motifs": len(records),
            "anchors_per_motif": config["anchor_count"],
            "chromosomes": config["chromosomes"],
            "threshold_set_id": config["threshold_set_id"],
            "threshold_application_rule": config["threshold_application_rule"],
            "feature_schema_version": 3,
            "distance_band_schema_id": (
                "exclusive_overlap_adjacent_0_5_gap_6_20_gap_21_50_"
                "gap_51_100_gap_101_150_v1"
            ),
            "payload_files_copied": False,
            "payload_checksum_verification_at_finalization":
                arguments.verify_checksums,
            "database": database.name,
            "outputs": output_inventory(staging),
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, final)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    print(f"I: finalized autosomal context maxima: {final}", file=sys.stderr)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser(
        "prepare", help="write immutable motif, chromosome, and scan-file plans"
    )
    prepare_parser.add_argument("--run-root", type=Path, required=True)
    prepare_parser.add_argument("--scan-package", type=Path, required=True)
    prepare_parser.add_argument("--evidence-package", type=Path, required=True)
    prepare_parser.add_argument("--threshold-registry", type=Path, required=True)
    prepare_parser.add_argument("--runtime-prefix", type=Path, required=True)
    prepare_parser.add_argument("--source", type=Path, required=True)
    prepare_parser.add_argument("--duckdb", default="duckdb")
    prepare_parser.add_argument("--run-id", required=True)
    prepare_parser.add_argument("--applied-threshold-set-id", required=True)
    prepare_parser.add_argument("--target-motif", default="MA0861.2")
    prepare_parser.add_argument("--flank", type=int, default=150)
    prepare_parser.add_argument(
        "--maximum-source-score-floor", type=float, default=-1.0,
        help=(
            "reject any non-target scan retained above this score; protects "
            "negative-reference and threshold-sensitivity analyses (default: -1)"
        ),
    )
    prepare_parser.add_argument("--motifs-per-batch", type=int, default=3)
    prepare_parser.add_argument("--scratch-root", type=Path,
                                default=Path("/scratch") / os.environ.get("USER", "sm718"))
    prepare_parser.add_argument("--threads", type=int, default=4)
    prepare_parser.add_argument("--memory-limit", default="28GB")
    prepare_parser.add_argument("--max-temp-size", default="100GB")
    prepare_parser.add_argument("--minimum-free-run-gb", type=float, default=20)
    prepare_parser.add_argument("--minimum-free-scratch-gb", type=float, default=20)
    prepare_parser.set_defaults(function=prepare)

    status_parser = subparsers.add_parser(
        "status", help="report exact completed motif counts"
    )
    status_parser.add_argument("--run-root", type=Path, required=True)
    status_parser.set_defaults(function=status)

    finalize_parser = subparsers.add_parser(
        "finalize", help="validate tasks and publish the compact query catalog"
    )
    finalize_parser.add_argument("--run-root", type=Path, required=True)
    finalize_parser.add_argument("--duckdb", default="duckdb")
    finalize_parser.add_argument("--verify-checksums", action="store_true")
    finalize_parser.set_defaults(function=finalize)
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        for name in ("flank", "motifs_per_batch", "threads"):
            if hasattr(arguments, name) and getattr(arguments, name) <= 0:
                raise GenomeContextError(f"--{name.replace('_', '-')} must be positive")
        if (hasattr(arguments, "maximum_source_score_floor")
                and not math.isfinite(arguments.maximum_source_score_floor)):
            raise GenomeContextError("--maximum-source-score-floor must be finite")
        arguments.function(arguments)
        return 0
    except (GenomeContextError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Manage restart-safe isoform-by-distance TP73 cofactor enrichment."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from typing import Any


DEFAULT_CHROMOSOMES = tuple(str(value) for value in range(1, 23))
DISTANCE_BANDS = (
    "overlap",
    "adjacent_0_5",
    "gap_6_20",
    "gap_21_50",
    "gap_51_100",
    "gap_101_150",
)
FINAL_SCHEMA_VERSION = 3
PRESPECIFIED_COFACTOR_NAMES = (
    "POU2F2",
    "SP1",
    "PATZ1",
    "REST",
    "E2F1",
)
SCIENTIFIC_FILES = (
    "scripts/build_tp73_distance_cofactor_counts.py",
    "scripts/manage_tp73_distance_cofactor_enrichment.py",
)
CURRENT_PHASE = "initializing"
CURRENT_MOTIF = "none"
CURRENT_CHROM = "none"
STARTED = time.monotonic()


class DistanceEnrichmentError(RuntimeError):
    """Raised when a planned run or generated result violates its contract."""


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat()


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


def sql_string_list(values: tuple[str, ...]) -> str:
    return "(" + ",".join(sql_string(value) for value in values) + ")"


def parse_chromosomes(value: str) -> tuple[str, ...]:
    chromosomes = tuple(part.strip() for part in value.split(",") if part.strip())
    if not chromosomes:
        raise DistanceEnrichmentError("at least one chromosome is required")
    if len(set(chromosomes)) != len(chromosomes):
        raise DistanceEnrichmentError("chromosomes must not be repeated")
    if any(not chromosome.isdigit() or int(chromosome) <= 0
           for chromosome in chromosomes):
        raise DistanceEnrichmentError(
            "--chromosomes accepts comma-separated positive integers"
        )
    return chromosomes


def configured_chromosomes(config: dict[str, Any]) -> tuple[str, ...]:
    value = config.get("chromosomes")
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise DistanceEnrichmentError("run configuration has invalid chromosomes")
    return parse_chromosomes(",".join(value))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if reader.fieldnames is None:
        raise DistanceEnrichmentError(f"TSV has no header: {path}")
    return rows


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise DistanceEnrichmentError(f"JSON object expected: {path}")
    return value


def write_json_new(path: Path, value: dict[str, Any]) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def run_duckdb(executable: str, database: str | Path, sql: str) -> None:
    process = subprocess.run(
        [executable, "-batch", str(database)],
        input=sql,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise DistanceEnrichmentError(
            process.stderr.strip() or "DuckDB command failed"
        )


def run_json(executable: str, database: str | Path, sql: str) -> list[dict[str, Any]]:
    process = subprocess.run(
        [executable, "-json", str(database), "-c", sql],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise DistanceEnrichmentError(
            process.stderr.strip() or "DuckDB JSON query failed"
        )
    value = json.loads(process.stdout or "[]")
    if not isinstance(value, list):
        raise DistanceEnrichmentError("DuckDB JSON result is not an array")
    return value


def source_identity(source: Path, source_commit: str) -> dict[str, str]:
    if not re.fullmatch(r"[0-9a-f]{40}", source_commit):
        raise DistanceEnrichmentError("source commit must be 40 lowercase hex digits")
    git_process = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        text=True,
        capture_output=True,
        check=False,
    )
    if (git_process.returncode != 0
            or git_process.stdout.strip() != source_commit):
        raise DistanceEnrichmentError(
            "source commit does not match the source checkout HEAD"
        )
    result: dict[str, str] = {}
    for relative in SCIENTIFIC_FILES:
        path = source / relative
        if not path.is_file():
            raise DistanceEnrichmentError(f"scientific source is missing: {path}")
        result[relative] = sha256(path)
    return result


def runtime_source_identity(
    source: Path, commit: str, dirty: bool
) -> dict[str, Any]:
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise DistanceEnrichmentError("finalizer source commit must be 40 hex digits")
    hashes: dict[str, str] = {}
    for relative in SCIENTIFIC_FILES:
        path = source / relative
        if not path.is_file():
            raise DistanceEnrichmentError(f"finalizer source is missing: {path}")
        hashes[relative] = sha256(path)
    return {
        "source": str(source),
        "commit": commit,
        "dirty": dirty,
        "scientific_source_sha256": hashes,
    }


def verify_source(config: dict[str, Any]) -> Path:
    source = Path(config["source"])
    for relative, expected in config["scientific_source_sha256"].items():
        path = source / relative
        if not path.is_file() or sha256(path) != expected:
            raise DistanceEnrichmentError(f"scientific source changed: {path}")
    return source


def verify_plan(run_root: Path, config: dict[str, Any]) -> None:
    for name, key in (
        ("tasks.tsv", "tasks_sha256"),
        ("scan_files.tsv", "scan_files_sha256"),
    ):
        path = run_root / "plan" / name
        if not path.is_file() or sha256(path) != config[key]:
            raise DistanceEnrichmentError(f"immutable run plan changed: {path}")


def config_path(run_root: Path) -> Path:
    return run_root / "plan" / "run_config.json"


def task_directory(run_root: Path, row: dict[str, str]) -> Path:
    return run_root / "tasks" / (
        f"task-{int(row['task_index']):06d}-{row['motif_id']}"
    )


def chromosome_directory(
    run_root: Path, row: dict[str, str], chrom: str
) -> Path:
    return run_root / "checkpoints" / (
        f"task-{int(row['task_index']):06d}-{row['motif_id']}"
    ) / f"chrom-{chrom}"


def prepare(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    if run_root.exists():
        raise DistanceEnrichmentError(f"run root already exists: {run_root}")
    run_root.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(
        prefix=f".{run_root.name}-prepare-", dir=run_root.parent
    ))
    try:
        prepare_plan(arguments, staging)
        os.replace(staging, run_root)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def prepare_plan(arguments: argparse.Namespace, run_root: Path) -> None:
    source = arguments.source.expanduser().resolve()
    scan_package = arguments.scan_package.expanduser().resolve()
    anchors = arguments.anchor_evidence.expanduser().resolve()
    thresholds = arguments.thresholds.expanduser().resolve()
    catalog = arguments.jaspar_catalog.expanduser().resolve()
    chromosomes = parse_chromosomes(arguments.chromosomes)
    for path in (anchors, thresholds, catalog / "manifest.json"):
        if not path.is_file():
            raise DistanceEnrichmentError(f"input is missing: {path}")
    if shutil.which(arguments.duckdb) is None:
        raise DistanceEnrichmentError(f"DuckDB not found: {arguments.duckdb}")

    scan_inventory = (
        scan_package / "tables" / "jaspar2026" / "scan_file_inventory" /
        "part-000000.parquet"
    )
    motif_metadata = (
        scan_package / "tables" / "jaspar2026" / "motif_metadata" /
        "part-000000.parquet"
    )
    for path in (scan_package / "manifest.json", scan_inventory, motif_metadata,
                 catalog / "jaspar_metadata.duckdb"):
        if not path.is_file():
            raise DistanceEnrichmentError(f"package input is missing: {path}")

    scientific_hashes = source_identity(source, arguments.source_commit)
    for name in ("plan", "input", "tasks", "checkpoints", "logs", "final"):
        (run_root / name).mkdir()
    task_file = run_root / "plan" / "tasks.tsv"
    scan_file = run_root / "plan" / "scan_files.tsv"
    chrom_values = ",".join(sql_string(value) for value in chromosomes)
    expected_files = 2 * len(chromosomes)
    tax_filter = (
        "TRUE" if arguments.tax_group == "all"
        else f"jm.tax_group = {sql_string(arguments.tax_group)}"
    )
    sql = f"""
ATTACH {sql_string(catalog / 'jaspar_metadata.duckdb')} AS jaspar (READ_ONLY);
CREATE TEMP TABLE selected AS
WITH inventory AS (
  SELECT motif_id,
         count(*)::BIGINT AS files,
         count(DISTINCT chrom)::BIGINT AS chromosomes,
         count(DISTINCT strand)::BIGINT AS strands,
         min(TRY_CAST(minimum_score AS DOUBLE)) AS minimum_score_min,
         max(TRY_CAST(minimum_score AS DOUBLE)) AS minimum_score_max
  FROM read_parquet({sql_string(scan_inventory)})
  WHERE CAST(chrom AS VARCHAR) IN ({chrom_values})
    AND motif_id <> {sql_string(arguments.anchor_motif)}
  GROUP BY motif_id
), threshold AS (
  SELECT *
  FROM read_parquet({sql_string(thresholds)})
  WHERE threshold_set_id = {sql_string(arguments.threshold_set_id)}
)
SELECT
  (row_number() OVER (ORDER BY mm.motif_id) - 1)::BIGINT AS task_index,
  mm.motif_id::VARCHAR AS motif_id,
  mm.motif_name::VARCHAR AS motif_name,
  mm.motif_length::BIGINT AS motif_length,
  t.recommended_threshold::DOUBLE AS positive_threshold,
  t.source_minimum_score::DOUBLE AS threshold_source_floor,
  t.calibration_status::VARCHAR AS threshold_status,
  jm.tax_group::VARCHAR AS tax_group,
  jm.includes_homo_sapiens::BOOLEAN AS includes_homo_sapiens,
  jm.species_count::BIGINT AS source_species_count,
  i.files, i.chromosomes, i.strands,
  i.minimum_score_min, i.minimum_score_max
FROM read_parquet({sql_string(motif_metadata)}) mm
JOIN inventory i USING (motif_id)
JOIN threshold t USING (motif_id)
JOIN jaspar.jaspar_matrix jm ON jm.matrix_id = mm.motif_id
WHERE {tax_filter}
  AND mm.motif_id <> {sql_string(arguments.anchor_motif)}
ORDER BY mm.motif_id;

SELECT CASE WHEN EXISTS (
  SELECT 1 FROM selected
  WHERE files <> {expected_files}
     OR chromosomes <> {len(chromosomes)} OR strands <> 2
     OR minimum_score_min <> -1 OR minimum_score_max <> -1
     OR threshold_source_floor <> -1
     OR positive_threshold < -1 OR NOT isfinite(positive_threshold)
) THEN error('selected motif has incomplete files or incompatible score floors') END;

COPY selected TO {sql_string(task_file)}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
COPY (
  SELECT s.task_index, s.motif_id, i.task_id::VARCHAR AS scan_task_id,
         CAST(i.chrom AS VARCHAR) AS chrom,
         i.strand::VARCHAR AS strand,
         {sql_string(str(scan_package) + '/task_data/task_id=')} ||
             CAST(i.task_id AS VARCHAR) || '/' || i.output_relative_path
             AS absolute_path,
         i.bytes::BIGINT AS bytes,
         i.sha256::VARCHAR AS sha256,
         i.emitted_hits::BIGINT AS emitted_hits,
         TRY_CAST(i.minimum_score AS DOUBLE) AS minimum_score
  FROM selected s
  JOIN read_parquet({sql_string(scan_inventory)}) i USING (motif_id)
  WHERE CAST(i.chrom AS VARCHAR) IN ({chrom_values})
  ORDER BY s.task_index, CAST(i.chrom AS INTEGER), i.strand
) TO {sql_string(scan_file)}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
"""
    run_duckdb(arguments.duckdb, ":memory:", sql)
    tasks = read_tsv(task_file)
    scans = read_tsv(scan_file)
    if not tasks or len(scans) != expected_files * len(tasks):
        raise DistanceEnrichmentError("prepared task/scan cardinality is wrong")
    expected_scan_keys = {
        (task["motif_id"], chrom, strand)
        for task in tasks for chrom in chromosomes for strand in ("+", "-")
    }
    observed_scan_keys = {
        (row["motif_id"], row["chrom"], row["strand"]) for row in scans
    }
    if len(scans) != len(observed_scan_keys) or observed_scan_keys != expected_scan_keys:
        raise DistanceEnrichmentError("prepared scan keys are incomplete or repeated")
    for scan in scans:
        path = Path(scan["absolute_path"])
        if not path.is_file():
            raise DistanceEnrichmentError(f"planned scan payload is missing: {path}")
        if path.stat().st_size != int(scan["bytes"]):
            raise DistanceEnrichmentError(f"planned scan payload size differs: {path}")

    scan_manifest = scan_package / "manifest.json"
    catalog_manifest = catalog / "manifest.json"
    config = {
        "schema_version": 1,
        "run_id": arguments.run_id,
        "analysis": "tp73_cutandrun_cofactor_enrichment_by_isoform_and_distance",
        "tax_group": arguments.tax_group,
        "taxonomic_scope_rule": (
            "all JASPAR taxonomic groups" if arguments.tax_group == "all"
            else f"JASPAR tax_group={arguments.tax_group}"
        ),
        "source_species_semantics": (
            "organisms contributing to JASPAR matrix construction; "
            "not a binding-species restriction"
        ),
        "chromosomes": list(chromosomes),
        "anchor_motif": arguments.anchor_motif,
        "source_score_floor": -1.0,
        "negative_reference_rule": "no score>=-1 locus in exclusive band",
        "positive_rule": "best score in exclusive band >= motif threshold",
        "intermediate_rule": "-1 <= best score < motif threshold",
        "distance_bands": list(DISTANCE_BANDS),
        "distance_metric": "signed_interval_edge_distance",
        "context_flank_bp": arguments.flank,
        "block_size_bp": arguments.block_size,
        "tp73_score_breaks": [-5, -1, 0, 1, 2, 5, 10, 15, "+Inf"],
        "series": ["saos2", "skmel29_2"],
        "isoforms": ["TA", "DN"],
        "excluded_series": ["skmel29_1"],
        "task_count": len(tasks),
        "human_source_task_count": sum(
            row["includes_homo_sapiens"].lower() == "true" for row in tasks
        ),
        "source_species_unspecified_task_count": sum(
            int(row["source_species_count"]) == 0 for row in tasks
        ),
        "source": str(source),
        "source_commit": arguments.source_commit,
        "scientific_source_sha256": scientific_hashes,
        "scan_package": str(scan_package),
        "scan_manifest_sha256": sha256(scan_manifest),
        "anchor_evidence": str(anchors),
        "anchor_evidence_sha256": sha256(anchors),
        "thresholds": str(thresholds),
        "thresholds_sha256": sha256(thresholds),
        "threshold_set_id": arguments.threshold_set_id,
        "jaspar_catalog": str(catalog),
        "jaspar_catalog_manifest_sha256": sha256(catalog_manifest),
        "jaspar_catalog_database_sha256": sha256(
            catalog / "jaspar_metadata.duckdb"
        ),
        "tasks_sha256": sha256(task_file),
        "scan_files_sha256": sha256(scan_file),
        "created_at_utc": utc_now(),
    }
    write_json_new(config_path(run_root), config)
    print(len(tasks))


def prepare_anchors(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(config_path(run_root))
    verify_source(config)
    verify_plan(run_root, config)
    chromosomes = configured_chromosomes(config)
    source = Path(config["anchor_evidence"])
    if sha256(source) != config["anchor_evidence_sha256"]:
        raise DistanceEnrichmentError("anchor evidence checksum changed")
    final = run_root / "input" / "anchors"
    if final.exists():
        marker = load_json(final / "complete.json")
        if marker.get("source_sha256") != config["anchor_evidence_sha256"]:
            raise DistanceEnrichmentError("existing anchor split has another source")
        inventory = final / "anchor_files.tsv"
        if (not inventory.is_file()
                or sha256(inventory) != marker.get("inventory_sha256")):
            raise DistanceEnrichmentError("existing anchor inventory is invalid")
        print(f"I: Reusing anchor partitions: {final}", file=sys.stderr)
        return
    staging = Path(tempfile.mkdtemp(prefix=".anchors-", dir=final.parent))
    try:
        rows: list[dict[str, Any]] = []
        for chrom in chromosomes:
            output = staging / f"chrom-{chrom}.parquet"
            sql = f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;
COPY (
  SELECT chrom, anchor_start, anchor_end, anchor_score,
         supported_tp73_saos2_TA,
         supported_negative_control_saos2_TA,
         supported_tp73_saos2_DN,
         supported_negative_control_saos2_DN,
         supported_tp73_skmel29_2_TA,
         supported_negative_control_skmel29_2_TA,
         supported_tp73_skmel29_2_DN,
         supported_negative_control_skmel29_2_DN
  FROM read_parquet({sql_string(source)})
  WHERE CAST(chrom AS VARCHAR) = {sql_string(chrom)}
  ORDER BY anchor_start, anchor_end
) TO {sql_string(output)}
  (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 131072);
"""
            run_duckdb(arguments.duckdb, ":memory:", sql)
            validation = run_json(
                arguments.duckdb,
                ":memory:",
                f"SELECT count(*)::BIGINT AS rows, "
                f"count(*)-count(DISTINCT (anchor_start,anchor_end)) AS duplicates "
                f"FROM read_parquet({sql_string(output)});",
            )
            if len(validation) != 1 or int(validation[0]["rows"]) == 0 \
                    or int(validation[0]["duplicates"]) != 0:
                raise DistanceEnrichmentError(f"invalid anchor partition: {chrom}")
            rows.append({
                "chrom": chrom,
                "relative_path": output.name,
                "rows": int(validation[0]["rows"]),
                "bytes": output.stat().st_size,
                "sha256": sha256(output),
            })
        inventory = staging / "anchor_files.tsv"
        write_tsv(
            inventory, rows, ["chrom", "relative_path", "rows", "bytes", "sha256"]
        )
        write_json_new(staging / "complete.json", {
            "schema_version": 1,
            "source_sha256": config["anchor_evidence_sha256"],
            "chromosomes": len(rows),
            "anchors": sum(row["rows"] for row in rows),
            "inventory_sha256": sha256(inventory),
            "completed_at_utc": utc_now(),
        })
        os.replace(staging, final)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    print(f"I: Prepared {len(chromosomes)} anchor partitions: {final}", file=sys.stderr)


def progress(_signal_number: int, _frame: object) -> None:
    print(
        "I: progress signal=SIGUSR1 "
        f"phase={CURRENT_PHASE} motif={CURRENT_MOTIF} chrom={CURRENT_CHROM} "
        f"elapsed_seconds={time.monotonic() - STARTED:.1f}",
        file=sys.stderr,
        flush=True,
    )


def set_phase(phase: str, motif: str | None = None, chrom: str | None = None) -> None:
    global CURRENT_PHASE, CURRENT_MOTIF, CURRENT_CHROM
    CURRENT_PHASE = phase
    if motif is not None:
        CURRENT_MOTIF = motif
    if chrom is not None:
        CURRENT_CHROM = chrom
    print(
        f"I: phase={phase} motif={CURRENT_MOTIF} chrom={CURRENT_CHROM}",
        file=sys.stderr,
        flush=True,
    )


def validate_existing_task(
    directory: Path,
    row: dict[str, str],
    config: dict[str, Any] | None = None,
    verify_checksum: bool = True,
) -> bool:
    if not directory.exists():
        return False
    marker = load_json(directory / "complete.json")
    if marker.get("motif_id") != row["motif_id"] \
            or marker.get("task_index") != int(row["task_index"]):
        raise DistanceEnrichmentError(f"existing task identity differs: {directory}")
    if config is not None and (
        marker.get("run_id") != config["run_id"]
        or marker.get("source_commit") != config["source_commit"]
    ):
        raise DistanceEnrichmentError(f"existing task run identity differs: {directory}")
    for name in ("block_components.parquet", "class_counts.parquet"):
        path = directory / name
        record = marker.get("files", {}).get(name, {})
        if (not path.is_file() or path.stat().st_size != record.get("bytes")
                or (verify_checksum and sha256(path) != record.get("sha256"))):
            raise DistanceEnrichmentError(f"existing task file is invalid: {path}")
    return True


def validate_existing_chromosome(
    directory: Path,
    row: dict[str, str],
    chrom: str,
    config: dict[str, Any],
    verify_checksum: bool = True,
) -> bool:
    if not directory.exists():
        return False
    marker = load_json(directory / "complete.json")
    if (marker.get("run_id") != config["run_id"]
            or marker.get("source_commit") != config["source_commit"]
            or marker.get("motif_id") != row["motif_id"]
            or marker.get("task_index") != int(row["task_index"])
            or marker.get("chrom") != chrom
            or marker.get("positive_threshold") != float(row["positive_threshold"])):
        raise DistanceEnrichmentError(
            f"existing chromosome checkpoint identity differs: {directory}"
        )
    for name in ("block_components.parquet", "class_counts.parquet"):
        path = directory / name
        record = marker.get("files", {}).get(name, {})
        if (not path.is_file() or path.stat().st_size != record.get("bytes")
                or (verify_checksum and sha256(path) != record.get("sha256"))):
            raise DistanceEnrichmentError(
                f"existing chromosome checkpoint is invalid: {path}"
            )
    return True


def run_task(arguments: argparse.Namespace) -> None:
    global STARTED
    STARTED = time.monotonic()
    signal.signal(signal.SIGUSR1, progress)
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(config_path(run_root))
    source = verify_source(config)
    verify_plan(run_root, config)
    chromosomes = configured_chromosomes(config)
    tasks = read_tsv(run_root / "plan" / "tasks.tsv")
    task_index = arguments.task_index
    if task_index is None:
        task_index = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0")) + arguments.task_offset
    if task_index < 0 or task_index >= len(tasks):
        raise DistanceEnrichmentError(f"task index is out of range: {task_index}")
    row = tasks[task_index]
    final = task_directory(run_root, row)
    if validate_existing_task(final, row, config):
        print(f"I: Reusing completed motif task: {row['motif_id']}", file=sys.stderr)
        return

    anchor_root = run_root / "input" / "anchors"
    anchor_marker = load_json(anchor_root / "complete.json")
    if anchor_marker.get("source_sha256") != config["anchor_evidence_sha256"]:
        raise DistanceEnrichmentError("anchor split identity differs")
    anchor_inventory = anchor_root / "anchor_files.tsv"
    if sha256(anchor_inventory) != anchor_marker.get("inventory_sha256"):
        raise DistanceEnrichmentError("anchor split inventory changed")
    anchor_rows = {row["chrom"]: row for row in read_tsv(
        anchor_root / "anchor_files.tsv"
    )}
    scan_rows = [scan for scan in read_tsv(run_root / "plan" / "scan_files.tsv")
                 if scan["motif_id"] == row["motif_id"]]
    expected_scan_files = 2 * len(chromosomes)
    if len(scan_rows) != expected_scan_files:
        raise DistanceEnrichmentError(
            f"motif scan plan does not contain {expected_scan_files} files"
        )
    scan_by_key = {(scan["chrom"], scan["strand"]): scan for scan in scan_rows}

    scratch_parent = arguments.scratch.expanduser().resolve()
    scratch_parent.mkdir(parents=True, exist_ok=True)
    required = max(
        int(anchor_rows[chrom]["bytes"])
        + sum(int(scan_by_key[(chrom, strand)]["bytes"])
              for strand in ("+", "-"))
        for chrom in chromosomes
    )
    free = shutil.disk_usage(scratch_parent).free
    if free < required + 10 * 1024**3:
        raise DistanceEnrichmentError(
            f"scratch space is tight: need {(required + 10 * 1024**3)/1024**3:.1f} GiB, "
            f"have {free/1024**3:.1f} GiB"
        )
    scratch = Path(tempfile.mkdtemp(
        prefix=f"tp73-distance-{row['motif_id']}-", dir=scratch_parent
    ))
    try:
        block_parts: list[Path] = []
        class_parts: list[Path] = []
        kernel = source / "scripts" / "build_tp73_distance_cofactor_counts.py"
        for chrom in chromosomes:
            checkpoint = chromosome_directory(run_root, row, chrom)
            if validate_existing_chromosome(checkpoint, row, chrom, config):
                print(
                    f"I: Reusing chromosome checkpoint motif={row['motif_id']} "
                    f"chrom={chrom}",
                    file=sys.stderr,
                    flush=True,
                )
                block_parts.append(checkpoint / "block_components.parquet")
                class_parts.append(checkpoint / "class_counts.parquet")
                continue

            set_phase("staging", row["motif_id"], chrom)
            anchor = anchor_rows[chrom]
            source_anchor = anchor_root / anchor["relative_path"]
            target_anchor = scratch / f"anchors-chrom-{chrom}.parquet"
            shutil.copy2(source_anchor, target_anchor)
            if (target_anchor.stat().st_size != int(anchor["bytes"])
                    or sha256(target_anchor) != anchor["sha256"]):
                raise DistanceEnrichmentError("staged anchor identity differs")
            staged_scans: dict[str, Path] = {}
            for strand, label in (("+", "plus"), ("-", "minus")):
                scan = scan_by_key[(chrom, strand)]
                source_scan = Path(scan["absolute_path"])
                target_scan = scratch / f"hits-chrom-{chrom}-{label}.parquet"
                shutil.copy2(source_scan, target_scan)
                if (target_scan.stat().st_size != int(scan["bytes"])
                        or sha256(target_scan) != scan["sha256"]):
                    raise DistanceEnrichmentError("staged scan identity differs")
                staged_scans[strand] = target_scan

            set_phase("distance_counts", row["motif_id"], chrom)
            block = scratch / f"block-chrom-{chrom}.parquet"
            classes = scratch / f"class-chrom-{chrom}.parquet"
            command = [
                sys.executable, str(kernel),
                "--anchors", str(target_anchor),
                "--plus-hits", str(staged_scans["+"]),
                "--minus-hits", str(staged_scans["-"]),
                "--motif-id", row["motif_id"],
                "--motif-name", row["motif_name"],
                "--chrom", chrom,
                "--positive-threshold", row["positive_threshold"],
                "--block-output", str(block),
                "--class-output", str(classes),
                "--block-size", str(config["block_size_bp"]),
                "--flank", str(config["context_flank_bp"]),
                "--threads", str(arguments.threads),
                "--memory-limit", arguments.memory_limit,
                "--max-temp-size", arguments.max_temp_size,
                "--temp-directory", str(scratch / "spill"),
                "--duckdb", arguments.duckdb,
            ]
            process = subprocess.run(command, check=False)
            if process.returncode != 0:
                raise DistanceEnrichmentError(
                    f"distance kernel failed for {row['motif_id']} chr{chrom}"
                )
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            checkpoint_staging = Path(tempfile.mkdtemp(
                prefix=f".chrom-{chrom}-", dir=checkpoint.parent
            ))
            try:
                checkpoint_files: dict[str, dict[str, Any]] = {}
                for name, source_file in (
                    ("block_components.parquet", block),
                    ("class_counts.parquet", classes),
                ):
                    destination = checkpoint_staging / name
                    shutil.copy2(source_file, destination)
                    checkpoint_files[name] = {
                        "bytes": destination.stat().st_size,
                        "sha256": sha256(destination),
                    }
                write_json_new(checkpoint_staging / "complete.json", {
                    "schema_version": 1,
                    "run_id": config["run_id"],
                    "source_commit": config["source_commit"],
                    "task_index": int(row["task_index"]),
                    "motif_id": row["motif_id"],
                    "chrom": chrom,
                    "positive_threshold": float(row["positive_threshold"]),
                    "files": checkpoint_files,
                    "completed_at_utc": utc_now(),
                })
                os.replace(checkpoint_staging, checkpoint)
            finally:
                if checkpoint_staging.exists():
                    shutil.rmtree(checkpoint_staging)
            block_parts.append(checkpoint / "block_components.parquet")
            class_parts.append(checkpoint / "class_counts.parquet")
            for path in (target_anchor, *staged_scans.values(), block, classes):
                path.unlink()

        set_phase("consolidating", row["motif_id"], "all")
        block_output = scratch / "block_components.parquet"
        class_output = scratch / "class_counts.parquet"
        run_duckdb(arguments.duckdb, ":memory:", f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;
COPY (SELECT * FROM read_parquet({sql_path_list(block_parts)})
      ORDER BY CAST(chrom AS INTEGER), block_index, series_id, isoform,
               distance_band_order)
TO {sql_string(block_output)} (FORMAT PARQUET, COMPRESSION ZSTD);
COPY (SELECT * FROM read_parquet({sql_path_list(class_parts)})
      ORDER BY CAST(chrom AS INTEGER), distance_band_order)
TO {sql_string(class_output)} (FORMAT PARQUET, COMPRESSION ZSTD);
""")
        validation = run_json(arguments.duckdb, ":memory:", f"""
SELECT
  (SELECT count(*) FROM read_parquet({sql_string(class_output)}))::BIGINT
    AS class_rows,
  (SELECT count(DISTINCT chrom) FROM read_parquet({sql_string(class_output)}))
    ::BIGINT AS class_chromosomes,
  (SELECT count(DISTINCT distance_band)
   FROM read_parquet({sql_string(class_output)}))::BIGINT AS class_bands,
  (SELECT count(*) FROM read_parquet({sql_string(block_output)}))::BIGINT
    AS block_rows,
  (SELECT count(*) FROM read_parquet({sql_string(class_output)})
   WHERE anchors_total <> anchors_positive + anchors_intermediate
                          + anchors_negative)::BIGINT AS invalid_classes;
""")
        expected = validation[0]
        if (int(expected["class_rows"]) != len(chromosomes) * len(DISTANCE_BANDS)
                or int(expected["class_chromosomes"]) != len(chromosomes)
                or int(expected["class_bands"]) != 6
                or int(expected["block_rows"]) == 0
                or int(expected["invalid_classes"]) != 0):
            raise DistanceEnrichmentError(f"task validation failed: {expected}")

        attempt = run_root / "tasks" / (
            f".attempt-{row['motif_id']}-job-"
            f"{os.environ.get('SLURM_JOB_ID', 'manual')}-pid-{os.getpid()}"
        )
        attempt.mkdir()
        files: dict[str, dict[str, Any]] = {}
        for name, source_file in (
            ("block_components.parquet", block_output),
            ("class_counts.parquet", class_output),
        ):
            destination = attempt / name
            shutil.copy2(source_file, destination)
            files[name] = {
                "bytes": destination.stat().st_size,
                "sha256": sha256(destination),
            }
        write_json_new(attempt / "complete.json", {
            "schema_version": 1,
            "task_index": int(row["task_index"]),
            "motif_id": row["motif_id"],
            "positive_threshold": float(row["positive_threshold"]),
            "run_id": config["run_id"],
            "source_commit": config["source_commit"],
            "validation": expected,
            "files": files,
            "completed_at_utc": utc_now(),
        })
        os.replace(attempt, final)
    finally:
        if scratch.exists():
            shutil.rmtree(scratch)
    set_phase("complete", row["motif_id"], "all")


def status(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(config_path(run_root))
    verify_plan(run_root, config)
    tasks = read_tsv(run_root / "plan" / "tasks.tsv")
    complete = sum(validate_existing_task(
                       task_directory(run_root, row), row, config,
                       verify_checksum=False,
                   )
                   for row in tasks)
    chromosomes = configured_chromosomes(config)
    chromosome_complete = sum(
        validate_existing_chromosome(
            chromosome_directory(run_root, row, chrom), row, chrom, config,
            verify_checksum=False,
        )
        for row in tasks for chrom in chromosomes
    )
    print(f"planned\t{len(tasks)}")
    print(f"complete\t{complete}")
    print(f"pending\t{len(tasks) - complete}")
    print(f"chromosome_checkpoints_complete\t{chromosome_complete}")
    print(f"chromosome_checkpoints_planned\t{len(tasks) * len(chromosomes)}")


def bh_adjust(rows: list[dict[str, Any]], family_size: int) -> None:
    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in rows:
        if row["evaluation_status"] == "ok" and math.isfinite(row["p_value"]):
            groups.setdefault((row["isoform"], row["distance_band"]), []).append(row)
    for group in groups.values():
        ordered = sorted(group, key=lambda row: row["p_value"])
        next_value = 1.0
        for rank in range(len(ordered), 0, -1):
            row = ordered[rank - 1]
            next_value = min(next_value, row["p_value"] * family_size / rank)
            row["q_value_bh_tax_group"] = min(1.0, next_value)


def bh_adjust_isoform_contrast(
    rows: list[dict[str, Any]], family_size: int
) -> None:
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        if (row["evaluation_status"] == "ok"
                and math.isfinite(row["p_value"])):
            groups.setdefault(row["distance_band"], []).append(row)
    for group in groups.values():
        ordered = sorted(group, key=lambda row: row["p_value"])
        next_value = 1.0
        for rank in range(len(ordered), 0, -1):
            row = ordered[rank - 1]
            next_value = min(next_value, row["p_value"] * family_size / rank)
            row["q_value_bh_tax_group"] = min(1.0, next_value)


def highlight_candidate_sql(source_table: str) -> str:
    panel = sql_string_list(PRESPECIFIED_COFACTOR_NAMES)
    return f"""
SELECT 'TA_association'::VARCHAR AS selection_criterion,
       CASE WHEN c.ta_adjusted_odds_ratio > 1
            THEN 'enriched' ELSE 'depleted' END::VARCHAR AS selection_direction,
       c.ta_adjusted_odds_ratio::DOUBLE AS selection_estimate,
       abs(ln(c.ta_adjusted_odds_ratio))::DOUBLE AS selection_score,
       c.*
FROM {source_table} c
WHERE c.ta_evaluation_status='ok' AND c.ta_q_value_bh_tax_group <= 0.05
  AND c.ta_adjusted_odds_ratio <> 1
UNION ALL
SELECT 'DN_association'::VARCHAR,
       CASE WHEN c.dn_adjusted_odds_ratio > 1
            THEN 'enriched' ELSE 'depleted' END::VARCHAR,
       c.dn_adjusted_odds_ratio::DOUBLE,
       abs(ln(c.dn_adjusted_odds_ratio))::DOUBLE,
       c.*
FROM {source_table} c
WHERE c.dn_evaluation_status='ok' AND c.dn_q_value_bh_tax_group <= 0.05
  AND c.dn_adjusted_odds_ratio <> 1
UNION ALL
SELECT 'isoform_difference'::VARCHAR,
       CASE WHEN c.ta_vs_dn_odds_ratio_ratio > 1
            THEN 'TA_larger_OR' ELSE 'DN_larger_OR' END::VARCHAR,
       c.ta_vs_dn_odds_ratio_ratio::DOUBLE,
       abs(c.ta_vs_dn_log_odds_difference)::DOUBLE,
       c.*
FROM {source_table} c
WHERE c.evaluation_status='ok' AND c.q_value_bh_tax_group <= 0.05
  AND c.ta_vs_dn_odds_ratio_ratio <> 1
UNION ALL
SELECT 'opposite_direction'::VARCHAR,
       c.supported_direction_pattern::VARCHAR,
       c.ta_vs_dn_odds_ratio_ratio::DOUBLE,
       abs(c.ta_vs_dn_log_odds_difference)::DOUBLE,
       c.*
FROM {source_table} c
WHERE c.supported_direction_pattern <> 'not_supported'
UNION ALL
SELECT 'prespecified_panel'::VARCHAR, 'tracked'::VARCHAR,
       c.ta_vs_dn_odds_ratio_ratio::DOUBLE,
       NULL::DOUBLE,
       c.*
FROM {source_table} c
WHERE c.motif_name IN {panel}
"""


def finalize(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(config_path(run_root))
    verify_source(config)
    verify_plan(run_root, config)
    catalog_db = Path(config["jaspar_catalog"]) / "jaspar_metadata.duckdb"
    if (not catalog_db.is_file()
            or sha256(catalog_db) != config["jaspar_catalog_database_sha256"]):
        raise DistanceEnrichmentError("JASPAR metadata catalog changed")
    tasks = read_tsv(run_root / "plan" / "tasks.tsv")
    block_files: list[Path] = []
    class_files: list[Path] = []
    missing: list[str] = []
    for row in tasks:
        directory = task_directory(run_root, row)
        if not validate_existing_task(directory, row, config):
            missing.append(row["motif_id"])
            continue
        marker = load_json(directory / "complete.json")
        for name, collection in (
            ("block_components.parquet", block_files),
            ("class_counts.parquet", class_files),
        ):
            path = directory / name
            if sha256(path) != marker["files"][name]["sha256"]:
                raise DistanceEnrichmentError(f"task checksum differs: {path}")
            collection.append(path)
    if missing:
        raise DistanceEnrichmentError(
            f"cannot finalize; {len(missing)} tasks missing, first={missing[:5]}"
        )
    final = run_root / "final" / arguments.final_name
    if final.exists():
        manifest = load_json(final / "manifest.json")
        if manifest.get("run_id") != config["run_id"]:
            raise DistanceEnrichmentError("existing final output has another identity")
        if manifest.get("schema_version") != FINAL_SCHEMA_VERSION:
            raise DistanceEnrichmentError(
                "existing final output uses another schema; choose a new --final-name"
            )
        print(f"I: Reusing final output: {final}", file=sys.stderr)
        return
    finalizer_source = Path(__file__).resolve().parent.parent
    finalizer_commit = arguments.finalizer_source_commit or config["source_commit"]
    finalizer_identity = runtime_source_identity(
        finalizer_source, finalizer_commit, arguments.finalizer_source_dirty
    )
    if (arguments.finalizer_source_commit is None
            and finalizer_identity["scientific_source_sha256"]
            != config["scientific_source_sha256"]):
        raise DistanceEnrichmentError(
            "finalizer source differs from the task source; declare its commit with "
            "--finalizer-source-commit"
        )
    staging = Path(tempfile.mkdtemp(prefix=".distance-final-", dir=final.parent))
    try:
        table_root = staging / "tables" / "jaspar2026"
        table_root.mkdir(parents=True)
        database = staging / "tp73_distance_cofactor_enrichment.duckdb"
        raw_tsv = staging / "raw_estimate.tsv"
        raw_contrast_tsv = staging / "raw_isoform_contrast.tsv"
        run_duckdb(arguments.duckdb, database, f"""
SET preserve_insertion_order=false;
ATTACH {sql_string(catalog_db)} AS jaspar (READ_ONLY);
CREATE TABLE cofactor_distance_block_component AS
SELECT * FROM read_parquet({sql_path_list(block_files)});
CREATE TABLE cofactor_distance_class_count_chrom AS
SELECT * FROM read_parquet({sql_path_list(class_files)});
CREATE TABLE jaspar_matrix AS SELECT * FROM jaspar.jaspar_matrix;
CREATE TABLE jaspar_matrix_species AS SELECT * FROM jaspar.jaspar_matrix_species;
CREATE TABLE jaspar_motif_set_matrix AS
SELECT * FROM jaspar.jaspar_motif_set_matrix;
CREATE TABLE cofactor_distance_class_count AS
SELECT motif_id, max(motif_name) AS motif_name, distance_band,
       max(distance_band_order) AS distance_band_order,
       max(source_score_floor) AS source_score_floor,
       max(positive_threshold) AS positive_threshold,
       sum(anchors_total)::BIGINT AS anchors_total,
       sum(anchors_source_present)::BIGINT AS anchors_source_present,
       sum(anchors_positive)::BIGINT AS anchors_positive,
       sum(anchors_intermediate)::BIGINT AS anchors_intermediate,
       sum(anchors_negative)::BIGINT AS anchors_negative
FROM cofactor_distance_class_count_chrom
GROUP BY motif_id, distance_band;

COPY (
WITH total AS (
  SELECT motif_id, max(motif_name) AS motif_name, isoform, distance_band,
         max(distance_band_order) AS distance_band_order,
         max(source_score_floor) AS source_score_floor,
         max(positive_threshold) AS positive_threshold,
         sum(mh_numerator)::DOUBLE AS numerator,
         sum(mh_denominator)::DOUBLE AS denominator,
         sum(positive_anti_discordant)::BIGINT AS positive_anti_discordant,
         sum(positive_control_discordant)::BIGINT AS positive_control_discordant,
         sum(negative_anti_discordant)::BIGINT AS negative_anti_discordant,
         sum(negative_control_discordant)::BIGINT AS negative_control_discordant,
         count(DISTINCT (chrom, block_index))::BIGINT AS genomic_blocks
  FROM cofactor_distance_block_component
  GROUP BY motif_id, isoform, distance_band
), block AS (
  SELECT motif_id, isoform, distance_band, chrom, block_index,
         sum(mh_numerator)::DOUBLE AS block_numerator,
         sum(mh_denominator)::DOUBLE AS block_denominator
  FROM cofactor_distance_block_component
  GROUP BY motif_id, isoform, distance_band, chrom, block_index
), leave_one_out AS (
  SELECT t.motif_id, t.isoform, t.distance_band,
         CASE WHEN t.numerator - b.block_numerator > 0
                AND t.denominator - b.block_denominator > 0
              THEN ln((t.numerator - b.block_numerator)
                      / (t.denominator - b.block_denominator))
              ELSE NULL END AS leave_one_out_log_or
  FROM total t JOIN block b USING (motif_id, isoform, distance_band)
), jackknife AS (
  SELECT motif_id, isoform, distance_band,
         count(leave_one_out_log_or)::BIGINT AS jackknife_blocks,
         sqrt((count(leave_one_out_log_or) - 1)
              * var_pop(leave_one_out_log_or))::DOUBLE AS jackknife_se
  FROM leave_one_out
  GROUP BY motif_id, isoform, distance_band
), series AS (
  SELECT motif_id, isoform, distance_band, series_id,
         CASE WHEN sum(mh_numerator) > 0 AND sum(mh_denominator) > 0
              THEN sum(mh_numerator) / sum(mh_denominator) ELSE NULL END
              AS series_odds_ratio
  FROM cofactor_distance_block_component
  GROUP BY motif_id, isoform, distance_band, series_id
), series_wide AS (
  SELECT motif_id, isoform, distance_band,
         max(series_odds_ratio) FILTER (WHERE series_id='saos2') AS saos2_or,
         max(series_odds_ratio) FILTER (WHERE series_id='skmel29_2') AS skmel29_2_or
  FROM series GROUP BY motif_id, isoform, distance_band
), species_label AS (
  SELECT matrix_id,
         string_agg(species_scientific_name, '::' ORDER BY species_ordinal)
           AS source_species
  FROM jaspar_matrix_species GROUP BY matrix_id
)
SELECT t.*, c.anchors_total, c.anchors_source_present, c.anchors_positive,
       c.anchors_intermediate, c.anchors_negative,
       c.anchors_positive / c.anchors_total::DOUBLE AS positive_anchor_fraction,
       c.anchors_negative / c.anchors_total::DOUBLE AS negative_anchor_fraction,
       CASE WHEN t.numerator > 0 AND t.denominator > 0
            THEN t.numerator / t.denominator ELSE NULL END AS adjusted_odds_ratio,
       j.jackknife_blocks, j.jackknife_se,
       s.saos2_or, s.skmel29_2_or,
       m.collection, m.tax_group, m.class AS jaspar_tf_class,
       m.family AS jaspar_tf_family, m.includes_homo_sapiens,
       m.includes_mus_musculus, m.includes_rattus_norvegicus,
       m.species_count AS source_species_count,
       coalesce(sl.source_species, 'source_species_unspecified') AS source_species
FROM total t
JOIN cofactor_distance_class_count c
  USING (motif_id, distance_band)
JOIN jackknife j USING (motif_id, isoform, distance_band)
JOIN series_wide s USING (motif_id, isoform, distance_band)
JOIN jaspar_matrix m ON m.matrix_id=t.motif_id
LEFT JOIN species_label sl ON sl.matrix_id=t.motif_id
ORDER BY t.isoform, t.distance_band_order, t.motif_id
) TO {sql_string(raw_tsv)}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');

COPY (
WITH total AS (
  SELECT motif_id, isoform, distance_band,
         max(distance_band_order) AS distance_band_order,
         sum(mh_numerator)::DOUBLE AS numerator,
         sum(mh_denominator)::DOUBLE AS denominator
  FROM cofactor_distance_block_component
  GROUP BY motif_id, isoform, distance_band
), paired_total AS (
  SELECT motif_id, distance_band, max(distance_band_order) AS distance_band_order,
         max(numerator) FILTER (WHERE isoform='TA') AS ta_numerator,
         max(denominator) FILTER (WHERE isoform='TA') AS ta_denominator,
         max(numerator) FILTER (WHERE isoform='DN') AS dn_numerator,
         max(denominator) FILTER (WHERE isoform='DN') AS dn_denominator
  FROM total
  GROUP BY motif_id, distance_band
  HAVING count(DISTINCT isoform) = 2
), block AS (
  SELECT motif_id, distance_band, chrom, block_index,
         coalesce(sum(mh_numerator) FILTER (WHERE isoform='TA'), 0)::DOUBLE
           AS ta_block_numerator,
         coalesce(sum(mh_denominator) FILTER (WHERE isoform='TA'), 0)::DOUBLE
           AS ta_block_denominator,
         coalesce(sum(mh_numerator) FILTER (WHERE isoform='DN'), 0)::DOUBLE
           AS dn_block_numerator,
         coalesce(sum(mh_denominator) FILTER (WHERE isoform='DN'), 0)::DOUBLE
           AS dn_block_denominator
  FROM cofactor_distance_block_component
  GROUP BY motif_id, distance_band, chrom, block_index
), leave_one_out AS (
  SELECT t.motif_id, t.distance_band,
         CASE WHEN t.ta_numerator - b.ta_block_numerator > 0
                AND t.ta_denominator - b.ta_block_denominator > 0
                AND t.dn_numerator - b.dn_block_numerator > 0
                AND t.dn_denominator - b.dn_block_denominator > 0
              THEN ln((t.ta_numerator - b.ta_block_numerator)
                      / (t.ta_denominator - b.ta_block_denominator))
                   - ln((t.dn_numerator - b.dn_block_numerator)
                        / (t.dn_denominator - b.dn_block_denominator))
              ELSE NULL END AS leave_one_out_log_or_difference
  FROM paired_total t JOIN block b USING (motif_id, distance_band)
), jackknife AS (
  SELECT motif_id, distance_band,
         count(*)::BIGINT AS genomic_blocks,
         count(leave_one_out_log_or_difference)::BIGINT AS jackknife_blocks,
         sqrt((count(leave_one_out_log_or_difference) - 1)
              * var_pop(leave_one_out_log_or_difference))::DOUBLE AS jackknife_se
  FROM leave_one_out
  GROUP BY motif_id, distance_band
), series AS (
  SELECT motif_id, isoform, distance_band, series_id,
         CASE WHEN sum(mh_numerator) > 0 AND sum(mh_denominator) > 0
              THEN sum(mh_numerator) / sum(mh_denominator) ELSE NULL END
              AS series_odds_ratio
  FROM cofactor_distance_block_component
  GROUP BY motif_id, isoform, distance_band, series_id
), series_contrast AS (
  SELECT motif_id, distance_band,
         max(series_odds_ratio) FILTER (
           WHERE isoform='TA' AND series_id='saos2'
         ) / nullif(max(series_odds_ratio) FILTER (
           WHERE isoform='DN' AND series_id='saos2'
         ), 0) AS saos2_ta_vs_dn_odds_ratio_ratio,
         max(series_odds_ratio) FILTER (
           WHERE isoform='TA' AND series_id='skmel29_2'
         ) / nullif(max(series_odds_ratio) FILTER (
           WHERE isoform='DN' AND series_id='skmel29_2'
         ), 0) AS skmel29_2_ta_vs_dn_odds_ratio_ratio
  FROM series
  GROUP BY motif_id, distance_band
)
SELECT t.motif_id, t.distance_band, t.distance_band_order,
       CASE WHEN t.ta_numerator > 0 AND t.ta_denominator > 0
            THEN t.ta_numerator / t.ta_denominator ELSE NULL END
         AS ta_adjusted_odds_ratio,
       CASE WHEN t.dn_numerator > 0 AND t.dn_denominator > 0
            THEN t.dn_numerator / t.dn_denominator ELSE NULL END
         AS dn_adjusted_odds_ratio,
       CASE WHEN t.ta_numerator > 0 AND t.ta_denominator > 0
              AND t.dn_numerator > 0 AND t.dn_denominator > 0
            THEN ln(t.ta_numerator / t.ta_denominator)
                 - ln(t.dn_numerator / t.dn_denominator)
            ELSE NULL END AS ta_vs_dn_log_odds_difference,
       CASE WHEN t.ta_numerator > 0 AND t.ta_denominator > 0
              AND t.dn_numerator > 0 AND t.dn_denominator > 0
            THEN (t.ta_numerator / t.ta_denominator)
                 / (t.dn_numerator / t.dn_denominator)
            ELSE NULL END AS ta_vs_dn_odds_ratio_ratio,
       j.genomic_blocks, j.jackknife_blocks, j.jackknife_se,
       s.saos2_ta_vs_dn_odds_ratio_ratio,
       s.skmel29_2_ta_vs_dn_odds_ratio_ratio
FROM paired_total t
JOIN jackknife j USING (motif_id, distance_band)
JOIN series_contrast s USING (motif_id, distance_band)
ORDER BY t.distance_band_order, t.motif_id
) TO {sql_string(raw_contrast_tsv)}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
""")
        raw_rows = read_tsv(raw_tsv)
        result_rows: list[dict[str, Any]] = []
        numeric_float = {
            "source_score_floor", "positive_threshold", "numerator", "denominator",
            "positive_anchor_fraction", "negative_anchor_fraction",
            "adjusted_odds_ratio", "jackknife_se", "saos2_or", "skmel29_2_or",
        }
        numeric_int = {
            "distance_band_order", "positive_anti_discordant",
            "positive_control_discordant", "negative_anti_discordant",
            "negative_control_discordant", "genomic_blocks", "anchors_total",
            "anchors_source_present", "anchors_positive", "anchors_intermediate",
            "anchors_negative", "jackknife_blocks", "source_species_count",
        }
        for raw in raw_rows:
            row: dict[str, Any] = dict(raw)
            for name in numeric_float:
                row[name] = None if raw[name] == "NA" else float(raw[name])
            for name in numeric_int:
                row[name] = int(raw[name])
            row["includes_homo_sapiens"] = (
                raw["includes_homo_sapiens"].lower() == "true"
            )
            row["includes_mus_musculus"] = (
                raw["includes_mus_musculus"].lower() == "true"
            )
            row["includes_rattus_norvegicus"] = (
                raw["includes_rattus_norvegicus"].lower() == "true"
            )
            odds_ratio = row["adjusted_odds_ratio"]
            standard_error = row["jackknife_se"]
            estimable = (
                odds_ratio is not None and odds_ratio > 0
                and standard_error is not None and standard_error > 0
                and row["jackknife_blocks"] == row["genomic_blocks"]
                and row["jackknife_blocks"] >= 20
            )
            if estimable:
                log_or = math.log(odds_ratio)
                z_value = log_or / standard_error
                row["adjusted_log_odds"] = log_or
                row["z_value"] = z_value
                row["p_value"] = math.erfc(abs(z_value) / math.sqrt(2))
                row["confidence_interval_95_lower"] = math.exp(
                    log_or - 1.959963984540054 * standard_error
                )
                row["confidence_interval_95_upper"] = math.exp(
                    log_or + 1.959963984540054 * standard_error
                )
                row["evaluation_status"] = "ok"
            else:
                row["adjusted_log_odds"] = None
                row["z_value"] = None
                row["p_value"] = None
                row["confidence_interval_95_lower"] = None
                row["confidence_interval_95_upper"] = None
                row["evaluation_status"] = "not_estimable"
            row["q_value_bh_tax_group"] = None
            series_values = (row["saos2_or"], row["skmel29_2_or"])
            if all(value is not None and value > 0 for value in series_values):
                row["series_direction_consistent"] = (
                    all(value > 1 for value in series_values)
                    or all(value < 1 for value in series_values)
                )
            else:
                row["series_direction_consistent"] = None
            row["class_support_flag"] = (
                "negative_below_1pct" if row["negative_anchor_fraction"] < 0.01
                else "positive_below_1pct" if row["positive_anchor_fraction"] < 0.01
                else "adequate"
            )
            result_rows.append(row)
        bh_adjust(result_rows, len(tasks))
        result_json = staging / "result.jsonl"
        with result_json.open("x", encoding="utf-8") as handle:
            for row in result_rows:
                handle.write(json.dumps(row, sort_keys=True) + "\n")

        result_by_key = {
            (row["motif_id"], row["isoform"], row["distance_band"]): row
            for row in result_rows
        }
        contrast_rows: list[dict[str, Any]] = []
        contrast_float_fields = {
            "ta_adjusted_odds_ratio", "dn_adjusted_odds_ratio",
            "ta_vs_dn_log_odds_difference", "ta_vs_dn_odds_ratio_ratio",
            "jackknife_se", "saos2_ta_vs_dn_odds_ratio_ratio",
            "skmel29_2_ta_vs_dn_odds_ratio_ratio",
        }
        for raw in read_tsv(raw_contrast_tsv):
            motif_id = raw["motif_id"]
            distance_band = raw["distance_band"]
            ta = result_by_key[(motif_id, "TA", distance_band)]
            dn = result_by_key[(motif_id, "DN", distance_band)]
            row = {
                "motif_id": motif_id,
                "motif_name": ta["motif_name"],
                "distance_band": distance_band,
                "distance_band_order": int(raw["distance_band_order"]),
                "source_score_floor": ta["source_score_floor"],
                "positive_threshold": ta["positive_threshold"],
                "anchors_total": ta["anchors_total"],
                "anchors_source_present": ta["anchors_source_present"],
                "anchors_positive": ta["anchors_positive"],
                "anchors_intermediate": ta["anchors_intermediate"],
                "anchors_negative": ta["anchors_negative"],
                "positive_anchor_fraction": ta["positive_anchor_fraction"],
                "negative_anchor_fraction": ta["negative_anchor_fraction"],
                "class_support_flag": ta["class_support_flag"],
                "collection": ta["collection"],
                "tax_group": ta["tax_group"],
                "jaspar_tf_class": ta["jaspar_tf_class"],
                "jaspar_tf_family": ta["jaspar_tf_family"],
                "includes_homo_sapiens": ta["includes_homo_sapiens"],
                "includes_mus_musculus": ta["includes_mus_musculus"],
                "includes_rattus_norvegicus": ta["includes_rattus_norvegicus"],
                "source_species_count": ta["source_species_count"],
                "source_species": ta["source_species"],
                "genomic_blocks": int(raw["genomic_blocks"]),
                "jackknife_blocks": int(raw["jackknife_blocks"]),
            }
            for name in contrast_float_fields:
                row[name] = None if raw[name] == "NA" else float(raw[name])
            for isoform_name, source in (("ta", ta), ("dn", dn)):
                for name in (
                    "confidence_interval_95_lower",
                    "confidence_interval_95_upper",
                    "q_value_bh_tax_group",
                    "evaluation_status",
                    "saos2_or",
                    "skmel29_2_or",
                    "series_direction_consistent",
                ):
                    row[f"{isoform_name}_{name}"] = source[name]
            log_difference = row["ta_vs_dn_log_odds_difference"]
            ratio = row["ta_vs_dn_odds_ratio_ratio"]
            standard_error = row["jackknife_se"]
            estimable = (
                log_difference is not None and ratio is not None and ratio > 0
                and standard_error is not None and standard_error > 0
                and row["jackknife_blocks"] == row["genomic_blocks"]
                and row["jackknife_blocks"] >= 20
            )
            if estimable:
                z_value = log_difference / standard_error
                row["z_value"] = z_value
                row["p_value"] = math.erfc(abs(z_value) / math.sqrt(2))
                row["confidence_interval_95_lower"] = math.exp(
                    log_difference - 1.959963984540054 * standard_error
                )
                row["confidence_interval_95_upper"] = math.exp(
                    log_difference + 1.959963984540054 * standard_error
                )
                row["evaluation_status"] = "ok"
            else:
                row["z_value"] = None
                row["p_value"] = None
                row["confidence_interval_95_lower"] = None
                row["confidence_interval_95_upper"] = None
                row["evaluation_status"] = "not_estimable"
            row["q_value_bh_tax_group"] = None
            series_ratios = (
                row["saos2_ta_vs_dn_odds_ratio_ratio"],
                row["skmel29_2_ta_vs_dn_odds_ratio_ratio"],
            )
            if all(value is not None and value > 0 for value in series_ratios):
                row["series_contrast_direction_consistent"] = (
                    all(value > 1 for value in series_ratios)
                    or all(value < 1 for value in series_ratios)
                )
            else:
                row["series_contrast_direction_consistent"] = None
            contrast_rows.append(row)
        bh_adjust_isoform_contrast(contrast_rows, len(tasks))
        for row in contrast_rows:
            ta_odds_ratio = row["ta_adjusted_odds_ratio"]
            dn_odds_ratio = row["dn_adjusted_odds_ratio"]
            if ta_odds_ratio is None or dn_odds_ratio is None:
                row["point_direction_pattern"] = "not_estimable"
            elif ta_odds_ratio > 1 and dn_odds_ratio > 1:
                row["point_direction_pattern"] = "both_enriched"
            elif ta_odds_ratio < 1 and dn_odds_ratio < 1:
                row["point_direction_pattern"] = "both_depleted"
            elif ta_odds_ratio > 1 and dn_odds_ratio < 1:
                row["point_direction_pattern"] = "TA_enriched_DN_depleted"
            elif ta_odds_ratio < 1 and dn_odds_ratio > 1:
                row["point_direction_pattern"] = "DN_enriched_TA_depleted"
            else:
                row["point_direction_pattern"] = "on_null_boundary"

            contrast_q = row["q_value_bh_tax_group"]
            ta_q = row["ta_q_value_bh_tax_group"]
            dn_q = row["dn_q_value_bh_tax_group"]
            all_tests_supported = (
                row["evaluation_status"] == "ok"
                and row["ta_evaluation_status"] == "ok"
                and row["dn_evaluation_status"] == "ok"
                and contrast_q is not None and contrast_q <= 0.05
                and ta_q is not None and ta_q <= 0.05
                and dn_q is not None and dn_q <= 0.05
            )
            if (all_tests_supported
                    and row["ta_confidence_interval_95_lower"] > 1
                    and row["dn_confidence_interval_95_upper"] < 1):
                row["supported_direction_pattern"] = "TA_enriched_DN_depleted"
            elif (all_tests_supported
                  and row["dn_confidence_interval_95_lower"] > 1
                  and row["ta_confidence_interval_95_upper"] < 1):
                row["supported_direction_pattern"] = "DN_enriched_TA_depleted"
            else:
                row["supported_direction_pattern"] = "not_supported"
        contrast_json = staging / "isoform_contrast.jsonl"
        with contrast_json.open("x", encoding="utf-8") as handle:
            for row in contrast_rows:
                handle.write(json.dumps(row, sort_keys=True) + "\n")

        paths = {
            "cofactor_distance_block_component":
                table_root / "cofactor_distance_block_component" / "part-000000.parquet",
            "cofactor_distance_class_count":
                table_root / "cofactor_distance_class_count" / "part-000000.parquet",
            "cofactor_distance_enrichment":
                table_root / "cofactor_distance_enrichment" / "part-000000.parquet",
            "cofactor_distance_isoform_contrast":
                table_root / "cofactor_distance_isoform_contrast" /
                "part-000000.parquet",
            "cofactor_distance_isoform_comparison":
                table_root / "cofactor_distance_isoform_comparison" /
                "part-000000.parquet",
            "cofactor_distance_isoform_comparison_human_source":
                table_root / "cofactor_distance_isoform_comparison_human_source" /
                "part-000000.parquet",
            "cofactor_distance_highlight_20":
                table_root / "cofactor_distance_highlight_20" /
                "part-000000.parquet",
            "cofactor_distance_highlight_20_human_source":
                table_root / "cofactor_distance_highlight_20_human_source" /
                "part-000000.parquet",
            "cofactor_distance_top_bottom_20":
                table_root / "cofactor_distance_top_bottom_20" / "part-000000.parquet",
            "cofactor_distance_top_bottom_20_human_source":
                table_root / "cofactor_distance_top_bottom_20_human_source" /
                "part-000000.parquet",
            "jaspar_matrix": table_root / "jaspar_matrix" / "part-000000.parquet",
            "jaspar_matrix_species":
                table_root / "jaspar_matrix_species" / "part-000000.parquet",
            "jaspar_motif_set_matrix":
                table_root / "jaspar_motif_set_matrix" / "part-000000.parquet",
        }
        for path in paths.values():
            path.parent.mkdir(parents=True, exist_ok=True)
        run_duckdb(arguments.duckdb, database, f"""
CREATE TABLE cofactor_distance_enrichment AS
SELECT * FROM read_json_auto({sql_string(result_json)}, format='newline_delimited');
CREATE TABLE cofactor_distance_isoform_contrast AS
SELECT * FROM read_json_auto(
  {sql_string(contrast_json)}, format='newline_delimited'
);
CREATE TABLE cofactor_distance_isoform_comparison AS
SELECT c.*,
       CASE WHEN ta_evaluation_status='ok' THEN row_number() OVER (
         PARTITION BY distance_band
         ORDER BY CASE WHEN ta_evaluation_status='ok' THEN 0 ELSE 1 END,
                  ta_adjusted_odds_ratio DESC NULLS LAST, motif_id
       ) END::BIGINT AS ta_enrichment_rank,
       CASE WHEN ta_evaluation_status='ok' THEN row_number() OVER (
         PARTITION BY distance_band
         ORDER BY CASE WHEN ta_evaluation_status='ok' THEN 0 ELSE 1 END,
                  ta_adjusted_odds_ratio ASC NULLS LAST, motif_id
       ) END::BIGINT AS ta_depletion_rank,
       CASE WHEN dn_evaluation_status='ok' THEN row_number() OVER (
         PARTITION BY distance_band
         ORDER BY CASE WHEN dn_evaluation_status='ok' THEN 0 ELSE 1 END,
                  dn_adjusted_odds_ratio DESC NULLS LAST, motif_id
       ) END::BIGINT AS dn_enrichment_rank,
       CASE WHEN dn_evaluation_status='ok' THEN row_number() OVER (
         PARTITION BY distance_band
         ORDER BY CASE WHEN dn_evaluation_status='ok' THEN 0 ELSE 1 END,
                  dn_adjusted_odds_ratio ASC NULLS LAST, motif_id
       ) END::BIGINT AS dn_depletion_rank,
       CASE WHEN evaluation_status='ok' THEN row_number() OVER (
         PARTITION BY distance_band
         ORDER BY CASE WHEN evaluation_status='ok' THEN 0 ELSE 1 END,
                  abs(ta_vs_dn_log_odds_difference) DESC NULLS LAST, motif_id
       ) END::BIGINT AS isoform_difference_absolute_rank
FROM cofactor_distance_isoform_contrast c;
CREATE TABLE cofactor_distance_isoform_comparison_human_source AS
SELECT c.*,
       CASE WHEN ta_evaluation_status='ok' THEN row_number() OVER (
         PARTITION BY distance_band
         ORDER BY CASE WHEN ta_evaluation_status='ok' THEN 0 ELSE 1 END,
                  ta_adjusted_odds_ratio DESC NULLS LAST, motif_id
       ) END::BIGINT AS ta_enrichment_rank,
       CASE WHEN ta_evaluation_status='ok' THEN row_number() OVER (
         PARTITION BY distance_band
         ORDER BY CASE WHEN ta_evaluation_status='ok' THEN 0 ELSE 1 END,
                  ta_adjusted_odds_ratio ASC NULLS LAST, motif_id
       ) END::BIGINT AS ta_depletion_rank,
       CASE WHEN dn_evaluation_status='ok' THEN row_number() OVER (
         PARTITION BY distance_band
         ORDER BY CASE WHEN dn_evaluation_status='ok' THEN 0 ELSE 1 END,
                  dn_adjusted_odds_ratio DESC NULLS LAST, motif_id
       ) END::BIGINT AS dn_enrichment_rank,
       CASE WHEN dn_evaluation_status='ok' THEN row_number() OVER (
         PARTITION BY distance_band
         ORDER BY CASE WHEN dn_evaluation_status='ok' THEN 0 ELSE 1 END,
                  dn_adjusted_odds_ratio ASC NULLS LAST, motif_id
       ) END::BIGINT AS dn_depletion_rank,
       CASE WHEN evaluation_status='ok' THEN row_number() OVER (
         PARTITION BY distance_band
         ORDER BY CASE WHEN evaluation_status='ok' THEN 0 ELSE 1 END,
                  abs(ta_vs_dn_log_odds_difference) DESC NULLS LAST, motif_id
       ) END::BIGINT AS isoform_difference_absolute_rank
FROM cofactor_distance_isoform_contrast c
WHERE includes_homo_sapiens;
CREATE TEMP VIEW cofactor_distance_highlight_candidate AS
{highlight_candidate_sql('cofactor_distance_isoform_comparison')};
CREATE TABLE cofactor_distance_highlight_20 AS
WITH ranked AS (
  SELECT CASE WHEN selection_criterion='prespecified_panel' THEN NULL
              ELSE row_number() OVER (
                PARTITION BY selection_criterion, selection_direction,
                             distance_band
                ORDER BY selection_score DESC NULLS LAST, motif_id
              ) END::BIGINT AS selection_rank,
         *
  FROM cofactor_distance_highlight_candidate
)
SELECT * FROM ranked
WHERE selection_criterion='prespecified_panel' OR selection_rank <= 20
ORDER BY distance_band_order, selection_criterion, selection_direction,
         selection_rank NULLS LAST, motif_name, motif_id;
CREATE TEMP VIEW cofactor_distance_highlight_candidate_human_source AS
{highlight_candidate_sql('cofactor_distance_isoform_comparison_human_source')};
CREATE TABLE cofactor_distance_highlight_20_human_source AS
WITH ranked AS (
  SELECT CASE WHEN selection_criterion='prespecified_panel' THEN NULL
              ELSE row_number() OVER (
                PARTITION BY selection_criterion, selection_direction,
                             distance_band
                ORDER BY selection_score DESC NULLS LAST, motif_id
              ) END::BIGINT AS selection_rank,
         *
  FROM cofactor_distance_highlight_candidate_human_source
)
SELECT * FROM ranked
WHERE selection_criterion='prespecified_panel' OR selection_rank <= 20
ORDER BY distance_band_order, selection_criterion, selection_direction,
         selection_rank NULLS LAST, motif_name, motif_id;
CREATE TABLE cofactor_distance_top_bottom_20 AS
WITH ranked AS (
  SELECT 'enriched'::VARCHAR AS direction,
         row_number() OVER (
           PARTITION BY isoform, distance_band
           ORDER BY adjusted_odds_ratio DESC, motif_id
         ) AS rank, *
  FROM cofactor_distance_enrichment WHERE evaluation_status='ok'
  UNION ALL
  SELECT 'depleted'::VARCHAR AS direction,
         row_number() OVER (
           PARTITION BY isoform, distance_band
           ORDER BY adjusted_odds_ratio ASC, motif_id
         ) AS rank, *
  FROM cofactor_distance_enrichment WHERE evaluation_status='ok'
), selected AS (
  SELECT * FROM ranked WHERE rank <= 20
)
SELECT s.*,
       c.ta_adjusted_odds_ratio, c.ta_confidence_interval_95_lower,
       c.ta_confidence_interval_95_upper, c.ta_q_value_bh_tax_group,
       c.dn_adjusted_odds_ratio, c.dn_confidence_interval_95_lower,
       c.dn_confidence_interval_95_upper, c.dn_q_value_bh_tax_group,
       c.ta_vs_dn_odds_ratio_ratio, c.ta_vs_dn_log_odds_difference,
       c.confidence_interval_95_lower AS ta_vs_dn_confidence_interval_95_lower,
       c.confidence_interval_95_upper AS ta_vs_dn_confidence_interval_95_upper,
       c.p_value AS ta_vs_dn_p_value,
       c.q_value_bh_tax_group AS ta_vs_dn_q_value_bh_tax_group,
       c.evaluation_status AS ta_vs_dn_evaluation_status,
       c.series_contrast_direction_consistent
FROM selected s
JOIN cofactor_distance_isoform_contrast c USING (motif_id, distance_band)
ORDER BY s.isoform, s.distance_band_order,
         CASE s.direction WHEN 'enriched' THEN 1 ELSE 2 END, s.rank;
CREATE TABLE cofactor_distance_top_bottom_20_human_source AS
WITH ranked AS (
  SELECT 'enriched'::VARCHAR AS direction,
         row_number() OVER (
           PARTITION BY isoform, distance_band
           ORDER BY adjusted_odds_ratio DESC, motif_id
         ) AS rank, *
  FROM cofactor_distance_enrichment
  WHERE evaluation_status='ok' AND includes_homo_sapiens
  UNION ALL
  SELECT 'depleted'::VARCHAR AS direction,
         row_number() OVER (
           PARTITION BY isoform, distance_band
           ORDER BY adjusted_odds_ratio ASC, motif_id
         ) AS rank, *
  FROM cofactor_distance_enrichment
  WHERE evaluation_status='ok' AND includes_homo_sapiens
), selected AS (
  SELECT * FROM ranked WHERE rank <= 20
)
SELECT s.*,
       c.ta_adjusted_odds_ratio, c.ta_confidence_interval_95_lower,
       c.ta_confidence_interval_95_upper, c.ta_q_value_bh_tax_group,
       c.dn_adjusted_odds_ratio, c.dn_confidence_interval_95_lower,
       c.dn_confidence_interval_95_upper, c.dn_q_value_bh_tax_group,
       c.ta_vs_dn_odds_ratio_ratio, c.ta_vs_dn_log_odds_difference,
       c.confidence_interval_95_lower AS ta_vs_dn_confidence_interval_95_lower,
       c.confidence_interval_95_upper AS ta_vs_dn_confidence_interval_95_upper,
       c.p_value AS ta_vs_dn_p_value,
       c.q_value_bh_tax_group AS ta_vs_dn_q_value_bh_tax_group,
       c.evaluation_status AS ta_vs_dn_evaluation_status,
       c.series_contrast_direction_consistent
FROM selected s
JOIN cofactor_distance_isoform_contrast c USING (motif_id, distance_band)
ORDER BY s.isoform, s.distance_band_order,
         CASE s.direction WHEN 'enriched' THEN 1 ELSE 2 END, s.rank;
CREATE UNIQUE INDEX distance_result_idx
  ON cofactor_distance_enrichment(motif_id, isoform, distance_band);
CREATE UNIQUE INDEX distance_isoform_contrast_idx
  ON cofactor_distance_isoform_contrast(motif_id, distance_band);
CREATE UNIQUE INDEX distance_isoform_comparison_idx
  ON cofactor_distance_isoform_comparison(motif_id, distance_band);
CREATE UNIQUE INDEX distance_isoform_comparison_human_idx
  ON cofactor_distance_isoform_comparison_human_source(motif_id, distance_band);
COPY cofactor_distance_block_component
  TO {sql_string(paths['cofactor_distance_block_component'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_distance_class_count
  TO {sql_string(paths['cofactor_distance_class_count'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_distance_enrichment
  TO {sql_string(paths['cofactor_distance_enrichment'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_distance_isoform_contrast
  TO {sql_string(paths['cofactor_distance_isoform_contrast'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_distance_isoform_comparison
  TO {sql_string(paths['cofactor_distance_isoform_comparison'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_distance_isoform_comparison_human_source
  TO {sql_string(paths['cofactor_distance_isoform_comparison_human_source'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_distance_highlight_20
  TO {sql_string(paths['cofactor_distance_highlight_20'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_distance_highlight_20_human_source
  TO {sql_string(paths['cofactor_distance_highlight_20_human_source'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_distance_top_bottom_20
  TO {sql_string(paths['cofactor_distance_top_bottom_20'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_distance_top_bottom_20_human_source
  TO {sql_string(paths['cofactor_distance_top_bottom_20_human_source'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY jaspar_matrix TO {sql_string(paths['jaspar_matrix'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY jaspar_matrix_species TO {sql_string(paths['jaspar_matrix_species'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY jaspar_motif_set_matrix TO {sql_string(paths['jaspar_motif_set_matrix'])}
  (FORMAT PARQUET, COMPRESSION ZSTD);
COPY cofactor_distance_enrichment
  TO {sql_string(staging / 'cofactor_distance_enrichment.tsv')}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
COPY cofactor_distance_isoform_contrast
  TO {sql_string(staging / 'cofactor_distance_isoform_contrast.tsv')}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
COPY cofactor_distance_isoform_comparison
  TO {sql_string(staging / 'cofactor_distance_isoform_comparison.tsv')}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
COPY cofactor_distance_isoform_comparison_human_source
  TO {sql_string(staging / 'cofactor_distance_isoform_comparison_human_source.tsv')}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
COPY cofactor_distance_highlight_20
  TO {sql_string(staging / 'cofactor_distance_highlight_20.tsv')}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
COPY cofactor_distance_highlight_20_human_source
  TO {sql_string(staging / 'cofactor_distance_highlight_20_human_source.tsv')}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
COPY cofactor_distance_top_bottom_20
  TO {sql_string(staging / 'cofactor_distance_top_bottom_20.tsv')}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
COPY cofactor_distance_top_bottom_20_human_source
  TO {sql_string(staging / 'cofactor_distance_top_bottom_20_human_source.tsv')}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
CHECKPOINT;
""")
        result_json.unlink()
        contrast_json.unlink()
        raw_tsv.unlink()
        raw_contrast_tsv.unlink()
        counts = run_json(arguments.duckdb, database, """
SELECT (SELECT count(*) FROM cofactor_distance_enrichment)::BIGINT AS results,
       (SELECT count(*) FROM cofactor_distance_top_bottom_20)::BIGINT AS rankings,
       (SELECT count(*) FROM cofactor_distance_enrichment
        WHERE evaluation_status='ok')::BIGINT AS estimable,
       (SELECT count(*) FROM cofactor_distance_isoform_contrast)::BIGINT
         AS isoform_contrasts,
       (SELECT count(*) FROM cofactor_distance_isoform_contrast
        WHERE evaluation_status='ok')::BIGINT AS estimable_isoform_contrasts,
       (SELECT count(*) FROM cofactor_distance_isoform_comparison)::BIGINT
         AS isoform_comparisons,
       (SELECT count(*) FROM cofactor_distance_isoform_comparison_human_source)
         ::BIGINT AS human_source_isoform_comparisons,
       (SELECT count(*) FROM cofactor_distance_highlight_20)::BIGINT
         AS highlights,
       (SELECT count(*) FROM cofactor_distance_highlight_20_human_source)::BIGINT
         AS human_source_highlights,
       (SELECT count(*) FROM cofactor_distance_highlight_20
        WHERE selection_criterion='prespecified_panel')::BIGINT
         AS prespecified_panel_rows;
""")[0]
        expected_results = len(tasks) * len(DISTANCE_BANDS) * 2
        if int(counts["results"]) != expected_results:
            raise DistanceEnrichmentError("final result cardinality differs")
        expected_contrasts = len(tasks) * len(DISTANCE_BANDS)
        if int(counts["isoform_contrasts"]) != expected_contrasts:
            raise DistanceEnrichmentError("isoform-contrast cardinality differs")
        if int(counts["isoform_comparisons"]) != expected_contrasts:
            raise DistanceEnrichmentError("isoform-comparison cardinality differs")
        expected_human_comparisons = (
            int(config["human_source_task_count"]) * len(DISTANCE_BANDS)
        )
        if int(counts["human_source_isoform_comparisons"]) \
                != expected_human_comparisons:
            raise DistanceEnrichmentError(
                "human-source isoform-comparison cardinality differs"
            )
        manifest = {
            "schema_version": FINAL_SCHEMA_VERSION,
            "run_id": config["run_id"],
            "analysis": config["analysis"],
            "tax_group": config["tax_group"],
            "taxonomic_scope_rule": config["taxonomic_scope_rule"],
            "source_species_semantics": config["source_species_semantics"],
            "chromosomes": config["chromosomes"],
            "source_commit": config["source_commit"],
            "task_source_commit": config["source_commit"],
            "finalizer_source": finalizer_identity["source"],
            "finalizer_source_commit": finalizer_identity["commit"],
            "finalizer_source_dirty": finalizer_identity["dirty"],
            "finalizer_scientific_source_sha256":
                finalizer_identity["scientific_source_sha256"],
            "task_count": len(tasks),
            "human_source_task_count": config["human_source_task_count"],
            "source_species_unspecified_task_count":
                config["source_species_unspecified_task_count"],
            "result_rows": int(counts["results"]),
            "estimable_rows": int(counts["estimable"]),
            "ranking_rows": int(counts["rankings"]),
            "isoform_contrast_rows": int(counts["isoform_contrasts"]),
            "estimable_isoform_contrast_rows":
                int(counts["estimable_isoform_contrasts"]),
            "isoform_comparison_rows": int(counts["isoform_comparisons"]),
            "human_source_isoform_comparison_rows":
                int(counts["human_source_isoform_comparisons"]),
            "highlight_rows": int(counts["highlights"]),
            "human_source_highlight_rows":
                int(counts["human_source_highlights"]),
            "prespecified_panel_rows": int(counts["prespecified_panel_rows"]),
            "prespecified_cofactor_names": list(PRESPECIFIED_COFACTOR_NAMES),
            "multiple_testing_family_size_per_isoform_band": len(tasks),
            "multiple_testing_family_size_per_isoform_contrast_band": len(tasks),
            "uncertainty": "5Mb genomic-block delete-one-cluster jackknife",
            "isoform_contrast_uncertainty":
                "paired 5Mb genomic-block delete-one-cluster jackknife",
            "common_effect": "Mantel-Haenszel across TP73 score strata and series",
            "distance_bands": list(DISTANCE_BANDS),
            "excluded_series": config["excluded_series"],
            "jaspar_catalog_manifest_sha256":
                config["jaspar_catalog_manifest_sha256"],
            "database": database.name,
            "database_sha256": sha256(database),
            "tables": {
                name: {
                    "path": str(path.relative_to(staging)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for name, path in paths.items()
            },
            "completed_at_utc": utc_now(),
        }
        write_json_new(staging / "manifest.json", manifest)
        os.replace(staging, final)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    print(f"I: Finalized distance enrichment: {final}", file=sys.stderr)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare", help="write immutable plans")
    prepare_parser.add_argument("--run-root", type=Path, required=True)
    prepare_parser.add_argument("--scan-package", type=Path, required=True)
    prepare_parser.add_argument("--anchor-evidence", type=Path, required=True)
    prepare_parser.add_argument("--thresholds", type=Path, required=True)
    prepare_parser.add_argument("--threshold-set-id", required=True)
    prepare_parser.add_argument("--jaspar-catalog", type=Path, required=True)
    prepare_parser.add_argument("--source", type=Path, required=True)
    prepare_parser.add_argument("--source-commit", required=True)
    prepare_parser.add_argument("--run-id", required=True)
    prepare_parser.add_argument("--tax-group", default="vertebrates")
    prepare_parser.add_argument(
        "--chromosomes", default=",".join(DEFAULT_CHROMOSOMES),
        help="comma-separated positive chromosome numbers (default: 1-22)",
    )
    prepare_parser.add_argument("--anchor-motif", default="MA0861.2")
    prepare_parser.add_argument("--flank", type=int, default=150)
    prepare_parser.add_argument("--block-size", type=int, default=5_000_000)
    prepare_parser.add_argument("--duckdb", default="duckdb")
    prepare_parser.set_defaults(function=prepare)

    anchor_parser = subparsers.add_parser(
        "prepare-anchors", help="split the shared anchor evidence by chromosome"
    )
    anchor_parser.add_argument("--run-root", type=Path, required=True)
    anchor_parser.add_argument("--threads", type=int, default=2)
    anchor_parser.add_argument("--memory-limit", default="16GB")
    anchor_parser.add_argument("--duckdb", default="duckdb")
    anchor_parser.set_defaults(function=prepare_anchors)

    task_parser = subparsers.add_parser("run-task", help="run one motif task")
    task_parser.add_argument("--run-root", type=Path, required=True)
    task_parser.add_argument("--task-index", type=int)
    task_parser.add_argument("--task-offset", type=int, default=0)
    task_parser.add_argument("--scratch", type=Path, required=True)
    task_parser.add_argument("--threads", type=int, default=2)
    task_parser.add_argument("--memory-limit", default="24GB")
    task_parser.add_argument("--max-temp-size", default="200GB")
    task_parser.add_argument("--duckdb", default="duckdb")
    task_parser.set_defaults(function=run_task)

    status_parser = subparsers.add_parser("status", help="report completion counts")
    status_parser.add_argument("--run-root", type=Path, required=True)
    status_parser.set_defaults(function=status)

    finalize_parser = subparsers.add_parser(
        "finalize", help="validate and publish DuckDB/Parquet rankings"
    )
    finalize_parser.add_argument("--run-root", type=Path, required=True)
    finalize_parser.add_argument("--final-name", default="distance_enrichment")
    finalize_parser.add_argument(
        "--finalizer-source-commit",
        help="commit of a finalizer newer than the task source (40 lowercase hex digits)",
    )
    finalize_parser.add_argument(
        "--finalizer-source-dirty", action="store_true",
        help="record that finalizer scientific files differ from the declared commit",
    )
    finalize_parser.add_argument("--duckdb", default="duckdb")
    finalize_parser.set_defaults(function=finalize)
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        arguments.function(arguments)
    except (DistanceEnrichmentError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

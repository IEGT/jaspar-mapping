#!/usr/bin/env python3

"""Plan, run, inspect, and finalize genome-wide H3K4me3 anchor signal."""

from __future__ import annotations

import argparse
import csv
import errno
import hashlib
import io
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class GenomeSignalError(RuntimeError):
    pass


WINDOWS = (
    ("central_150", 0, 150),
    ("central_500", 0, 500),
    ("central_1000", 0, 1000),
    ("flank_150_1000", 150, 1000),
)
PRIMARY_WINDOW = "flank_150_1000"
PSEUDOCOUNT = 1.0
CONDITIONS = ("GFP", "TA", "DN")
SCIENTIFIC_SOURCE_FILES = (
    "scripts/build_h3k4me3_anchor_signal.py",
    "scripts/export_bigwig_chrom_bedgraph.py",
    "scripts/manage_h3k4me3_genome_signal.py",
    "scripts/run_h3k4me3_genome_signal_finalize.sh",
    "scripts/submit_h3k4me3_genome_signal_slurm.sh",
)

STARTED = time.monotonic()
CURRENT_PHASE = "startup"
CURRENT_CHROM = "unresolved"
CURRENT_ATTEMPT: Path | None = None
CURRENT_SCRATCH: Path | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_path_list(paths: list[Path]) -> str:
    return "[" + ",".join(sql_string(path) for path in paths) + "]"


def safe_label(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", value).strip("-") or "region"


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise GenomeSignalError(f"JSON file is missing: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise GenomeSignalError(f"JSON value is not an object: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise GenomeSignalError(f"TSV file is missing: {path}")
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        rows = list(reader)
    if reader.fieldnames is None or not rows:
        raise GenomeSignalError(f"TSV is empty or lacks a header: {path}")
    return rows


def immutable_write(path: Path, content: str) -> None:
    encoded = content.encode("utf-8")
    if path.exists():
        if path.read_bytes() != encoded:
            raise GenomeSignalError(f"immutable file differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        temporary.write_bytes(encoded)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def write_rows(path: Path, rows: list[dict[str, Any]],
               fields: tuple[str, ...]) -> None:
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    immutable_write(path, output.getvalue())


def query_json(duckdb: Path, database: Path | str,
               query: str) -> list[dict[str, Any]]:
    database_arguments = [] if str(database) == ":memory:" else ["-readonly"]
    process = subprocess.run(
        [str(duckdb), "-light-mode", *database_arguments,
         "-json", str(database), "-c", query],
        text=True, capture_output=True, check=False,
    )
    if process.returncode != 0:
        raise GenomeSignalError(process.stderr.strip() or "DuckDB query failed")
    value = json.loads(process.stdout or "[]")
    if not isinstance(value, list):
        raise GenomeSignalError("DuckDB result is not a row array")
    return value


def run(command: list[str], *, cwd: Path | None = None) -> None:
    displayed = list(command)
    if "-c" in displayed:
        index = displayed.index("-c") + 1
        if index < len(displayed) and len(displayed[index]) > 200:
            displayed[index] = f"<SQL:{len(displayed[index])} chars>"
    print("I: running: " + " ".join(displayed), file=sys.stderr, flush=True)
    process = subprocess.run(command, cwd=cwd, check=False)
    if process.returncode != 0:
        raise GenomeSignalError(
            f"command failed with exit code {process.returncode}: {command[0]}"
        )


def git_identity(source: Path) -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        text=True, capture_output=True, check=False,
    )
    if commit.returncode != 0:
        raise GenomeSignalError(commit.stderr.strip() or "cannot read Git commit")
    dirty = False
    for arguments in (
        ["diff", "--quiet", "--ignore-submodules", "--"],
        ["diff", "--cached", "--quiet", "--ignore-submodules", "--"],
    ):
        result = subprocess.run(
            ["git", "-C", str(source), *arguments], check=False
        )
        if result.returncode not in (0, 1):
            raise GenomeSignalError("cannot inspect Git state")
        dirty = dirty or result.returncode == 1
    return commit.stdout.strip(), dirty


def scientific_hashes(source: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for relative in SCIENTIFIC_SOURCE_FILES:
        path = source / relative
        if not path.is_file():
            raise GenomeSignalError(f"scientific source is missing: {path}")
        result[relative] = sha256(path)
    return result


def verify_scientific_hashes(config: dict[str, Any]) -> None:
    source = Path(str(config["source"]))
    expected = config.get("scientific_source_file_sha256")
    if not isinstance(expected, dict) or not expected:
        raise GenomeSignalError("run plan lacks scientific source hashes")
    for relative, digest in expected.items():
        path = source / relative
        if not path.is_file() or sha256(path) != digest:
            raise GenomeSignalError(f"scientific source changed: {path}")


def verify_run_inputs(run_root: Path, config: dict[str, Any]) -> None:
    checks = (
        (run_root / "plan" / "chromosome_tasks.tsv",
         config["task_plan_sha256"], "task plan"),
        (Path(str(config["evidence_package"])) / "manifest.json",
         config["evidence_manifest_sha256"], "evidence manifest"),
        (Path(str(config["evidence_package"])) /
         "chromosome_file_inventory.tsv",
         config["evidence_inventory_sha256"], "evidence inventory"),
        (Path(str(config["track_manifest"])),
         config["track_manifest_sha256"], "track manifest"),
        (Path(str(config["track_file_inventory"])),
         config["track_file_inventory_sha256"], "track file inventory"),
    )
    for path, expected, label in checks:
        if not path.is_file() or sha256(path) != expected:
            raise GenomeSignalError(f"{label} changed after planning: {path}")


def tree_bytes(path: Path | None) -> int:
    if path is None or not path.exists():
        return 0
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def progress(_signal_number: int | None = None,
             _frame: object | None = None) -> None:
    print(
        "I: progress "
        f"chrom={CURRENT_CHROM} phase={CURRENT_PHASE} "
        f"elapsed_seconds={time.monotonic() - STARTED:.1f} "
        f"durable_bytes={tree_bytes(CURRENT_ATTEMPT)} "
        f"scratch_bytes={tree_bytes(CURRENT_SCRATCH)}",
        file=sys.stderr, flush=True,
    )


def set_phase(name: str) -> None:
    global CURRENT_PHASE
    CURRENT_PHASE = name
    progress()


def check_free_space(path: Path, minimum_gb: float, label: str) -> None:
    free = shutil.disk_usage(path).free
    print(
        f"I: {label} free space: {free / 1024**3:.1f} GiB "
        f"(required {minimum_gb:.1f} GiB)",
        file=sys.stderr, flush=True,
    )
    if free < int(minimum_gb * 1024**3):
        raise GenomeSignalError(f"insufficient {label} space")


def contract_sha256(contract: dict[str, Any]) -> str:
    encoded = json.dumps(contract, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def is_mitochondrial(row: dict[str, str]) -> bool:
    normalized = row["chrom"].removeprefix("chr").upper()
    return (row["analysis_partition"] == "mitochondrial_bystander_control"
            or normalized in {"M", "MT", "25"})


def track_contract(
    path: Path, root: Path, inventory_path: Path,
) -> tuple[list[str], list[str]]:
    rows = read_tsv(path)
    included_signal = [
        row for row in rows
        if row.get("analysis_included", "").lower() == "true"
        and row.get("channel") in {"h3k4me3", "input"}
    ]
    included_keys = [
        (row["series_id"], row["condition"], row["channel"])
        for row in included_signal
    ]
    if len(included_keys) != len(set(included_keys)):
        raise GenomeSignalError("track manifest duplicates an included signal key")
    included_series = sorted({row["series_id"] for row in included_signal})
    if not included_series:
        raise GenomeSignalError("track manifest has no included signal series")
    for series_id in included_series:
        observed = {
            (row["condition"], row["channel"])
            for row in included_signal if row["series_id"] == series_id
        }
        expected = {
            (condition, channel) for condition in CONDITIONS
            for channel in ("h3k4me3", "input")
        }
        if observed != expected:
            raise GenomeSignalError(
                f"included signal series is not factorial: {series_id}"
            )
    previous: dict[tuple[str, str, str], dict[str, str]] = {}
    if inventory_path.is_file():
        for row in read_tsv(inventory_path):
            key = (row["series_id"], row["condition"], row["channel"])
            if key in previous:
                raise GenomeSignalError(
                    f"track file inventory duplicates key: {key}"
                )
            previous[key] = row
    inventory: list[dict[str, Any]] = []
    for row in included_signal:
        source = (root / row["filename"]).resolve()
        try:
            source.relative_to(root)
        except ValueError as error:
            raise GenomeSignalError(
                f"track path escapes its root: {row['filename']}"
            ) from error
        if not source.is_file():
            raise GenomeSignalError(f"included signal track is missing: {source}")
        stat = source.stat()
        key = (row["series_id"], row["condition"], row["channel"])
        cached = previous.get(key)
        digest = (
            cached["sha256"]
            if cached is not None
            and cached.get("resolved_source") == str(source)
            and int(cached.get("bytes", -1)) == stat.st_size
            and int(cached.get("mtime_ns", -1)) == stat.st_mtime_ns
            and re.fullmatch(r"[0-9a-f]{64}", cached.get("sha256", ""))
            else sha256(source)
        )
        inventory.append({
            "series_id": row["series_id"],
            "condition": row["condition"],
            "channel": row["channel"],
            "filename": row["filename"],
            "resolved_source": str(source),
            "bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "sha256": digest,
        })
    excluded_series = sorted({
        row["series_id"] for row in rows
        if row.get("analysis_included", "").lower() != "true"
    })
    inventory.sort(key=lambda row: (
        row["series_id"], CONDITIONS.index(row["condition"]), row["channel"]
    ))
    write_rows(
        inventory_path, inventory,
        (
            "series_id", "condition", "channel", "filename",
            "resolved_source", "bytes", "mtime_ns", "sha256",
        ),
    )
    return included_series, excluded_series


def evidence_rows(package: Path) -> list[dict[str, str]]:
    manifest = load_json(package / "manifest.json")
    if manifest.get("state") != "complete":
        raise GenomeSignalError("TP73 evidence package is not complete")
    rows = read_tsv(package / "chromosome_file_inventory.tsv")
    required = {
        "sequence_order", "chrom", "chrom_length", "analysis_partition",
        "primary_inference", "relative_path", "bytes", "sha256",
    }
    if required - set(rows[0]):
        raise GenomeSignalError("TP73 evidence inventory lacks required columns")
    seen: set[str] = set()
    selected: list[dict[str, str]] = []
    for row in rows:
        chrom = row["chrom"]
        if chrom in seen:
            raise GenomeSignalError(f"duplicate evidence chromosome: {chrom}")
        seen.add(chrom)
        relative = Path(row["relative_path"])
        source = (package / relative).resolve()
        try:
            source.relative_to(package)
        except ValueError as error:
            raise GenomeSignalError(
                f"evidence path escapes its package: {relative}"
            ) from error
        if not source.is_file() or source.stat().st_size != int(row["bytes"]):
            raise GenomeSignalError(f"evidence file differs: {source}")
        row = dict(row)
        row["evidence_path"] = str(source)
        if not is_mitochondrial(row):
            selected.append(row)
    if sum(row["analysis_partition"] == "autosome" for row in selected) != 22:
        raise GenomeSignalError("evidence package does not contain 22 autosomes")
    if sum(row["analysis_partition"] == "sex_chromosome" for row in selected) != 2:
        raise GenomeSignalError("evidence package does not contain X and Y")
    if len(selected) != 24:
        raise GenomeSignalError(
            "evidence package contains an unexpected nuclear chromosome"
        )
    return sorted(selected, key=lambda row: int(row["sequence_order"]))


def prepare(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    source = arguments.source.expanduser().resolve()
    evidence_package = arguments.evidence_package.expanduser().resolve()
    track_root = arguments.track_root.expanduser().resolve()
    runtime_prefix = arguments.runtime_prefix.expanduser().resolve()
    track_manifest = (
        arguments.track_manifest.expanduser().resolve()
        if arguments.track_manifest else
        source / "config" / "h3k4me3_cutandrun_tracks_v1.tsv"
    )
    for path, label in (
        (source, "source"), (evidence_package, "evidence package"),
        (track_root, "track root"), (runtime_prefix, "runtime prefix"),
    ):
        if not path.is_dir():
            raise GenomeSignalError(f"{label} is missing: {path}")
    if not track_manifest.is_file():
        raise GenomeSignalError(f"track manifest is missing: {track_manifest}")
    duckdb = runtime_prefix / "duckdb" / "bin" / "duckdb"
    bigwig_python = runtime_prefix / "duckdb" / "bin" / "python3"
    if not duckdb.is_file() or not bigwig_python.is_file():
        raise GenomeSignalError("runtime lacks DuckDB CLI or pyBigWig Python")
    for name in ("plan", "logs", "chromosomes", "final"):
        (run_root / name).mkdir(parents=True, exist_ok=True)
    track_inventory_path = run_root / "plan" / "track_file_inventory.tsv"
    included_series, excluded_series = track_contract(
        track_manifest, track_root, track_inventory_path
    )
    rows = evidence_rows(evidence_package)

    content = io.StringIO()
    fields = (
        "task_index", "sequence_order", "chrom", "chrom_length",
        "analysis_partition", "primary_inference", "evidence_path",
        "evidence_bytes", "evidence_sha256",
    )
    writer = csv.DictWriter(content, fieldnames=fields, delimiter="\t",
                            lineterminator="\n")
    writer.writeheader()
    partitions: dict[str, int] = {}
    for task_index, row in enumerate(rows):
        partition = row["analysis_partition"]
        partitions[partition] = partitions.get(partition, 0) + 1
        writer.writerow({
            "task_index": task_index,
            "sequence_order": row["sequence_order"],
            "chrom": row["chrom"],
            "chrom_length": row["chrom_length"],
            "analysis_partition": partition,
            "primary_inference": row["primary_inference"],
            "evidence_path": row["evidence_path"],
            "evidence_bytes": row["bytes"],
            "evidence_sha256": row["sha256"],
        })
    task_path = run_root / "plan" / "chromosome_tasks.tsv"
    immutable_write(task_path, content.getvalue())
    commit, dirty = git_identity(source)
    config = {
        "schema_version": 1,
        "run_id": arguments.run_id,
        "analysis": "whole_genome_h3k4me3_gfp_referenced_tp73_anchor_signal",
        "source": str(source),
        "source_commit": commit,
        "source_dirty": dirty,
        "scientific_source_file_sha256": scientific_hashes(source),
        "evidence_package": str(evidence_package),
        "evidence_manifest_sha256": sha256(evidence_package / "manifest.json"),
        "evidence_inventory_sha256": sha256(
            evidence_package / "chromosome_file_inventory.tsv"
        ),
        "track_root": str(track_root),
        "track_manifest": str(track_manifest),
        "track_manifest_sha256": sha256(track_manifest),
        "track_file_inventory": str(track_inventory_path),
        "track_file_inventory_sha256": sha256(track_inventory_path),
        "runtime_prefix": str(runtime_prefix),
        "scratch_root": str(arguments.scratch_root),
        "task_count": len(rows),
        "partition_counts": partitions,
        "primary_inference_partition": "autosome",
        "sex_chromosome_policy": "separate_sensitivity",
        "mitochondrial_policy": "excluded_from_histone_mark_inference",
        "included_series": included_series,
        "excluded_series": excluded_series,
        "minimum_anchor_score": -1.0,
        "windows": [
            {"window_name": "motif_span", "geometry": "anchor_span"},
            *[
                {"window_name": name, "inner_bp": inner, "outer_bp": outer}
                for name, inner, outer in WINDOWS
            ],
        ],
        "primary_window": PRIMARY_WINDOW,
        "normalization": "log2((h3k4me3_area+alpha)/(input_area+alpha))",
        "pseudocount_alpha": PSEUDOCOUNT,
        "primary_changes": ["TA-GFP", "DN-GFP"],
        "empty_chromosome_track_policy": "explicit_zero_signal",
        "threads": arguments.threads,
        "memory_limit": arguments.memory_limit,
        "minimum_free_run_gb": arguments.minimum_free_run_gb,
        "minimum_free_scratch_gb": arguments.minimum_free_scratch_gb,
        "task_plan_sha256": sha256(task_path),
    }
    immutable_write(
        run_root / "plan" / "run_config.json",
        json.dumps(config, indent=2, sort_keys=True) + "\n",
    )
    print(len(rows))


def planned_tasks(run_root: Path) -> list[dict[str, str]]:
    rows = read_tsv(run_root / "plan" / "chromosome_tasks.tsv")
    if [int(row["task_index"]) for row in rows] != list(range(len(rows))):
        raise GenomeSignalError("task indices are not contiguous")
    if len({row["chrom"] for row in rows}) != len(rows):
        raise GenomeSignalError("task plan contains duplicate chromosomes")
    return rows


def chromosome_root(run_root: Path, row: dict[str, str]) -> Path:
    return run_root / "chromosomes" / f"chrom-{safe_label(row['chrom'])}"


def build_contract(config: dict[str, Any], row: dict[str, str]) -> dict[str, Any]:
    return {
        "format_version": 1,
        "analysis": "h3k4me3_chromosome_anchor_signal",
        "run_id": config["run_id"],
        "source_commit": config["source_commit"],
        "evidence_package": config["evidence_package"],
        "evidence_manifest_sha256": config["evidence_manifest_sha256"],
        "evidence_path": row["evidence_path"],
        "evidence_bytes": int(row["evidence_bytes"]),
        "evidence_sha256": row["evidence_sha256"],
        "track_manifest": config["track_manifest"],
        "track_manifest_sha256": config["track_manifest_sha256"],
        "track_file_inventory": config["track_file_inventory"],
        "track_file_inventory_sha256": config["track_file_inventory_sha256"],
        "track_root": config["track_root"],
        "chrom": row["chrom"],
        "chrom_length": int(row["chrom_length"]),
        "analysis_partition": row["analysis_partition"],
        "primary_inference": row["primary_inference"] == "true",
        "included_series": config["included_series"],
        "excluded_series": config["excluded_series"],
        "windows": config["windows"],
        "primary_window": config["primary_window"],
        "normalization": config["normalization"],
        "pseudocount_alpha": config["pseudocount_alpha"],
        "empty_chromosome_track_policy": config["empty_chromosome_track_policy"],
    }


def canonicalize_provenance(attempt: Path, final: Path,
                            scratch_evidence: Path,
                            durable_evidence: Path) -> None:
    for path in sorted(item for item in attempt.rglob("*")
                       if item.is_file() and item.suffix in {".json", ".tsv"}):
        content = path.read_text(encoding="utf-8")
        updated = content.replace(str(attempt), str(final))
        updated = updated.replace(str(scratch_evidence), str(durable_evidence))
        if updated != content:
            temporary = path.with_name(f".{path.name}.canonicalizing")
            temporary.write_text(updated, encoding="utf-8")
            os.replace(temporary, path)


def build_change_table(duckdb: Path, signal_path: Path, output: Path,
                       threads: int, memory_limit: str,
                       temp_directory: Path) -> None:
    run([str(duckdb), "-light-mode", ":memory:", "-c", f"""
SET threads={threads};
SET memory_limit={sql_string(memory_limit)};
SET temp_directory={sql_string(temp_directory)};
SET preserve_insertion_order=false;
COPY (
  WITH primary_signal AS (
    SELECT * FROM read_parquet({sql_string(signal_path)}, hive_partitioning=false)
    WHERE window_name = {sql_string(PRIMARY_WINDOW)}
  ), pivoted AS (
    SELECT chrom, anchor_start, anchor_end, anchor_score,
           series_id, cell_line, replicate, window_name,
           min(inner_bp)::BIGINT AS inner_bp,
           max(outer_bp)::BIGINT AS outer_bp,
           max(segment_count)::INTEGER AS segment_count,
           max(effective_window_bp)::BIGINT AS effective_window_bp,
           max(h3k4me3_area) FILTER (condition = 'GFP') AS gfp_h3k4me3_area,
           max(input_area) FILTER (condition = 'GFP') AS gfp_input_area,
           max(h3k4me3_max) FILTER (condition = 'GFP') AS gfp_h3k4me3_max,
           max(input_max) FILTER (condition = 'GFP') AS gfp_input_max,
           max(h3k4me3_area) FILTER (condition = 'TA') AS ta_h3k4me3_area,
           max(input_area) FILTER (condition = 'TA') AS ta_input_area,
           max(h3k4me3_max) FILTER (condition = 'TA') AS ta_h3k4me3_max,
           max(input_max) FILTER (condition = 'TA') AS ta_input_max,
           max(h3k4me3_area) FILTER (condition = 'DN') AS dn_h3k4me3_area,
           max(input_area) FILTER (condition = 'DN') AS dn_input_area,
           max(h3k4me3_max) FILTER (condition = 'DN') AS dn_h3k4me3_max,
           max(input_max) FILTER (condition = 'DN') AS dn_input_max
    FROM primary_signal
    GROUP BY chrom, anchor_start, anchor_end, anchor_score,
             series_id, cell_line, replicate, window_name
    HAVING count(*) = 3
  ), normalized AS (
    SELECT *,
      log2((gfp_h3k4me3_area + {PSEUDOCOUNT}) /
           (gfp_input_area + {PSEUDOCOUNT})) AS gfp_log2_h3k4me3_input_ratio,
      log2((ta_h3k4me3_area + {PSEUDOCOUNT}) /
           (ta_input_area + {PSEUDOCOUNT})) AS ta_log2_h3k4me3_input_ratio,
      log2((dn_h3k4me3_area + {PSEUDOCOUNT}) /
           (dn_input_area + {PSEUDOCOUNT})) AS dn_log2_h3k4me3_input_ratio
    FROM pivoted
  )
  SELECT *,
    ta_log2_h3k4me3_input_ratio - gfp_log2_h3k4me3_input_ratio
      AS delta_ta_vs_gfp,
    dn_log2_h3k4me3_input_ratio - gfp_log2_h3k4me3_input_ratio
      AS delta_dn_vs_gfp,
    greatest(gfp_h3k4me3_area, ta_h3k4me3_area, dn_h3k4me3_area) > 0
      AS has_any_h3k4me3_signal,
    greatest(gfp_input_area, ta_input_area, dn_input_area) > 0
      AS has_any_input_signal,
    greatest(gfp_h3k4me3_area, ta_h3k4me3_area, dn_h3k4me3_area) = 0
      AS uninformative_all_h3k4me3_zero
  FROM normalized
  ORDER BY chrom, anchor_start, anchor_end, series_id
) TO {sql_string(output)}
  (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 131072);
"""])


def validate_outputs(duckdb: Path, row: dict[str, str], evidence: Path,
                     signal_path: Path, change_path: Path,
                     included_series: list[str]) -> dict[str, int]:
    expected_signal_factor = len(included_series) * len(CONDITIONS) * (
        len(WINDOWS) + 1
    )
    expected_change_factor = len(included_series)
    expected_series = ",".join(sql_string(value) for value in included_series)
    values = query_json(duckdb, ":memory:", f"""
WITH evidence AS (
  SELECT chrom, anchor_start, anchor_end, anchor_score
  FROM read_parquet({sql_string(evidence)}, hive_partitioning=false)
), signal AS (
  SELECT * FROM read_parquet({sql_string(signal_path)}, hive_partitioning=false)
), change AS (
  SELECT * FROM read_parquet({sql_string(change_path)}, hive_partitioning=false)
)
SELECT
  (SELECT count(*) FROM evidence)::BIGINT AS anchors,
  (SELECT count(*) - count(DISTINCT (chrom, anchor_start, anchor_end))
     FROM evidence)::BIGINT AS duplicate_evidence_keys,
  (SELECT count(*) FROM signal)::BIGINT AS signal_rows,
  (SELECT count(*) - count(DISTINCT (
      chrom, anchor_start, anchor_end, series_id, condition, window_name))
     FROM signal)::BIGINT AS duplicate_signal_keys,
  (SELECT count(*) FROM change)::BIGINT AS change_rows,
  (SELECT count(*) - count(DISTINCT (
      chrom, anchor_start, anchor_end, series_id))
     FROM change)::BIGINT AS duplicate_change_keys,
  (SELECT count(*) FROM signal
    WHERE chrom <> {sql_string(row['chrom'])}
       OR h3k4me3_area IS NULL OR input_area IS NULL
       OR h3k4me3_area < 0 OR input_area < 0
       OR NOT isfinite(h3k4me3_area) OR NOT isfinite(input_area)
       OR series_id NOT IN ({expected_series}))::BIGINT AS invalid_signal_rows,
  (SELECT count(*) FROM change
    WHERE chrom <> {sql_string(row['chrom'])}
       OR NOT isfinite(delta_ta_vs_gfp)
       OR NOT isfinite(delta_dn_vs_gfp)
       OR window_name <> {sql_string(PRIMARY_WINDOW)}
       OR series_id NOT IN ({expected_series}))::BIGINT AS invalid_change_rows,
  (SELECT count(*) FROM (
     SELECT chrom, anchor_start, anchor_end FROM evidence
     EXCEPT
     SELECT chrom, anchor_start, anchor_end FROM signal
   ))::BIGINT AS evidence_missing_from_signal,
  (SELECT count(*) FROM (
     SELECT chrom, anchor_start, anchor_end FROM signal
     EXCEPT
     SELECT chrom, anchor_start, anchor_end FROM evidence
   ))::BIGINT AS signal_missing_from_evidence;
""")
    if len(values) != 1:
        raise GenomeSignalError("signal validation returned no summary")
    result = {key: int(value) for key, value in values[0].items()}
    if (result["anchors"] <= 0
            or result["duplicate_evidence_keys"] != 0
            or result["signal_rows"] != result["anchors"] * expected_signal_factor
            or result["duplicate_signal_keys"] != 0
            or result["change_rows"] != result["anchors"] * expected_change_factor
            or result["duplicate_change_keys"] != 0
            or result["invalid_signal_rows"] != 0
            or result["invalid_change_rows"] != 0
            or result["evidence_missing_from_signal"] != 0
            or result["signal_missing_from_evidence"] != 0):
        raise GenomeSignalError(f"chromosome signal validation failed: {result}")
    return result


def output_inventory(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path.name == "manifest.json":
            continue
        rows.append({
            "relative_path": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    return rows


def run_task(arguments: argparse.Namespace) -> None:
    global CURRENT_ATTEMPT, CURRENT_CHROM, CURRENT_SCRATCH
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(run_root / "plan" / "run_config.json")
    verify_run_inputs(run_root, config)
    rows = planned_tasks(run_root)
    task_index = arguments.task_index
    if task_index is None:
        value = os.environ.get("SLURM_ARRAY_TASK_ID", "")
        if not value.isdigit():
            raise GenomeSignalError(
                "--task-index or SLURM_ARRAY_TASK_ID is required"
            )
        task_index = int(value)
    if task_index < 0 or task_index >= len(rows):
        raise GenomeSignalError(f"task index is outside the plan: {task_index}")
    verify_scientific_hashes(config)
    row = rows[task_index]
    CURRENT_CHROM = row["chrom"]
    signal.signal(signal.SIGUSR1, progress)

    evidence = Path(row["evidence_path"])
    if (not evidence.is_file()
            or evidence.stat().st_size != int(row["evidence_bytes"])
            or sha256(evidence) != row["evidence_sha256"]):
        raise GenomeSignalError(f"planned TP73 evidence differs: {evidence}")
    contract = build_contract(config, row)
    fingerprint = contract_sha256(contract)
    root = chromosome_root(run_root, row)
    final = root / "final"
    if final.is_dir():
        manifest = load_json(final / "manifest.json")
        if (manifest.get("state") == "complete"
                and manifest.get("contract_sha256") == fingerprint):
            print(f"I: reusing completed chromosome package: {final}",
                  file=sys.stderr)
            return
        raise GenomeSignalError("final chromosome package has another contract")

    root.mkdir(parents=True, exist_ok=True)
    (root / "attempts").mkdir(exist_ok=True)
    check_free_space(root, float(config["minimum_free_run_gb"]), "durable")
    job_id = os.environ.get("SLURM_JOB_ID", "manual")
    restart = os.environ.get("SLURM_RESTART_COUNT", "0")
    attempt = root / "attempts" / f"job-{job_id}-restart-{restart}"
    if attempt.exists():
        attempt = attempt.with_name(f"{attempt.name}-pid-{os.getpid()}")
    attempt.mkdir()
    CURRENT_ATTEMPT = attempt

    scratch_root = Path(str(config["scratch_root"])).expanduser().resolve()
    scratch_root.mkdir(parents=True, exist_ok=True)
    check_free_space(
        scratch_root, float(config["minimum_free_scratch_gb"]), "scratch"
    )
    scratch = scratch_root / f"h3-signal-{job_id}-{restart}-{os.getpid()}"
    scratch.mkdir()
    CURRENT_SCRATCH = scratch
    scratch_evidence = scratch / f"tp73_evidence_chr{safe_label(row['chrom'])}.parquet"

    set_phase("copy-evidence-to-scratch")
    shutil.copy2(evidence, scratch_evidence)
    signal_path = attempt / "h3k4me3_anchor_signal.parquet"
    set_phase("extract-windowed-h3k4me3-input")
    runtime = Path(str(config["runtime_prefix"]))
    duckdb = runtime / "duckdb" / "bin" / "duckdb"
    bigwig_python = runtime / "duckdb" / "bin" / "python3"
    command = [
        str(bigwig_python),
        str(Path(str(config["source"])) /
            "scripts" / "build_h3k4me3_anchor_signal.py"),
        "--tp73-evidence-input", str(scratch_evidence),
        "--track-manifest", str(config["track_manifest"]),
        "--track-file-inventory", str(config["track_file_inventory"]),
        "--track-root", str(config["track_root"]),
        "--signal-output", str(signal_path),
        "--chrom", row["chrom"],
        "--chrom-length", row["chrom_length"],
        "--minimum-anchor-score", str(config["minimum_anchor_score"]),
        "--threads", str(config["threads"]),
        "--memory-limit", str(config["memory_limit"]),
        "--duckdb", str(duckdb),
        "--rscript", str(bigwig_python),
        "--bigwig-exporter", str(
            Path(str(config["source"])) /
            "scripts" / "export_bigwig_chrom_bedgraph.py"
        ),
        "--tmpdir", str(scratch),
    ]
    for name, inner, outer in WINDOWS:
        command.extend(["--window", f"{name}:{inner}:{outer}"])
    run(command)

    set_phase("derive-gfp-referenced-change")
    change_path = attempt / "h3k4me3_anchor_change.parquet"
    (scratch / "duckdb_tmp").mkdir(exist_ok=True)
    build_change_table(
        duckdb, signal_path, change_path, int(config["threads"]),
        str(config["memory_limit"]), scratch / "duckdb_tmp",
    )
    set_phase("canonicalize-provenance")
    canonicalize_provenance(attempt, final, scratch_evidence, evidence)
    set_phase("validate-and-publish")
    validation = validate_outputs(
        duckdb, row, evidence, signal_path, change_path,
        list(config["included_series"]),
    )
    manifest = {
        **contract,
        "contract_sha256": fingerprint,
        "state": "complete",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_restart_count": int(os.environ.get("SLURM_RESTART_COUNT", "0")),
        "validation": validation,
        "outputs": output_inventory(attempt),
    }
    (attempt / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(attempt, final)
    CURRENT_ATTEMPT = final
    set_phase("complete")


def validate_child(row: dict[str, str], root: Path,
                   verify_checksum: bool = False) -> dict[str, Any]:
    manifest = load_json(root / "final" / "manifest.json")
    if (manifest.get("state") != "complete"
            or str(manifest.get("chrom")) != row["chrom"]
            or int(manifest.get("chrom_length", -1)) != int(row["chrom_length"])
            or manifest.get("analysis_partition") != row["analysis_partition"]):
        raise GenomeSignalError(f"chromosome package identity differs: {root}")
    records = {
        str(record.get("relative_path")): record
        for record in manifest.get("outputs", []) if isinstance(record, dict)
    }
    for relative in (
        "h3k4me3_anchor_signal.parquet", "h3k4me3_anchor_change.parquet"
    ):
        path = root / "final" / relative
        record = records.get(relative)
        if (record is None or not path.is_file()
                or path.stat().st_size != int(record.get("bytes", -1))):
            raise GenomeSignalError(f"chromosome output differs: {path}")
        if verify_checksum and sha256(path) != record.get("sha256"):
            raise GenomeSignalError(f"chromosome checksum differs: {path}")
    validation = manifest.get("validation", {})
    anchors = int(validation.get("anchors", 0))
    if (anchors <= 0
            or int(validation.get("signal_rows", -1)) <= 0
            or int(validation.get("change_rows", -1)) <= 0):
        raise GenomeSignalError(f"chromosome validation is incomplete: {root}")
    return manifest


def status(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    rows = planned_tasks(run_root)
    complete = invalid = 0
    partitions: dict[str, int] = {}
    for row in rows:
        root = chromosome_root(run_root, row)
        if not (root / "final" / "manifest.json").is_file():
            continue
        try:
            validate_child(row, root)
            complete += 1
            key = row["analysis_partition"]
            partitions[key] = partitions.get(key, 0) + 1
        except (GenomeSignalError, OSError, ValueError, json.JSONDecodeError):
            invalid += 1
    print("key\tvalue")
    print(f"planned\t{len(rows)}")
    print(f"complete\t{complete}")
    print(f"pending\t{len(rows) - complete - invalid}")
    print(f"invalid_complete_markers\t{invalid}")
    for key, value in sorted(partitions.items()):
        print(f"complete_{key}\t{value}")


def link_or_copy(source: Path, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(source, target)
    except OSError as error:
        if error.errno != errno.EXDEV:
            raise
        shutil.copy2(source, target)


def finalize(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(run_root / "plan" / "run_config.json")
    verify_run_inputs(run_root, config)
    rows = planned_tasks(run_root)
    if int(config.get("task_count", -1)) != len(rows):
        raise GenomeSignalError("run configuration task count differs")
    verify_scientific_hashes(config)
    manifests = [
        validate_child(row, chromosome_root(run_root, row), verify_checksum=True)
        for row in rows
    ]
    final = run_root / "final" / "genome_h3k4me3_signal"
    if final.exists():
        manifest = load_json(final / "manifest.json")
        if (manifest.get("state") == "complete"
                and manifest.get("run_id") == config.get("run_id")):
            print(f"I: reusing finalized H3K4me3 signal: {final}",
                  file=sys.stderr)
            return
        raise GenomeSignalError(f"final output has another identity: {final}")

    final.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".h3-genome-signal.", dir=final.parent))
    try:
        inventory: list[dict[str, Any]] = []
        dataset_paths: dict[str, list[Path]] = {
            "h3k4me3_anchor_signal": [], "h3k4me3_anchor_change": [],
        }
        partition_paths: dict[tuple[str, str], list[Path]] = {}
        for row in rows:
            partition = row["analysis_partition"]
            for dataset, filename in (
                ("h3k4me3_anchor_signal", "h3k4me3_anchor_signal.parquet"),
                ("h3k4me3_anchor_change", "h3k4me3_anchor_change.parquet"),
            ):
                source = chromosome_root(run_root, row) / "final" / filename
                target = (
                    staging / "tables" / dataset /
                    f"analysis_partition={partition}" /
                    f"chrom={safe_label(row['chrom'])}" / "part-000000.parquet"
                )
                link_or_copy(source, target)
                relative = target.relative_to(staging)
                dataset_paths[dataset].append(relative)
                partition_paths.setdefault((dataset, partition), []).append(relative)
                inventory.append({
                    "dataset": dataset,
                    "sequence_order": int(row["sequence_order"]),
                    "chrom": row["chrom"],
                    "chrom_length": int(row["chrom_length"]),
                    "analysis_partition": partition,
                    "primary_inference": row["primary_inference"],
                    "relative_path": str(relative),
                    "bytes": source.stat().st_size,
                    "sha256": sha256(source),
                })
        inventory_path = staging / "chromosome_file_inventory.tsv"
        with inventory_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=tuple(inventory[0]), delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(inventory)

        provenance = staging / "provenance"
        provenance.mkdir()
        shutil.copy2(run_root / "plan" / "run_config.json", provenance)
        shutil.copy2(run_root / "plan" / "chromosome_tasks.tsv", provenance)
        shutil.copy2(
            run_root / "plan" / "track_file_inventory.tsv", provenance
        )
        shutil.copy2(
            Path(str(config["track_manifest"])),
            provenance / "track_manifest.tsv",
        )

        duckdb = arguments.duckdb.expanduser().resolve()
        summary = staging / "h3k4me3_change_summary_by_chromosome.parquet"
        change_paths = dataset_paths["h3k4me3_anchor_change"]
        run([str(duckdb), "-light-mode", ":memory:", "-c", f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET temp_directory={sql_string(arguments.temp_directory)};
COPY (
  SELECT chrom, series_id, cell_line, count(*)::BIGINT AS anchors,
         count(*) FILTER (WHERE has_any_h3k4me3_signal)::BIGINT
           AS anchors_with_h3k4me3_signal,
         count(*) FILTER (WHERE uninformative_all_h3k4me3_zero)::BIGINT
           AS uninformative_all_h3k4me3_zero,
         avg(delta_ta_vs_gfp)::DOUBLE AS mean_delta_ta_vs_gfp,
         median(delta_ta_vs_gfp)::DOUBLE AS median_delta_ta_vs_gfp,
         avg(delta_dn_vs_gfp)::DOUBLE AS mean_delta_dn_vs_gfp,
         median(delta_dn_vs_gfp)::DOUBLE AS median_delta_dn_vs_gfp
  FROM read_parquet({sql_path_list(change_paths)}, hive_partitioning=false)
  GROUP BY chrom, series_id, cell_line
  ORDER BY chrom, series_id
) TO {sql_string(summary)} (FORMAT PARQUET, COMPRESSION ZSTD);
"""], cwd=staging)

        database = staging / "genome_h3k4me3_signal.duckdb"
        statements = [
            "CREATE VIEW chromosome_file_inventory AS SELECT * FROM "
            "read_csv_auto('chromosome_file_inventory.tsv', delim='\\t', "
            "header=true);",
            "CREATE VIEW h3k4me3_change_summary_by_chromosome AS SELECT * FROM "
            "read_parquet('h3k4me3_change_summary_by_chromosome.parquet');",
        ]
        for dataset, paths in dataset_paths.items():
            statements.append(
                f"CREATE VIEW {dataset} AS SELECT * FROM "
                f"read_parquet({sql_path_list(paths)}, hive_partitioning=false);"
            )
        for (dataset, partition), paths in sorted(partition_paths.items()):
            statements.append(
                f"CREATE VIEW {dataset}_{safe_label(partition)} AS SELECT * FROM "
                f"read_parquet({sql_path_list(paths)}, hive_partitioning=false);"
            )
        run([str(duckdb), str(database), "-c", "\n".join(statements)], cwd=staging)

        totals = {
            "chromosomes": len(rows),
            "anchors": sum(int(item["validation"]["anchors"])
                           for item in manifests),
            "signal_rows": sum(int(item["validation"]["signal_rows"])
                               for item in manifests),
            "change_rows": sum(int(item["validation"]["change_rows"])
                               for item in manifests),
        }
        expected_signal_factor = (
            len(config["included_series"]) * len(CONDITIONS) * (len(WINDOWS) + 1)
        )
        expected_change_factor = len(config["included_series"])
        if (totals["chromosomes"] != 24 or totals["anchors"] <= 0
                or totals["signal_rows"] !=
                totals["anchors"] * expected_signal_factor
                or totals["change_rows"] !=
                totals["anchors"] * expected_change_factor):
            raise GenomeSignalError(f"genome signal validation failed: {totals}")
        manifest = {
            "schema_version": 1,
            "state": "complete",
            "run_id": config["run_id"],
            "analysis": config["analysis"],
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "source_commit": config["source_commit"],
            "partition_counts": config["partition_counts"],
            "primary_inference_partition": config["primary_inference_partition"],
            "sex_chromosome_policy": config["sex_chromosome_policy"],
            "mitochondrial_policy": config["mitochondrial_policy"],
            "primary_window": config["primary_window"],
            "normalization": config["normalization"],
            "pseudocount_alpha": config["pseudocount_alpha"],
            "validation": totals,
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
    print(f"I: finalized genome-wide H3K4me3 signal: {final}", file=sys.stderr)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser(
        "prepare", help="write an immutable nuclear-chromosome task plan"
    )
    prepare_parser.add_argument("--run-root", type=Path, required=True)
    prepare_parser.add_argument("--evidence-package", type=Path, required=True)
    prepare_parser.add_argument("--track-root", type=Path, required=True)
    prepare_parser.add_argument("--runtime-prefix", type=Path, required=True)
    prepare_parser.add_argument("--source", type=Path, required=True)
    prepare_parser.add_argument("--track-manifest", type=Path)
    prepare_parser.add_argument("--scratch-root", type=Path, required=True)
    prepare_parser.add_argument("--run-id", required=True)
    prepare_parser.add_argument("--threads", type=int, default=4)
    prepare_parser.add_argument("--memory-limit", default="28GB")
    prepare_parser.add_argument("--minimum-free-run-gb", type=float, default=5)
    prepare_parser.add_argument("--minimum-free-scratch-gb", type=float, default=10)
    prepare_parser.set_defaults(function=prepare)

    task_parser = subparsers.add_parser(
        "run-task", help="execute one restart-safe chromosome measurement"
    )
    task_parser.add_argument("--run-root", type=Path, required=True)
    task_parser.add_argument("--task-index", type=int)
    task_parser.set_defaults(function=run_task)

    status_parser = subparsers.add_parser(
        "status", help="report exact chromosome completion counts"
    )
    status_parser.add_argument("--run-root", type=Path, required=True)
    status_parser.set_defaults(function=status)

    finalize_parser = subparsers.add_parser(
        "finalize", help="validate and publish the genome signal catalog"
    )
    finalize_parser.add_argument("--run-root", type=Path, required=True)
    finalize_parser.add_argument("--duckdb", type=Path, required=True)
    finalize_parser.add_argument("--threads", type=int, default=4)
    finalize_parser.add_argument("--memory-limit", default="28GB")
    finalize_parser.add_argument("--temp-directory", type=Path, required=True)
    finalize_parser.set_defaults(function=finalize)
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        if hasattr(arguments, "threads") and arguments.threads <= 0:
            raise GenomeSignalError("--threads must be positive")
        arguments.function(arguments)
        return 0
    except (GenomeSignalError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        if arguments.command == "run-task":
            progress()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

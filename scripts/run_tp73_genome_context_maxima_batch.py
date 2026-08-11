#!/usr/bin/env python3

"""Build and atomically publish one bounded batch of autosomal motif contexts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class ContextBatchError(RuntimeError):
    pass


STARTED = time.monotonic()
CURRENT_PHASE = "startup"
CURRENT_MOTIF = "none"
CURRENT_CHROM = "none"
CURRENT_CHILD_PID: int | None = None
CURRENT_BATCH = -1


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


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise ContextBatchError(f"JSON file is missing: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ContextBatchError(f"JSON value is not an object: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ContextBatchError(f"TSV file is missing: {path}")
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        rows = list(reader)
    if reader.fieldnames is None or not rows:
        raise ContextBatchError(f"TSV is empty or lacks a header: {path}")
    return rows


def progress(_signal_number: int | None = None,
             _frame: object | None = None) -> None:
    print(
        "I: progress signal=SIGUSR1 "
        f"phase={CURRENT_PHASE} batch={CURRENT_BATCH} motif={CURRENT_MOTIF} "
        f"chrom={CURRENT_CHROM} elapsed_seconds={time.monotonic() - STARTED:.1f} "
        f"child_pid={CURRENT_CHILD_PID if CURRENT_CHILD_PID is not None else 'none'}",
        file=sys.stderr, flush=True,
    )


def set_phase(name: str, *, motif: str | None = None,
              chrom: str | None = None) -> None:
    global CURRENT_PHASE, CURRENT_MOTIF, CURRENT_CHROM
    CURRENT_PHASE = name
    if motif is not None:
        CURRENT_MOTIF = motif
    if chrom is not None:
        CURRENT_CHROM = chrom
    print(
        f"I: phase={CURRENT_PHASE} batch={CURRENT_BATCH} "
        f"motif={CURRENT_MOTIF} chrom={CURRENT_CHROM}",
        file=sys.stderr, flush=True,
    )


def run_child(command: list[str]) -> None:
    global CURRENT_CHILD_PID
    displayed = list(command)
    if "-c" in displayed:
        sql_index = displayed.index("-c") + 1
        if sql_index < len(displayed) and len(displayed[sql_index]) > 200:
            displayed[sql_index] = f"<SQL:{len(displayed[sql_index])} chars>"
    print("I: running: " + " ".join(displayed), file=sys.stderr, flush=True)
    process = subprocess.Popen(command)
    CURRENT_CHILD_PID = process.pid
    try:
        while True:
            try:
                return_code = process.wait()
                break
            except InterruptedError:
                continue
    finally:
        CURRENT_CHILD_PID = None
    if return_code != 0:
        raise ContextBatchError(
            f"command failed with exit code {return_code}: {command[0]}"
        )


def run_json(duckdb: Path, query: str) -> list[dict[str, Any]]:
    process = subprocess.run(
        [str(duckdb), "-light-mode", "-json", ":memory:", "-c", query],
        text=True, capture_output=True, check=False,
    )
    if process.returncode != 0:
        raise ContextBatchError(process.stderr.strip() or "DuckDB query failed")
    value = json.loads(process.stdout or "[]")
    if not isinstance(value, list):
        raise ContextBatchError("DuckDB JSON result is not a row array")
    return value


def check_space(path: Path, minimum_gb: float, label: str) -> None:
    free = shutil.disk_usage(path).free
    required = int(minimum_gb * 1024**3)
    print(
        f"I: {label} free space {free / 1024**3:.1f} GiB; "
        f"required {minimum_gb:.1f} GiB",
        file=sys.stderr, flush=True,
    )
    if free < required:
        raise ContextBatchError(f"insufficient {label} space")


def verify_plan(run_root: Path, config: dict[str, Any]) -> None:
    paths = {
        "task_plan_sha256": run_root / "plan" / "calibration_tasks.tsv",
        "anchor_file_plan_sha256": run_root / "plan" / "anchor_files.tsv",
        "scan_file_plan_sha256": run_root / "plan" / "scan_files.tsv",
        "threshold_registry_sha256": Path(str(config["threshold_registry"])),
    }
    for key, path in paths.items():
        if not path.is_file() or sha256(path) != config.get(key):
            raise ContextBatchError(f"immutable planned input changed: {path}")
    source = Path(str(config["source"]))
    for relative, digest in config["scientific_source_file_sha256"].items():
        path = source / relative
        if not path.is_file() or sha256(path) != digest:
            raise ContextBatchError(f"scientific source changed: {path}")


def task_directory(run_root: Path, row: dict[str, str]) -> Path:
    return run_root / "tasks" / (
        f"task-{int(row['task_index']):06d}-{row['motif_id']}"
    )


def validate_existing(directory: Path, row: dict[str, str],
                      config: dict[str, Any]) -> bool:
    if not directory.exists():
        return False
    marker = load_json(directory / "complete.json")
    record = marker.get("files", {}).get("cofactor_maxima.parquet")
    maxima = directory / "cofactor_maxima.parquet"
    if (marker.get("schema_version") != 1
            or marker.get("task_index") != int(row["task_index"])
            or marker.get("motif_id") != row["motif_id"]
            or marker.get("threshold_registry_sha256")
               != config["threshold_registry_sha256"]
            or not isinstance(record, dict)
            or not maxima.is_file()
            or maxima.stat().st_size != int(record.get("bytes", -1))):
        raise ContextBatchError(f"existing task has another identity: {directory}")
    return True


def selected_scan_files(scan_plan: Path, motif_ids: set[str],
                        scan_package: Path) -> dict[tuple[str, str, str], Path]:
    result: dict[tuple[str, str, str], Path] = {}
    with scan_plan.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["motif_id"] not in motif_ids:
                continue
            relative = Path(row["relative_path"])
            if relative.is_absolute() or ".." in relative.parts:
                raise ContextBatchError(f"unsafe scan path: {relative}")
            path = scan_package / relative
            if not path.is_file() or path.stat().st_size != int(row["bytes"]):
                raise ContextBatchError(f"scan file is missing or size changed: {path}")
            key = (row["motif_id"], row["chrom"], row["strand"])
            if key in result:
                raise ContextBatchError(f"duplicate scan file in plan: {key}")
            result[key] = path
    expected = len(motif_ids) * 22 * 2
    if len(result) != expected:
        raise ContextBatchError(
            f"selected scan inventory has {len(result)} files; expected {expected}"
        )
    return result


def aggregate_motif(
    duckdb: Path,
    partials: list[Path],
    output: Path,
    temp_directory: Path,
    threads: int,
    memory_limit: str,
) -> None:
    run_child([
        str(duckdb), "-light-mode", "-batch", ":memory:", "-c", f"""
SET threads={threads};
SET memory_limit={sql_string(memory_limit)};
SET preserve_insertion_order=false;
SET temp_directory={sql_string(temp_directory)};
COPY (
  SELECT * FROM read_parquet({sql_path_list(partials)}, hive_partitioning=false)
  ORDER BY CAST(chrom AS INTEGER), anchor_start, anchor_end
) TO {sql_string(output)}
  (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 131072);
""",
    ])


def validate_output(duckdb: Path, output: Path, row: dict[str, str],
                    expected_anchors: int, threshold_set_id: str) -> dict[str, int]:
    values = run_json(duckdb, f"""
SELECT count(*)::BIGINT AS rows,
       count(DISTINCT (chrom, anchor_start, anchor_end))::BIGINT AS anchors,
       count(DISTINCT chrom)::BIGINT AS chromosomes,
       count(DISTINCT motif_id)::BIGINT AS motifs,
       count(*) - count(DISTINCT (chrom, anchor_start, anchor_end, motif_id))::BIGINT
           AS duplicate_keys,
       count(*) FILTER (WHERE motif_id <> {sql_string(row['motif_id'])})::BIGINT
           AS wrong_motif_rows,
       count(*) FILTER (WHERE context_score IS NOT NULL
                         AND context_score < source_score_floor)::BIGINT
           AS scores_below_source_floor,
       count(*) FILTER (WHERE n_neighbor_loci_above_threshold IS NULL
                         OR has_neighbor_locus_above_threshold IS NULL
                         OR anchor_locus_id IS NULL)::BIGINT AS missing_count_rows,
       count(*) FILTER (WHERE threshold_set_id <>
                         {sql_string(threshold_set_id)})::BIGINT
           AS wrong_threshold_rows,
       count(DISTINCT source_score_floor)::BIGINT AS source_floors,
       count(DISTINCT recommended_threshold)::BIGINT AS applied_thresholds
FROM read_parquet({sql_string(output)}, hive_partitioning=false);
""")
    if len(values) != 1:
        raise ContextBatchError("context validation returned no summary")
    result = {key: int(value) for key, value in values[0].items()}
    if (result["rows"] != expected_anchors
            or result["anchors"] != expected_anchors
            or result["chromosomes"] != 22 or result["motifs"] != 1
            or result["duplicate_keys"] != 0 or result["wrong_motif_rows"] != 0
            or result["scores_below_source_floor"] != 0
            or result["missing_count_rows"] != 0
            or result["wrong_threshold_rows"] != 0
            or result["source_floors"] != 1
            or result["applied_thresholds"] != 1):
        raise ContextBatchError(f"context output validation failed: {result}")
    return result


def publish_task(
    run_root: Path,
    row: dict[str, str],
    output: Path,
    validation: dict[str, int],
    config: dict[str, Any],
    batch_index: int,
) -> None:
    final = task_directory(run_root, row)
    attempt = run_root / "tasks" / (
        f".attempt-task-{int(row['task_index']):06d}-{row['motif_id']}-"
        f"job-{os.environ.get('SLURM_JOB_ID', 'manual')}-"
        f"restart-{os.environ.get('SLURM_RESTART_COUNT', '0')}-pid-{os.getpid()}"
    )
    if attempt.exists():
        raise ContextBatchError(f"durable attempt already exists: {attempt}")
    attempt.mkdir(parents=True)
    staged_output = attempt / "cofactor_maxima.parquet"
    shutil.copy2(output, staged_output)
    run_config = {
        "schema_version": 1,
        "analysis": config["analysis"],
        "analysis_scope": config["analysis_scope"],
        "task_index": int(row["task_index"]),
        "batch_index": batch_index,
        "motif_id": row["motif_id"],
        "motif_name": row["motif_name"],
        "motif_length": int(row["motif_length"]),
        "scan_minimum_score": float(row["scan_minimum_score"]),
        "source_recommended_threshold": (
            None if row["source_recommended_threshold"] == "NA"
            else float(row["source_recommended_threshold"])
        ),
        "applied_context_threshold": float(row["applied_context_threshold"]),
        "threshold_application_status": row["threshold_application_status"],
        "threshold_set_id": config["threshold_set_id"],
        "threshold_registry_sha256": config["threshold_registry_sha256"],
        "context_flank_bp": config["context_flank_bp"],
        "context_distance_metric": config["context_distance_metric"],
        "chromosomes": config["chromosomes"],
        "source_commit": config["source_commit"],
        "scan_manifest_sha256": config["scan_manifest_sha256"],
        "evidence_manifest_sha256": config["evidence_manifest_sha256"],
        "anchor_file_plan_sha256": config["anchor_file_plan_sha256"],
        "scan_file_plan_sha256": config["scan_file_plan_sha256"],
        "validation": validation,
    }
    run_config_path = attempt / "run_config.json"
    run_config_path.write_text(
        json.dumps(run_config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    marker = {
        **run_config,
        "state": "complete",
        "analysis_partition": "autosome",
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "slurm_array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "slurm_restart_count": int(os.environ.get("SLURM_RESTART_COUNT", "0")),
        "files": {
            "cofactor_maxima.parquet": {
                "bytes": staged_output.stat().st_size,
                "sha256": sha256(staged_output),
            },
            "run_config.json": {
                "bytes": run_config_path.stat().st_size,
                "sha256": sha256(run_config_path),
            },
        },
    }
    (attempt / "complete.json").write_text(
        json.dumps(marker, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if final.exists():
        raise ContextBatchError(f"task completed concurrently: {final}")
    os.replace(attempt, final)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Stage exact autosomal TP73 anchors once, process a bounded motif "
            "batch chromosome by chromosome, and publish each completed motif "
            "atomically. SIGUSR1 reports progress without stopping the child."
        )
    )
    result.add_argument("--run-root", type=Path, required=True)
    result.add_argument("--batch-index", type=int)
    result.add_argument(
        "--batch-offset", type=int, default=0,
        help="add this offset to SLURM_ARRAY_TASK_ID for chunked arrays",
    )
    return result


def main() -> int:
    global CURRENT_BATCH
    arguments = parser().parse_args()
    scratch: Path | None = None
    try:
        signal.signal(signal.SIGUSR1, progress)
        run_root = arguments.run_root.expanduser().resolve()
        config = load_json(run_root / "plan" / "run_config.json")
        verify_plan(run_root, config)
        rows = read_tsv(run_root / "plan" / "calibration_tasks.tsv")
        batch_index = arguments.batch_index
        if batch_index is None:
            value = os.environ.get("SLURM_ARRAY_TASK_ID", "")
            if not value.isdigit():
                raise ContextBatchError(
                    "--batch-index or SLURM_ARRAY_TASK_ID is required"
                )
            batch_index = int(value) + arguments.batch_offset
        elif arguments.batch_offset != 0:
            raise ContextBatchError(
                "--batch-offset cannot be combined with --batch-index"
            )
        if batch_index < 0 or batch_index >= int(config["batch_count"]):
            raise ContextBatchError(f"batch index is outside the plan: {batch_index}")
        CURRENT_BATCH = batch_index
        first = batch_index * int(config["motifs_per_batch"])
        selected = rows[first:first + int(config["motifs_per_batch"])]
        active = [
            row for row in selected
            if not validate_existing(task_directory(run_root, row), row, config)
        ]
        if not active:
            print(f"I: reusing completed context batch {batch_index}", file=sys.stderr)
            return 0

        source = Path(str(config["source"]))
        runtime = Path(str(config["runtime_prefix"]))
        duckdb = runtime / "duckdb" / "bin" / "duckdb"
        if not duckdb.is_file():
            raise ContextBatchError(f"DuckDB runtime is missing: {duckdb}")
        run_root.mkdir(parents=True, exist_ok=True)
        check_space(run_root, float(config["minimum_free_run_gb"]), "durable")
        scratch_root = Path(os.environ.get(
            "SLURM_TMPDIR", str(config["scratch_root"])
        )).expanduser().resolve()
        scratch_root.mkdir(parents=True, exist_ok=True)
        check_space(
            scratch_root, float(config["minimum_free_scratch_gb"]), "scratch"
        )
        scratch = scratch_root / (
            f"tp73-context-batch-{batch_index}-"
            f"job-{os.environ.get('SLURM_JOB_ID', 'manual')}-"
            f"restart-{os.environ.get('SLURM_RESTART_COUNT', '0')}-pid-{os.getpid()}"
        )
        scratch.mkdir()
        (scratch / "duckdb-tmp").mkdir()

        set_phase("stage-threshold-registry")
        registry = Path(str(config["threshold_registry"]))
        staged_registry = scratch / "applied_context_thresholds.parquet"
        shutil.copy2(registry, staged_registry)
        if sha256(staged_registry) != config["threshold_registry_sha256"]:
            raise ContextBatchError("staged threshold registry checksum differs")

        set_phase("stage-autosome-anchors")
        anchors: dict[str, Path] = {}
        for anchor_row in read_tsv(run_root / "plan" / "anchor_files.tsv"):
            source_anchor = Path(anchor_row["path"])
            if (not source_anchor.is_file()
                    or source_anchor.stat().st_size != int(anchor_row["bytes"])):
                raise ContextBatchError(f"anchor source differs: {source_anchor}")
            staged = scratch / f"anchor-chr{anchor_row['chrom']}.parquet"
            shutil.copy2(source_anchor, staged)
            if sha256(staged) != anchor_row["sha256"]:
                raise ContextBatchError(f"staged anchor checksum differs: {staged}")
            anchors[anchor_row["chrom"]] = staged
        if set(anchors) != set(config["chromosomes"]):
            raise ContextBatchError("staged anchor chromosome set differs")

        motif_ids = {row["motif_id"] for row in active}
        scan_files = selected_scan_files(
            run_root / "plan" / "scan_files.tsv", motif_ids,
            Path(str(config["scan_package"])),
        )
        for row in active:
            motif_id = row["motif_id"]
            motif_scratch = scratch / f"motif-{motif_id}"
            partial_root = motif_scratch / "partials"
            partial_root.mkdir(parents=True)
            partials: list[Path] = []
            for chrom in config["chromosomes"]:
                set_phase("build-chromosome-context", motif=motif_id, chrom=chrom)
                output = partial_root / f"chr{chrom}.parquet"
                run_child([
                    sys.executable,
                    str(source / "scripts" / "build_sparse_context_maxima.py"),
                    "--anchor-parquet", str(anchors[chrom]),
                    "--cofactor", motif_id,
                    str(scan_files[(motif_id, chrom, "+")]),
                    str(scan_files[(motif_id, chrom, "-")]),
                    "--cofactor-span", motif_id, row["motif_length"],
                    "--threshold-parquet", str(staged_registry),
                    "--threshold-set-id", str(config["threshold_set_id"]),
                    "--threshold-role", str(config["threshold_role"]),
                    "--calibration-stratum-id",
                    str(config["calibration_stratum_id"]),
                    "--output", str(output),
                    "--flank", str(config["context_flank_bp"]),
                    "--source-score-floor", row["scan_minimum_score"],
                    "--duckdb", str(duckdb),
                    "--threads", str(config["threads"]),
                    "--memory-limit", str(config["memory_limit"]),
                    "--max-temp-size", str(config["max_temp_size"]),
                    "--temp-directory", str(scratch / "duckdb-tmp"),
                ])
                partials.append(output)

            set_phase("aggregate-autosomes", motif=motif_id, chrom="all")
            aggregate = motif_scratch / "cofactor_maxima.parquet"
            aggregate_motif(
                duckdb, partials, aggregate, scratch / "duckdb-tmp",
                int(config["threads"]), str(config["memory_limit"]),
            )
            set_phase("validate-autosomes", motif=motif_id, chrom="all")
            validation = validate_output(
                duckdb, aggregate, row, int(config["anchor_count"]),
                str(config["threshold_set_id"]),
            )
            set_phase("publish", motif=motif_id, chrom="all")
            publish_task(
                run_root, row, aggregate, validation, config, batch_index
            )
            shutil.rmtree(motif_scratch)
            print(
                f"I: published context task {row['task_index']} ({motif_id})",
                file=sys.stderr, flush=True,
            )
        set_phase("complete", motif="none", chrom="none")
        return 0
    except (ContextBatchError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        progress()
        return 1
    finally:
        if scratch is not None and scratch.exists():
            shutil.rmtree(scratch)


if __name__ == "__main__":
    raise SystemExit(main())

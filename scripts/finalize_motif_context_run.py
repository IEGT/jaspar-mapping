#!/usr/bin/env python3

"""Validate and catalog a complete batched motif-context run.

The finalizer follows the immutable task plan. It never discovers completed
tasks from a package glob, never rewrites task packages, and never removes
orphan staging data. A completion manifest is published only after every
planned package has passed provenance and schema checks.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any


class ContextFinalizationError(RuntimeError):
    pass


SAFE_IDENTIFIER = re.compile(r"^[A-Za-z0-9._-]+$")
SHA256 = re.compile(r"^[0-9a-f]{64}$")
ANCHOR_MOTIF = "MA0861.2"


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_read(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ContextFinalizationError(f"cannot read JSON {path}: {error}") from error


def json_write(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def write_json_lines(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True) + "\n")


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def parse_task_plan(path: Path) -> tuple[list[dict[str, Any]], str]:
    if not path.is_file():
        raise ContextFinalizationError(f"context task plan not found: {path}")
    plan_sha256 = sha256_file(path)
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {
            "task_index", "chrom", "cofactor_motif_ids", "output_tier",
            "builder_source_commit",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ContextFinalizationError(
                "task plan lacks required columns: "
                + ", ".join(sorted(required))
            )
        rows = list(reader)
    if not rows:
        raise ContextFinalizationError("context task plan contains no tasks")

    tasks: list[dict[str, Any]] = []
    for expected_index, row in enumerate(rows):
        try:
            task_index = int(row["task_index"])
        except (TypeError, ValueError) as error:
            raise ContextFinalizationError("task_index must be an integer") from error
        if task_index != expected_index:
            raise ContextFinalizationError(
                f"task indices must be contiguous from zero; expected "
                f"{expected_index}, found {task_index}"
            )
        chrom = row["chrom"]
        output_tier = row["output_tier"]
        source_commit = row["builder_source_commit"]
        task_kind = row.get("task_kind") or "cofactor_context"
        if (not SAFE_IDENTIFIER.fullmatch(chrom)
                or not SAFE_IDENTIFIER.fullmatch(source_commit)):
            raise ContextFinalizationError(f"unsafe task identity at index {task_index}")
        if output_tier not in {"selected", "summary", "band"}:
            raise ContextFinalizationError(
                f"unsupported output tier at task {task_index}: {output_tier}"
            )
        if task_kind == "anchor_annotation":
            if row["cofactor_motif_ids"] != "none" or output_tier != "summary":
                raise ContextFinalizationError(
                    f"anchor-annotation task {task_index} must have no "
                    "cofactors and use summary tier"
                )
            cofactors: list[str] = []
        elif task_kind == "cofactor_context":
            cofactors = row["cofactor_motif_ids"].split(",")
            if (not cofactors or any(
                not SAFE_IDENTIFIER.fullmatch(motif) or motif == ANCHOR_MOTIF
                for motif in cofactors
            ) or len(cofactors) != len(set(cofactors))):
                raise ContextFinalizationError(
                    f"invalid cofactor list at task {task_index}"
                )
        else:
            raise ContextFinalizationError(
                f"unsupported task kind at task {task_index}: {task_kind}"
            )

        schema_text = row.get("context_schema_version") or "5"
        try:
            schema_version = int(schema_text)
        except ValueError as error:
            raise ContextFinalizationError(
                f"invalid context schema at task {task_index}: {schema_text}"
            ) from error
        if schema_version not in {5, 6, 7}:
            raise ContextFinalizationError(
                f"unsupported context schema at task {task_index}: {schema_version}"
            )
        if schema_version == 5:
            if output_tier != "band":
                raise ContextFinalizationError(
                    "legacy schema-5 finalization is restricted to annotation-free "
                    "band runs because their GTF content was not pinned"
                )
            gtf_size_bytes = None
            gtf_sha256 = None
        else:
            gtf_size_text = row.get("gtf_size_bytes") or ""
            gtf_sha256_text = row.get("gtf_sha256") or ""
            if output_tier == "band":
                if gtf_size_text != "0" or gtf_sha256_text != "none":
                    raise ContextFinalizationError(
                        f"band task {task_index} carries unexpected GTF identity"
                    )
                gtf_size_bytes = None
                gtf_sha256 = None
            else:
                try:
                    gtf_size_bytes = int(gtf_size_text)
                except ValueError as error:
                    raise ContextFinalizationError(
                        f"invalid GTF size at task {task_index}"
                    ) from error
                if gtf_size_bytes < 0 or not SHA256.fullmatch(gtf_sha256_text):
                    raise ContextFinalizationError(
                        f"invalid GTF identity at task {task_index}"
                    )
                gtf_sha256 = gtf_sha256_text

        if schema_version >= 7:
            annotation_release = row.get("annotation_release") or ""
            promoter_definition_id = row.get("promoter_definition_id") or ""
            if (not SAFE_IDENTIFIER.fullmatch(annotation_release)
                    or not SAFE_IDENTIFIER.fullmatch(promoter_definition_id)):
                raise ContextFinalizationError(
                    f"invalid annotation identity at task {task_index}"
                )
            try:
                promoter_upstream_bp = int(row.get("promoter_upstream_bp") or "")
                promoter_downstream_bp = int(row.get("promoter_downstream_bp") or "")
            except ValueError as error:
                raise ContextFinalizationError(
                    f"invalid promoter extent at task {task_index}"
                ) from error
            if promoter_upstream_bp < 0 or promoter_downstream_bp < 0:
                raise ContextFinalizationError(
                    f"negative promoter extent at task {task_index}"
                )
        else:
            annotation_release = None
            promoter_definition_id = None
            promoter_upstream_bp = None
            promoter_downstream_bp = None

        tasks.append({
            "task_index": task_index,
            "task_kind": task_kind,
            "chrom": chrom,
            "cofactor_motif_ids": cofactors,
            "output_tier": output_tier,
            "builder_source_commit": source_commit,
            "context_schema_version": schema_version,
            "gtf_size_bytes": gtf_size_bytes,
            "gtf_sha256": gtf_sha256,
            "annotation_release": annotation_release,
            "promoter_definition_id": promoter_definition_id,
            "promoter_upstream_bp": promoter_upstream_bp,
            "promoter_downstream_bp": promoter_downstream_bp,
        })

    for field in (
        "task_kind", "output_tier", "builder_source_commit",
        "context_schema_version",
    ):
        values = {task[field] for task in tasks}
        if len(values) != 1:
            raise ContextFinalizationError(
                f"task plan mixes {field} values: {sorted(values)}"
            )
    gtf_identities = {
        (task["gtf_size_bytes"], task["gtf_sha256"]) for task in tasks
    }
    if len(gtf_identities) != 1:
        raise ContextFinalizationError("task plan mixes GTF content identities")
    annotation_identities = {
        (task["annotation_release"], task["promoter_definition_id"],
         task["promoter_upstream_bp"], task["promoter_downstream_bp"])
        for task in tasks
    }
    if len(annotation_identities) != 1:
        raise ContextFinalizationError(
            "task plan mixes annotation or promoter definitions"
        )
    return tasks, plan_sha256


def query_package_config(duckdb: str, package: Path) -> dict[str, Any]:
    database = package / "context.duckdb"
    query = """
SELECT c.*,
       (SELECT COUNT(*) FROM (
            SELECT * FROM anchor_motif_band_feature LIMIT 0
        ))::BIGINT AS feature_schema_probe
FROM motif_context_run_config c;
"""
    process = subprocess.run(
        [duckdb, "-light-mode", "-readonly", "-json", str(database), "-c", query],
        cwd=package,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ContextFinalizationError(
            process.stderr.strip() or f"cannot query context package: {package}"
        )
    try:
        rows = json.loads(process.stdout or "[]")
    except json.JSONDecodeError as error:
        raise ContextFinalizationError(
            f"context package returned invalid JSON: {package}"
        ) from error
    if not isinstance(rows, list) or len(rows) != 1:
        raise ContextFinalizationError(
            f"context package must contain one run-config row: {package}"
        )
    return rows[0]


def require_equal(config: dict[str, Any], field: str, expected: Any,
                  package: Path) -> None:
    observed = config.get(field)
    if observed != expected:
        raise ContextFinalizationError(
            f"{package}: {field} is {observed!r}, expected {expected!r}"
        )


def validate_package_config(config: dict[str, Any], task: dict[str, Any],
                            package: Path) -> None:
    expected = {
        "schema_version": task["context_schema_version"],
        "builder_source_commit": task["builder_source_commit"],
        "input_uniqueness": "validated_scan_inventory",
        "genome_id": "homo_sapiens_grch38_ensembl113_primary",
        "motif_set_id": "jaspar2026_core_nonredundant",
        "anchor_motif_id": ANCHOR_MOTIF,
        "score_mode": "log2_relative_risk",
        "pseudocount": 1.0,
        "background_model_id": "uniform_acgt_v1",
        "pseudocount_scheme": "additive_per_base",
        "anchor_minimum_score": -1.0,
        "partner_minimum_score": 0.0,
        "anchor_selection_mode": "local_peak",
        "anchor_local_peak_flank_bp": 150,
        "capture_flank_bp": 150,
        "context_flank_bp": 150,
        "tandem_flank_bp": 20,
        "cofactor_pair_flank_bp": 150,
        "output_tier": task["output_tier"],
        "cofactor_pair_scope": "at_least_one_member_is_a_tp73_context_locus",
        "cofactor_motif_locus_scope": "tp73_context_loci_plus_their_pair_partners",
        "cofactor_locus_pair_feature_scope": "tp73_context_loci_only",
        "feature_schema_probe": 0,
    }
    for field, value in expected.items():
        require_equal(config, field, value, package)

    if task["context_schema_version"] >= 6:
        require_equal(config, "gtf_sha256", task["gtf_sha256"], package)
        require_equal(config, "gtf_size_bytes", task["gtf_size_bytes"], package)
        if task["output_tier"] == "band":
            require_equal(config, "gtf_source", None, package)
        elif not config.get("gtf_source"):
            raise ContextFinalizationError(f"{package}: GTF source path is absent")
    if task["context_schema_version"] >= 7:
        for field in (
            "annotation_release", "promoter_definition_id",
            "promoter_upstream_bp", "promoter_downstream_bp",
        ):
            require_equal(config, field, task[field], package)
        require_equal(
            config, "tss_coordinate_rule",
            "one_base_bed_start_plus_end_minus_1", package,
        )
        require_equal(
            config, "nearest_tss_rule",
            "all_physical_tss_at_minimum_center_distance", package,
        )
        require_equal(
            config, "promoter_membership_rule",
            "anchor_overlaps_resolved_promoter_interval", package,
        )


def dataset_name(relative: Path) -> str:
    parts = relative.parts
    if len(parts) >= 3 and parts[:2] == ("tables", "jaspar2026"):
        return parts[2]
    if relative.name == "context.duckdb":
        return "duckdb_catalog"
    if relative.name == "input_manifest.json":
        return "input_manifest"
    return "package_metadata"


def package_files(run_root: Path, package: Path,
                  task: dict[str, Any]) -> list[dict[str, Any]]:
    required = {
        "context.duckdb",
        "input_manifest.json",
        "tables/jaspar2026/context_run_config.parquet",
    }
    if task["output_tier"] == "band":
        required.add(
            "tables/jaspar2026/anchor_motif_band_feature/data.parquet"
        )
    else:
        required.update({
            "tables/jaspar2026/transcription_start_site.parquet",
            "tables/jaspar2026/transcript_tss.parquet",
            "tables/jaspar2026/promoter.parquet",
            "tables/jaspar2026/promoter_gene.parquet",
            "tables/jaspar2026/tp73_anchor_nearest_tss.parquet",
            "tables/jaspar2026/tp73_anchor_promoter.parquet",
        })
    observed: set[str] = set()
    rows: list[dict[str, Any]] = []
    for root, directories, filenames in os.walk(package):
        directories.sort()
        filenames.sort()
        root_path = Path(root)
        for filename in filenames:
            path = root_path / filename
            relative = path.relative_to(package)
            relative_text = relative.as_posix()
            observed.add(relative_text)
            if filename.endswith(".wal"):
                raise ContextFinalizationError(
                    f"unclosed DuckDB write-ahead log remains: {path}"
                )
            stat = path.stat()
            rows.append({
                "task_index": task["task_index"],
                "chrom": task["chrom"],
                "dataset": dataset_name(relative),
                "package_relative_path": relative_text,
                "run_relative_path": path.relative_to(run_root).as_posix(),
                "absolute_path": str(path),
                "bytes": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
                "is_parquet": path.suffix == ".parquet",
            })
    missing = sorted(required - observed)
    if missing:
        raise ContextFinalizationError(
            f"{package}: required package files are missing: {missing}"
        )
    if task["output_tier"] != "band" and not any(
        path.startswith("tables/jaspar2026/tp73_context_anchor/")
        and path.endswith(".parquet") for path in observed
    ):
        raise ContextFinalizationError(
            f"{package}: TP73 context-anchor Parquet payload is missing"
        )
    return rows


def validate_task_package(duckdb: str, run_root: Path,
                          task: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, Any]]:
    task_index = task["task_index"]
    chrom = task["chrom"]
    input_root = run_root / "inputs" / f"task-{task_index}-chrom-{chrom}"
    package = run_root / "packages" / f"chrom-{chrom}" / f"task-{task_index}"
    expected_manifest = input_root / "input_manifest.json"
    package_manifest = package / "input_manifest.json"
    database = package / "context.duckdb"
    for required in (expected_manifest, package_manifest, database):
        if not required.is_file():
            raise ContextFinalizationError(
                f"task {task_index} is incomplete; missing {required}"
            )
    expected_bytes = expected_manifest.read_bytes()
    if package_manifest.read_bytes() != expected_bytes:
        raise ContextFinalizationError(
            f"task {task_index} package input manifest differs from staged input"
        )
    input_contract = json_read(package_manifest)
    expected_motifs = [ANCHOR_MOTIF, *task["cofactor_motif_ids"]]
    if input_contract.get("motifs") != expected_motifs:
        raise ContextFinalizationError(
            f"task {task_index} input motifs differ from the immutable plan"
        )
    if input_contract.get("chromosomes") != [chrom]:
        raise ContextFinalizationError(
            f"task {task_index} input chromosome differs from the immutable plan"
        )
    source_manifest_sha256 = input_contract.get("source_manifest_sha256")
    if not isinstance(source_manifest_sha256, str) or not SHA256.fullmatch(
        source_manifest_sha256
    ):
        raise ContextFinalizationError(
            f"task {task_index} input manifest lacks source-package identity"
        )

    config = query_package_config(duckdb, package)
    validate_package_config(config, task, package)
    files = package_files(run_root, package, task)
    package_row = {
        **task,
        "cofactor_count": len(task["cofactor_motif_ids"]),
        "package_relative_path": package.relative_to(run_root).as_posix(),
        "package_absolute_path": str(package),
        "input_manifest_sha256": hashlib.sha256(expected_bytes).hexdigest(),
        "source_scan_manifest_sha256": source_manifest_sha256,
        "file_count": len(files),
        "parquet_file_count": sum(row["is_parquet"] for row in files),
        "package_bytes": sum(row["bytes"] for row in files),
    }
    return package_row, files, config


def orphan_staging_paths(run_root: Path, excluded: Path | None = None) -> list[str]:
    staging_root = run_root / "staging"
    if not staging_root.is_dir():
        return []
    result: list[str] = []
    for entry in sorted(staging_root.iterdir(), key=lambda path: path.name):
        if excluded is not None and entry == excluded:
            continue
        if entry.name.startswith("task-") and entry.is_dir():
            for child in sorted(entry.iterdir(), key=lambda path: path.name):
                result.append(child.relative_to(run_root).as_posix())
        elif entry.name.startswith(("finalize-", ".finalize-")):
            result.append(entry.relative_to(run_root).as_posix())
    return result


def build_catalog(duckdb: str, destination: Path,
                  run_row: dict[str, Any], task_rows: list[dict[str, Any]],
                  file_rows: list[dict[str, Any]]) -> None:
    source_root = destination.parent
    run_json = source_root / "context_run.jsonl"
    task_json = source_root / "context_task_inventory.jsonl"
    file_json = source_root / "context_file_inventory.jsonl"
    write_json_lines(run_json, [run_row])
    write_json_lines(task_json, task_rows)
    write_json_lines(file_json, file_rows)
    sql = f"""
CREATE TABLE context_run AS
SELECT * FROM read_json_auto({sql_string(run_json)}, format='newline_delimited');
CREATE TABLE context_task_inventory AS
SELECT * FROM read_json_auto({sql_string(task_json)}, format='newline_delimited');
CREATE TABLE context_file_inventory AS
SELECT * FROM read_json_auto({sql_string(file_json)}, format='newline_delimited');
CREATE VIEW context_feature_file_inventory AS
SELECT f.*, t.task_kind, t.cofactor_motif_ids, t.output_tier,
       t.context_schema_version, t.builder_source_commit
FROM context_file_inventory f
JOIN context_task_inventory t USING (task_index, chrom);
CREATE MACRO anchor_motif_band_feature_files(file_paths) AS TABLE
SELECT * FROM read_parquet(file_paths, hive_partitioning=true);
"""
    process = subprocess.run(
        [duckdb, "-bail", str(destination)],
        input=sql,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ContextFinalizationError(
            process.stderr.strip() or "could not build context DuckDB catalog"
        )


def publish_directory_no_replace(source: Path, destination: Path,
                                 lock_path: Path) -> bool:
    """Publish once among cooperating retries without replacing any path."""
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a", encoding="utf-8") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if os.path.lexists(destination):
            return False
        os.rename(source, destination)
        return True


def validate_existing_final(duckdb: str, output: Path, plan_sha256: str,
                            expected_tasks: int) -> bool:
    manifest_path = output / "manifest.json"
    database = output / "context.duckdb"
    if not manifest_path.is_file() or not database.is_file():
        return False
    manifest = json_read(manifest_path)
    if (manifest.get("state") != "complete"
            or manifest.get("task_plan_sha256") != plan_sha256
            or manifest.get("task_count") != expected_tasks):
        return False
    query = (
        "SELECT COUNT(*)::BIGINT AS tasks FROM context_task_inventory;"
    )
    process = subprocess.run(
        [duckdb, "-readonly", "-json", str(database), "-c", query],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        return False
    try:
        rows = json.loads(process.stdout or "[]")
    except json.JSONDecodeError:
        return False
    return rows == [{"tasks": expected_tasks}]


def finalize(arguments: argparse.Namespace) -> None:
    if shutil.which(arguments.duckdb) is None:
        raise ContextFinalizationError(
            f"DuckDB executable not found: {arguments.duckdb}"
        )
    run_root = arguments.run_root.expanduser().resolve()
    if not run_root.is_dir():
        raise ContextFinalizationError(f"context run root not found: {run_root}")
    task_file = arguments.task_file.expanduser().resolve()
    tasks, plan_sha256 = parse_task_plan(task_file)
    output = (
        arguments.output.expanduser().resolve()
        if arguments.output is not None else run_root / "final"
    )
    if output.exists():
        if validate_existing_final(
            arguments.duckdb, output, plan_sha256, len(tasks)
        ):
            print(f"I: Reusing finalized context catalog: {output}", file=sys.stderr)
            return
        raise ContextFinalizationError(
            f"existing final context catalog is incomplete or incompatible: {output}"
        )

    task_rows: list[dict[str, Any]] = []
    file_rows: list[dict[str, Any]] = []
    configs: list[dict[str, Any]] = []
    errors: list[str] = []
    for offset, task in enumerate(tasks, start=1):
        try:
            task_row, task_files, config = validate_task_package(
                arguments.duckdb, run_root, task
            )
            task_rows.append(task_row)
            file_rows.extend(task_files)
            configs.append(config)
        except ContextFinalizationError as error:
            errors.append(str(error))
        if offset % 100 == 0 or offset == len(tasks):
            print(
                f"I: Context finalization inspected {offset}/{len(tasks)} tasks; "
                f"failures={len(errors)}",
                file=sys.stderr,
            )
    if errors:
        preview = "\n".join(f"  - {error}" for error in errors[:20])
        remainder = len(errors) - min(len(errors), 20)
        if remainder:
            preview += f"\n  - ... and {remainder} more"
        raise ContextFinalizationError(
            f"cannot finalize: {len(errors)} of {len(tasks)} tasks failed:\n{preview}"
        )

    scan_manifest_ids = {
        row["source_scan_manifest_sha256"] for row in task_rows
    }
    if len(scan_manifest_ids) != 1:
        raise ContextFinalizationError(
            "task inputs originate from different finalized scan manifests"
        )
    gtf_sources = {
        (config.get("gtf_source"), config.get("gtf_size_bytes"),
         config.get("gtf_sha256"))
        for config in configs if config.get("gtf_source") is not None
    }
    if len(gtf_sources) > 1:
        raise ContextFinalizationError("task packages report different GTF sources")
    if gtf_sources:
        source_text, expected_size, expected_sha256 = next(iter(gtf_sources))
        source = Path(source_text)
        if (not source.is_file() or source.stat().st_size != expected_size
                or sha256_file(source) != expected_sha256):
            raise ContextFinalizationError(
                "current GTF bytes do not match finalized task provenance"
            )

    staging_root = run_root / "staging"
    staging_root.mkdir(parents=True, exist_ok=True)
    attempt = Path(tempfile.mkdtemp(prefix=".finalize-", dir=staging_root))
    staged_output = attempt / "final"
    staged_output.mkdir()
    try:
        orphans = orphan_staging_paths(run_root, excluded=attempt)
        completed_at = utc_now()
        run_row = {
            "state": "complete",
            "completed_at_utc": completed_at,
            "task_plan_sha256": plan_sha256,
            "task_count": len(task_rows),
            "file_count": len(file_rows),
            "parquet_file_count": sum(row["is_parquet"] for row in file_rows),
            "package_bytes": sum(row["package_bytes"] for row in task_rows),
            "output_tier": task_rows[0]["output_tier"],
            "task_kind": task_rows[0]["task_kind"],
            "context_schema_version": task_rows[0]["context_schema_version"],
            "builder_source_commit": task_rows[0]["builder_source_commit"],
            "source_scan_manifest_sha256": next(iter(scan_manifest_ids)),
            "gtf_size_bytes": task_rows[0]["gtf_size_bytes"],
            "gtf_sha256": task_rows[0]["gtf_sha256"],
            "annotation_release": task_rows[0]["annotation_release"],
            "promoter_definition_id": task_rows[0]["promoter_definition_id"],
            "promoter_upstream_bp": task_rows[0]["promoter_upstream_bp"],
            "promoter_downstream_bp": task_rows[0]["promoter_downstream_bp"],
            "orphan_staging_count": len(orphans),
            "finalization_validation_mode": (
                "exact_task_plan_input_manifest_run_config_schema_probe_file_inventory"
            ),
        }
        database = staged_output / "context.duckdb"
        build_catalog(
            arguments.duckdb, database, run_row, task_rows, file_rows
        )
        json_write(staged_output / "orphan_staging_paths.json", {
            "policy": "report_only_no_automatic_deletion",
            "paths": orphans,
        })
        manifest = {
            "format_version": 1,
            **run_row,
            "database": database.name,
            "task_inventory_table": "context_task_inventory",
            "file_inventory_table": "context_file_inventory",
            "feature_macro": "anchor_motif_band_feature_files",
        }
        json_write(staged_output / "manifest.json", manifest)
        if not publish_directory_no_replace(
            staged_output, output, run_root / "locks" / "finalizer.lock"
        ):
            if not validate_existing_final(
                arguments.duckdb, output, plan_sha256, len(tasks)
            ):
                raise ContextFinalizationError(
                    "another finalizer published an incompatible catalog"
                )
        print(
            f"I: Finalized {len(task_rows)} context tasks and {len(file_rows)} "
            f"files: {output}",
            file=sys.stderr,
        )
        if orphans:
            print(
                f"W: Reported {len(orphans)} orphan staging paths; none were removed.",
                file=sys.stderr,
            )
    finally:
        shutil.rmtree(attempt, ignore_errors=True)


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Require every immutable motif-context task, then publish an "
            "exact-file DuckDB inventory and completion manifest."
        )
    )
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument(
        "--task-file", required=True, type=Path,
        help="immutable context_tasks.tsv used by the Slurm workers",
    )
    parser.add_argument(
        "--output", type=Path,
        help="new finalized catalog directory (default: RUN_ROOT/final)",
    )
    parser.add_argument("--duckdb", default="duckdb")
    return parser


def main() -> int:
    arguments = argument_parser().parse_args()
    try:
        finalize(arguments)
    except ContextFinalizationError as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

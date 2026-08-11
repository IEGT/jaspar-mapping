#!/usr/bin/env python3

"""Plan, run, inspect, and finalize whole-genome TP73 CUT&RUN evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class GenomeEvidenceError(RuntimeError):
    pass


SCIENTIFIC_SOURCE_FILES = (
    "scripts/build_h3k4me3_anchor_signal.py",
    "scripts/build_tp73_anchor_evidence.py",
    "scripts/export_bigwig_chrom_bedgraph.py",
    "scripts/manage_tp73_genome_evidence.py",
    "scripts/run_tp73_chromosome_production.py",
)


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


def run_json(command: list[str], *, cwd: Path | None = None) -> list[dict[str, Any]]:
    process = subprocess.run(
        command, text=True, capture_output=True, check=False, cwd=cwd
    )
    if process.returncode != 0:
        raise GenomeEvidenceError(process.stderr.strip() or "command failed")
    value = json.loads(process.stdout or "[]")
    if not isinstance(value, list):
        raise GenomeEvidenceError("JSON result is not a row array")
    return value


def run(command: list[str], *, cwd: Path | None = None) -> None:
    displayed = list(command)
    if "-c" in displayed:
        sql_index = displayed.index("-c") + 1
        if sql_index < len(displayed) and len(displayed[sql_index]) > 200:
            displayed[sql_index] = f"<SQL:{len(displayed[sql_index])} chars>"
    print("I: running: " + " ".join(displayed), file=sys.stderr, flush=True)
    process = subprocess.run(command, check=False, cwd=cwd)
    if process.returncode != 0:
        raise GenomeEvidenceError(
            f"command failed with exit code {process.returncode}: {command[0]}"
        )


def immutable_write(path: Path, content: str) -> None:
    encoded = content.encode("utf-8")
    if path.exists():
        if path.read_bytes() != encoded:
            raise GenomeEvidenceError(f"immutable file differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        temporary.write_bytes(encoded)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise GenomeEvidenceError(f"JSON file is missing: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise GenomeEvidenceError(f"JSON value is not an object: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise GenomeEvidenceError(f"TSV file is missing: {path}")
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        rows = list(reader)
    if reader.fieldnames is None or not rows:
        raise GenomeEvidenceError(f"TSV is empty or lacks a header: {path}")
    return rows


def package_database(package: Path) -> Path:
    manifest = load_json(package / "manifest.json")
    database_value = manifest.get("database")
    if not database_value:
        raise GenomeEvidenceError("scan manifest does not identify its database")
    database = Path(str(database_value))
    if not database.is_absolute():
        database = package / database
    database = database.resolve()
    if not database.is_file():
        raise GenomeEvidenceError(f"scan database is missing: {database}")
    return database


def classify_chromosome(chrom: str) -> tuple[str, bool]:
    if re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2])", chrom):
        return "autosome", True
    if chrom in {"X", "Y"}:
        return "sex_chromosome", False
    if chrom in {"MT", "M", "chrM", "chrMT", "25", "chr25"}:
        return "mitochondrial_bystander_control", False
    return "other_primary_assembly", False


def is_mitochondrial_chromosome(chrom: str) -> bool:
    return chrom.removeprefix("chr") in {"M", "MT", "25"}


def safe_label(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", value).strip("-") or "region"


def git_identity(source: Path) -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        text=True, capture_output=True, check=False,
    )
    if commit.returncode != 0:
        raise GenomeEvidenceError(commit.stderr.strip() or "cannot read Git commit")
    dirty = False
    for arguments in (
        ["diff", "--quiet", "--ignore-submodules", "--"],
        ["diff", "--cached", "--quiet", "--ignore-submodules", "--"],
    ):
        result = subprocess.run(
            ["git", "-C", str(source), *arguments], check=False
        )
        if result.returncode not in (0, 1):
            raise GenomeEvidenceError("cannot inspect Git state")
        dirty = dirty or result.returncode == 1
    return commit.stdout.strip(), dirty


def scientific_hashes(source: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for relative in SCIENTIFIC_SOURCE_FILES:
        path = source / relative
        if not path.is_file():
            raise GenomeEvidenceError(f"scientific source is missing: {path}")
        result[relative] = sha256(path)
    return result


def verify_scientific_hashes(config: dict[str, Any]) -> None:
    source = Path(str(config["source"]))
    expected = config.get("scientific_source_file_sha256")
    if not isinstance(expected, dict) or not expected:
        raise GenomeEvidenceError("run plan lacks scientific source hashes")
    for relative, digest in expected.items():
        path = source / relative
        if not path.is_file() or sha256(path) != digest:
            raise GenomeEvidenceError(f"scientific source changed: {path}")


def prepare(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    source = arguments.source.expanduser().resolve()
    annotation_run = arguments.annotation_run.expanduser().resolve()
    scan_package = arguments.scan_package.expanduser().resolve()
    track_root = arguments.track_root.expanduser().resolve()
    runtime_prefix = arguments.runtime_prefix.expanduser().resolve()
    track_manifest = (
        arguments.track_manifest.expanduser().resolve()
        if arguments.track_manifest is not None
        else source / "config" / "h3k4me3_cutandrun_tracks_v1.tsv"
    )
    for path, label in (
        (source, "source"), (annotation_run, "annotation run"),
        (scan_package, "scan package"), (track_root, "track root"),
        (runtime_prefix, "runtime prefix"),
    ):
        if not path.is_dir():
            raise GenomeEvidenceError(f"{label} is missing: {path}")
    for path, label in (
        (track_manifest, "track manifest"),
        (annotation_run / "final" / "manifest.json", "annotation manifest"),
        (scan_package / "manifest.json", "scan manifest"),
    ):
        if not path.is_file():
            raise GenomeEvidenceError(f"{label} is missing: {path}")
    duckdb = shutil.which(arguments.duckdb) or arguments.duckdb
    if not Path(duckdb).is_file():
        raise GenomeEvidenceError(f"DuckDB is unavailable: {arguments.duckdb}")

    scan_database = package_database(scan_package)
    motif_summary = run_json([
        duckdb, "-readonly", "-json", str(scan_database), "-c",
        "SELECT count(*)::BIGINT AS motifs, "
        "count(DISTINCT motif_id)::BIGINT AS distinct_motifs "
        "FROM motif_metadata;",
    ], cwd=scan_package)
    if (len(motif_summary) != 1
            or int(motif_summary[0]["motifs"]) <= 0
            or motif_summary[0]["motifs"] != motif_summary[0]["distinct_motifs"]):
        raise GenomeEvidenceError("scan motif catalog is empty or non-unique")
    motif_count = int(motif_summary[0]["motifs"])
    regions = run_json([
        duckdb, "-readonly", "-json", str(scan_database), "-c",
        "SELECT s.sequence_order::INTEGER AS sequence_order, "
        "CAST(i.chrom AS VARCHAR) AS chrom, s.length::BIGINT AS chrom_length, "
        "s.sequence_sha256::VARCHAR AS sequence_sha256, "
        "count(*)::BIGINT AS files, "
        "count(DISTINCT i.motif_id)::BIGINT AS motifs, "
        "count(DISTINCT i.strand)::BIGINT AS strands, "
        "count(*) FILTER (WHERE i.strand = '+')::BIGINT AS plus_files, "
        "count(*) FILTER (WHERE i.strand = '-')::BIGINT AS minus_files, "
        "count(*) - count(DISTINCT (i.motif_id, i.strand))::BIGINT "
        "AS duplicate_motif_strands "
        "FROM scan_file_inventory i "
        "LEFT JOIN sequence_region s "
        "ON CAST(s.chrom AS VARCHAR) = CAST(i.chrom AS VARCHAR) "
        "GROUP BY s.sequence_order, i.chrom, s.length, s.sequence_sha256 "
        "ORDER BY s.sequence_order;",
    ], cwd=scan_package)
    if not regions:
        raise GenomeEvidenceError("scan package has no finalized hit regions")
    invalid_regions = [
        str(row["chrom"]) for row in regions
        if row.get("sequence_order") is None
        or row.get("chrom_length") is None
        or row.get("sequence_sha256") is None
        or int(row["files"]) != motif_count * 2
        or int(row["motifs"]) != motif_count
        or int(row["strands"]) != 2
        or int(row["plus_files"]) != motif_count
        or int(row["minus_files"]) != motif_count
        or int(row["duplicate_motif_strands"]) != 0
    ]
    if invalid_regions:
        raise GenomeEvidenceError(
            "scan hit inventory is not motif-by-strand complete for: "
            + ", ".join(invalid_regions)
        )
    chromosomes = [str(row["chrom"]) for row in regions]
    if len(chromosomes) != len(set(chromosomes)):
        raise GenomeEvidenceError("scan sequence regions contain duplicates")
    annotation_database = annotation_run / "final" / "context.duckdb"
    annotation_rows = run_json([
        duckdb, "-readonly", "-json", str(annotation_database), "-c",
        "SELECT CAST(chrom AS VARCHAR) AS chrom, count(*)::BIGINT AS files "
        "FROM context_file_inventory WHERE dataset = 'tp73_context_anchor' "
        "AND is_parquet "
        "GROUP BY chrom;",
    ])
    observed = {str(row["chrom"]): int(row["files"]) for row in annotation_rows}
    annotation_chromosomes: dict[str, str] = {}
    invalid: list[str] = []
    for chrom in chromosomes:
        candidates = [chrom] if chrom in observed else []
        if is_mitochondrial_chromosome(chrom):
            candidates = sorted(
                candidate for candidate in observed
                if is_mitochondrial_chromosome(candidate)
            )
        if len(candidates) != 1 or observed.get(candidates[0]) != 1:
            invalid.append(chrom)
        else:
            annotation_chromosomes[chrom] = candidates[0]
    if invalid:
        raise GenomeEvidenceError(
            "annotation needs exactly one TP73 anchor file for: "
            + ", ".join(invalid)
        )

    fields = (
        "task_index", "sequence_order", "chrom", "annotation_chrom",
        "chrom_length",
        "sequence_sha256", "analysis_partition", "primary_inference",
    )
    content = io.StringIO(newline="")
    writer = csv.DictWriter(content, fieldnames=fields, delimiter="\t",
                            lineterminator="\n")
    writer.writeheader()
    partitions: dict[str, int] = {}
    for task_index, row in enumerate(regions):
        chrom = str(row["chrom"])
        partition, primary = classify_chromosome(chrom)
        partitions[partition] = partitions.get(partition, 0) + 1
        writer.writerow({
            "task_index": task_index,
            "sequence_order": int(row["sequence_order"]),
            "chrom": chrom,
            "annotation_chrom": annotation_chromosomes[chrom],
            "chrom_length": int(row["chrom_length"]),
            "sequence_sha256": str(row["sequence_sha256"]),
            "analysis_partition": partition,
            "primary_inference": str(primary).lower(),
        })
    task_path = run_root / "plan" / "chromosome_tasks.tsv"
    immutable_write(task_path, content.getvalue())
    commit, dirty = git_identity(source)
    config = {
        "schema_version": 1,
        "run_id": arguments.run_id,
        "analysis": "whole_genome_tp73_cutandrun_evidence",
        "source": str(source),
        "source_commit": commit,
        "source_dirty": dirty,
        "scientific_source_file_sha256": scientific_hashes(source),
        "annotation_run": str(annotation_run),
        "annotation_manifest_sha256": sha256(
            annotation_run / "final" / "manifest.json"
        ),
        "scan_package": str(scan_package),
        "scan_manifest_sha256": sha256(scan_package / "manifest.json"),
        "scan_region_selection": (
            "distinct sequence regions in scan_file_inventory; each requires "
            "one file per catalog motif and strand"
        ),
        "scan_motif_count": motif_count,
        "track_root": str(track_root),
        "track_manifest": str(track_manifest),
        "track_manifest_sha256": sha256(track_manifest),
        "runtime_prefix": str(runtime_prefix),
        "scratch_root": str(arguments.scratch_root),
        "task_count": len(regions),
        "partition_counts": partitions,
        "primary_inference_partition": "autosome",
        "sex_chromosome_policy": "separate_sensitivity",
        "mitochondrial_policy": (
            "bystander_control_excluded_from_primary_inference_due_to_copy_number"
        ),
        "minimum_anchor_score": -1,
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
    for name in ("logs", "chromosomes", "final"):
        (run_root / name).mkdir(parents=True, exist_ok=True)
    print(len(regions))


def planned_tasks(run_root: Path) -> list[dict[str, str]]:
    rows = read_tsv(run_root / "plan" / "chromosome_tasks.tsv")
    if [int(row["task_index"]) for row in rows] != list(range(len(rows))):
        raise GenomeEvidenceError("task indices are not contiguous")
    if len({row["chrom"] for row in rows}) != len(rows):
        raise GenomeEvidenceError("task plan contains duplicate chromosomes")
    return rows


def chromosome_root(run_root: Path, row: dict[str, str]) -> Path:
    return run_root / "chromosomes" / f"chrom-{safe_label(row['chrom'])}"


def validate_child(row: dict[str, str], root: Path,
                   verify_checksum: bool = False) -> dict[str, Any]:
    manifest_path = root / "final" / "manifest.json"
    manifest = load_json(manifest_path)
    if (manifest.get("state") != "complete"
            or str(manifest.get("chrom")) != row["chrom"]
            or int(manifest.get("chrom_length", -1)) != int(row["chrom_length"])
            or manifest.get("analysis_partition") != row["analysis_partition"]):
        raise GenomeEvidenceError(f"chromosome package identity differs: {root}")
    records = {
        str(record.get("relative_path")): record
        for record in manifest.get("outputs", []) if isinstance(record, dict)
    }
    relative = "tp73_anchor_evidence.parquet"
    if relative not in records:
        raise GenomeEvidenceError(f"chromosome package lacks evidence: {root}")
    evidence = root / "final" / relative
    record = records[relative]
    if not evidence.is_file() or evidence.stat().st_size != int(record.get("bytes", -1)):
        raise GenomeEvidenceError(f"chromosome evidence size differs: {evidence}")
    if verify_checksum and sha256(evidence) != record.get("sha256"):
        raise GenomeEvidenceError(f"chromosome evidence checksum differs: {evidence}")
    return manifest


def run_task(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(run_root / "plan" / "run_config.json")
    rows = planned_tasks(run_root)
    task_index = arguments.task_index
    if task_index is None:
        value = os.environ.get("SLURM_ARRAY_TASK_ID", "")
        if not value.isdigit():
            raise GenomeEvidenceError(
                "--task-index or SLURM_ARRAY_TASK_ID is required"
            )
        task_index = int(value)
    if task_index < 0 or task_index >= len(rows):
        raise GenomeEvidenceError(f"task index is outside the plan: {task_index}")
    verify_scientific_hashes(config)
    row = rows[task_index]
    source = Path(str(config["source"]))
    runtime = Path(str(config["runtime_prefix"]))
    command = [
        sys.executable,
        str(source / "scripts" / "run_tp73_chromosome_production.py"),
        "--run-root", str(chromosome_root(run_root, row)),
        "--annotation-run", str(config["annotation_run"]),
        "--track-root", str(config["track_root"]),
        "--track-manifest", str(config["track_manifest"]),
        "--source", str(source),
        "--source-commit", str(config["source_commit"]),
        "--chrom", row["chrom"],
        "--anchor-chrom", row["annotation_chrom"],
        "--chrom-length", row["chrom_length"],
        "--analysis-partition", row["analysis_partition"],
        "--duckdb", str(runtime / "duckdb" / "bin" / "duckdb"),
        "--bigwig-python", str(runtime / "duckdb" / "bin" / "python3"),
        "--scratch-root", str(config["scratch_root"]),
        "--minimum-anchor-score", str(config["minimum_anchor_score"]),
        "--threads", str(config["threads"]),
        "--memory-limit", str(config["memory_limit"]),
        "--minimum-free-run-gb", str(config["minimum_free_run_gb"]),
        "--minimum-free-scratch-gb", str(config["minimum_free_scratch_gb"]),
    ]
    if row["primary_inference"] == "true":
        command.append("--primary-inference")
    print(
        f"I: executing chromosome task {task_index}/{len(rows) - 1}: "
        f"{row['chrom']} ({row['analysis_partition']})",
        file=sys.stderr, flush=True,
    )
    os.execv(sys.executable, command)


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
        except (GenomeEvidenceError, OSError, ValueError, json.JSONDecodeError):
            invalid += 1
    print("key\tvalue")
    print(f"planned\t{len(rows)}")
    print(f"complete\t{complete}")
    print(f"pending\t{len(rows) - complete - invalid}")
    print(f"invalid_complete_markers\t{invalid}")
    for key, value in sorted(partitions.items()):
        print(f"complete_{key}\t{value}")


def output_inventory(root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path.name == "manifest.json" or "chromosomes" in path.parts:
            continue
        rows.append({
            "relative_path": str(path.relative_to(root)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    return rows


def finalize(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(run_root / "plan" / "run_config.json")
    rows = planned_tasks(run_root)
    if int(config.get("task_count", -1)) != len(rows):
        raise GenomeEvidenceError("run configuration task count differs")
    verify_scientific_hashes(config)
    for row in rows:
        validate_child(row, chromosome_root(run_root, row), verify_checksum=True)
    final = run_root / "final" / "genome_evidence"
    if final.exists():
        manifest = load_json(final / "manifest.json")
        if (manifest.get("state") == "complete"
                and manifest.get("run_id") == config.get("run_id")):
            print(f"I: reusing finalized genome evidence: {final}", file=sys.stderr)
            return
        raise GenomeEvidenceError(f"final output has another identity: {final}")

    final.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".tp73-genome-evidence.", dir=final.parent))
    try:
        partition_paths: dict[str, list[Path]] = {}
        inventory: list[dict[str, Any]] = []
        for row in rows:
            source = (
                chromosome_root(run_root, row) / "final" /
                "tp73_anchor_evidence.parquet"
            )
            partition = row["analysis_partition"]
            target = (
                staging / "chromosomes" / f"analysis_partition={partition}" /
                f"chrom={safe_label(row['chrom'])}" / "tp73_anchor_evidence.parquet"
            )
            target.parent.mkdir(parents=True, exist_ok=True)
            os.link(source, target)
            partition_paths.setdefault(partition, []).append(target)
            inventory.append({
                "sequence_order": int(row["sequence_order"]),
                "chrom": row["chrom"],
                "annotation_chrom": row["annotation_chrom"],
                "chrom_length": int(row["chrom_length"]),
                "analysis_partition": partition,
                "primary_inference": row["primary_inference"],
                "relative_path": str(target.relative_to(staging)),
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

        tables = staging / "tables"
        tables.mkdir()
        aggregate_paths: dict[str, Path] = {}
        duckdb = str(arguments.duckdb)
        for partition, paths in sorted(partition_paths.items()):
            output = tables / f"tp73_anchor_evidence_{partition}.parquet"
            aggregate_paths[partition] = output
            run([duckdb, "-light-mode", ":memory:", "-c", f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;
SET temp_directory={sql_string(arguments.temp_directory)};
COPY (
  SELECT * FROM read_parquet({sql_path_list(paths)}, hive_partitioning=false)
  ORDER BY chrom, anchor_start, anchor_end
) TO {sql_string(output)} (FORMAT PARQUET, COMPRESSION ZSTD);
"""])

        autosomes = aggregate_paths.get("autosome")
        if autosomes is None:
            raise GenomeEvidenceError("primary autosome evidence is absent")
        description = run_json([
            duckdb, "-light-mode", "-json", ":memory:", "-c",
            f"DESCRIBE SELECT * FROM read_parquet({sql_string(autosomes)});",
        ])
        columns = {str(row["column_name"]) for row in description}
        support_columns = sorted(
            name for name in columns if name.startswith("supported_")
        )
        if len(support_columns) != 12:
            raise GenomeEvidenceError(
                "expected 12 included TP73/control support columns, found "
                f"{len(support_columns)}"
            )
        if any("skmel29_1" in name for name in support_columns):
            raise GenomeEvidenceError("excluded skmel29_1 entered genome evidence")
        depth_columns = [
            name.replace("supported_", "depth_", 1)
            for name in support_columns
        ]
        missing_depth = sorted(set(depth_columns) - columns)
        if missing_depth:
            raise GenomeEvidenceError(
                "support columns lack matched depth columns: "
                + ", ".join(missing_depth)
            )
        summary_selects: list[str] = []
        all_paths = [path for paths in partition_paths.values() for path in paths]
        for support in support_columns:
            depth = support.replace("supported_", "depth_", 1)
            summary_selects.append(f"""
SELECT chrom, {sql_string(support)}::VARCHAR AS support_column,
       {sql_string(depth)}::VARCHAR AS depth_column,
       count(*)::BIGINT AS anchors,
       count(*) FILTER (WHERE {support})::BIGINT AS supported_anchors,
       avg({support}::INTEGER)::DOUBLE AS supported_fraction,
       avg({depth}) FILTER (WHERE {depth} > 0)::DOUBLE AS mean_positive_depth,
       quantile_cont({depth}, 0.5) FILTER (WHERE {depth} > 0)::DOUBLE
           AS median_positive_depth,
       max({depth})::DOUBLE AS maximum_depth
FROM evidence GROUP BY chrom
""".strip())
        coverage_summary = tables / "tp73_coverage_summary_by_chromosome.parquet"
        summary_sql = "\nUNION ALL\n".join(summary_selects)
        run([duckdb, "-light-mode", ":memory:", "-c", f"""
SET threads={arguments.threads};
COPY (
  WITH evidence AS (
    SELECT * FROM read_parquet({sql_path_list(all_paths)}, hive_partitioning=false)
  )
  SELECT * FROM ({summary_sql})
  ORDER BY chrom, support_column
) TO {sql_string(coverage_summary)} (FORMAT PARQUET, COMPRESSION ZSTD);
"""])

        partition_chromosomes: dict[str, list[str]] = {}
        for partition, path in aggregate_paths.items():
            observed_rows = run_json([
                duckdb, "-light-mode", "-json", ":memory:", "-c",
                "SELECT DISTINCT CAST(chrom AS VARCHAR) AS chrom "
                f"FROM read_parquet({sql_string(path)}) ORDER BY chrom;",
            ])
            partition_chromosomes[partition] = [
                str(item["chrom"]) for item in observed_rows
            ]
            expected_chromosomes = sorted(
                row["chrom"] for row in rows
                if row["analysis_partition"] == partition
            )
            if sorted(partition_chromosomes[partition]) != expected_chromosomes:
                raise GenomeEvidenceError(
                    f"{partition} chromosome membership differs: "
                    f"{partition_chromosomes[partition]} != {expected_chromosomes}"
                )

        mismatch_predicate = " OR ".join(
            f"{support} IS NULL OR {depth} IS NULL OR {depth} < 0 "
            f"OR {support} <> ({depth} > 0)"
            for support, depth in zip(support_columns, depth_columns)
        )
        validation = run_json([
            duckdb, "-light-mode", "-json", ":memory:", "-c", f"""
WITH all_evidence AS (
  SELECT * FROM read_parquet({sql_path_list(all_paths)}, hive_partitioning=false)
), primary_evidence AS (
  SELECT * FROM read_parquet({sql_string(autosomes)}, hive_partitioning=false)
)
SELECT (SELECT count(*) FROM all_evidence)::BIGINT AS anchors,
       (SELECT count(DISTINCT chrom) FROM all_evidence)::BIGINT AS chromosomes,
       (SELECT count(*) - count(DISTINCT (chrom, anchor_start, anchor_end))
          FROM all_evidence)::BIGINT AS duplicate_keys,
       (SELECT count(*) FROM all_evidence WHERE anchor_score < -1)::BIGINT
           AS below_score_floor,
       (SELECT count(*) FROM all_evidence WHERE {mismatch_predicate})::BIGINT
           AS support_depth_mismatches,
       (SELECT count(*) FROM primary_evidence)::BIGINT AS primary_anchors,
       (SELECT count(DISTINCT chrom) FROM primary_evidence)::BIGINT
           AS primary_chromosomes;
""",
        ])
        values = {key: int(value) for key, value in validation[0].items()}
        if (values["anchors"] <= 0 or values["chromosomes"] != len(rows)
                or values["duplicate_keys"] != 0
                or values["below_score_floor"] != 0
                or values["support_depth_mismatches"] != 0
                or values["primary_anchors"] <= 0
                or values["primary_chromosomes"] != 22):
            raise GenomeEvidenceError(f"genome evidence validation failed: {values}")

        database = staging / "tp73_genome_evidence.duckdb"
        view_statements = [
            "CREATE VIEW chromosome_file_inventory AS SELECT * FROM "
            "read_csv_auto('chromosome_file_inventory.tsv', delim='\\t', header=true);",
            "CREATE VIEW tp73_coverage_summary_by_chromosome AS SELECT * FROM "
            "read_parquet('tables/tp73_coverage_summary_by_chromosome.parquet');",
        ]
        for partition, path in aggregate_paths.items():
            view_statements.append(
                f"CREATE VIEW tp73_anchor_evidence_{partition} AS SELECT * FROM "
                f"read_parquet('tables/{path.name}');"
            )
        run([duckdb, str(database), "-c", "\n".join(view_statements)], cwd=staging)

        manifest = {
            "schema_version": 1,
            "state": "complete",
            "run_id": config["run_id"],
            "analysis": config["analysis"],
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "source_commit": config["source_commit"],
            "chromosomes": len(rows),
            "partition_counts": config["partition_counts"],
            "partition_chromosomes": partition_chromosomes,
            "primary_inference_partition": "autosome",
            "mitochondrial_policy": config["mitochondrial_policy"],
            "validation": values,
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
    print(f"I: finalized TP73 genome evidence: {final}", file=sys.stderr)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser(
        "prepare", help="write an immutable chromosome task plan"
    )
    prepare_parser.add_argument("--run-root", type=Path, required=True)
    prepare_parser.add_argument("--annotation-run", type=Path, required=True)
    prepare_parser.add_argument("--scan-package", type=Path, required=True)
    prepare_parser.add_argument("--track-root", type=Path, required=True)
    prepare_parser.add_argument("--runtime-prefix", type=Path, required=True)
    prepare_parser.add_argument("--source", type=Path, required=True)
    prepare_parser.add_argument("--track-manifest", type=Path)
    prepare_parser.add_argument("--duckdb", default="duckdb")
    prepare_parser.add_argument("--scratch-root", type=Path,
                                default=Path("/scratch") / os.environ.get("USER", "sm718"))
    prepare_parser.add_argument("--run-id", required=True)
    prepare_parser.add_argument("--threads", type=int, default=4)
    prepare_parser.add_argument("--memory-limit", default="28GB")
    prepare_parser.add_argument("--minimum-free-run-gb", type=float, default=5)
    prepare_parser.add_argument("--minimum-free-scratch-gb", type=float, default=10)
    prepare_parser.set_defaults(function=prepare)

    task_parser = subparsers.add_parser(
        "run-task", help="execute one planned chromosome"
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
        "finalize", help="validate and aggregate chromosome evidence"
    )
    finalize_parser.add_argument("--run-root", type=Path, required=True)
    finalize_parser.add_argument("--duckdb", required=True)
    finalize_parser.add_argument("--threads", type=int, default=4)
    finalize_parser.add_argument("--memory-limit", default="28GB")
    finalize_parser.add_argument("--temp-directory", type=Path, required=True)
    finalize_parser.set_defaults(function=finalize)
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        if hasattr(arguments, "threads") and arguments.threads <= 0:
            raise GenomeEvidenceError("--threads must be positive")
        arguments.function(arguments)
        return 0
    except (GenomeEvidenceError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

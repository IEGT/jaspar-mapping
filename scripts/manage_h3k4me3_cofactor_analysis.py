#!/usr/bin/env python3

"""Plan, run, inspect, and finalize all-motif H3K4me3 cofactor analysis."""

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
import signal
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class AnalysisError(RuntimeError):
    pass


AUTOSOMES = tuple(str(value) for value in range(1, 23))
DISTANCE_BAND_COLUMNS = {
    "all_150": "context_score",
    "overlap": "context_score_overlap",
    "adjacent_0_5": "context_score_adjacent_0_5",
    "gap_6_20": "context_score_gap_6_20",
    "gap_21_50": "context_score_gap_21_50",
    "gap_51_100": "context_score_gap_51_100",
    "gap_101_150": "context_score_gap_101_150",
}
OUTPUTS = {
    "intensity_effect": "intensity_effect.tsv",
    "isoform_contrast": "isoform_contrast.tsv",
    "tp73_interaction": "tp73_interaction.tsv",
    "series_summary": "series_summary.tsv",
    "binding_state_summary": "binding_state_summary.tsv",
    "occurrence_summary": "occurrence_summary.tsv",
    "context_stratified_intensity_effect":
        "context_stratified_intensity_effect.tsv",
    "gene_relation_stratified_intensity_effect":
        "gene_relation_stratified_intensity_effect.tsv",
    "gene_relation_stratified_tp73_occupancy":
        "gene_relation_stratified_tp73_occupancy.tsv",
    "score_gradient": "score_gradient.tsv",
    "run_config": "run_config.tsv",
}
BH_PARTITIONS = {
    "intensity_effect": (
        "series_id", "isoform", "negative_reference_threshold", "distance_band",
    ),
    "isoform_contrast": (
        "series_id", "negative_reference_threshold", "distance_band",
    ),
    "tp73_interaction": (
        "series_id", "isoform", "negative_reference_threshold", "contrast",
        "distance_band",
    ),
    "context_stratified_intensity_effect": (
        "series_id", "isoform", "negative_reference_threshold",
        "genomic_context_class", "distance_band",
    ),
    "gene_relation_stratified_intensity_effect": (
        "series_id", "isoform", "negative_reference_threshold",
        "gene_relation_class", "distance_band",
    ),
    "gene_relation_stratified_tp73_occupancy": (
        "negative_reference_threshold", "gene_relation_class", "distance_band",
    ),
    "score_gradient": (
        "series_id", "isoform", "score_clamp_reference", "distance_band",
    ),
}
SCIENTIFIC_SOURCE_FILES = (
    "scripts/analyze_h3k4me3_cofactor_change.R",
    "scripts/manage_h3k4me3_cofactor_analysis.py",
    "scripts/run_h3k4me3_cofactor_analysis_finalize.sh",
    "scripts/submit_h3k4me3_cofactor_analysis_slurm.sh",
)

STARTED = time.monotonic()
CURRENT_PHASE = "startup"
CURRENT_BATCH = "unresolved"
CURRENT_MOTIF = "none"
CURRENT_CHILD = "none"


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
    return re.sub(r"[^A-Za-z0-9._-]+", "-", value).strip("-") or "item"


def load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise AnalysisError(f"JSON file is missing: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise AnalysisError(f"JSON value is not an object: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise AnalysisError(f"TSV file is missing: {path}")
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        rows = list(reader)
    if reader.fieldnames is None or not rows:
        raise AnalysisError(f"TSV is empty or lacks a header: {path}")
    return rows


def immutable_write(path: Path, content: str) -> None:
    encoded = content.encode("utf-8")
    if path.exists():
        if path.read_bytes() != encoded:
            raise AnalysisError(f"immutable file differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        temporary.write_bytes(encoded)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def query_json(duckdb: Path, database: Path | str,
               query: str, cwd: Path | None = None) -> list[dict[str, Any]]:
    mode = [] if str(database) == ":memory:" else ["-readonly"]
    process = subprocess.run(
        [str(duckdb), "-light-mode", *mode, "-json", str(database), "-c", query],
        cwd=cwd, text=True, capture_output=True, check=False,
    )
    if process.returncode != 0:
        raise AnalysisError(process.stderr.strip() or "DuckDB query failed")
    value = json.loads(process.stdout or "[]")
    if not isinstance(value, list):
        raise AnalysisError("DuckDB query did not return a row array")
    return value


def run_sql(duckdb: Path, database: Path | str, query: str,
            cwd: Path | None = None) -> None:
    process = subprocess.run(
        [str(duckdb), "-light-mode", "-bail", str(database)],
        cwd=cwd, input=query, text=True, capture_output=True, check=False,
    )
    if process.returncode != 0:
        raise AnalysisError(process.stderr.strip() or "DuckDB command failed")


def git_identity(source: Path) -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        text=True, capture_output=True, check=False,
    )
    if commit.returncode != 0:
        raise AnalysisError(commit.stderr.strip() or "cannot read Git commit")
    dirty = False
    for arguments in (
        ["diff", "--quiet", "--ignore-submodules", "--"],
        ["diff", "--cached", "--quiet", "--ignore-submodules", "--"],
    ):
        result = subprocess.run(
            ["git", "-C", str(source), *arguments], check=False
        )
        if result.returncode not in (0, 1):
            raise AnalysisError("cannot inspect Git state")
        dirty = dirty or result.returncode == 1
    return commit.stdout.strip(), dirty


def scientific_hashes(source: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for relative in SCIENTIFIC_SOURCE_FILES:
        path = source / relative
        if not path.is_file():
            raise AnalysisError(f"scientific source is missing: {path}")
        result[relative] = sha256(path)
    return result


def verify_scientific_hashes(config: dict[str, Any]) -> None:
    source = Path(str(config["source"]))
    expected = config.get("scientific_source_file_sha256")
    if not isinstance(expected, dict) or not expected:
        raise AnalysisError("run plan lacks scientific source hashes")
    for relative, digest in expected.items():
        path = source / relative
        if not path.is_file() or sha256(path) != digest:
            raise AnalysisError(f"scientific source changed: {path}")


def progress(_number: int | None = None, _frame: object | None = None) -> None:
    print(
        "I: progress signal=SIGUSR1 "
        f"phase={CURRENT_PHASE} batch={CURRENT_BATCH} motif={CURRENT_MOTIF} "
        f"elapsed_s={time.monotonic() - STARTED:.1f} child_pid={CURRENT_CHILD}",
        file=sys.stderr, flush=True,
    )


def set_phase(value: str) -> None:
    global CURRENT_PHASE
    CURRENT_PHASE = value


def run_process(command: list[str], cwd: Path | None = None) -> None:
    global CURRENT_CHILD
    process = subprocess.Popen(command, cwd=cwd)
    CURRENT_CHILD = str(process.pid)
    try:
        while True:
            try:
                status = process.wait()
                break
            except InterruptedError:
                continue
    finally:
        CURRENT_CHILD = "none"
    if status != 0:
        raise AnalysisError(
            f"command failed with exit code {status}: {command[0]}"
        )


def check_free_space(path: Path, minimum_gb: float, label: str) -> None:
    free = shutil.disk_usage(path).free
    print(
        f"I: {label} free space: {free / 1024**3:.1f} GiB "
        f"(required {minimum_gb:.1f} GiB)",
        file=sys.stderr, flush=True,
    )
    if free < minimum_gb * 1024**3:
        raise AnalysisError(f"insufficient {label} free space")


def resolve_package_file(package: Path, relative: str) -> Path:
    relative_path = Path(relative)
    if relative_path.is_absolute() or ".." in relative_path.parts:
        raise AnalysisError(f"unsafe package-relative path: {relative}")
    result = (package / relative_path).resolve()
    try:
        result.relative_to(package)
    except ValueError as error:
        raise AnalysisError(f"package path escapes root: {relative}") from error
    return result


def write_rows(path: Path, rows: list[dict[str, Any]],
               fields: tuple[str, ...]) -> None:
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    immutable_write(path, output.getvalue())


def prepare(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    source = arguments.source.expanduser().resolve()
    h3_package = arguments.h3_package.expanduser().resolve()
    evidence_package = arguments.evidence_package.expanduser().resolve()
    context_package = arguments.context_maxima_package.expanduser().resolve()
    annotation_catalog = arguments.annotation_catalog.expanduser().resolve()
    runtime = arguments.runtime_prefix.expanduser().resolve()
    distance_bands = tuple(
        item.strip() for item in arguments.distance_bands.split(",")
    )
    if (not distance_bands or any(not item for item in distance_bands)
            or len(set(distance_bands)) != len(distance_bands)
            or set(distance_bands) - set(DISTANCE_BAND_COLUMNS)):
        raise AnalysisError(
            "--distance-bands contains an empty, duplicate, or unknown band"
        )
    for path, label in (
        (source, "source"), (h3_package, "H3K4me3 package"),
        (evidence_package, "evidence package"),
        (context_package, "context-maxima package"),
        (annotation_catalog, "annotation catalog"), (runtime, "runtime"),
    ):
        if not path.is_dir():
            raise AnalysisError(f"{label} is missing: {path}")
    duckdb = runtime / "duckdb" / "bin" / "duckdb"
    if not duckdb.is_file():
        raise AnalysisError(f"runtime DuckDB is missing: {duckdb}")
    for name in ("plan", "logs", "tasks", "final"):
        (run_root / name).mkdir(parents=True, exist_ok=True)

    h3_manifest_path = h3_package / "manifest.json"
    h3_manifest = load_json(h3_manifest_path)
    if (h3_manifest.get("state") != "complete"
            or h3_manifest.get("primary_inference_partition") != "autosome"):
        raise AnalysisError("H3K4me3 package is not a complete autosome contract")
    h3_inventory_path = h3_package / "chromosome_file_inventory.tsv"
    h3_rows = [
        row for row in read_tsv(h3_inventory_path)
        if row["dataset"] == "h3k4me3_anchor_change"
        and row["analysis_partition"] == "autosome"
    ]
    if ({row["chrom"] for row in h3_rows} != set(AUTOSOMES)
            or len(h3_rows) != len(AUTOSOMES)):
        raise AnalysisError("H3K4me3 inventory does not contain 22 autosomes")

    annotation_manifest_path = annotation_catalog / "manifest.json"
    annotation_manifest = load_json(annotation_manifest_path)
    annotation_database = annotation_catalog / str(
        annotation_manifest.get("database", "context.duckdb")
    )
    if (annotation_manifest.get("state") != "complete"
            or int(annotation_manifest.get("context_schema_version", -1)) != 9
            or annotation_manifest.get("task_kind") != "anchor_annotation"
            or not annotation_database.is_file()):
        raise AnalysisError("annotation catalog is not completed schema 9")
    annotation_rows = query_json(duckdb, annotation_database, """
SELECT CAST(chrom AS VARCHAR) AS chrom, absolute_path, bytes
FROM context_file_inventory
WHERE dataset = 'tp73_context_anchor'
ORDER BY try_cast(chrom AS INTEGER), chrom;
""")
    annotation_rows = [
        row for row in annotation_rows if str(row["chrom"]) in AUTOSOMES
    ]
    if ({str(row["chrom"]) for row in annotation_rows} != set(AUTOSOMES)
            or len(annotation_rows) != len(AUTOSOMES)):
        raise AnalysisError("annotation catalog lacks one anchor file per autosome")

    evidence_manifest_path = evidence_package / "manifest.json"
    evidence_manifest = load_json(evidence_manifest_path)
    evidence = (
        evidence_package / "tables" / "tp73_anchor_evidence_autosome.parquet"
    )
    if (evidence_manifest.get("state") != "complete"
            or evidence_manifest.get("primary_inference_partition") != "autosome"
            or not evidence.is_file()):
        raise AnalysisError("TP73 evidence package lacks its autosome aggregate")

    fixed_rows: list[dict[str, Any]] = []
    for row in sorted(h3_rows, key=lambda item: int(item["sequence_order"])):
        path = resolve_package_file(h3_package, row["relative_path"])
        if not path.is_file() or path.stat().st_size != int(row["bytes"]):
            raise AnalysisError(f"H3K4me3 change file differs: {path}")
        fixed_rows.append({
            "kind": "h3k4me3_change", "chrom": row["chrom"],
            "path": str(path), "bytes": int(row["bytes"]),
            "sha256": row["sha256"],
        })
    for row in sorted(annotation_rows, key=lambda item: int(item["chrom"])):
        path = Path(str(row["absolute_path"])).resolve()
        if not path.is_file() or path.stat().st_size != int(row["bytes"]):
            raise AnalysisError(f"schema-9 annotation file differs: {path}")
        fixed_rows.append({
            "kind": "schema9_annotation", "chrom": row["chrom"],
            "path": str(path), "bytes": int(row["bytes"]),
            "sha256": sha256(path),
        })
    fixed_rows.append({
        "kind": "tp73_evidence", "chrom": "autosomes",
        "path": str(evidence.resolve()), "bytes": evidence.stat().st_size,
        "sha256": sha256(evidence),
    })
    fixed_path = run_root / "plan" / "fixed_inputs.tsv"
    write_rows(
        fixed_path, fixed_rows, ("kind", "chrom", "path", "bytes", "sha256")
    )

    context_manifest_path = context_package / "manifest.json"
    context_manifest = load_json(context_manifest_path)
    context_database = context_package / str(
        context_manifest.get("database", "tp73_genome_context_maxima.duckdb")
    )
    if (context_manifest.get("state") != "complete"
            or not context_database.is_file()):
        raise AnalysisError("context-maxima package is incomplete")
    task_rows = query_json(duckdb, context_database, """
SELECT i.task_index::BIGINT AS source_task_index,
       i.motif_id::VARCHAR AS motif_id,
       i.motif_name::VARCHAR AS factor_name,
       i.absolute_path::VARCHAR AS maxima_path,
       i.bytes::BIGINT AS maxima_bytes,
       i.sha256::VARCHAR AS maxima_sha256,
       i.scan_minimum_score::DOUBLE AS source_score_floor,
       i.applied_context_threshold::DOUBLE AS positive_threshold,
       coalesce(t.threshold_set_id, 'unknown')::VARCHAR AS threshold_set_id,
       coalesce(t.calibration_status, 'unknown')::VARCHAR AS calibration_status
FROM context_maxima_file_inventory i
LEFT JOIN motif_score_threshold t USING (motif_id)
ORDER BY i.task_index;
""")
    if not task_rows or len({str(row["motif_id"]) for row in task_rows}) != len(task_rows):
        raise AnalysisError("context-maxima inventory is empty or duplicates motifs")
    planned: list[dict[str, Any]] = []
    for task_index, row in enumerate(task_rows):
        motif_id = str(row["motif_id"])
        maxima = Path(str(row["maxima_path"])).resolve()
        if (motif_id == arguments.target_motif
                or not maxima.is_file()
                or maxima.stat().st_size != int(row["maxima_bytes"])
                or re.fullmatch(r"[0-9a-f]{64}", str(row["maxima_sha256"])) is None):
            raise AnalysisError(f"invalid context-maxima task: {motif_id}")
        source_floor = float(row["source_score_floor"])
        source_positive_threshold = float(row["positive_threshold"])
        positive_threshold = (
            arguments.fixed_positive_threshold
            if arguments.fixed_positive_threshold is not None
            else source_positive_threshold
        )
        if positive_threshold < source_floor:
            raise AnalysisError(
                f"positive threshold {positive_threshold:g} is below the "
                f"source floor {source_floor:g} for {motif_id}"
            )
        planned.append({
            "task_index": task_index,
            "source_task_index": int(row["source_task_index"]),
            "motif_id": motif_id,
            "factor_name": str(row["factor_name"]),
            "positive_threshold": format(positive_threshold, ".17g"),
            "source_applied_context_threshold": format(
                source_positive_threshold, ".17g"
            ),
            "source_score_floor": format(source_floor, ".17g"),
            "threshold_set_id": (
                "fixed_score_zero_v1"
                if arguments.fixed_positive_threshold is not None
                else str(row["threshold_set_id"])
            ),
            "calibration_status": (
                "fixed_common_score_zero_positive"
                if arguments.fixed_positive_threshold is not None
                else str(row["calibration_status"])
            ),
            "maxima_path": str(maxima),
            "maxima_bytes": int(row["maxima_bytes"]),
            "maxima_sha256": str(row["maxima_sha256"]),
        })
    if any(band != "all_150" for band in distance_bands):
        schema_rows = query_json(
            duckdb, ":memory:",
            "DESCRIBE SELECT * FROM read_parquet("
            f"{sql_string(planned[0]['maxima_path'])});",
        )
        columns = {str(row["column_name"]) for row in schema_rows}
        required_columns = {
            DISTANCE_BAND_COLUMNS[band] for band in distance_bands
        } | {"distance_band_schema_id"}
        missing_columns = sorted(required_columns - columns)
        if missing_columns:
            raise AnalysisError(
                "context-maxima package lacks distance-band columns: "
                + ", ".join(missing_columns)
            )
    task_path = run_root / "plan" / "cofactor_tasks.tsv"
    task_fields = (
        "task_index", "source_task_index", "motif_id", "factor_name",
        "positive_threshold", "source_applied_context_threshold",
        "source_score_floor", "threshold_set_id", "calibration_status",
        "maxima_path", "maxima_bytes", "maxima_sha256",
    )
    write_rows(task_path, planned, task_fields)

    commit, dirty = git_identity(source)
    config = {
        "schema_version": 1,
        "run_id": arguments.run_id,
        "analysis": "whole_autosome_h3k4me3_cofactor_change",
        "analysis_role": (
            "gfp_baseline_adjusted_sensitivity"
            if arguments.adjust_gfp_baseline
            else "exploratory_whole_autosome_chr1_threshold_transfer"
        ),
        "source": str(source),
        "source_commit": commit,
        "source_dirty": dirty,
        "scientific_source_file_sha256": scientific_hashes(source),
        "h3k4me3_package": str(h3_package),
        "h3k4me3_manifest_sha256": sha256(h3_manifest_path),
        "h3k4me3_inventory_sha256": sha256(h3_inventory_path),
        "evidence_package": str(evidence_package),
        "evidence_manifest_sha256": sha256(evidence_manifest_path),
        "context_maxima_package": str(context_package),
        "context_maxima_manifest_sha256": sha256(context_manifest_path),
        "annotation_catalog": str(annotation_catalog),
        "annotation_manifest_sha256": sha256(annotation_manifest_path),
        "runtime_prefix": str(runtime),
        "scratch_root": str(arguments.scratch_root),
        "target_motif": arguments.target_motif,
        "analysis_partition": "autosome",
        "chromosomes": list(AUTOSOMES),
        "primary_window": "flank_150_1000",
        "distance_bands": list(distance_bands),
        "positive_threshold_policy": (
            f"fixed_{arguments.fixed_positive_threshold:g}"
            if arguments.fixed_positive_threshold is not None
            else "context_package_applied_threshold"
        ),
        "negative_references": [-1.0, 0.0],
        "block_size_bp": arguments.block_size,
        "spline_df": arguments.spline_df,
        "adjust_gfp_baseline": arguments.adjust_gfp_baseline,
        "minimum_class_fraction": arguments.minimum_class_fraction,
        "minimum_class_count": arguments.minimum_class_count,
        "minimum_interaction_cell_count":
            arguments.minimum_interaction_cell_count,
        "task_count": len(planned),
        "motifs_per_batch": arguments.motifs_per_batch,
        "batch_count": math.ceil(len(planned) / arguments.motifs_per_batch),
        "minimum_free_run_gb": arguments.minimum_free_run_gb,
        "minimum_free_scratch_gb": arguments.minimum_free_scratch_gb,
        "fixed_input_plan_sha256": sha256(fixed_path),
        "task_plan_sha256": sha256(task_path),
        "multiple_testing_policy":
            "task_q_is_single_motif_diagnostic;final_q_is_all_planned_motifs",
        "tp73_confirmation_role": "secondary_descriptive_post_treatment",
        "strict_intergenic_definition":
            "no_transcript_promoter_or_downstream_region_overlap",
        "gene_relation_precedence":
            "promoter_then_downstream_then_gene_body_then_intergenic",
        "future_covariate_policy":
            "add_pretreatment_anchor_annotations_without_changing_outcome_grain",
    }
    immutable_write(
        run_root / "plan" / "run_config.json",
        json.dumps(config, indent=2, sort_keys=True) + "\n",
    )
    print(config["batch_count"])


def planned_tasks(run_root: Path) -> list[dict[str, str]]:
    rows = read_tsv(run_root / "plan" / "cofactor_tasks.tsv")
    if [int(row["task_index"]) for row in rows] != list(range(len(rows))):
        raise AnalysisError("task indices are not contiguous")
    if len({row["motif_id"] for row in rows}) != len(rows):
        raise AnalysisError("task plan duplicates motifs")
    return rows


def verify_plan(run_root: Path, config: dict[str, Any]) -> None:
    checks = (
        (run_root / "plan" / "fixed_inputs.tsv",
         config["fixed_input_plan_sha256"], "fixed-input plan"),
        (run_root / "plan" / "cofactor_tasks.tsv",
         config["task_plan_sha256"], "task plan"),
        (Path(config["h3k4me3_package"]) / "manifest.json",
         config["h3k4me3_manifest_sha256"], "H3K4me3 manifest"),
        (Path(config["evidence_package"]) / "manifest.json",
         config["evidence_manifest_sha256"], "evidence manifest"),
        (Path(config["context_maxima_package"]) / "manifest.json",
         config["context_maxima_manifest_sha256"], "context-maxima manifest"),
        (Path(config["annotation_catalog"]) / "manifest.json",
         config["annotation_manifest_sha256"], "annotation manifest"),
    )
    for path, digest, label in checks:
        if not path.is_file() or sha256(path) != digest:
            raise AnalysisError(f"{label} changed after planning: {path}")


def verify_file_row(row: dict[str, str], checksum: bool = True) -> Path:
    path = Path(row["path"]).resolve()
    if not path.is_file() or path.stat().st_size != int(row["bytes"]):
        raise AnalysisError(f"planned input differs: {path}")
    if checksum and sha256(path) != row["sha256"]:
        raise AnalysisError(f"planned input checksum differs: {path}")
    return path


def preflight(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(run_root / "plan" / "run_config.json")
    verify_plan(run_root, config)
    verify_scientific_hashes(config)
    fixed = read_tsv(run_root / "plan" / "fixed_inputs.tsv")
    for row in fixed:
        verify_file_row(row)
    changes = [Path(row["path"]) for row in fixed
               if row["kind"] == "h3k4me3_change"]
    annotations = [Path(row["path"]) for row in fixed
                   if row["kind"] == "schema9_annotation"]
    evidence_rows = [row for row in fixed if row["kind"] == "tp73_evidence"]
    if len(changes) != 22 or len(annotations) != 22 or len(evidence_rows) != 1:
        raise AnalysisError("fixed-input plan does not contain the nuclear contract")
    evidence = Path(evidence_rows[0]["path"])
    duckdb = Path(config["runtime_prefix"]) / "duckdb" / "bin" / "duckdb"
    values = query_json(duckdb, ":memory:", f"""
WITH h AS (
  SELECT chrom, anchor_start, anchor_end, series_id
  FROM read_parquet({sql_path_list(changes)}, hive_partitioning=false)
), e AS (
  SELECT chrom, anchor_start, anchor_end
  FROM read_parquet({sql_string(evidence)}, hive_partitioning=false)
), a AS (
  SELECT CAST(annotation.chrom AS VARCHAR) AS chrom,
         annotation.start AS anchor_start,
         annotation."end" AS anchor_end,
         count(DISTINCT primary_genomic_context) AS context_values,
         count(DISTINCT strict_intergenic) AS intergenic_values,
         count(DISTINCT overlaps_any_downstream_region) AS downstream_values,
         count(DISTINCT gene_relation_class) AS gene_relation_values,
         min(gene_relation_class) AS gene_relation_class,
         bool_or(strict_intergenic) AS strict_intergenic,
         bool_or(overlaps_any_promoter) AS overlaps_any_promoter,
         bool_or(overlaps_any_downstream_region)
           AS overlaps_any_downstream_region,
         bool_or(in_any_transcript) AS in_any_transcript,
         count(DISTINCT nearest_tss_distance_bp) AS tss_values,
         count(DISTINCT nearest_cds_distance_bp) AS cds_values,
         count(DISTINCT nearest_tss_genomic_distance_bp) AS tss_genomic_values,
         count(DISTINCT nearest_cds_genomic_distance_bp) AS cds_genomic_values,
         count(DISTINCT nearest_tss_has_mixed_strands) AS tss_mixed_values,
         count(DISTINCT nearest_cds_has_mixed_strands) AS cds_mixed_values,
         count(DISTINCT nearest_tss_relation) AS tss_relation_values,
         count(DISTINCT nearest_cds_relation) AS cds_relation_values
  -- tp73_context_anchor is partitioned by genome_id/chrom. DuckDB omits those
  -- columns from the Parquet payload and recovers them from the directory path.
  FROM read_parquet(
    {sql_path_list(annotations)}, hive_partitioning=true
  ) AS annotation
  GROUP BY CAST(annotation.chrom AS VARCHAR), annotation.start,
           annotation."end"
)
SELECT
  (SELECT count(*) FROM e)::BIGINT AS evidence_anchors,
  (SELECT count(DISTINCT (chrom, anchor_start, anchor_end)) FROM h)::BIGINT
    AS h3_anchors,
  (SELECT count(DISTINCT series_id) FROM h)::BIGINT AS h3_series,
  (SELECT count(*) FROM a)::BIGINT AS annotation_anchors,
  (SELECT count(*) FROM h)::BIGINT AS h3_rows,
  (SELECT count(*) - count(DISTINCT (chrom, anchor_start, anchor_end)) FROM e)
    ::BIGINT AS duplicate_evidence,
  (SELECT count(*) FROM a WHERE context_values <> 1 OR intergenic_values <> 1
      OR downstream_values <> 1 OR gene_relation_values <> 1
      OR tss_values > 1 OR cds_values > 1
      OR tss_genomic_values > 1 OR cds_genomic_values > 1
      OR tss_mixed_values > 1 OR cds_mixed_values > 1
      OR tss_relation_values > 1 OR cds_relation_values > 1
      OR gene_relation_class <> CASE
           WHEN overlaps_any_promoter THEN 'promoter'
           WHEN overlaps_any_downstream_region THEN 'downstream'
           WHEN in_any_transcript THEN 'gene_body'
           ELSE 'intergenic' END
      OR strict_intergenic <> (gene_relation_class = 'intergenic'))::BIGINT
    AS annotation_conflicts,
  (SELECT count(*) FROM (SELECT * FROM e EXCEPT
                         SELECT chrom, anchor_start, anchor_end FROM h))::BIGINT
    AS evidence_missing_h3,
  (SELECT count(*) FROM (SELECT * FROM e EXCEPT
                         SELECT chrom, anchor_start, anchor_end FROM a))::BIGINT
    AS evidence_missing_annotation,
  (SELECT count(*) FROM (SELECT chrom, anchor_start, anchor_end FROM a EXCEPT
                         SELECT * FROM e))::BIGINT AS annotation_missing_evidence;
""")
    if len(values) != 1:
        raise AnalysisError("preflight query returned no summary")
    result = {key: int(value) for key, value in values[0].items()}
    if (result["evidence_anchors"] <= 0
            or result["h3_anchors"] != result["evidence_anchors"]
            or result["annotation_anchors"] != result["evidence_anchors"]
            or result["h3_series"] != 2
            or result["h3_rows"] != result["evidence_anchors"] * 2
            or any(result[key] != 0 for key in (
                "duplicate_evidence", "annotation_conflicts",
                "evidence_missing_h3", "evidence_missing_annotation",
                "annotation_missing_evidence",
            ))):
        raise AnalysisError(f"whole-autosome preflight failed: {result}")
    marker = {
        "schema_version": 1,
        "state": "complete",
        "run_id": config["run_id"],
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "fixed_input_plan_sha256": config["fixed_input_plan_sha256"],
        "task_plan_sha256": config["task_plan_sha256"],
        "validation": result,
    }
    marker_path = run_root / "plan" / "preflight.json"
    if marker_path.exists():
        existing = load_json(marker_path)
        for key in (
            "schema_version", "state", "run_id", "fixed_input_plan_sha256",
            "task_plan_sha256", "validation",
        ):
            if existing.get(key) != marker.get(key):
                raise AnalysisError("existing preflight marker has another contract")
        print(json.dumps(result, sort_keys=True))
        return
    immutable_write(marker_path, json.dumps(marker, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


def task_directory(run_root: Path, row: dict[str, str]) -> Path:
    return run_root / "tasks" / (
        f"task-{int(row['task_index']):06d}-{safe_label(row['motif_id'])}"
    )


def output_inventory(root: Path) -> dict[str, dict[str, Any]]:
    return {
        path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
        for path in sorted(root.iterdir()) if path.is_file()
        and path.name != "complete.json"
    }


def validate_task(row: dict[str, str], directory: Path,
                  checksums: bool = False) -> dict[str, Any]:
    marker = load_json(directory / "complete.json")
    if (marker.get("state") != "complete"
            or marker.get("task_index") != int(row["task_index"])
            or marker.get("motif_id") != row["motif_id"]
            or marker.get("source_maxima_sha256") != row["maxima_sha256"]):
        raise AnalysisError(f"task identity differs: {directory}")
    records = marker.get("files")
    if not isinstance(records, dict) or set(records) != set(OUTPUTS.values()):
        raise AnalysisError(f"task output inventory differs: {directory}")
    for name, record in records.items():
        path = directory / name
        if (not path.is_file()
                or path.stat().st_size != int(record.get("bytes", -1))):
            raise AnalysisError(f"task output differs: {path}")
        if checksums and sha256(path) != record.get("sha256"):
            raise AnalysisError(f"task output checksum differs: {path}")
    return marker


def validate_result(
    prefix: Path, motif_id: str, expected_distance_bands: tuple[str, ...]
) -> None:
    band_count = len(expected_distance_bands)
    expected_counts = {
        "intensity_effect": 8 * band_count,
        "isoform_contrast": 4 * band_count,
        "tp73_interaction": 24 * band_count,
        "series_summary": 4 * band_count,
        "gene_relation_stratified_intensity_effect": 32 * band_count,
        "gene_relation_stratified_tp73_occupancy": 8 * band_count,
        "score_gradient": 8 * band_count,
        "run_config": 1,
    }
    for dataset, name in OUTPUTS.items():
        path = Path(f"{prefix}_{dataset}.tsv")
        rows = read_tsv(path)
        if dataset in expected_counts and len(rows) != expected_counts[dataset]:
            raise AnalysisError(
                f"{dataset} row count differs for {motif_id}: {len(rows)}"
            )
        if "motif_id" in rows[0] and any(
            row["motif_id"] != motif_id for row in rows
        ):
            raise AnalysisError(f"{dataset} contains another motif")
        if "distance_band" in rows[0] and {
            row["distance_band"] for row in rows
        } != set(expected_distance_bands):
            raise AnalysisError(f"{dataset} lacks a requested distance band")
        if dataset == "run_config":
            if (rows[0].get("schema_version") != "6"
                    or set(rows[0].get("chromosomes", "").split(",")) !=
                    set(AUTOSOMES)
                    or tuple(rows[0].get("distance_bands", "").split(",")) !=
                    expected_distance_bands):
                raise AnalysisError("evaluator run configuration is incomplete")


def canonicalize_run_config(prefix: Path, row: dict[str, str]) -> None:
    path = Path(f"{prefix}_run_config.tsv")
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        records = list(reader)
        fields = list(reader.fieldnames or ())
    if len(records) != 1 or not fields:
        raise AnalysisError("evaluator run configuration is not one TSV row")
    record = records[0]
    record.update({
        "change": "provenance/fixed_inputs.tsv#kind=h3k4me3_change",
        "tp73_evidence": "provenance/fixed_inputs.tsv#kind=tp73_evidence",
        "cofactor_maxima": (
            "provenance/cofactor_tasks.tsv#motif_id=" + row["motif_id"]
        ),
        "annotation": "provenance/fixed_inputs.tsv#kind=schema9_annotation",
        "thresholds": (
            "provenance/cofactor_tasks.tsv#motif_id=" + row["motif_id"]
        ),
        "input_reference_semantics": "package_provenance_selector",
        "execution_inputs_staged_to_scratch": "true",
    })
    for field in ("input_reference_semantics", "execution_inputs_staged_to_scratch"):
        if field not in fields:
            fields.append(field)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerow(record)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def run_batch(arguments: argparse.Namespace) -> None:
    global CURRENT_BATCH, CURRENT_MOTIF
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(run_root / "plan" / "run_config.json")
    verify_plan(run_root, config)
    verify_scientific_hashes(config)
    preflight_marker = load_json(run_root / "plan" / "preflight.json")
    if (preflight_marker.get("state") != "complete"
            or preflight_marker.get("run_id") != config["run_id"]
            or preflight_marker.get("task_plan_sha256") !=
                config["task_plan_sha256"]):
        raise AnalysisError("whole-autosome preflight is absent or stale")
    tasks = planned_tasks(run_root)
    batch_index = arguments.batch_index
    if batch_index is None:
        value = os.environ.get("SLURM_ARRAY_TASK_ID", "")
        if not value.isdigit():
            raise AnalysisError("--batch-index or SLURM_ARRAY_TASK_ID is required")
        batch_index = int(value) + arguments.batch_offset
    batch_size = int(config["motifs_per_batch"])
    first = batch_index * batch_size
    last = min(first + batch_size, len(tasks))
    if first < 0 or first >= len(tasks):
        raise AnalysisError(f"batch index is outside the task plan: {batch_index}")
    CURRENT_BATCH = str(batch_index)
    signal.signal(signal.SIGUSR1, progress)
    batch_tasks = tasks[first:last]
    completed = 0
    for row in batch_tasks:
        final = task_directory(run_root, row)
        if final.exists():
            validate_task(row, final)
            completed += 1
    if completed == len(batch_tasks):
        set_phase("complete_reused")
        print(
            f"I: reusing completed H3K4me3 batch {batch_index} "
            f"({len(batch_tasks)} motifs)", file=sys.stderr,
        )
        return

    check_free_space(run_root, float(config["minimum_free_run_gb"]), "durable")
    scratch_root = Path(config["scratch_root"]).expanduser().resolve()
    scratch_root.mkdir(parents=True, exist_ok=True)
    check_free_space(
        scratch_root, float(config["minimum_free_scratch_gb"]), "scratch"
    )
    scratch = scratch_root / (
        f"h3-cofactor-{os.environ.get('SLURM_JOB_ID', 'manual')}-"
        f"batch-{batch_index}-restart-{os.environ.get('SLURM_RESTART_COUNT', '0')}-"
        f"pid-{os.getpid()}"
    )
    scratch.mkdir()
    try:
        fixed = read_tsv(run_root / "plan" / "fixed_inputs.tsv")
        staged_change: list[Path] = []
        staged_annotation: list[Path] = []
        staged_evidence: Path | None = None
        set_phase("staging_fixed_inputs")
        for row in fixed:
            source = verify_file_row(row, checksum=False)
            directory = scratch / row["kind"]
            if row["kind"] == "schema9_annotation":
                # The partitioned annotation payload omits chrom physically.
                # Preserve it as a Hive path while staging to node-local scratch.
                directory /= f"chrom={safe_label(row['chrom'])}"
                target = directory / "data.parquet"
            else:
                target = directory / f"chrom-{safe_label(row['chrom'])}.parquet"
            directory.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            if sha256(target) != row["sha256"]:
                raise AnalysisError(f"staged input checksum differs: {target}")
            if row["kind"] == "h3k4me3_change":
                staged_change.append(target)
            elif row["kind"] == "schema9_annotation":
                staged_annotation.append(target)
            elif row["kind"] == "tp73_evidence":
                staged_evidence = target
        if (len(staged_change) != 22 or len(staged_annotation) != 22
                or staged_evidence is None):
            raise AnalysisError("staged fixed input set is incomplete")

        runtime = Path(config["runtime_prefix"])
        duckdb = runtime / "duckdb" / "bin" / "duckdb"
        rscript = shutil.which(arguments.rscript)
        if rscript is None:
            raise AnalysisError(f"Rscript is unavailable: {arguments.rscript}")
        for row in batch_tasks:
            CURRENT_MOTIF = row["motif_id"]
            final = task_directory(run_root, row)
            if final.exists():
                validate_task(row, final)
                print(f"I: reusing completed H3K4me3 task {row['task_index']} "
                      f"({row['motif_id']})", file=sys.stderr)
                continue
            maxima_source = Path(row["maxima_path"])
            if (not maxima_source.is_file()
                    or maxima_source.stat().st_size != int(row["maxima_bytes"])):
                raise AnalysisError(f"cofactor maxima differs: {maxima_source}")
            maxima = scratch / "cofactor-maxima.parquet"
            if maxima.exists():
                maxima.unlink()
            set_phase("staging_motif_maxima")
            shutil.copy2(maxima_source, maxima)
            if sha256(maxima) != row["maxima_sha256"]:
                raise AnalysisError(f"cofactor maxima checksum differs: {maxima}")
            threshold = scratch / "threshold.tsv"
            with threshold.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(
                    stream,
                    fieldnames=(
                        "motif_id", "positive_threshold", "factor_name",
                        "positive_threshold_source", "selection_semantics",
                    ),
                    delimiter="\t", lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow({
                    "motif_id": row["motif_id"],
                    "positive_threshold": row["positive_threshold"],
                    "factor_name": row["factor_name"],
                    "positive_threshold_source": row["threshold_set_id"],
                    "selection_semantics": row["calibration_status"],
                })
            prefix = scratch / "result"
            for dataset in OUTPUTS:
                output = Path(f"{prefix}_{dataset}.tsv")
                if output.exists():
                    output.unlink()
            command = [
                rscript,
                str(Path(config["source"]) /
                    "scripts" / "analyze_h3k4me3_cofactor_change.R"),
            ]
            for path in staged_change:
                command.extend(["--change", str(path)])
            command.extend(["--tp73-evidence", str(staged_evidence)])
            command.extend(["--cofactor-maxima", str(maxima)])
            for path in staged_annotation:
                command.extend(["--annotation", str(path)])
            command.extend([
                "--thresholds", str(threshold),
                "--output-prefix", str(prefix),
                "--negative-references", "-1,0",
                "--distance-bands", ",".join(config["distance_bands"]),
                "--block-size", str(config["block_size_bp"]),
                "--spline-df", str(config["spline_df"]),
                "--minimum-class-fraction",
                str(config["minimum_class_fraction"]),
                "--minimum-class-count", str(config["minimum_class_count"]),
                "--minimum-interaction-cell-count",
                str(config["minimum_interaction_cell_count"]),
                "--duckdb", str(duckdb),
                "--analysis-role", str(config["analysis_role"]),
            ])
            if config.get("adjust_gfp_baseline", False):
                command.append("--adjust-gfp-baseline")
            set_phase("evaluating_motif")
            run_process(command, cwd=Path(config["source"]))
            canonicalize_run_config(prefix, row)
            set_phase("validating_motif")
            validate_result(
                prefix, row["motif_id"], tuple(config["distance_bands"])
            )

            attempt = run_root / "tasks" / (
                f".attempt-task-{int(row['task_index']):06d}-"
                f"{safe_label(row['motif_id'])}-job-"
                f"{os.environ.get('SLURM_JOB_ID', 'manual')}-"
                f"restart-{os.environ.get('SLURM_RESTART_COUNT', '0')}-"
                f"pid-{os.getpid()}"
            )
            attempt.mkdir()
            for dataset, name in OUTPUTS.items():
                shutil.copy2(Path(f"{prefix}_{dataset}.tsv"), attempt / name)
            marker = {
                "schema_version": 1,
                "state": "complete",
                "task_index": int(row["task_index"]),
                "motif_id": row["motif_id"],
                "source_maxima_sha256": row["maxima_sha256"],
                "positive_threshold": float(row["positive_threshold"]),
                "source_commit": config["source_commit"],
                "completed_at_utc": datetime.now(timezone.utc).isoformat(),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
                "slurm_restart_count":
                    int(os.environ.get("SLURM_RESTART_COUNT", "0")),
                "files": output_inventory(attempt),
            }
            (attempt / "complete.json").write_text(
                json.dumps(marker, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            set_phase("publishing_motif")
            try:
                os.rename(attempt, final)
            except OSError:
                if not final.exists():
                    raise
                validate_task(row, final)
                shutil.rmtree(attempt)
        set_phase("complete")
        CURRENT_MOTIF = "none"
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def status(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    rows = planned_tasks(run_root)
    complete = invalid = 0
    for row in rows:
        directory = task_directory(run_root, row)
        if not (directory / "complete.json").is_file():
            continue
        try:
            validate_task(row, directory)
            complete += 1
        except (AnalysisError, OSError, ValueError, json.JSONDecodeError):
            invalid += 1
    config = load_json(run_root / "plan" / "run_config.json")
    print("key\tvalue")
    print(f"planned_motifs\t{len(rows)}")
    print(f"planned_batches\t{config['batch_count']}")
    print(f"complete_motifs\t{complete}")
    print(f"pending_motifs\t{len(rows) - complete - invalid}")
    print(f"invalid_complete_markers\t{invalid}")


def bh_sql(source_sql: str, partitions: tuple[str, ...]) -> str:
    partition = ", ".join(partitions)
    return f"""
WITH raw_source AS ({source_sql}), source AS (
  SELECT * REPLACE (
    try_cast(p_value AS DOUBLE) AS p_value,
    try_cast(q_value_bh AS DOUBLE) AS q_value_bh
  )
  FROM raw_source
), ranked AS (
  SELECT * EXCLUDE (q_value_bh), q_value_bh AS q_value_bh_task,
         row_number() OVER (
           PARTITION BY {partition}
           ORDER BY p_value ASC NULLS LAST, motif_id ASC
         )::DOUBLE AS bh_rank,
         count(p_value) OVER (PARTITION BY {partition})::DOUBLE AS bh_count
  FROM source
), adjusted AS (
  SELECT * EXCLUDE (bh_rank, bh_count),
         CASE WHEN p_value IS NULL THEN NULL ELSE least(
           1.0,
           min(p_value * bh_count / bh_rank) OVER (
             PARTITION BY {partition}
             ORDER BY p_value DESC NULLS LAST, motif_id DESC
             ROWS BETWEEN UNBOUNDED PRECEDING AND CURRENT ROW
           )
         ) END AS q_value_bh_all_motifs
  FROM ranked
)
SELECT * FROM adjusted
"""


def finalize(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_json(run_root / "plan" / "run_config.json")
    verify_plan(run_root, config)
    verify_scientific_hashes(config)
    rows = planned_tasks(run_root)
    manifests = [
        validate_task(row, task_directory(run_root, row), checksums=True)
        for row in rows
    ]
    final = run_root / "final" / "h3k4me3_cofactor_analysis"
    if final.exists():
        manifest = load_json(final / "manifest.json")
        if (manifest.get("state") == "complete"
                and manifest.get("run_id") == config["run_id"]):
            print(f"I: reusing finalized H3K4me3 cofactor analysis: {final}",
                  file=sys.stderr)
            return
        raise AnalysisError(f"final output has another identity: {final}")
    final.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".h3-cofactor.", dir=final.parent))
    try:
        inventory_rows: list[dict[str, Any]] = []
        paths_by_dataset: dict[str, list[Path]] = {
            dataset: [] for dataset in OUTPUTS
        }
        for row, marker in zip(rows, manifests):
            directory = task_directory(run_root, row)
            for dataset, name in OUTPUTS.items():
                path = directory / name
                paths_by_dataset[dataset].append(path)
                record = marker["files"][name]
                inventory_rows.append({
                    "task_index": int(row["task_index"]),
                    "motif_id": row["motif_id"],
                    "dataset": dataset,
                    "absolute_path": str(path.resolve()),
                    "bytes": int(record["bytes"]),
                    "sha256": record["sha256"],
                })
        inventory_path = staging / "task_output_inventory.tsv"
        with inventory_path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=tuple(inventory_rows[0]), delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(inventory_rows)

        provenance = staging / "provenance"
        provenance.mkdir()
        for name in (
            "run_config.json", "fixed_inputs.tsv", "cofactor_tasks.tsv",
            "preflight.json",
        ):
            shutil.copy2(run_root / "plan" / name, provenance)

        duckdb = arguments.duckdb.expanduser().resolve()
        table_root = staging / "tables"
        table_root.mkdir()
        for dataset, paths in paths_by_dataset.items():
            source_sql = (
                "SELECT * FROM read_csv_auto("
                f"{sql_path_list(paths)}, delim='\\t', header=true, "
                "nullstr='NA', union_by_name=true, sample_size=-1)"
            )
            selected = (
                bh_sql(source_sql, BH_PARTITIONS[dataset])
                if dataset in BH_PARTITIONS else source_sql
            )
            output = table_root / f"{dataset}.parquet"
            run_sql(duckdb, ":memory:", f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET temp_directory={sql_string(arguments.temp_directory)};
COPY ({selected}) TO {sql_string(output)}
  (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 131072);
""")

        database = staging / "h3k4me3_cofactor_analysis.duckdb"
        statements = [
            "CREATE VIEW task_output_inventory AS SELECT * FROM "
            "read_csv_auto('task_output_inventory.tsv', delim='\\t', header=true);"
        ]
        for dataset in OUTPUTS:
            statements.append(
                f"CREATE VIEW {dataset} AS SELECT * FROM "
                f"read_parquet('tables/{dataset}.parquet');"
            )
        run_sql(duckdb, database, "\n".join(statements), cwd=staging)

        validation = query_json(duckdb, database, """
SELECT
  (SELECT count(DISTINCT motif_id) FROM intensity_effect)::BIGINT AS motifs,
  (SELECT count(*) FROM intensity_effect)::BIGINT AS intensity_rows,
  (SELECT count(*) FROM isoform_contrast)::BIGINT AS isoform_contrast_rows,
  (SELECT count(*) FROM tp73_interaction)::BIGINT AS interaction_rows,
  (SELECT count(*) FROM gene_relation_stratified_intensity_effect)::BIGINT
    AS gene_relation_rows,
  (SELECT count(*) FROM gene_relation_stratified_tp73_occupancy)::BIGINT
    AS gene_relation_occupancy_rows,
  (SELECT count(*) FROM score_gradient)::BIGINT AS score_rows,
  (SELECT count(*) FROM run_config)::BIGINT AS run_config_rows,
  (SELECT count(*) FROM intensity_effect
    WHERE p_value IS NOT NULL AND q_value_bh_all_motifs IS NULL)::BIGINT
    AS missing_global_q,
  (SELECT count(*) FROM isoform_contrast
    WHERE p_value IS NOT NULL AND q_value_bh_all_motifs IS NULL)::BIGINT
    AS missing_isoform_contrast_q,
  (SELECT count(*) FROM gene_relation_stratified_intensity_effect
    WHERE p_value IS NOT NULL AND q_value_bh_all_motifs IS NULL)::BIGINT
    AS missing_gene_relation_q,
  (SELECT count(*) FROM gene_relation_stratified_tp73_occupancy
    WHERE p_value IS NOT NULL AND q_value_bh_all_motifs IS NULL)::BIGINT
    AS missing_gene_relation_occupancy_q;
""", cwd=staging)
        if len(validation) != 1:
            raise AnalysisError("final validation returned no summary")
        values = {key: int(value) for key, value in validation[0].items()}
        task_count = int(config["task_count"])
        distance_band_count = len(config["distance_bands"])
        if (values["motifs"] != task_count
                or values["intensity_rows"] != task_count * 8 * distance_band_count
                or values["isoform_contrast_rows"] !=
                    task_count * 4 * distance_band_count
                or values["interaction_rows"] != task_count * 24 * distance_band_count
                or values["gene_relation_rows"] !=
                    task_count * 32 * distance_band_count
                or values["gene_relation_occupancy_rows"] !=
                    task_count * 8 * distance_band_count
                or values["score_rows"] != task_count * 8 * distance_band_count
                or values["run_config_rows"] != task_count
                or values["missing_global_q"] != 0
                or values["missing_isoform_contrast_q"] != 0
                or values["missing_gene_relation_q"] != 0
                or values["missing_gene_relation_occupancy_q"] != 0):
            raise AnalysisError(f"final analysis validation failed: {values}")
        manifest = {
            "schema_version": 1,
            "state": "complete",
            "run_id": config["run_id"],
            "analysis": config["analysis"],
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "source_commit": config["source_commit"],
            "analysis_partition": config["analysis_partition"],
            "primary_window": config["primary_window"],
            "distance_bands": config["distance_bands"],
            "positive_threshold_policy": config["positive_threshold_policy"],
            "motifs": task_count,
            "multiple_testing_scope":
                "all_planned_non_TP73_JASPAR_motifs_by_declared_result_family",
            "validation": values,
            "database": database.name,
            "outputs": {
                str(path.relative_to(staging)): {
                    "bytes": path.stat().st_size, "sha256": sha256(path)
                }
                for path in sorted(item for item in staging.rglob("*")
                                   if item.is_file())
            },
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, final)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    print(f"I: finalized all-motif H3K4me3 cofactor analysis: {final}",
          file=sys.stderr)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)

    prepare_parser = commands.add_parser(
        "prepare", help="write immutable fixed-input and cofactor task plans"
    )
    prepare_parser.add_argument("--run-root", type=Path, required=True)
    prepare_parser.add_argument("--h3-package", type=Path, required=True)
    prepare_parser.add_argument("--evidence-package", type=Path, required=True)
    prepare_parser.add_argument("--context-maxima-package", type=Path, required=True)
    prepare_parser.add_argument("--annotation-catalog", type=Path, required=True)
    prepare_parser.add_argument("--runtime-prefix", type=Path, required=True)
    prepare_parser.add_argument("--source", type=Path, required=True)
    prepare_parser.add_argument("--run-id", required=True)
    prepare_parser.add_argument("--target-motif", default="MA0861.2")
    prepare_parser.add_argument(
        "--distance-bands", default="all_150",
        help=(
            "comma-separated all_150/overlap/adjacent_0_5/gap_6_20/"
            "gap_21_50/gap_51_100/gap_101_150 selections"
        ),
    )
    prepare_parser.add_argument(
        "--fixed-positive-threshold", type=float,
        help="override every motif's positive threshold (use 0 for common-score analysis)",
    )
    prepare_parser.add_argument("--motifs-per-batch", type=int, default=8)
    prepare_parser.add_argument("--scratch-root", type=Path, required=True)
    prepare_parser.add_argument("--block-size", type=int, default=5_000_000)
    prepare_parser.add_argument("--spline-df", type=int, default=3)
    prepare_parser.add_argument(
        "--adjust-gfp-baseline", action="store_true",
        help="fit the GFP-baseline-adjusted sensitivity model",
    )
    prepare_parser.add_argument("--minimum-class-fraction", type=float,
                                default=0.005)
    prepare_parser.add_argument("--minimum-class-count", type=int, default=100)
    prepare_parser.add_argument("--minimum-interaction-cell-count", type=int,
                                default=100)
    prepare_parser.add_argument("--minimum-free-run-gb", type=float, default=5)
    prepare_parser.add_argument("--minimum-free-scratch-gb", type=float,
                                default=10)
    prepare_parser.set_defaults(function=prepare)

    preflight_parser = commands.add_parser(
        "preflight", help="validate whole-autosome anchor-key equivalence"
    )
    preflight_parser.add_argument("--run-root", type=Path, required=True)
    preflight_parser.set_defaults(function=preflight)

    batch_parser = commands.add_parser(
        "run-batch", help="evaluate one restart-safe bounded motif batch"
    )
    batch_parser.add_argument("--run-root", type=Path, required=True)
    batch_parser.add_argument("--batch-index", type=int)
    batch_parser.add_argument("--batch-offset", type=int, default=0)
    batch_parser.add_argument("--rscript", default="Rscript")
    batch_parser.set_defaults(function=run_batch)

    status_parser = commands.add_parser(
        "status", help="report exact completed motif counts"
    )
    status_parser.add_argument("--run-root", type=Path, required=True)
    status_parser.set_defaults(function=status)

    finalize_parser = commands.add_parser(
        "finalize", help="publish combined all-motif tables and DuckDB catalog"
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
        for name in (
            "motifs_per_batch", "block_size", "spline_df", "minimum_class_count",
            "minimum_interaction_cell_count", "threads",
        ):
            if hasattr(arguments, name) and getattr(arguments, name) <= 0:
                raise AnalysisError(f"--{name.replace('_', '-')} must be positive")
        if hasattr(arguments, "minimum_class_fraction") and not (
            0 < arguments.minimum_class_fraction < 0.5
        ):
            raise AnalysisError("--minimum-class-fraction must be in (0,0.5)")
        if (hasattr(arguments, "fixed_positive_threshold")
                and arguments.fixed_positive_threshold is not None
                and not math.isfinite(arguments.fixed_positive_threshold)):
            raise AnalysisError("--fixed-positive-threshold must be finite")
        arguments.function(arguments)
        return 0
    except (AnalysisError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        if arguments.command == "run-batch":
            progress()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

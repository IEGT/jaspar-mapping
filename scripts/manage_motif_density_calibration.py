#!/usr/bin/env python3

"""Prepare, run, inspect, and finalize motif-density threshold calibration."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from build_density_capped_thresholds import (
    ThresholdError,
    apply_negative_sensitivity,
    apply_overrides,
    read_distribution,
    read_informative_thresholds,
)
from manage_genome_scan import (
    ScanError,
    absolute_file,
    check_free_space,
    file_stat_identity,
    inspect_fasta,
    json_read,
    json_write_new,
    parse_fasta_index,
    parse_jaspar,
    safe_label,
    sha256_file,
)
from stage_fasta_region import StageError, stage_fasta_region


class DensityCalibrationError(RuntimeError):
    pass


COMMIT = re.compile(r"^[0-9a-f]{40}$")
PLAN_SCHEMA_VERSION = 1


def immutable_write(path: Path, content: str) -> None:
    encoded = content.encode("ascii")
    if path.exists():
        if path.read_bytes() != encoded:
            raise DensityCalibrationError(f"existing immutable file differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("xb") as stream:
        stream.write(encoded)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def input_identity(path: Path) -> dict[str, Any]:
    stat = file_stat_identity(path)
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "bytes": stat["file_bytes"],
        "mtime_ns": stat["file_mtime_ns"],
    }


def validate_identity(identity: dict[str, Any], label: str) -> Path:
    path = absolute_file(identity["path"])
    stat = path.stat()
    if stat.st_size != identity["bytes"] or stat.st_mtime_ns != identity["mtime_ns"]:
        raise DensityCalibrationError(f"planned {label} size or mtime changed: {path}")
    if sha256_file(path) != identity["sha256"]:
        raise DensityCalibrationError(f"planned {label} checksum changed: {path}")
    return path


def load_plan(run_root: Path) -> tuple[dict[str, Any], Path, str]:
    path = run_root / "plan" / "density_plan.json"
    plan = json_read(path)
    if plan.get("schema_version") != PLAN_SCHEMA_VERSION:
        raise DensityCalibrationError(f"unsupported density plan schema: {path}")
    return plan, path, sha256_file(path)


def validate_plan_inputs(plan: dict[str, Any]) -> None:
    genome = absolute_file(plan["genome"]["path"])
    metadata = genome.stat()
    if (metadata.st_size != plan["genome"]["bytes"]
            or metadata.st_mtime_ns != plan["genome"]["mtime_ns"]):
        raise DensityCalibrationError("planned genome size or mtime changed")
    validate_identity(plan["fasta_index"], "FASTA index")
    validate_identity(plan["jaspar"], "JASPAR source")
    for label, identity in plan["threshold_inputs"].items():
        if identity is not None:
            validate_identity(identity, label.replace("_", " "))


def prepare(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan_path = run_root / "plan" / "density_plan.json"
    if plan_path.exists():
        raise DensityCalibrationError(f"density plan already exists: {plan_path}")
    if not COMMIT.fullmatch(arguments.source_commit):
        raise DensityCalibrationError("--source-commit must be a full lowercase Git SHA")
    if arguments.motif_batch_size < 1:
        raise DensityCalibrationError("--motif-batch-size must be positive")
    if arguments.minimum_spacing_bp <= 0:
        raise DensityCalibrationError("--minimum-spacing-bp must be positive")
    if (not math.isfinite(arguments.pseudocount) or arguments.pseudocount < 0
            or not math.isfinite(arguments.default_minimum_score)
            or abs(arguments.default_minimum_score
                   - round(arguments.default_minimum_score)) > 1e-9):
        raise DensityCalibrationError(
            "pseudocount must be non-negative and the default score an integer"
        )

    genome = absolute_file(arguments.genome)
    fasta_index = absolute_file(arguments.fasta_index)
    jaspar = absolute_file(arguments.jaspar)
    informative = absolute_file(arguments.informative_thresholds)
    negative = absolute_file(arguments.negative_sensitivity) \
        if arguments.negative_sensitivity else None
    overrides = absolute_file(arguments.override_thresholds) \
        if arguments.override_thresholds else None

    motifs = parse_jaspar(jaspar)
    motif_ids = [str(row["motif_id"]) for row in motifs]
    threshold_rows = read_informative_thresholds(informative)
    apply_negative_sensitivity(threshold_rows, negative)
    apply_overrides(threshold_rows, overrides)
    if set(threshold_rows) != set(motif_ids):
        raise DensityCalibrationError(
            "informative thresholds plus overrides do not exactly cover JASPAR; "
            f"missing={sorted(set(motif_ids) - set(threshold_rows))[:10]}, "
            f"extra={sorted(set(threshold_rows) - set(motif_ids))[:10]}"
        )

    fasta_sha256, fasta_regions, sequence_set_sha256 = inspect_fasta(genome)
    index_regions = parse_fasta_index(fasta_index)
    fasta_by_name = {row["chrom"]: row for row in fasta_regions}
    index_by_name = {row["chrom"]: row for row in index_regions}
    if arguments.chrom not in fasta_by_name or arguments.chrom not in index_by_name:
        raise DensityCalibrationError(f"chromosome is absent from FASTA/index: {arguments.chrom}")
    if [row["chrom"] for row in fasta_regions] != [row["chrom"] for row in index_regions]:
        raise DensityCalibrationError("FASTA and index sequence order differ")
    region = dict(index_by_name[arguments.chrom])
    observed_region = fasta_by_name[arguments.chrom]
    if region["length"] != observed_region["length"]:
        raise DensityCalibrationError("FASTA and index chromosome lengths differ")
    region["sequence_sha256"] = observed_region["sequence_sha256"]

    run_root.mkdir(parents=True, exist_ok=True)
    check_free_space(run_root, arguments.minimum_free_gib)
    motif_list_dir = run_root / "plan" / "motif_lists"
    motif_list_dir.mkdir(parents=True, exist_ok=False)
    for directory in ("tasks", "staging", "logs"):
        (run_root / directory).mkdir(exist_ok=True)

    batches: list[dict[str, Any]] = []
    for offset in range(0, len(motif_ids), arguments.motif_batch_size):
        selected = motif_ids[offset:offset + arguments.motif_batch_size]
        batch_id = f"batch-{offset // arguments.motif_batch_size:04d}"
        relative = Path("plan") / "motif_lists" / f"{batch_id}.txt"
        content = "".join(f"{motif_id}\n" for motif_id in selected)
        immutable_write(run_root / relative, content)
        batches.append({
            "batch_index": len(batches),
            "batch_id": batch_id,
            "motif_ids": selected,
            "motif_count": len(selected),
            "motif_list": str(relative),
            "motif_list_sha256": sha256_file(run_root / relative),
        })

    genome_stat = file_stat_identity(genome)
    threshold_inputs = {
        "informative_thresholds": input_identity(informative),
        "negative_sensitivity": input_identity(negative) if negative else None,
        "override_thresholds": input_identity(overrides) if overrides else None,
    }
    plan = {
        "schema_version": PLAN_SCHEMA_VERSION,
        "run_id": arguments.run_id,
        "source_commit": arguments.source_commit,
        "minimum_free_gib": arguments.minimum_free_gib,
        "threshold_set_id": arguments.threshold_set_id,
        "genome_id": arguments.genome_id,
        "motif_set_id": arguments.motif_set_id,
        "genome": {
            "path": str(genome),
            "sha256": fasta_sha256,
            "bytes": genome_stat["file_bytes"],
            "mtime_ns": genome_stat["file_mtime_ns"],
            "sequence_set_sha256": sequence_set_sha256,
        },
        "fasta_index": input_identity(fasta_index),
        "jaspar": input_identity(jaspar),
        "threshold_inputs": threshold_inputs,
        "chromosome": region,
        "motif_count": len(motif_ids),
        "score_mode": arguments.score_mode,
        "pseudocount": arguments.pseudocount,
        "distribution_bin_width": 1,
        "orientation_aggregation": "max_score_per_alignment_span",
        "default_minimum_score": arguments.default_minimum_score,
        "minimum_spacing_bp": arguments.minimum_spacing_bp,
        "batches": batches,
    }
    json_write_new(plan_path, plan)
    print(f"I: Wrote immutable density plan: {plan_path}", file=sys.stderr)
    print(
        f"I: Planned {len(batches)} batches for {len(motif_ids)} motifs on "
        f"chromosome {arguments.chrom}.",
        file=sys.stderr,
    )
    print(f"BATCH_COUNT={len(batches)}")


def batch_directory(run_root: Path, batch: dict[str, Any]) -> Path:
    return run_root / "tasks" / batch["batch_id"]


def validate_complete_batch(
    run_root: Path, plan: dict[str, Any], plan_sha256: str,
    batch: dict[str, Any], *, verify_hashes: bool,
) -> dict[str, Any] | None:
    directory = batch_directory(run_root, batch)
    marker = directory / "complete.json"
    if not marker.is_file():
        return None
    value = json_read(marker)
    if (value.get("batch_id") != batch["batch_id"]
            or value.get("plan_sha256") != plan_sha256
            or value.get("motif_ids") != batch["motif_ids"]):
        raise DensityCalibrationError(f"invalid complete marker: {marker}")
    files = value.get("files")
    if not isinstance(files, list) or len(files) != batch["motif_count"]:
        raise DensityCalibrationError(f"incomplete file inventory: {marker}")
    by_motif = {str(row.get("motif_id")): row for row in files}
    if set(by_motif) != set(batch["motif_ids"]):
        raise DensityCalibrationError(f"batch motif inventory differs: {marker}")
    for motif_id, row in by_motif.items():
        relative = Path(str(row["relative_path"]))
        if relative.is_absolute() or ".." in relative.parts:
            raise DensityCalibrationError(f"unsafe distribution path for {motif_id}")
        path = directory / relative
        if not path.is_file() or path.stat().st_size != int(row["bytes"]):
            raise DensityCalibrationError(f"missing distribution for {motif_id}: {path}")
        if verify_hashes and sha256_file(path) != row["sha256"]:
            raise DensityCalibrationError(f"distribution checksum mismatch: {path}")
    return value


def validate_scanner(
    scanner: Path, plan: dict[str, Any], allow_mismatch: bool,
) -> dict[str, Any]:
    process = subprocess.run(
        [str(scanner), "--version-json"], text=True, capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise DensityCalibrationError(
            process.stderr.strip() or "scanner --version-json failed"
        )
    try:
        build = json.loads(process.stdout)
    except json.JSONDecodeError as error:
        raise DensityCalibrationError(
            "scanner --version-json returned invalid JSON"
        ) from error
    required = {
        "program", "source_commit", "source_dirty", "build_flags",
        "lto_enabled", "ndebug", "parquet_enabled",
    }
    if (not isinstance(build, dict) or not required.issubset(build)
            or build.get("program") != "pssm_scan"):
        raise DensityCalibrationError("scanner build information is incomplete")
    matches = (
        build.get("source_commit") == plan["source_commit"]
        and not build.get("source_dirty", True)
    )
    if not matches and not allow_mismatch:
        raise DensityCalibrationError(
            "scanner build is dirty or does not match the density plan"
        )
    return build


def distribution_motif_id(path: Path) -> str:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        first = next(reader, None)
    if first is None or not first.get("MotifID"):
        raise DensityCalibrationError(f"distribution has no motif row: {path}")
    return first["MotifID"]


def run_batch(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan, _, plan_sha256 = load_plan(run_root)
    if not arguments.planned_inputs_verified:
        validate_plan_inputs(plan)
    if arguments.batch_index < 0 or arguments.batch_index >= len(plan["batches"]):
        raise DensityCalibrationError("batch index is outside the density plan")
    batch = plan["batches"][arguments.batch_index]
    if validate_complete_batch(
        run_root, plan, plan_sha256, batch, verify_hashes=True
    ) is not None:
        print(f"I: Reusing completed density batch {batch['batch_id']}", file=sys.stderr)
        return

    scanner = absolute_file(arguments.scanner)
    scanner_build = validate_scanner(
        scanner, plan, arguments.allow_scanner_provenance_mismatch
    )
    staged_fasta = absolute_file(arguments.staged_fasta)
    staged_index = absolute_file(arguments.staged_fasta_index)
    staged_jaspar = absolute_file(arguments.staged_jaspar)
    if sha256_file(staged_jaspar) != plan["jaspar"]["sha256"]:
        raise DensityCalibrationError("staged JASPAR checksum differs from plan")
    staged_metadata = json_read(absolute_file(arguments.staged_metadata))
    region = plan["chromosome"]
    if (staged_metadata.get("sequence") != region["chrom"]
            or staged_metadata.get("length") != region["length"]
            or staged_metadata.get("sequence_sha256") != region["sequence_sha256"]):
        raise DensityCalibrationError("staged chromosome metadata differs from plan")

    attempt_parent = run_root / "staging" / batch["batch_id"]
    attempt_parent.mkdir(parents=True, exist_ok=True)
    attempt = attempt_parent / (
        f"attempt-job-{os.environ.get('SLURM_JOB_ID', 'local')}-"
        f"restart-{os.environ.get('SLURM_RESTART_COUNT', '0')}-pid-{os.getpid()}"
    )
    attempt.mkdir(exist_ok=False)
    work_parent = arguments.task_work_directory.expanduser().resolve()
    work_parent.mkdir(parents=True, exist_ok=True)
    work = work_parent / f"{batch['batch_id']}-pid-{os.getpid()}"
    work.mkdir(exist_ok=False)
    output = work / "distributions"
    output.mkdir()

    command = [
        str(scanner), "--genome", str(staged_fasta),
        "--fasta-index", str(staged_index),
        "--pssm", str(staged_jaspar),
        "--motif-list", str(run_root / batch["motif_list"]),
        "--chr", region["chrom"], "--strand", "both",
        "--score-mode", plan["score_mode"],
        "--pseudocount", str(plan["pseudocount"]),
        "--score-distribution", "--collapse-orientations",
        "--distribution-bin-width", "1", "--skip-N",
        "--outdir", str(output),
    ]
    active: subprocess.Popen[str] | None = None

    def forward_progress(_signal_number: int, _frame: Any) -> None:
        if active is not None and active.poll() is None:
            try:
                active.send_signal(signal.SIGUSR1)
            except (ProcessLookupError, OSError):
                pass

    if hasattr(signal, "SIGUSR1"):
        signal.signal(signal.SIGUSR1, forward_progress)
    try:
        print(f"I: Starting density batch {batch['batch_id']}", file=sys.stderr)
        active = subprocess.Popen(command, text=True)
        return_code = active.wait()
        if return_code != 0:
            raise DensityCalibrationError(
                f"scanner failed for {batch['batch_id']} with status {return_code}"
            )

        observed_paths = [Path(entry.path) for entry in os.scandir(output)
                          if entry.is_file()]
        if len(observed_paths) != batch["motif_count"]:
            raise DensityCalibrationError(
                f"{batch['batch_id']} produced {len(observed_paths)} files; "
                f"expected {batch['motif_count']}"
            )
        observed: dict[str, Path] = {}
        for path in observed_paths:
            motif_id = distribution_motif_id(path)
            if motif_id in observed:
                raise DensityCalibrationError(f"duplicate distribution for {motif_id}")
            observed[motif_id] = path
        if set(observed) != set(batch["motif_ids"]):
            raise DensityCalibrationError(
                f"distribution motif set differs for {batch['batch_id']}"
            )

        durable = attempt / "distributions"
        durable.mkdir()
        files: list[dict[str, Any]] = []
        for motif_id in batch["motif_ids"]:
            source = observed[motif_id]
            destination = durable / source.name
            with source.open("rb") as input_stream, destination.open("xb") as output_stream:
                shutil.copyfileobj(input_stream, output_stream, length=8 * 1024 * 1024)
                output_stream.flush()
                os.fsync(output_stream.fileno())
            digest = sha256_file(destination)
            read_distribution(
                motif_id, {"path": str(destination), "sha256": digest},
                region["chrom"], plan["score_mode"], plan["pseudocount"],
            )
            files.append({
                "motif_id": motif_id,
                "relative_path": str(Path("distributions") / destination.name),
                "bytes": destination.stat().st_size,
                "sha256": digest,
            })
        json_write_new(attempt / "complete.json", {
            "schema_version": 1,
            "batch_id": batch["batch_id"],
            "batch_index": batch["batch_index"],
            "motif_ids": batch["motif_ids"],
            "plan_sha256": plan_sha256,
            "scanner_sha256": sha256_file(scanner),
            "scanner_build": scanner_build,
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", "local"),
            "slurm_restart_count": int(os.environ.get("SLURM_RESTART_COUNT", "0")),
            "completed_at_epoch": int(time.time()),
            "files": files,
        })
        final = batch_directory(run_root, batch)
        if final.exists():
            raise DensityCalibrationError(f"batch appeared concurrently: {final}")
        os.replace(attempt, final)
        print(f"I: Promoted density batch {batch['batch_id']}", file=sys.stderr)
    except Exception as error:
        if not (attempt / "failure.json").exists():
            json_write_new(attempt / "failure.json", {
                "batch_id": batch["batch_id"],
                "error": str(error),
                "command": command,
            })
        raise


def run(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan, _, plan_sha256 = load_plan(run_root)
    if not arguments.allow_local and not os.environ.get("SLURM_JOB_ID"):
        raise DensityCalibrationError(
            "run requires Slurm; use --allow-local only for testing"
        )
    validate_plan_inputs(plan)
    check_free_space(run_root, plan["minimum_free_gib"])
    scanner = absolute_file(arguments.scanner)
    validate_scanner(scanner, plan, arguments.allow_scanner_provenance_mismatch)

    scratch = arguments.scratch_directory.expanduser().resolve()
    scratch.mkdir(parents=True, exist_ok=False)
    region = plan["chromosome"]
    required_gib = arguments.minimum_scratch_free_gib + region["length"] / 1024**3
    check_free_space(scratch, required_gib, "Scratch filesystem")
    staged_fasta = scratch / f"{safe_label(region['chrom'])}.fa"
    stage_result = stage_fasta_region(
        absolute_file(plan["genome"]["path"]),
        absolute_file(plan["fasta_index"]["path"]),
        region["chrom"], staged_fasta,
        expected_length=region["length"],
        expected_sha256=region["sequence_sha256"],
    )
    staged_metadata = scratch / "staged_sequence.json"
    json_write_new(staged_metadata, stage_result)
    staged_jaspar = scratch / Path(plan["jaspar"]["path"]).name
    shutil.copy2(plan["jaspar"]["path"], staged_jaspar)
    if sha256_file(staged_jaspar) != plan["jaspar"]["sha256"]:
        raise DensityCalibrationError("JASPAR changed while staging to scratch")
    print(
        f"I: Staged chromosome {region['chrom']} once at {staged_fasta} and "
        f"JASPAR at {staged_jaspar}; all density batches use these local copies.",
        file=sys.stderr,
    )

    pending = [
        batch for batch in plan["batches"]
        if validate_complete_batch(
            run_root, plan, plan_sha256, batch, verify_hashes=True
        ) is None
    ]
    print(
        f"I: {len(pending)} of {len(plan['batches'])} density batches remain.",
        file=sys.stderr,
    )
    active: dict[subprocess.Popen[str], dict[str, Any]] = {}
    next_batch = 0
    failed: str | None = None

    def forward_progress(_signal_number: int, _frame: Any) -> None:
        for process in active:
            if process.poll() is None:
                try:
                    process.send_signal(signal.SIGUSR1)
                    return
                except (ProcessLookupError, OSError):
                    continue
        print(
            f"I: density progress phase=batch_dispatch complete="
            f"{next_batch - len(active)}/{len(pending)} active={len(active)}",
            file=sys.stderr,
        )

    if hasattr(signal, "SIGUSR1"):
        signal.signal(signal.SIGUSR1, forward_progress)
    while active or (next_batch < len(pending) and failed is None):
        while (failed is None and next_batch < len(pending)
               and len(active) < arguments.batch_workers):
            batch = pending[next_batch]
            command = [
                sys.executable, str(Path(__file__).resolve()), "run-batch",
                "--run-root", str(run_root),
                "--batch-index", str(batch["batch_index"]),
                "--scanner", str(scanner),
                "--staged-fasta", str(staged_fasta),
                "--staged-fasta-index", stage_result["staged_fasta_index"],
                "--staged-jaspar", str(staged_jaspar),
                "--staged-metadata", str(staged_metadata),
                "--task-work-directory", str(scratch / "task-work"),
                "--planned-inputs-verified",
            ]
            if arguments.allow_scanner_provenance_mismatch:
                command.append("--allow-scanner-provenance-mismatch")
            process = subprocess.Popen(command, text=True)
            active[process] = batch
            next_batch += 1
        if not active:
            break
        time.sleep(0.1)
        for process, batch in list(active.items()):
            return_code = process.poll()
            if return_code is None:
                continue
            del active[process]
            if return_code != 0 and failed is None:
                failed = batch["batch_id"]
    if failed is not None:
        raise DensityCalibrationError(
            f"density worker stopped after failed batch {failed}; "
            "completed batches remain reusable"
        )
    print("I: All density batches are complete.", file=sys.stderr)


def status(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan, _, plan_sha256 = load_plan(run_root)
    complete = 0
    invalid = 0
    for batch in plan["batches"]:
        try:
            if validate_complete_batch(
                run_root, plan, plan_sha256, batch, verify_hashes=False
            ) is not None:
                complete += 1
        except DensityCalibrationError:
            invalid += 1
    print("key\tvalue")
    print(f"batches_total\t{len(plan['batches'])}")
    print(f"batches_complete\t{complete}")
    print(f"batches_pending\t{len(plan['batches']) - complete - invalid}")
    print(f"batches_invalid\t{invalid}")
    print(f"motifs_total\t{plan['motif_count']}")
    print(f"finalized\t{str((run_root / 'final' / 'manifest.json').is_file()).lower()}")


def finalize(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan, plan_path, plan_sha256 = load_plan(run_root)
    validate_plan_inputs(plan)
    final = run_root / "final"
    manifest_path = final / "manifest.json"
    if manifest_path.is_file():
        manifest = json_read(manifest_path)
        threshold_tsv = final / manifest["threshold_tsv"]
        threshold_json = final / manifest["threshold_json"]
        if (sha256_file(threshold_tsv) != manifest["threshold_tsv_sha256"]
                or sha256_file(threshold_json) != manifest["threshold_json_sha256"]):
            raise DensityCalibrationError("finalized density registry checksum differs")
        print(f"I: Reusing finalized density registry: {final}", file=sys.stderr)
        return
    if final.exists():
        raise DensityCalibrationError(
            f"partial final directory exists without completion marker: {final}"
        )

    complete_rows: list[tuple[dict[str, Any], dict[str, Any]]] = []
    for batch in plan["batches"]:
        marker = validate_complete_batch(
            run_root, plan, plan_sha256, batch, verify_hashes=True
        )
        if marker is None:
            raise DensityCalibrationError(f"density batch is incomplete: {batch['batch_id']}")
        complete_rows.append((batch, marker))

    attempt = run_root / (
        f"final-attempt-job-{os.environ.get('SLURM_JOB_ID', 'local')}-"
        f"pid-{os.getpid()}"
    )
    attempt.mkdir(exist_ok=False)
    distribution_manifest = attempt / "distribution_manifest.tsv"
    lines = ["motif_id\tpath\tsha256"]
    for batch, marker in complete_rows:
        by_motif = {row["motif_id"]: row for row in marker["files"]}
        for motif_id in batch["motif_ids"]:
            row = by_motif[motif_id]
            relative = Path("..") / "tasks" / batch["batch_id"] / row["relative_path"]
            lines.append(f"{motif_id}\t{relative}\t{row['sha256']}")
    immutable_write(distribution_manifest, "\n".join(lines) + "\n")

    threshold_tsv = attempt / "motif_scan_thresholds.tsv"
    threshold_json = attempt / "motif_scan_thresholds.json"
    command = [
        sys.executable,
        str(Path(__file__).with_name("build_density_capped_thresholds.py")),
        "--informative-thresholds",
        plan["threshold_inputs"]["informative_thresholds"]["path"],
        "--distribution-manifest", str(distribution_manifest),
        "--output-tsv", str(threshold_tsv),
        "--output-json", str(threshold_json),
        "--threshold-set-id", plan["threshold_set_id"],
        "--genome-id", plan["genome_id"],
        "--motif-set-id", plan["motif_set_id"],
        "--genome-fasta-sha256", plan["genome"]["sha256"],
        "--density-chrom-sequence-sha256",
        plan["chromosome"]["sequence_sha256"],
        "--jaspar-sha256", plan["jaspar"]["sha256"],
        "--source-commit", plan["source_commit"],
        "--chrom", plan["chromosome"]["chrom"],
        "--default-minimum-score", str(plan["default_minimum_score"]),
        "--minimum-spacing-bp", str(plan["minimum_spacing_bp"]),
        "--score-mode", plan["score_mode"],
        "--pseudocount", str(plan["pseudocount"]),
    ]
    if plan["threshold_inputs"]["negative_sensitivity"] is not None:
        command.extend([
            "--negative-sensitivity",
            plan["threshold_inputs"]["negative_sensitivity"]["path"],
        ])
    if plan["threshold_inputs"]["override_thresholds"] is not None:
        command.extend([
            "--override-thresholds",
            plan["threshold_inputs"]["override_thresholds"]["path"],
        ])
    process = subprocess.run(command, text=True, check=False)
    if process.returncode != 0:
        raise DensityCalibrationError("density-threshold finalizer failed")
    manifest = {
        "schema_version": 1,
        "run_id": plan["run_id"],
        "threshold_set_id": plan["threshold_set_id"],
        "source_commit": plan["source_commit"],
        "plan": str(Path("..") / "plan" / plan_path.name),
        "plan_sha256": plan_sha256,
        "distribution_manifest": distribution_manifest.name,
        "distribution_manifest_sha256": sha256_file(distribution_manifest),
        "threshold_tsv": threshold_tsv.name,
        "threshold_tsv_sha256": sha256_file(threshold_tsv),
        "threshold_json": threshold_json.name,
        "threshold_json_sha256": sha256_file(threshold_json),
        "motif_count": plan["motif_count"],
    }
    json_write_new(attempt / "manifest.json", manifest)
    os.replace(attempt, final)
    print(f"I: Finalized density registry: {final}", file=sys.stderr)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare", help="write an immutable run plan")
    prepare_parser.add_argument("--run-root", required=True, type=Path)
    prepare_parser.add_argument("--run-id", required=True)
    prepare_parser.add_argument("--source-commit", required=True)
    prepare_parser.add_argument("--genome", required=True, type=Path)
    prepare_parser.add_argument("--fasta-index", required=True, type=Path)
    prepare_parser.add_argument("--jaspar", required=True, type=Path)
    prepare_parser.add_argument("--informative-thresholds", required=True, type=Path)
    prepare_parser.add_argument("--negative-sensitivity", type=Path)
    prepare_parser.add_argument("--override-thresholds", type=Path)
    prepare_parser.add_argument("--threshold-set-id", required=True)
    prepare_parser.add_argument("--genome-id", required=True)
    prepare_parser.add_argument("--motif-set-id", required=True)
    prepare_parser.add_argument("--chrom", default="1")
    prepare_parser.add_argument("--motif-batch-size", type=int, default=64)
    prepare_parser.add_argument("--default-minimum-score", type=float, default=-1)
    prepare_parser.add_argument("--minimum-spacing-bp", type=float, default=200)
    prepare_parser.add_argument(
        "--score-mode", choices=("log2_relative_risk", "log_odds"),
        default="log2_relative_risk",
    )
    prepare_parser.add_argument("--pseudocount", type=float, default=1)
    prepare_parser.add_argument("--minimum-free-gib", type=float, default=500)
    prepare_parser.set_defaults(function=prepare)

    run_parser = subparsers.add_parser(
        "run", help="stage one chromosome once and execute remaining batches"
    )
    run_parser.add_argument("--run-root", required=True, type=Path)
    run_parser.add_argument("--scanner", required=True, type=Path)
    run_parser.add_argument("--scratch-directory", required=True, type=Path)
    run_parser.add_argument("--batch-workers", type=int, default=1)
    run_parser.add_argument("--minimum-scratch-free-gib", type=float, default=5)
    run_parser.add_argument("--allow-local", action="store_true")
    run_parser.add_argument("--allow-scanner-provenance-mismatch", action="store_true")
    run_parser.set_defaults(function=run)

    batch_parser = subparsers.add_parser("run-batch", help=argparse.SUPPRESS)
    batch_parser.add_argument("--run-root", required=True, type=Path)
    batch_parser.add_argument("--batch-index", required=True, type=int)
    batch_parser.add_argument("--scanner", required=True, type=Path)
    batch_parser.add_argument("--staged-fasta", required=True, type=Path)
    batch_parser.add_argument("--staged-fasta-index", required=True, type=Path)
    batch_parser.add_argument("--staged-jaspar", required=True, type=Path)
    batch_parser.add_argument("--staged-metadata", required=True, type=Path)
    batch_parser.add_argument("--task-work-directory", required=True, type=Path)
    batch_parser.add_argument("--planned-inputs-verified", action="store_true")
    batch_parser.add_argument("--allow-scanner-provenance-mismatch", action="store_true")
    batch_parser.set_defaults(function=run_batch)

    status_parser = subparsers.add_parser("status", help="report checkpoint progress")
    status_parser.add_argument("--run-root", required=True, type=Path)
    status_parser.set_defaults(function=status)

    finalize_parser = subparsers.add_parser(
        "finalize", help="build the exact manifest and density-capped registry"
    )
    finalize_parser.add_argument("--run-root", required=True, type=Path)
    finalize_parser.set_defaults(function=finalize)
    return parser


def main() -> int:
    try:
        arguments = build_parser().parse_args()
        if hasattr(arguments, "batch_workers") and arguments.batch_workers < 1:
            raise DensityCalibrationError("--batch-workers must be positive")
        arguments.function(arguments)
    except (
        OSError, ValueError, DensityCalibrationError, ScanError, StageError,
        ThresholdError,
    ) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

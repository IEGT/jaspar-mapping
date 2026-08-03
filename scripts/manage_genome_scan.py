#!/usr/bin/env python3

"""Plan, execute, inspect, and finalize sparse JASPAR genome scans.

Genomic hits are written directly by ``pssm_scan_parquet``. This coordinator
stores only small JSON control records before finalization; it never creates a
BED/TSV hit intermediate and never removes an existing path. Each Slurm task is
validated as a complete motif batch and atomically promoted from task-local
staging. Failed staging directories are deliberately retained for diagnosis.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
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
from typing import Any, BinaryIO, Iterable

from stage_fasta_region import StageError, stage_fasta_region


class ScanError(RuntimeError):
    pass


IDENTIFIER = re.compile(r"^[A-Za-z0-9._-]+$")
COMMIT = re.compile(r"^[0-9a-f]{40}$")
PLAN_SCHEMA_VERSION = 2
SUPPORTED_PLAN_SCHEMA_VERSIONS = {1, 2}


class OperationProgress:
    """One-line, signal-triggered progress for coordinator operations."""

    def __init__(self, phase: str = "starting") -> None:
        self.started = time.monotonic()
        self.phase = phase
        self.files_complete = 0
        self.files_total = 0
        self.bytes_complete = 0
        self.bytes_total = 0
        self.current_file_bytes = 0
        self.current_path = "none"
        self.current_task = "none"

    def install(self) -> None:
        if hasattr(signal, "SIGUSR1"):
            signal.signal(signal.SIGUSR1, self._handle_signal)

    def _handle_signal(self, _signal_number: int, _frame: Any) -> None:
        self.emit()

    def set_phase(self, phase: str, *, current_path: str | Path | None = None,
                  current_task: str | None = None) -> None:
        self.phase = phase
        if current_path is not None:
            self.current_path = str(current_path)
        if current_task is not None:
            self.current_task = current_task

    def emit(self) -> None:
        elapsed = max(0.0, time.monotonic() - self.started)
        observed_bytes = self.bytes_complete + self.current_file_bytes
        throughput = observed_bytes / elapsed if elapsed > 0 else 0.0
        remaining = max(0, self.bytes_total - observed_bytes)
        eta = remaining / throughput if throughput > 0 else math.nan
        eta_text = f"{eta:.1f}" if math.isfinite(eta) else "unknown"
        print(
            "I: manager progress request"
            f" phase={self.phase}"
            f" files_complete={self.files_complete}"
            f" files_total={self.files_total}"
            f" bytes_complete={observed_bytes}"
            f" bytes_total={self.bytes_total}"
            f" throughput_mib_s={throughput / 1024**2:.3f}"
            f" elapsed_seconds={elapsed:.3f}"
            f" eta_seconds={eta_text}"
            f" task={json.dumps(self.current_task)}"
            f" path={json.dumps(self.current_path)}",
            file=sys.stderr,
            flush=True,
        )


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def positive_integer(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("expected an integer greater than zero")
    return result


def nonnegative_float(value: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result < 0:
        raise argparse.ArgumentTypeError("expected a finite non-negative number")
    return result


def positive_float(value: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise argparse.ArgumentTypeError("expected a finite number greater than zero")
    return result


def finite_float(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise argparse.ArgumentTypeError("expected a finite number")
    return result


def duckdb_memory_limit(value: str) -> str:
    normalized = value.upper()
    if not re.fullmatch(r"[1-9][0-9]*(?:MB|GB)", normalized):
        raise argparse.ArgumentTypeError("expected a value such as 512MB or 12GB")
    return normalized


def identifier(value: str) -> str:
    if not IDENTIFIER.fullmatch(value):
        raise argparse.ArgumentTypeError(
            "use only letters, digits, '.', '_', or '-'"
        )
    return value


def absolute_file(value: str | Path) -> Path:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise ScanError(f"file not found: {path}")
    return path


def sha256_file(path: Path, progress: OperationProgress | None = None,
                max_read_mib_per_second: float | None = None) -> str:
    digest = hashlib.sha256()
    started = time.monotonic()
    bytes_read = 0
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
            bytes_read += len(chunk)
            if progress is not None:
                progress.current_file_bytes = bytes_read
            if max_read_mib_per_second is not None:
                expected_elapsed = bytes_read / (max_read_mib_per_second * 1024**2)
                delay = expected_elapsed - (time.monotonic() - started)
                if delay > 0:
                    time.sleep(delay)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def file_stat_identity(path: Path) -> dict[str, int]:
    metadata = path.stat()
    return {"file_bytes": metadata.st_size, "file_mtime_ns": metadata.st_mtime_ns}


def json_write_new(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def json_write_atomic(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def json_read(path: Path) -> Any:
    try:
        with path.open(encoding="utf-8") as stream:
            return json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise ScanError(f"cannot read JSON {path}: {error}") from error


def check_free_space(path: Path, minimum_gib: float,
                     label: str = "Output filesystem") -> None:
    free = shutil.disk_usage(path).free
    required = int(minimum_gib * 1024**3)
    print(
        f"I: {label} has {free / 1024**3:.1f} GiB free; "
        f"required preflight reserve is {minimum_gib:.1f} GiB.",
        file=sys.stderr,
    )
    if free < required:
        raise ScanError(
            f"insufficient free space: {free / 1024**3:.1f} GiB available, "
            f"{minimum_gib:.1f} GiB required"
        )


def scanner_build_info(scanner: Path) -> dict[str, Any]:
    process = subprocess.run(
        [str(scanner), "--version-json"],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ScanError(
            process.stderr.strip() or "scanner --version-json failed"
        )
    try:
        result = json.loads(process.stdout)
    except json.JSONDecodeError as error:
        raise ScanError("scanner --version-json returned invalid JSON") from error
    required = {
        "program", "version", "source_commit", "source_dirty", "compiler",
        "cplusplus", "build_flags",
        "lto_enabled", "ndebug", "parquet_enabled", "arrow_version",
        "parquet_version",
    }
    if not isinstance(result, dict) or not required.issubset(result):
        raise ScanError("scanner build information is incomplete")
    if result.get("program") != "pssm_scan":
        raise ScanError("scanner build information identifies a different program")
    if result.get("parquet_enabled") is not True:
        raise ScanError(
            "scanner lacks direct Parquet support; build pssm_scan_parquet"
        )
    return result


def duckdb_build_info(executable: str) -> str:
    process = subprocess.run(
        [executable, "--version"], text=True, capture_output=True, check=False
    )
    if process.returncode != 0:
        raise ScanError(process.stderr.strip() or "duckdb --version failed")
    return process.stdout.strip()


def parse_fasta_index(path: Path) -> list[dict[str, Any]]:
    regions: list[dict[str, Any]] = []
    seen: set[str] = set()
    with path.open(encoding="ascii") as stream:
        for order, line in enumerate(stream):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ScanError(f"invalid FASTA index line {order + 1}: {path}")
            name = fields[0]
            if name in seen:
                raise ScanError(f"duplicate sequence name in FASTA index: {name}")
            seen.add(name)
            try:
                length = int(fields[1])
            except ValueError as error:
                raise ScanError(f"invalid sequence length for {name}") from error
            regions.append(
                {"sequence_order": order, "chrom": name, "length": length}
            )
    if not regions:
        raise ScanError(f"empty FASTA index: {path}")
    return regions


def open_fasta_binary(path: Path) -> tuple[BinaryIO, bool]:
    if path.suffix == ".gz":
        return gzip.open(path, "rb"), True
    return path.open("rb"), False


def inspect_fasta(path: Path) -> tuple[str, list[dict[str, Any]], str]:
    raw_digest = hashlib.sha256()
    sequences: list[dict[str, Any]] = []
    name: str | None = None
    sequence_digest = hashlib.sha256()
    sequence_length = 0
    stream, compressed = open_fasta_binary(path)

    def finish_sequence() -> None:
        nonlocal name, sequence_digest, sequence_length
        if name is None:
            return
        sequences.append(
            {
                "chrom": name,
                "length": sequence_length,
                "sequence_sha256": sequence_digest.hexdigest(),
            }
        )

    with stream:
        for raw_line in stream:
            if not compressed:
                raw_digest.update(raw_line)
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(b">"):
                finish_sequence()
                header = line[1:].split(None, 1)[0]
                if not header:
                    raise ScanError(f"empty FASTA sequence name in {path}")
                name = header.decode("ascii")
                sequence_digest = hashlib.sha256()
                sequence_length = 0
                continue
            if name is None:
                raise ScanError(f"FASTA sequence precedes first header in {path}")
            bases = line.upper()
            sequence_digest.update(bases)
            sequence_length += len(bases)
    finish_sequence()
    if compressed:
        raw_sha256 = sha256_file(path)
    else:
        raw_sha256 = raw_digest.hexdigest()
    if not sequences:
        raise ScanError(f"no sequences in FASTA: {path}")
    set_digest = hashlib.sha256()
    for sequence in sequences:
        set_digest.update(
            (
                f"{sequence['chrom']}\t{sequence['length']}\t"
                f"{sequence['sequence_sha256']}\n"
            ).encode("ascii")
        )
    return raw_sha256, sequences, set_digest.hexdigest()


def parse_jaspar(path: Path) -> list[dict[str, Any]]:
    motifs: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    row_lengths: list[int] = []

    def finish() -> None:
        nonlocal current, row_lengths
        if current is None:
            return
        if len(row_lengths) != 4 or len(set(row_lengths)) != 1:
            raise ScanError(f"motif {current['motif_id']} has an invalid matrix")
        current["motif_length"] = row_lengths[0]
        motifs.append(current)

    with path.open(encoding="utf-8") as stream:
        for line_number, raw_line in enumerate(stream, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                fields = line[1:].split(None, 1)
                current = {
                    "motif_id": fields[0],
                    "motif_name": fields[1].strip() if len(fields) > 1 else fields[0],
                }
                row_lengths = []
                continue
            if current is None:
                raise ScanError(f"matrix row before motif header at line {line_number}")
            fields = line.replace("[", " ").replace("]", " ").split()
            if not fields or fields[0].upper() not in {"A", "C", "G", "T"}:
                raise ScanError(f"invalid JASPAR row at line {line_number}")
            try:
                [float(value) for value in fields[1:]]
            except ValueError as error:
                raise ScanError(f"invalid JASPAR count at line {line_number}") from error
            row_lengths.append(len(fields) - 1)
    finish()
    ids = [motif["motif_id"] for motif in motifs]
    if not motifs or len(ids) != len(set(ids)):
        raise ScanError("JASPAR file is empty or contains duplicate motif IDs")
    return motifs


def parse_chromosomes(values: Iterable[str]) -> list[str]:
    result: list[str] = []
    for value in values:
        for item in value.split(","):
            chrom = item.strip()
            if not chrom:
                raise ScanError("--chrom contains an empty value")
            if chrom not in result:
                result.append(chrom)
    return result


def safe_label(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "-", value).strip("-") or "sequence"


def optional_float(value: str, label: str) -> float | None:
    if value in {"", "NA"}:
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ScanError(f"{label} must be finite or NA")
    return result


def parse_boolean(value: str, label: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise ScanError(f"{label} must be true or false")


def read_motif_threshold_registry(
    threshold_path: Path,
    metadata_path: Path,
    motif_ids: set[str],
    arguments: argparse.Namespace,
    fasta_sha256: str,
    sequence_sha256_by_chrom: dict[str, str],
    jaspar_sha256: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    metadata = json_read(metadata_path)
    if metadata.get("threshold_tsv_sha256") != sha256_file(threshold_path):
        raise ScanError("motif-threshold TSV checksum disagrees with its metadata")
    expected_metadata = {
        "schema_version": 1,
        "genome_id": arguments.genome_id,
        "motif_set_id": arguments.motif_set_id,
        "score_mode": arguments.score_mode,
        "genome_fasta_sha256": fasta_sha256,
        "jaspar_sha256": jaspar_sha256,
        "source_commit": arguments.source_commit,
        "orientation_aggregation": "max_score_per_alignment_span",
        "distribution_bin_width": 1,
        "candidate_formula": "min(informative_threshold, default_minimum_score)",
        "final_formula": "max(candidate_minimum_score, density_threshold)",
    }
    for field, expected in expected_metadata.items():
        if metadata.get(field) != expected:
            raise ScanError(f"motif-threshold metadata {field} mismatch")
    if abs(float(metadata.get("pseudocount", math.nan)) - arguments.pseudocount) > 1e-9:
        raise ScanError("motif-threshold metadata pseudocount mismatch")
    density_chrom = str(metadata.get("density_chrom", ""))
    if density_chrom not in sequence_sha256_by_chrom:
        raise ScanError("motif-threshold density chromosome is absent from the genome")
    if (metadata.get("density_chrom_sequence_sha256")
            != sequence_sha256_by_chrom[density_chrom]):
        raise ScanError("motif-threshold density chromosome checksum mismatch")
    density_spacing = float(metadata.get("density_minimum_spacing_bp", math.nan))
    default_score = float(metadata.get("default_minimum_score", math.nan))
    if not math.isfinite(density_spacing) or density_spacing <= 0:
        raise ScanError("motif-threshold density spacing is invalid")
    if not math.isfinite(default_score):
        raise ScanError("motif-threshold default score is invalid")
    threshold_set_id = metadata.get("threshold_set_id")
    if not isinstance(threshold_set_id, str) or not IDENTIFIER.fullmatch(threshold_set_id):
        raise ScanError("motif-threshold set ID is invalid")
    for field in (
        "informative_thresholds_sha256", "distribution_manifest_sha256",
    ):
        digest = metadata.get(field)
        if (not isinstance(digest, str) or len(digest) != 64
                or any(character not in "0123456789abcdef" for character in digest)):
            raise ScanError(f"motif-threshold metadata {field} is invalid")
    for field in ("negative_sensitivity_sha256", "override_thresholds_sha256"):
        digest = metadata.get(field)
        if digest is not None and (
            not isinstance(digest, str) or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
        ):
            raise ScanError(f"motif-threshold metadata {field} is invalid")

    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    with threshold_path.open(encoding="ascii", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {
            "motif_id", "informative_threshold", "informative_source",
            "default_minimum_score", "candidate_minimum_score",
            "density_minimum_spacing_bp", "density_maximum_loci",
            "density_threshold", "final_minimum_score", "density_limited",
            "density_chrom", "valid_locus_starts", "skipped_locus_starts",
            "loci_at_candidate_threshold", "loci_at_final_threshold",
            "mean_spacing_bp_at_final_threshold", "distribution_sha256",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ScanError("motif-threshold registry has an incomplete schema")
        for source in reader:
            motif_id = source["motif_id"]
            if not motif_id or motif_id in seen:
                raise ScanError(f"blank or duplicate motif threshold: {motif_id!r}")
            seen.add(motif_id)
            final_score = optional_float(
                source["final_minimum_score"], f"final threshold for {motif_id}"
            )
            if final_score is None:
                raise ScanError(f"final threshold is absent for {motif_id}")
            informative_threshold = optional_float(
                source["informative_threshold"],
                f"informative threshold for {motif_id}",
            )
            row = {
                    "motif_id": motif_id,
                    "informative_threshold": informative_threshold,
                    "informative_source": source["informative_source"],
                    "default_minimum_score": float(source["default_minimum_score"]),
                    "candidate_minimum_score": float(source["candidate_minimum_score"]),
                    "density_minimum_spacing_bp": float(
                        source["density_minimum_spacing_bp"]
                    ),
                    "density_maximum_loci": int(source["density_maximum_loci"]),
                    "density_threshold": float(source["density_threshold"]),
                    "final_minimum_score": final_score,
                    "density_limited": parse_boolean(
                        source["density_limited"], f"density flag for {motif_id}"
                    ),
                    "density_chrom": source["density_chrom"],
                    "valid_locus_starts": int(source["valid_locus_starts"]),
                    "skipped_locus_starts": int(source["skipped_locus_starts"]),
                    "loci_at_candidate_threshold": int(
                        source["loci_at_candidate_threshold"]
                    ),
                    "loci_at_final_threshold": int(
                        source["loci_at_final_threshold"]
                    ),
                    "mean_spacing_bp_at_final_threshold": optional_float(
                        source["mean_spacing_bp_at_final_threshold"],
                        f"mean spacing for {motif_id}",
                    ),
                    "distribution_sha256": source["distribution_sha256"],
                }
            numeric_thresholds = (
                row["default_minimum_score"], row["candidate_minimum_score"],
                row["density_threshold"], row["final_minimum_score"],
            )
            if any(not math.isfinite(value) or abs(value - round(value)) > 1e-9
                   for value in numeric_thresholds):
                raise ScanError(f"thresholds are not on the integer grid for {motif_id}")
            if (informative_threshold is not None
                    and abs(informative_threshold - round(informative_threshold)) > 1e-9):
                raise ScanError(
                    f"informative threshold is not on the integer grid for {motif_id}"
                )
            expected_candidate = default_score if informative_threshold is None \
                else min(default_score, informative_threshold)
            if (row["informative_source"] == ""
                    or row["default_minimum_score"] != default_score
                    or row["candidate_minimum_score"] != expected_candidate
                    or row["density_minimum_spacing_bp"] != density_spacing
                    or row["density_chrom"] != density_chrom):
                raise ScanError(f"threshold policy fields are inconsistent for {motif_id}")
            expected_maximum = math.floor(
                row["valid_locus_starts"] / density_spacing
            )
            expected_final = max(
                row["candidate_minimum_score"], row["density_threshold"]
            )
            if (row["density_maximum_loci"] != expected_maximum
                    or row["final_minimum_score"] != expected_final
                    or row["density_limited"]
                        != (expected_final > row["candidate_minimum_score"])
                    or min(
                        row["valid_locus_starts"], row["skipped_locus_starts"],
                        row["loci_at_candidate_threshold"],
                        row["loci_at_final_threshold"],
                    ) < 0
                    or row["loci_at_final_threshold"] > expected_maximum
                    or row["loci_at_final_threshold"]
                        > row["loci_at_candidate_threshold"]):
                raise ScanError(f"density decision is inconsistent for {motif_id}")
            expected_spacing = (
                None if row["loci_at_final_threshold"] == 0
                else row["valid_locus_starts"] / row["loci_at_final_threshold"]
            )
            observed_spacing = row["mean_spacing_bp_at_final_threshold"]
            if ((expected_spacing is None) != (observed_spacing is None)
                    or (expected_spacing is not None
                        and not math.isclose(
                            expected_spacing, observed_spacing,
                            rel_tol=1e-9, abs_tol=1e-9,
                        ))):
                raise ScanError(f"mean density spacing is inconsistent for {motif_id}")
            digest = row["distribution_sha256"]
            if (len(digest) != 64
                    or any(character not in "0123456789abcdef" for character in digest)):
                raise ScanError(f"distribution checksum is invalid for {motif_id}")
            rows.append(row)
    if seen != motif_ids:
        raise ScanError(
            "motif-threshold registry does not exactly cover the JASPAR motif set; "
            f"missing={sorted(motif_ids - seen)[:10]}, extra={sorted(seen - motif_ids)[:10]}"
        )
    if metadata.get("motifs") != len(rows):
        raise ScanError("motif-threshold metadata row count mismatch")
    return rows, metadata


def load_plan(run_root: Path) -> tuple[dict[str, Any], Path, str]:
    path = run_root / "plan" / "scan_plan.json"
    plan = json_read(path)
    if plan.get("schema_version") not in SUPPORTED_PLAN_SCHEMA_VERSIONS:
        raise ScanError(f"unsupported scan-plan schema in {path}")
    return plan, path, sha256_file(path)


def prepare(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    run_root.mkdir(parents=True, exist_ok=True)
    check_free_space(run_root, arguments.minimum_free_gib)
    plan_path = run_root / "plan" / "scan_plan.json"
    if plan_path.exists():
        raise ScanError(f"scan plan already exists: {plan_path}")
    if not COMMIT.fullmatch(arguments.source_commit):
        raise ScanError("--source-commit must be a full lowercase Git SHA")

    genome = absolute_file(arguments.genome)
    fasta_index = absolute_file(arguments.fasta_index)
    jaspar = absolute_file(arguments.jaspar)
    genome_stat = file_stat_identity(genome)
    index_regions = parse_fasta_index(fasta_index)
    fasta_sha256, fasta_regions, sequence_set_sha256 = inspect_fasta(genome)
    fasta_by_name = {row["chrom"]: row for row in fasta_regions}
    if [row["chrom"] for row in index_regions] != [row["chrom"] for row in fasta_regions]:
        raise ScanError("FASTA and .fai sequence order differ")
    for row in index_regions:
        observed = fasta_by_name[row["chrom"]]
        if row["length"] != observed["length"]:
            raise ScanError(f"FASTA and .fai lengths differ for {row['chrom']}")
        row["sequence_sha256"] = observed["sequence_sha256"]

    requested = parse_chromosomes(arguments.chrom)
    indexed_names = {row["chrom"] for row in index_regions}
    missing = [chrom for chrom in requested if chrom not in indexed_names]
    if missing:
        raise ScanError("requested sequence region(s) absent from .fai: " + ", ".join(missing))
    selected_names = set(requested or indexed_names)
    for row in index_regions:
        row["included_in_scan"] = row["chrom"] in selected_names

    motifs = parse_jaspar(jaspar)
    jaspar_sha256 = sha256_file(jaspar)
    motif_by_id = {motif["motif_id"]: motif for motif in motifs}
    if arguments.special_motif not in motif_by_id:
        raise ScanError(f"special motif absent from JASPAR: {arguments.special_motif}")
    for motif in motifs:
        motif["motif_set_id"] = arguments.motif_set_id
        motif["jaspar_version"] = arguments.jaspar_version

    plan_dir = run_root / "plan"
    motif_list_dir = plan_dir / "motif_lists"
    motif_list_dir.mkdir(parents=True, exist_ok=False)
    batches: list[dict[str, Any]] = []

    def add_batch(batch_id: str, policy_id: str, motif_ids: list[str],
                  minimum_score: float) -> None:
        content = "".join(f"{motif_id}\n" for motif_id in motif_ids)
        relative = Path("plan") / "motif_lists" / f"{batch_id}.txt"
        path = run_root / relative
        with path.open("x", encoding="ascii") as stream:
            stream.write(content)
        batches.append(
            {
                "motif_batch_id": batch_id,
                "policy_id": policy_id,
                "motif_ids": motif_ids,
                "motif_count": len(motif_ids),
                "minimum_score": minimum_score,
                "motif_list": str(relative),
                "motif_list_sha256": sha256_text(content),
            }
        )

    threshold_registry: dict[str, Any] | None = None
    threshold_metadata: dict[str, Any] | None = None
    motif_thresholds: list[dict[str, Any]] = []
    if bool(arguments.motif_thresholds) != bool(arguments.motif_threshold_metadata):
        raise ScanError(
            "--motif-thresholds and --motif-threshold-metadata must be supplied together"
        )
    if arguments.motif_thresholds:
        threshold_path = absolute_file(arguments.motif_thresholds)
        threshold_metadata_path = absolute_file(arguments.motif_threshold_metadata)
        motif_thresholds, threshold_metadata = read_motif_threshold_registry(
            threshold_path, threshold_metadata_path, set(motif_by_id), arguments,
            fasta_sha256,
            {row["chrom"]: row["sequence_sha256"] for row in fasta_regions},
            jaspar_sha256,
        )
        threshold_registry = {
            "source_file": str(threshold_path),
            "source_sha256": sha256_file(threshold_path),
            "source_bytes": threshold_path.stat().st_size,
            "source_mtime_ns": threshold_path.stat().st_mtime_ns,
            "metadata_file": str(threshold_metadata_path),
            "metadata_sha256": sha256_file(threshold_metadata_path),
            "metadata_bytes": threshold_metadata_path.stat().st_size,
            "metadata_mtime_ns": threshold_metadata_path.stat().st_mtime_ns,
            "threshold_set_id": threshold_metadata["threshold_set_id"],
        }
        policy_id = "per_motif_density_capped"
        policies = [
            {
                "policy_id": policy_id,
                "selector_type": "per_motif_registry",
                "motif_id": None,
                "minimum_score": None,
                "precedence": 1,
                "threshold_set_id": threshold_metadata["threshold_set_id"],
                "threshold_registry_sha256": threshold_registry["source_sha256"],
                "threshold_metadata_sha256": threshold_registry["metadata_sha256"],
                "candidate_formula": threshold_metadata["candidate_formula"],
                "final_formula": threshold_metadata["final_formula"],
                "density_chrom": threshold_metadata["density_chrom"],
                "density_minimum_spacing_bp": threshold_metadata[
                    "density_minimum_spacing_bp"
                ],
                "orientation_aggregation": threshold_metadata[
                    "orientation_aggregation"
                ],
                "distribution_bin_width": threshold_metadata[
                    "distribution_bin_width"
                ],
                "informative_thresholds_sha256": threshold_metadata[
                    "informative_thresholds_sha256"
                ],
                "negative_sensitivity_sha256": threshold_metadata.get(
                    "negative_sensitivity_sha256"
                ),
                "override_thresholds_sha256": threshold_metadata.get(
                    "override_thresholds_sha256"
                ),
                "distribution_manifest_sha256": threshold_metadata[
                    "distribution_manifest_sha256"
                ],
                "rationale": (
                    "Per-motif informative-or-minus-one floor raised only as "
                    "needed to satisfy the chromosome-1 physical-locus density ceiling"
                ),
            }
        ]
        threshold_groups: dict[float, list[str]] = {}
        for row in motif_thresholds:
            threshold_groups.setdefault(row["final_minimum_score"], []).append(
                row["motif_id"]
            )
        for threshold_value in sorted(threshold_groups):
            motif_ids = sorted(threshold_groups[threshold_value])
            label = f"{threshold_value:.12g}".replace("-", "m").replace(".", "p")
            for offset in range(0, len(motif_ids), arguments.motif_batch_size):
                add_batch(
                    f"score-{label}-{offset // arguments.motif_batch_size:04d}",
                    policy_id,
                    motif_ids[offset : offset + arguments.motif_batch_size],
                    threshold_value,
                )
    else:
        add_batch(
            "special-0000", "tp73_calibrated", [arguments.special_motif],
            arguments.special_minimum_score,
        )
        other_ids = sorted(set(motif_by_id) - {arguments.special_motif})
        for offset in range(0, len(other_ids), arguments.motif_batch_size):
            add_batch(
                f"default-{offset // arguments.motif_batch_size:04d}",
                "default_uncalibrated",
                other_ids[offset : offset + arguments.motif_batch_size],
                arguments.default_minimum_score,
            )
        policies = [
            {
                "policy_id": "tp73_calibrated",
                "selector_type": "motif_id",
                "motif_id": arguments.special_motif,
                "minimum_score": arguments.special_minimum_score,
                "precedence": 1,
                "rationale": "TP73 chr1 CUT&RUN calibration; retain deeper anchors for derived floor sensitivity",
            },
            {
                "policy_id": "default_uncalibrated",
                "selector_type": "all_except",
                "motif_id": arguments.special_motif,
                "minimum_score": arguments.default_minimum_score,
                "precedence": 2,
                "rationale": "Conservative provisional floor for motifs without matched CUT&RUN calibration",
            },
        ]
        for motif_id in sorted(motif_by_id):
            minimum_score = (
                arguments.special_minimum_score
                if motif_id == arguments.special_motif
                else arguments.default_minimum_score
            )
            motif_thresholds.append(
                {
                    "motif_id": motif_id,
                    "informative_threshold": None,
                    "informative_source": "legacy_scan_policy",
                    "default_minimum_score": arguments.default_minimum_score,
                    "candidate_minimum_score": minimum_score,
                    "density_minimum_spacing_bp": None,
                    "density_maximum_loci": None,
                    "density_threshold": None,
                    "final_minimum_score": minimum_score,
                    "density_limited": None,
                    "density_chrom": None,
                    "valid_locus_starts": None,
                    "skipped_locus_starts": None,
                    "loci_at_candidate_threshold": None,
                    "loci_at_final_threshold": None,
                    "mean_spacing_bp_at_final_threshold": None,
                    "distribution_sha256": None,
                }
            )
    tasks: list[dict[str, Any]] = []
    selected_regions = [row for row in index_regions if row["included_in_scan"]]
    # Batch-major ordering distributes the first Slurm array wave across
    # chromosomes instead of making every concurrent task read chromosome 1.
    for batch in batches:
        for region in selected_regions:
            task_index = len(tasks)
            task_id = (
                f"seq-{region['sequence_order']:04d}-{safe_label(region['chrom'])}-"
                f"{batch['motif_batch_id']}"
            )
            tasks.append(
                {
                    "task_index": task_index,
                    "task_id": task_id,
                    "chrom": region["chrom"],
                    "sequence_length": region["length"],
                    "motif_batch_id": batch["motif_batch_id"],
                    "motif_list": batch["motif_list"],
                    "motif_list_sha256": batch["motif_list_sha256"],
                    "motif_ids": batch["motif_ids"],
                    "motif_count": batch["motif_count"],
                    "policy_id": batch["policy_id"],
                    "minimum_score": batch["minimum_score"],
                }
            )

    plan = {
        "schema_version": PLAN_SCHEMA_VERSION,
        "created_at_utc": utc_now(),
        "run_id": arguments.run_id,
        "source_commit": arguments.source_commit,
        "minimum_free_gib": arguments.minimum_free_gib,
        "genome": {
            "genome_id": arguments.genome_id,
            "taxon_id": arguments.taxon_id,
            "species_scientific_name": arguments.species,
            "assembly_name": arguments.assembly_name,
            "assembly_accession": arguments.assembly_accession,
            "ensembl_release": arguments.ensembl_release,
            "fasta_url": arguments.fasta_url,
            "fasta_file": str(genome),
            "fasta_sha256": fasta_sha256,
            "fasta_bytes": genome_stat["file_bytes"],
            "fasta_mtime_ns": genome_stat["file_mtime_ns"],
            "fasta_index_file": str(fasta_index),
            "fasta_index_sha256": sha256_file(fasta_index),
            "fasta_index_bytes": fasta_index.stat().st_size,
            "fasta_index_mtime_ns": fasta_index.stat().st_mtime_ns,
            "sequence_set_sha256": sequence_set_sha256,
        },
        "sequence_regions": index_regions,
        "motif_set": {
            "motif_set_id": arguments.motif_set_id,
            "motif_collection": "JASPAR CORE non-redundant",
            "jaspar_version": arguments.jaspar_version,
            "source_url": arguments.jaspar_url,
            "source_file": str(jaspar),
            "source_sha256": jaspar_sha256,
            "source_bytes": jaspar.stat().st_size,
            "source_mtime_ns": jaspar.stat().st_mtime_ns,
            "motif_count": len(motifs),
        },
        "motif_metadata": motifs,
        "scanner": {
            "score_mode": arguments.score_mode,
            "pseudocount": arguments.pseudocount,
            "background_model_id": "uniform_acgt_v1",
            "pseudocount_scheme": "additive_per_base",
            "minimum_pwm_relative_score": None,
            "maximum_pwm_relative_score": None,
            "coordinate_mode": "bed",
            "n_policy": "skip",
            "matched_sequence_policy": "omitted",
            "strand": "both",
        },
        "threshold_policies": policies,
        "threshold_registry": threshold_registry,
        "threshold_metadata": threshold_metadata,
        "motif_thresholds": motif_thresholds,
        "motif_batches": batches,
        "tasks": tasks,
    }
    json_write_new(plan_path, plan)
    print(f"I: Wrote immutable scan plan: {plan_path}", file=sys.stderr)
    print(
        f"I: Planned {len(tasks)} tasks across {len(selected_regions)} sequence "
        f"regions and {len(batches)} non-overlapping motif batches.",
        file=sys.stderr,
    )
    print(f"TASK_COUNT={len(tasks)}")


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_string_list(values: Iterable[str | Path]) -> str:
    return "[" + ",".join(sql_string(value) for value in values) + "]"


def run_duckdb_json(
    executable: str,
    query: str,
    memory_limit: str | None = None,
    temp_directory: Path | None = None,
) -> list[dict[str, Any]]:
    settings = []
    if memory_limit is not None:
        settings.append(f"SET memory_limit = {sql_string(memory_limit)};")
    if temp_directory is not None:
        settings.append(f"SET temp_directory = {sql_string(temp_directory)};")
    configured_query = "\n".join(settings + [query])
    process = subprocess.run(
        [executable, "-json", ":memory:", "-c", configured_query],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        detail = process.stderr.strip() or process.stdout.strip()
        raise ScanError(f"DuckDB validation failed: {detail}")
    try:
        value = json.loads(process.stdout or "[]")
    except json.JSONDecodeError as error:
        raise ScanError(f"DuckDB returned invalid JSON: {error}") from error
    if not isinstance(value, list):
        raise ScanError("DuckDB validation did not return a row list")
    return value


def task_for_arguments(plan: dict[str, Any], arguments: argparse.Namespace) -> dict[str, Any]:
    if arguments.task_id is not None:
        matches = [task for task in plan["tasks"] if task["task_id"] == arguments.task_id]
    else:
        matches = [
            task for task in plan["tasks"]
            if task["task_index"] == arguments.task_index
        ]
    if len(matches) != 1:
        requested = arguments.task_id if arguments.task_id is not None else arguments.task_index
        raise ScanError(f"scan plan has no unique task {requested!r}")
    return matches[0]


def validate_planned_inputs(plan: dict[str, Any]) -> None:
    genome = absolute_file(plan["genome"]["fasta_file"])
    genome_stat = file_stat_identity(genome)
    if genome_stat["file_bytes"] != plan["genome"].get("fasta_bytes") \
            or genome_stat["file_mtime_ns"] != plan["genome"].get("fasta_mtime_ns"):
        raise ScanError(
            "genome FASTA changed after planning (size or mtime differs); "
            "prepare a new immutable scan plan"
        )
    fasta_index = absolute_file(plan["genome"]["fasta_index_file"])
    if sha256_file(fasta_index) != plan["genome"]["fasta_index_sha256"]:
        raise ScanError("FASTA index checksum changed after planning")
    jaspar = absolute_file(plan["motif_set"]["source_file"])
    if sha256_file(jaspar) != plan["motif_set"]["source_sha256"]:
        raise ScanError("JASPAR source checksum changed after planning")
    threshold_registry = plan.get("threshold_registry")
    if threshold_registry is not None:
        threshold_file = absolute_file(threshold_registry["source_file"])
        threshold_metadata = absolute_file(threshold_registry["metadata_file"])
        if sha256_file(threshold_file) != threshold_registry["source_sha256"]:
            raise ScanError("motif-threshold registry checksum changed after planning")
        if sha256_file(threshold_metadata) != threshold_registry["metadata_sha256"]:
            raise ScanError("motif-threshold metadata checksum changed after planning")


def validate_task_outputs(
    duckdb: str,
    stats_path: Path,
    package: Path,
    task: dict[str, Any],
    plan: dict[str, Any],
    duckdb_memory: str,
    duckdb_temp: Path | None,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with stats_path.open(encoding="utf-8") as stream:
        for line_number, line in enumerate(stream, start=1):
            if not line.strip():
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as error:
                raise ScanError(f"invalid scanner statistics line {line_number}") from error
    expected_files = 2 * task["motif_count"]
    if len(records) != expected_files:
        raise ScanError(
            f"task produced {len(records)} statistics rows; expected {expected_files}"
        )
    expected_pairs = {
        (motif_id, strand)
        for motif_id in task["motif_ids"]
        for strand in ("+", "-")
    }
    observed_pairs = {(row["motif_id"], row["strand"]) for row in records}
    if observed_pairs != expected_pairs:
        raise ScanError("scanner statistics do not cover the planned motif/strand set")
    motif_lengths: dict[str, int] = {}
    for row in records:
        expected_text = {
            "motif_set_id": plan["motif_set"]["motif_set_id"],
            "genome_id": plan["genome"]["genome_id"],
            "chrom": task["chrom"],
            "score_mode": plan["scanner"]["score_mode"],
            "background_model_id": plan["scanner"]["background_model_id"],
            "pseudocount_scheme": plan["scanner"]["pseudocount_scheme"],
            "coordinate_mode": plan["scanner"]["coordinate_mode"],
            "n_policy": plan["scanner"]["n_policy"],
            "matched_sequence_policy": plan["scanner"]["matched_sequence_policy"],
        }
        for field, expected in expected_text.items():
            if row.get(field) != expected:
                raise ScanError(
                    f"scanner statistics {field} mismatch for {row.get('motif_id')}"
                )
        expected_numbers = {
            "pseudocount": plan["scanner"]["pseudocount"],
            "minimum_score": task["minimum_score"],
            "minimum_pwm_relative_score": plan["scanner"][
                "minimum_pwm_relative_score"
            ],
            "maximum_pwm_relative_score": plan["scanner"][
                "maximum_pwm_relative_score"
            ],
        }
        for field, expected in expected_numbers.items():
            observed = row.get(field)
            matches = observed is None if expected is None else (
                isinstance(observed, (int, float))
                and math.isclose(float(observed), float(expected), abs_tol=1e-12)
            )
            if not matches:
                raise ScanError(
                    f"scanner statistics {field} mismatch for {row.get('motif_id')}"
                )
        motif_lengths[row["motif_id"]] = row["motif_length"]
        expected_windows = max(
            0, task["sequence_length"] - row["motif_length"] + 1
        )
        if row["expected_windows"] != expected_windows:
            raise ScanError(f"expected-window mismatch for {row['motif_id']}")
        accounted = sum(
            row[field]
            for field in (
                "skipped_n_windows",
                "sentinel_score_windows",
                "below_minimum_score_windows",
                "pwm_filtered_windows",
                "emitted_hits",
            )
        )
        if accounted != expected_windows:
            raise ScanError(f"window accounting mismatch for {row['motif_id']}")
        output = Path(row["output_file"])
        if not output.is_file() or package not in output.parents:
            raise ScanError(f"scanner output lies outside task package: {output}")
        if row.get("bytes") != output.stat().st_size:
            raise ScanError(f"scanner byte count mismatch for {output}")
        if not isinstance(row.get("elapsed_seconds"), (int, float)) \
                or row["elapsed_seconds"] < 0:
            raise ScanError(f"invalid scanner elapsed time for {output}")

    output_paths = [Path(row["output_file"]).resolve() for row in records]
    query = f"""
WITH expected AS (
    SELECT output_file::VARCHAR AS filename, emitted_hits::BIGINT AS expected_rows,
           motif_length::BIGINT AS motif_length,
           minimum_score::DOUBLE AS minimum_score,
           genome_id::VARCHAR AS genome_id,
           motif_set_id::VARCHAR AS motif_set_id,
           motif_id::VARCHAR AS motif_id,
           chrom::VARCHAR AS chrom,
           strand::VARCHAR AS strand,
           score_mode::VARCHAR AS score_mode,
           pseudocount::DOUBLE AS pseudocount,
           background_model_id::VARCHAR AS background_model_id,
           pseudocount_scheme::VARCHAR AS pseudocount_scheme,
           minimum_pwm_relative_score::DOUBLE AS minimum_pwm_relative_score,
           maximum_pwm_relative_score::DOUBLE AS maximum_pwm_relative_score,
           n_policy::VARCHAR AS n_policy,
           matched_sequence_policy::VARCHAR AS matched_sequence_policy
    FROM read_json_auto({sql_string(stats_path)}, format='newline_delimited')
), observed AS (
    SELECT h.filename,
           COUNT(*)::BIGINT AS n_rows,
           COUNT(*) FILTER (
               WHERE h.score < e.minimum_score
                  OR h."end" - h.start <> e.motif_length
                  OR h.start < 0 OR h."end" > {task['sequence_length']}
                  OR NOT isfinite(h.score)
                  OR NOT isfinite(h.pwm_relative_score)
                  OR h.pwm_relative_score < 0 OR h.pwm_relative_score > 1
                  OR h.genome_id <> e.genome_id
                  OR h.motif_set_id <> e.motif_set_id
                  OR h.motif_id <> e.motif_id
                  OR CAST(h.chrom AS VARCHAR) <> e.chrom
                  OR CASE h.strand WHEN 'plus' THEN '+' WHEN 'minus' THEN '-'
                         ELSE h.strand END <> e.strand
                  OR h.score_mode <> e.score_mode
                  OR CAST(h.pseudocount AS DOUBLE) <> e.pseudocount
                  OR h.background_model_id <> e.background_model_id
                  OR h.pseudocount_scheme <> e.pseudocount_scheme
                  OR TRY_CAST(h.minimum_score AS DOUBLE) <> e.minimum_score
                  OR TRY_CAST(h.minimum_pwm_relative_score AS DOUBLE)
                         IS DISTINCT FROM e.minimum_pwm_relative_score
                  OR TRY_CAST(h.maximum_pwm_relative_score AS DOUBLE)
                         IS DISTINCT FROM e.maximum_pwm_relative_score
                  OR h.n_policy <> e.n_policy
                  OR h.matched_sequence <> e.matched_sequence_policy
           )::BIGINT AS invalid_rows
    FROM read_parquet({sql_string_list(output_paths)}, filename=true,
                      hive_partitioning=true) h
    JOIN expected e ON e.filename = h.filename
    GROUP BY h.filename
)
SELECT e.filename, e.expected_rows, COALESCE(o.n_rows, 0)::BIGINT AS n_rows,
       COALESCE(o.invalid_rows, 0)::BIGINT AS invalid_rows
FROM expected e
LEFT JOIN observed o USING (filename)
ORDER BY e.filename;
"""
    validation = run_duckdb_json(
        duckdb, query, memory_limit=duckdb_memory,
        temp_directory=duckdb_temp,
    )
    if len(validation) != expected_files:
        raise ScanError("DuckDB did not validate every planned Parquet file")
    for row in validation:
        if row["n_rows"] != row["expected_rows"] or row["invalid_rows"] != 0:
            raise ScanError(f"invalid sparse Parquet contents: {row['filename']}")

    for row in records:
        output = Path(row.pop("output_file"))
        row["output_relative_path"] = str(output.relative_to(package))
        metadata = output.stat()
        row["bytes"] = metadata.st_size
        row["mtime_ns"] = metadata.st_mtime_ns
        row["sha256"] = sha256_file(output)
        row["state"] = "complete"
    return records


def validate_promoted_task(
    task_root: Path,
    result: dict[str, Any],
    task: dict[str, Any],
    plan: dict[str, Any],
    plan_sha256: str,
    *,
    verify_checksums: bool = False,
    progress: OperationProgress | None = None,
    max_read_mib_per_second: float | None = None,
) -> None:
    if result.get("state") != "complete" \
            or result.get("plan_sha256") != plan_sha256 \
            or result.get("task") != task:
        raise ScanError(f"existing promoted task is inconsistent: {task_root}")
    inventory = result.get("inventory")
    if not isinstance(inventory, list) or len(inventory) != 2 * task["motif_count"]:
        raise ScanError(f"existing promoted task inventory is incomplete: {task_root}")
    expected_pairs = {
        (motif_id, strand)
        for motif_id in task["motif_ids"]
        for strand in ("+", "-")
    }
    observed_pairs = {(row.get("motif_id"), row.get("strand")) for row in inventory}
    if observed_pairs != expected_pairs:
        raise ScanError(f"existing promoted task has the wrong motif set: {task_root}")
    resolved_root = task_root.resolve()
    for row in inventory:
        output = (task_root / row.get("output_relative_path", "")).resolve()
        if resolved_root not in output.parents or not output.is_file():
            raise ScanError(f"promoted output is absent or outside its task: {output}")
        metadata = output.stat()
        if metadata.st_size != row.get("bytes"):
            raise ScanError(f"promoted output size mismatch: {output}")
        recorded_mtime = row.get("mtime_ns")
        if recorded_mtime is not None and metadata.st_mtime_ns != recorded_mtime:
            raise ScanError(f"promoted output modification time changed: {output}")
        if verify_checksums:
            if progress is not None:
                progress.current_path = str(output)
                progress.current_task = task["task_id"]
                progress.current_file_bytes = 0
            observed_sha256 = sha256_file(
                output, progress, max_read_mib_per_second
            )
            if observed_sha256 != row.get("sha256"):
                raise ScanError(f"promoted output checksum mismatch: {output}")
            if progress is not None:
                progress.files_complete += 1
                progress.bytes_complete += metadata.st_size
                progress.current_file_bytes = 0
        if row.get("run_id") != plan["run_id"] \
                or row.get("task_id") != task["task_id"] \
                or row.get("source_commit") != plan["source_commit"] \
                or row.get("state") != "complete":
            raise ScanError(f"promoted output provenance mismatch: {output}")


def validate_staged_input(plan: dict[str, Any], task: dict[str, Any],
                          genome: Path, fasta_index: Path) -> None:
    regions = parse_fasta_index(fasta_index)
    if regions != [{
        "sequence_order": 0,
        "chrom": task["chrom"],
        "length": task["sequence_length"],
    }]:
        raise ScanError(
            f"staged FASTA index is not exactly chromosome {task['chrom']}"
        )
    _, sequences, _ = inspect_fasta(genome)
    expected_region = next(
        row for row in plan["sequence_regions"] if row["chrom"] == task["chrom"]
    )
    if len(sequences) != 1 \
            or sequences[0]["chrom"] != task["chrom"] \
            or sequences[0]["length"] != task["sequence_length"] \
            or sequences[0]["sequence_sha256"] != expected_region["sequence_sha256"]:
        raise ScanError(f"staged FASTA identity mismatch for {task['chrom']}")


def validate_staged_input_metadata(plan: dict[str, Any], task: dict[str, Any],
                                   genome: Path, fasta_index: Path,
                                   metadata_path: Path) -> None:
    metadata = json_read(metadata_path)
    expected_region = next(
        row for row in plan["sequence_regions"] if row["chrom"] == task["chrom"]
    )
    current = genome.stat()
    current_index = fasta_index.stat()
    if metadata.get("sequence") != task["chrom"] \
            or metadata.get("length") != task["sequence_length"] \
            or metadata.get("sequence_sha256") != expected_region["sequence_sha256"] \
            or Path(str(metadata.get("staged_fasta", ""))).resolve() != genome \
            or Path(str(metadata.get("staged_fasta_index", ""))).resolve() != fasta_index \
            or metadata.get("staged_bytes") != current.st_size \
            or metadata.get("staged_mtime_ns") != current.st_mtime_ns \
            or metadata.get("staged_fasta_index_bytes") != current_index.st_size \
            or metadata.get("staged_fasta_index_mtime_ns") != current_index.st_mtime_ns \
            or metadata.get("staged_fasta_index_sha256") != sha256_file(fasta_index):
        raise ScanError(f"staged FASTA metadata mismatch for {task['chrom']}")


def copy_file_verified(source: Path, destination: Path,
                       expected_sha256: str) -> None:
    if destination.exists():
        raise ScanError(f"refusing to replace durable task output: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256()
    with source.open("rb") as input_stream, destination.open("xb") as output_stream:
        while chunk := input_stream.read(8 * 1024 * 1024):
            digest.update(chunk)
            output_stream.write(chunk)
        output_stream.flush()
        os.fsync(output_stream.fileno())
    if digest.hexdigest() != expected_sha256:
        raise ScanError(f"scratch output changed while publishing: {source}")


def publish_scratch_task(source_package: Path, durable_package: Path,
                         inventory: list[dict[str, Any]]) -> None:
    durable_package.mkdir(parents=True, exist_ok=False)
    for row in inventory:
        relative = Path(row["output_relative_path"])
        source = source_package / relative
        destination = durable_package / relative
        copy_file_verified(source, destination, row["sha256"])
        metadata = destination.stat()
        if metadata.st_size != row["bytes"]:
            raise ScanError(f"durable output size mismatch after copy: {destination}")
        row["mtime_ns"] = metadata.st_mtime_ns


def run_task(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan, plan_path, plan_sha256 = load_plan(run_root)
    check_free_space(run_root, plan["minimum_free_gib"])
    task = task_for_arguments(plan, arguments)
    if not getattr(arguments, "planned_inputs_verified", False):
        validate_planned_inputs(plan)
    scanner = absolute_file(arguments.scanner)
    if shutil.which(arguments.duckdb) is None:
        raise ScanError(f"DuckDB executable not found: {arguments.duckdb}")
    duckdb_temp: Path | None = None
    if arguments.duckdb_temp_directory is not None:
        duckdb_temp = arguments.duckdb_temp_directory.expanduser().resolve()
        duckdb_temp.mkdir(parents=True, exist_ok=True)
    if not arguments.allow_local and not os.environ.get("SLURM_JOB_ID"):
        raise ScanError("run-task requires Slurm; use --allow-local only for testing")
    motif_list = run_root / task["motif_list"]
    if sha256_file(motif_list) != task["motif_list_sha256"]:
        raise ScanError(f"motif list checksum mismatch: {motif_list}")

    package_root = run_root / "package"
    final_task = package_root / "task_data" / f"task_id={task['task_id']}"
    if final_task.exists():
        result = json_read(final_task / "task_result.json")
        validate_promoted_task(final_task, result, task, plan, plan_sha256)
        print(f"I: Reusing validated task {task['task_id']}", file=sys.stderr)
        return

    genome_override = getattr(arguments, "genome_override", None)
    fasta_index_override = getattr(arguments, "fasta_index_override", None)
    if (genome_override is None) != (fasta_index_override is None):
        raise ScanError("staged genome and FASTA index must be supplied together")
    genome_input = absolute_file(
        genome_override if genome_override is not None
        else plan["genome"]["fasta_file"]
    )
    fasta_index_input = absolute_file(
        fasta_index_override if fasta_index_override is not None
        else plan["genome"]["fasta_index_file"]
    )
    if genome_override is not None:
        staged_metadata = getattr(arguments, "staged_input_metadata", None)
        if staged_metadata is not None:
            validate_staged_input_metadata(
                plan, task, genome_input, fasta_index_input,
                absolute_file(staged_metadata),
            )
        elif not getattr(arguments, "staged_sequence_verified", False):
            validate_staged_input(plan, task, genome_input, fasta_index_input)
    jaspar_override = getattr(arguments, "jaspar_override", None)
    jaspar_input = absolute_file(
        jaspar_override if jaspar_override is not None
        else plan["motif_set"]["source_file"]
    )
    if (jaspar_override is not None
            and sha256_file(jaspar_input) != plan["motif_set"]["source_sha256"]):
        raise ScanError("staged JASPAR checksum differs from the scan plan")

    job_id = os.environ.get("SLURM_JOB_ID", "local")
    array_id = os.environ.get("SLURM_ARRAY_TASK_ID", str(task["task_index"]))
    restart_count = os.environ.get("SLURM_RESTART_COUNT", "0")
    attempt = f"job-{job_id}-array-{array_id}-restart-{restart_count}-pid-{os.getpid()}"
    staging = run_root / "staging" / task["task_id"] / attempt
    staging.mkdir(parents=True, exist_ok=False)
    task_work_directory = getattr(arguments, "task_work_directory", None)
    if task_work_directory is None:
        work_attempt = staging
    else:
        work_attempt = (
            Path(task_work_directory).expanduser().resolve()
            / task["task_id"] / attempt
        )
        work_attempt.mkdir(parents=True, exist_ok=False)
    task_package = work_attempt / "task_package"
    stats_path = work_attempt / "scan_file_stats.jsonl"
    command = [
        str(scanner),
        "--sparse-parquet",
        "--outdir", str(task_package),
        "--genome", str(genome_input),
        "--fasta-index", str(fasta_index_input),
        "--pssm", str(jaspar_input),
        "--motif-list", str(motif_list),
        "--motif-set-id", plan["motif_set"]["motif_set_id"],
        "--genome-id", plan["genome"]["genome_id"],
        "--chr", task["chrom"],
        "--strand", "both",
        "--coordinate-mode", "bed",
        "--score-mode", plan["scanner"]["score_mode"],
        "--pseudocount", str(plan["scanner"]["pseudocount"]),
        "--background-model-id", plan["scanner"]["background_model_id"],
        "--pseudocount-scheme", plan["scanner"]["pseudocount_scheme"],
        "--threshold", str(task["minimum_score"]),
        "--skip-N",
        "--scan-file-stats", str(stats_path),
    ]
    child: subprocess.Popen[str] | None = None
    task_progress = OperationProgress("task_prepare")
    task_progress.current_task = task["task_id"]
    task_progress.current_path = str(task_package)

    def forward_progress(_signal_number: int, _frame: Any) -> None:
        if child is not None and child.poll() is None:
            try:
                child.send_signal(signal.SIGUSR1)
                return
            except (ProcessLookupError, OSError):
                # The scanner may finish between poll() and signal delivery.
                pass
        task_progress.emit()

    if hasattr(signal, "SIGUSR1"):
        signal.signal(signal.SIGUSR1, forward_progress)
    try:
        print(f"I: Starting task {task['task_id']}", file=sys.stderr)
        scanner_sha256 = sha256_file(scanner)
        build_info = getattr(arguments, "scanner_build_info", None)
        if build_info is None:
            build_info = scanner_build_info(scanner)
        build_matches_plan = (
            build_info.get("source_commit") == plan["source_commit"]
            and not build_info.get("source_dirty", True)
        )
        if not build_matches_plan \
                and not arguments.allow_scanner_provenance_mismatch:
            raise ScanError(
                "scanner build is dirty or does not match the planned source "
                "commit; rebuild from the clean planned commit or pass the "
                "explicit provenance-mismatch override"
            )
        task_progress.set_phase("scanner")
        child = subprocess.Popen(command, text=True)
        return_code = child.wait()
        if return_code != 0:
            raise ScanError(f"scanner exited with status {return_code}")
        task_progress.set_phase("parquet_validation")
        inventory = validate_task_outputs(
            arguments.duckdb, stats_path, task_package, task, plan,
            arguments.duckdb_memory_limit, duckdb_temp,
        )
        for row in inventory:
            row.update(
                {
                    "run_id": plan["run_id"],
                    "task_id": task["task_id"],
                    "policy_id": task["policy_id"],
                    "source_commit": plan["source_commit"],
                    "scanner_sha256": scanner_sha256,
                    "slurm_job_id": job_id,
                    "slurm_array_task_id": array_id,
                    "slurm_restart_count": int(restart_count),
                }
            )
        durable_package = task_package
        if work_attempt != staging:
            task_progress.set_phase("durable_publication")
            durable_package = staging / "task_package"
            publish_scratch_task(task_package, durable_package, inventory)
        result = {
            "state": "complete",
            "completed_at_utc": utc_now(),
            "plan_sha256": plan_sha256,
            "plan_file": str(plan_path),
            "source_commit": plan["source_commit"],
            "scanner_file": str(scanner),
            "scanner_sha256": scanner_sha256,
            "scanner_build": build_info,
            "scanner_source_commit_matches_plan": (
                build_matches_plan
            ),
            "slurm_job_id": job_id,
            "slurm_array_task_id": array_id,
            "slurm_restart_count": int(restart_count),
            "task": task,
            "inventory": inventory,
        }
        json_write_new(durable_package / "task_result.json", result)
        final_task.parent.mkdir(parents=True, exist_ok=True)
        if final_task.exists():
            raise ScanError(f"task appeared concurrently: {final_task}")
        task_progress.set_phase("atomic_promotion")
        os.replace(durable_package, final_task)
        print(f"I: Validated and promoted task {task['task_id']}", file=sys.stderr)
    except Exception as error:
        failure_path = staging / "failure.json"
        if not failure_path.exists():
            json_write_new(
                failure_path,
                {
                    "state": "failed",
                    "failed_at_utc": utc_now(),
                    "task_id": task["task_id"],
                    "error": str(error),
                    "command": command,
                },
            )
        raise


def selected_sequence_regions(plan: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        region for region in plan["sequence_regions"]
        if region.get("included_in_scan", True)
    ]


def run_chromosome(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan, _, _ = load_plan(run_root)
    regions = selected_sequence_regions(plan)
    if arguments.chrom is not None:
        matches = [region for region in regions if region["chrom"] == arguments.chrom]
        if len(matches) != 1:
            raise ScanError(f"scan plan has no unique chromosome {arguments.chrom!r}")
        region = matches[0]
    else:
        if arguments.chrom_index < 0 or arguments.chrom_index >= len(regions):
            raise ScanError(
                f"chromosome index {arguments.chrom_index} is outside "
                f"0..{len(regions) - 1}"
            )
        region = regions[arguments.chrom_index]

    if not arguments.allow_local and not os.environ.get("SLURM_JOB_ID"):
        raise ScanError(
            "run-chromosome requires Slurm; use --allow-local only for testing"
        )
    scratch = arguments.scratch_directory.expanduser().resolve()
    scratch.mkdir(parents=True, exist_ok=False)
    required_scratch_gib = (
        arguments.minimum_scratch_free_gib
        + (region["length"] + 1024**2) / 1024**3
    )
    check_free_space(scratch, required_scratch_gib, "Scratch filesystem")
    check_free_space(run_root, plan["minimum_free_gib"])
    validate_planned_inputs(plan)

    reporter = OperationProgress("fasta_stage")
    reporter.install()
    reporter.bytes_total = region["length"]
    reporter.files_total = 1
    staged_fasta = scratch / f"{safe_label(region['chrom'])}.fa"
    reporter.current_path = str(staged_fasta)
    reporter.current_task = f"chromosome:{region['chrom']}"

    def stage_progress(copied: int, _total: int) -> None:
        reporter.current_file_bytes = copied

    stage_result = stage_fasta_region(
        absolute_file(plan["genome"]["fasta_file"]),
        absolute_file(plan["genome"]["fasta_index_file"]),
        region["chrom"],
        staged_fasta,
        expected_length=region["length"],
        expected_sha256=region["sequence_sha256"],
        progress_callback=stage_progress,
    )
    reporter.files_complete = 1
    reporter.bytes_complete = region["length"]
    reporter.current_file_bytes = 0
    json_write_new(scratch / "staged_sequence.json", stage_result)
    staged_jaspar = scratch / Path(plan["motif_set"]["source_file"]).name
    copy_file_verified(
        absolute_file(plan["motif_set"]["source_file"]), staged_jaspar,
        plan["motif_set"]["source_sha256"],
    )
    print(
        f"I: Staged and verified chromosome {region['chrom']} at {staged_fasta} "
        f"and JASPAR at {staged_jaspar}",
        file=sys.stderr,
    )

    scanner = absolute_file(arguments.scanner)
    build_info = scanner_build_info(scanner)
    build_matches_plan = (
        build_info.get("source_commit") == plan["source_commit"]
        and not build_info.get("source_dirty", True)
    )
    if not build_matches_plan \
            and not arguments.allow_scanner_provenance_mismatch:
        raise ScanError(
            "scanner build is dirty or does not match the planned source commit"
        )
    tasks = [task for task in plan["tasks"] if task["chrom"] == region["chrom"]]
    if not tasks:
        raise ScanError(f"scan plan has no tasks for chromosome {region['chrom']}")
    print(
        f"I: Chromosome worker will process {len(tasks)} motif batches for "
        f"{region['chrom']} from node-local FASTA.",
        file=sys.stderr,
    )
    active: dict[subprocess.Popen[str], tuple[int, dict[str, Any]]] = {}
    next_task = 0
    failed_task: str | None = None

    def worker_progress(_signal_number: int, _frame: Any) -> None:
        for process in active:
            if process.poll() is None:
                try:
                    process.send_signal(signal.SIGUSR1)
                    return
                except (ProcessLookupError, OSError):
                    # The batch may finish between poll() and signal delivery.
                    continue
        reporter.emit()

    if hasattr(signal, "SIGUSR1"):
        signal.signal(signal.SIGUSR1, worker_progress)

    while active or (next_task < len(tasks) and failed_task is None):
        while failed_task is None \
                and next_task < len(tasks) \
                and len(active) < arguments.batch_workers:
            task = tasks[next_task]
            batch_index = next_task + 1
            check_free_space(run_root, plan["minimum_free_gib"])
            check_free_space(
                scratch, arguments.minimum_scratch_free_gib, "Scratch filesystem"
            )
            reporter.set_phase(
                "motif_batch",
                current_path=staged_fasta,
                current_task=f"{task['task_id']}:{batch_index}/{len(tasks)}",
            )
            print(
                f"I: Chromosome {region['chrom']} batch {batch_index}/{len(tasks)}: "
                f"{task['task_id']}",
                file=sys.stderr,
            )
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "run-task",
                "--run-root", str(run_root),
                "--task-id", task["task_id"],
                "--scanner", str(scanner),
                "--duckdb", arguments.duckdb,
                "--duckdb-memory-limit", arguments.duckdb_memory_limit,
                "--duckdb-temp-directory", str(scratch / "duckdb" / task["task_id"]),
                "--genome-override", str(staged_fasta),
                "--fasta-index-override", str(stage_result["staged_fasta_index"]),
                "--jaspar-override", str(staged_jaspar),
                "--staged-input-metadata", str(scratch / "staged_sequence.json"),
                "--planned-inputs-verified",
            ]
            if arguments.scratch_task_output:
                command.extend(["--task-work-directory", str(scratch / "task-work")])
            if arguments.allow_local:
                command.append("--allow-local")
            if arguments.allow_scanner_provenance_mismatch:
                command.append("--allow-scanner-provenance-mismatch")
            process = subprocess.Popen(command, text=True)
            active[process] = (batch_index, task)
            next_task += 1

        if not active:
            break
        time.sleep(0.1)
        for process, (_batch_index, task) in list(active.items()):
            return_code = process.poll()
            if return_code is None:
                continue
            del active[process]
            if return_code != 0 and failed_task is None:
                failed_task = task["task_id"]

    if failed_task is not None:
        raise ScanError(
            f"chromosome worker stopped after failed batch {failed_task}; "
            "completed batches remain reusable"
        )

    reporter.set_phase("chromosome_complete", current_path=staged_fasta)
    print(
        f"I: Completed all motif batches for chromosome {region['chrom']}; "
        "node-local files are left for the cluster scratch lifecycle.",
        file=sys.stderr,
    )


def write_json_lines(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True, separators=(",", ":")))
            stream.write("\n")


def parquet_from_json(duckdb: str, source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    query = (
        f"COPY (SELECT * FROM read_json_auto({sql_string(source)}, "
        "format='newline_delimited')) "
        f"TO {sql_string(destination)} (FORMAT PARQUET, COMPRESSION ZSTD);"
    )
    process = subprocess.run(
        [duckdb, ":memory:", "-bail", "-c", query],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ScanError(process.stderr.strip() or "DuckDB Parquet export failed")


def validate_existing_metadata(
    duckdb: str,
    existing_root: Path,
    staged_root: Path,
    dataset_names: Iterable[str],
) -> None:
    for name in dataset_names:
        existing = existing_root / name / "part-000000.parquet"
        staged = staged_root / name / "part-000000.parquet"
        query = f"""
SELECT COUNT(*)::BIGINT AS differences
FROM (
            (SELECT * FROM read_parquet({sql_string(existing)})
     EXCEPT ALL
     SELECT * FROM read_parquet({sql_string(staged)}))
    UNION ALL
    (SELECT * FROM read_parquet({sql_string(staged)})
     EXCEPT ALL
     SELECT * FROM read_parquet({sql_string(existing)}))
);
"""
        rows = run_duckdb_json(duckdb, query)
        if len(rows) != 1 or rows[0].get("differences") != 0:
            raise ScanError(
                f"existing final metadata differs from staged {name}: {existing_root}"
            )


def build_catalog_database(duckdb: str, package: Path, database: Path,
                           progress: OperationProgress | None = None) -> None:
    if database.exists():
        raise ScanError(f"refusing to replace DuckDB catalog: {database}")
    schema = Path(__file__).resolve().parent.parent / "sql" / "genome_scan_schema.sql"
    if progress is not None:
        progress.set_phase("duckdb_catalog", current_path=database)
    process = subprocess.run(
        [duckdb, str(database), "-bail", "-f", str(schema)],
        cwd=package,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ScanError(process.stderr.strip() or "could not build DuckDB query index")
    threshold_table = (
        package / "tables" / "jaspar2026" / "scan_motif_threshold" /
        "part-000000.parquet"
    )
    if threshold_table.is_file():
        query = (
            "CREATE OR REPLACE TABLE scan_motif_threshold AS "
            f"SELECT * FROM read_parquet({sql_string(threshold_table)});"
        )
        replacement = subprocess.run(
            [duckdb, str(database), "-bail", "-c", query],
            cwd=package, text=True, capture_output=True, check=False,
        )
        if replacement.returncode != 0:
            raise ScanError(
                replacement.stderr.strip()
                or "could not install detailed motif-threshold table"
            )


def validate_existing_database(
    duckdb: str,
    database: Path,
    package: Path,
    plan: dict[str, Any],
    expected_files: int,
    expected_hits: int,
) -> None:
    query = f"""
SELECT
    (SELECT COUNT(*) FROM scan_run
     WHERE run_id = {sql_string(plan['run_id'])}
       AND source_commit = {sql_string(plan['source_commit'])}
       AND state = 'complete')::BIGINT AS matching_runs,
    (SELECT COUNT(*) FROM scan_file_inventory)::BIGINT AS file_count,
    (SELECT COALESCE(SUM(emitted_hits), 0)
     FROM scan_file_inventory)::BIGINT AS emitted_hits;
"""
    process = subprocess.run(
        [duckdb, "-readonly", "-json", str(database), "-c", query],
        cwd=package,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ScanError(
            process.stderr.strip() or f"existing DuckDB index is invalid: {database}"
        )
    try:
        rows = json.loads(process.stdout or "[]")
    except json.JSONDecodeError as error:
        raise ScanError(f"existing DuckDB index returned invalid JSON: {error}") from error
    if len(rows) != 1 or rows[0] != {
        "matching_runs": 1,
        "file_count": expected_files,
        "emitted_hits": expected_hits,
    }:
        raise ScanError(f"existing DuckDB index does not match the scan plan: {database}")
    if plan.get("schema_version", 1) >= 2:
        expected_threshold_set = (
            (plan.get("threshold_metadata") or {}).get("threshold_set_id")
            or f"legacy_{plan['run_id']}"
        )
        threshold_query = f"""
WITH inventory_threshold AS (
    SELECT motif_id, MIN(CAST(minimum_score AS DOUBLE)) AS minimum_score,
           COUNT(DISTINCT CAST(minimum_score AS DOUBLE)) AS score_count
    FROM scan_file_inventory
    GROUP BY motif_id
)
SELECT
    (SELECT COUNT(*) FROM scan_motif_threshold)::BIGINT AS threshold_rows,
    (SELECT COUNT(DISTINCT motif_id)
     FROM scan_motif_threshold)::BIGINT AS threshold_motifs,
    (SELECT COUNT(*) FROM scan_motif_threshold
     WHERE threshold_set_id <> {sql_string(expected_threshold_set)})::BIGINT
        AS wrong_threshold_set,
    (SELECT COUNT(*)
     FROM scan_motif_threshold t
     JOIN inventory_threshold i USING (motif_id)
     WHERE i.score_count <> 1
        OR t.final_minimum_score <> i.minimum_score)::BIGINT AS score_mismatches;
"""
        threshold_process = subprocess.run(
            [duckdb, "-readonly", "-json", str(database), "-c", threshold_query],
            cwd=package, text=True, capture_output=True, check=False,
        )
        if threshold_process.returncode != 0:
            raise ScanError(
                threshold_process.stderr.strip()
                or f"existing threshold catalog is invalid: {database}"
            )
        try:
            threshold_rows = json.loads(threshold_process.stdout or "[]")
        except json.JSONDecodeError as error:
            raise ScanError(
                f"existing threshold catalog returned invalid JSON: {error}"
            ) from error
        expected_motifs = int(plan["motif_set"]["motif_count"])
        expected = {
            "threshold_rows": expected_motifs,
            "threshold_motifs": expected_motifs,
            "wrong_threshold_set": 0,
            "score_mismatches": 0,
        }
        if len(threshold_rows) != 1 or threshold_rows[0] != expected:
            raise ScanError(
                f"existing motif-threshold catalog differs from the plan: {database}"
            )


def finalize(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan, _, plan_sha256 = load_plan(run_root)
    package = run_root / "package"
    manifest_path = package / "manifest.json"
    if shutil.which(arguments.duckdb) is None:
        raise ScanError(f"DuckDB executable not found: {arguments.duckdb}")
    if manifest_path.exists():
        manifest = json_read(manifest_path)
        if manifest.get("plan_sha256") != plan_sha256:
            raise ScanError("existing manifest belongs to a different scan plan")
        expected_files = manifest.get("file_count")
        expected_hits = manifest.get("emitted_hit_count")
        if not isinstance(expected_files, int) or not isinstance(expected_hits, int):
            raise ScanError("existing manifest lacks file or emitted-hit totals")
        database_name = manifest.get("database") or "jaspar_genome_scan.duckdb"
        database = package / str(database_name)
        if not database.is_file():
            raise ScanError(f"finalized DuckDB catalog is absent: {database}")
        validate_existing_database(
            arguments.duckdb,
            database,
            package,
            plan,
            expected_files,
            expected_hits,
        )
        print(f"I: Reusing finalized package: {package}", file=sys.stderr)
        return

    reporter = OperationProgress("task_inventory")
    reporter.install()
    reporter.files_total = 2 * sum(
        task["motif_count"] for task in plan["tasks"]
    )
    task_results: list[dict[str, Any]] = []
    missing: list[str] = []
    for task in plan["tasks"]:
        task_path = package / "task_data" / f"task_id={task['task_id']}"
        result_path = task_path / "task_result.json"
        if not result_path.is_file():
            missing.append(task["task_id"])
            continue
        result = json_read(result_path)
        validate_promoted_task(task_path, result, task, plan, plan_sha256)
        task_results.append(result)
        reporter.files_complete += len(result["inventory"])
        reporter.bytes_complete += sum(row["bytes"] for row in result["inventory"])
    if missing:
        raise ScanError(
            f"cannot finalize: {len(missing)} of {len(plan['tasks'])} tasks are missing"
        )

    inventory: list[dict[str, Any]] = []
    seen_files: set[tuple[str, str, str]] = set()
    scanner_sha256s = {result.get("scanner_sha256") for result in task_results}
    if None in scanner_sha256s or len(scanner_sha256s) != 1:
        raise ScanError("promoted tasks were produced by different scanner binaries")
    scanner_builds = {
        json.dumps(result.get("scanner_build"), sort_keys=True)
        for result in task_results if result.get("scanner_build") is not None
    }
    if len(scanner_builds) > 1:
        raise ScanError("promoted tasks report different scanner builds")
    scanner_build = (
        json.loads(next(iter(scanner_builds))) if scanner_builds else None
    )
    for result in task_results:
        for row in result["inventory"]:
            key = (row["motif_id"], row["chrom"], row["strand"])
            if key in seen_files:
                raise ScanError(f"duplicate promoted motif/chromosome/strand file: {key}")
            seen_files.add(key)
            inventory.append(row)

    expected_files = 2 * sum(task["motif_count"] for task in plan["tasks"])
    if len(inventory) != expected_files:
        raise ScanError(
            f"inventory has {len(inventory)} files; expected {expected_files}"
        )

    attempt = f"finalize-{utc_now().replace(':', '')}-pid-{os.getpid()}"
    staging = run_root / "staging" / attempt
    staging_package = staging / "package"
    staging_package.mkdir(parents=True, exist_ok=False)
    metadata_root = staging_package / "tables" / "jaspar2026"
    completed_at = utc_now()

    motif_set = dict(plan["motif_set"])
    motif_set.pop("source_file", None)
    motif_metadata = []
    for motif in plan["motif_metadata"]:
        row = dict(motif)
        row["source_sha256"] = plan["motif_set"]["source_sha256"]
        motif_metadata.append(row)
    genome = dict(plan["genome"])
    genome.pop("fasta_file", None)
    genome.pop("fasta_index_file", None)
    sequence_regions = [
        {"genome_id": plan["genome"]["genome_id"], **row}
        for row in plan["sequence_regions"]
    ]
    threshold_policies = [
        {"run_id": plan["run_id"], **row} for row in plan["threshold_policies"]
    ]
    threshold_metadata = plan.get("threshold_metadata") or {}
    threshold_registry = plan.get("threshold_registry") or {}
    planned_motif_thresholds = plan.get("motif_thresholds")
    if not planned_motif_thresholds:
        threshold_by_motif: dict[str, float] = {}
        for task in plan["tasks"]:
            for motif_id in task["motif_ids"]:
                observed = threshold_by_motif.setdefault(
                    motif_id, float(task["minimum_score"])
                )
                if observed != float(task["minimum_score"]):
                    raise ScanError(
                        f"legacy plan uses multiple minimum scores for {motif_id}"
                    )
        planned_motif_thresholds = [
            {
                "motif_id": motif_id,
                "informative_threshold": None,
                "informative_source": "legacy_scan_policy",
                "default_minimum_score": None,
                "candidate_minimum_score": minimum_score,
                "density_minimum_spacing_bp": None,
                "density_maximum_loci": None,
                "density_threshold": None,
                "final_minimum_score": minimum_score,
                "density_limited": None,
                "density_chrom": None,
                "valid_locus_starts": None,
                "skipped_locus_starts": None,
                "loci_at_candidate_threshold": None,
                "loci_at_final_threshold": None,
                "mean_spacing_bp_at_final_threshold": None,
                "distribution_sha256": None,
            }
            for motif_id, minimum_score in sorted(threshold_by_motif.items())
        ]
    motif_thresholds = [
        {
            "run_id": plan["run_id"],
            "genome_id": plan["genome"]["genome_id"],
            "motif_set_id": plan["motif_set"]["motif_set_id"],
            "threshold_set_id": threshold_metadata.get(
                "threshold_set_id", f"legacy_{plan['run_id']}"
            ),
            "threshold_registry_sha256": threshold_registry.get("source_sha256"),
            **row,
        }
        for row in planned_motif_thresholds
    ]
    scan_tasks = []
    for result in task_results:
        task = result["task"]
        scan_tasks.append(
            {
                "run_id": plan["run_id"],
                "task_id": task["task_id"],
                "task_index": task["task_index"],
                "chrom": task["chrom"],
                "motif_batch_id": task["motif_batch_id"],
                "policy_id": task["policy_id"],
                "minimum_score": task["minimum_score"],
                "motif_count": task["motif_count"],
                "motif_list_sha256": task["motif_list_sha256"],
                "state": "complete",
                "slurm_job_id": result["slurm_job_id"],
                "slurm_array_task_id": result["slurm_array_task_id"],
                "slurm_restart_count": result["slurm_restart_count"],
                "completed_at_utc": result["completed_at_utc"],
                "source_commit": result["source_commit"],
                "scanner_sha256": result["scanner_sha256"],
                "scanner_build_json": (
                    json.dumps(result.get("scanner_build"), sort_keys=True)
                    if result.get("scanner_build") is not None else None
                ),
                "scanner_source_commit_matches_plan": result.get(
                    "scanner_source_commit_matches_plan"
                ),
            }
        )
    scan_run = {
        "run_id": plan["run_id"],
        "state": "complete",
        "planned_at_utc": plan["created_at_utc"],
        "completed_at_utc": completed_at,
        "source_commit": plan["source_commit"],
        "genome_id": plan["genome"]["genome_id"],
        "motif_set_id": plan["motif_set"]["motif_set_id"],
        "motif_threshold_set_id": threshold_metadata.get("threshold_set_id"),
        "motif_threshold_registry_sha256": threshold_registry.get("source_sha256"),
        "motif_threshold_metadata_sha256": threshold_registry.get("metadata_sha256"),
        "score_mode": plan["scanner"]["score_mode"],
        "pseudocount": plan["scanner"]["pseudocount"],
        "background_model_id": plan["scanner"]["background_model_id"],
        "pseudocount_scheme": plan["scanner"]["pseudocount_scheme"],
        "minimum_pwm_relative_score": plan["scanner"][
            "minimum_pwm_relative_score"
        ],
        "maximum_pwm_relative_score": plan["scanner"][
            "maximum_pwm_relative_score"
        ],
        "coordinate_mode": plan["scanner"]["coordinate_mode"],
        "n_policy": plan["scanner"]["n_policy"],
        "matched_sequence_policy": plan["scanner"]["matched_sequence_policy"],
        "task_count": len(scan_tasks),
        "file_count": len(inventory),
        "emitted_hit_count": sum(row["emitted_hits"] for row in inventory),
        "parquet_bytes": sum(row["bytes"] for row in inventory),
        "scanner_sha256": next(iter(scanner_sha256s)),
        "scanner_build_json": (
            json.dumps(scanner_build, sort_keys=True)
            if scanner_build is not None else None
        ),
        "duckdb_version": duckdb_build_info(arguments.duckdb),
        "finalization_validation_mode": "exact_inventory_size_mtime_when_available",
    }
    datasets = {
        "motif_set": [motif_set],
        "genome": [genome],
        "sequence_region": sequence_regions,
        "motif_metadata": motif_metadata,
        "scan_run": [scan_run],
        "scan_threshold_policy": threshold_policies,
        "scan_motif_threshold": motif_thresholds,
        "scan_task": scan_tasks,
        "scan_file_inventory": inventory,
    }
    for name, rows in datasets.items():
        source = staging / f"{name}.jsonl"
        destination = metadata_root / name / "part-000000.parquet"
        write_json_lines(source, rows)
        parquet_from_json(arguments.duckdb, source, destination)

    database = staging_package / "jaspar_genome_scan.duckdb"
    build_catalog_database(arguments.duckdb, staging_package, database, reporter)

    manifest = {
        "schema_version": 2,
        "run_id": plan["run_id"],
        "state": "complete",
        "completed_at_utc": completed_at,
        "plan_sha256": plan_sha256,
        "source_commit": plan["source_commit"],
        "genome_id": plan["genome"]["genome_id"],
        "motif_set_id": plan["motif_set"]["motif_set_id"],
        "motif_threshold_set_id": threshold_metadata.get("threshold_set_id"),
        "motif_threshold_registry_sha256": threshold_registry.get("source_sha256"),
        "motif_threshold_metadata_sha256": threshold_registry.get("metadata_sha256"),
        "task_count": len(scan_tasks),
        "file_count": len(inventory),
        "emitted_hit_count": scan_run["emitted_hit_count"],
        "parquet_bytes": scan_run["parquet_bytes"],
        "scanner_sha256": next(iter(scanner_sha256s)),
        "scanner_build": scanner_build,
        "duckdb_version": scan_run["duckdb_version"],
        "finalization_validation_mode": scan_run["finalization_validation_mode"],
        "checksum_audit": "verification/checksum_audit.json",
        "database": database.name,
    }
    staged_manifest = staging_package / "manifest.json"
    json_write_new(staged_manifest, manifest)

    metadata_destination = package / "tables" / "jaspar2026"
    metadata_destination.parent.mkdir(parents=True, exist_ok=True)
    if metadata_destination.exists():
        validate_existing_metadata(
            arguments.duckdb, metadata_destination, metadata_root, datasets
        )
    else:
        os.replace(metadata_root, metadata_destination)
    final_database = package / database.name
    if final_database.exists():
        validate_existing_database(
            arguments.duckdb,
            final_database,
            package,
            plan,
            len(inventory),
            scan_run["emitted_hit_count"],
        )
    else:
        os.replace(database, final_database)
    os.replace(staged_manifest, manifest_path)
    print(f"I: Finalized sparse scan package: {package}", file=sys.stderr)


def rebuild_catalog(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan, _, _ = load_plan(run_root)
    package = run_root / "package"
    metadata = package / "tables" / "jaspar2026" / "scan_file_inventory"
    if not metadata.is_dir():
        raise ScanError(f"finalized scan metadata is absent: {metadata}")
    if shutil.which(arguments.duckdb) is None:
        raise ScanError(f"DuckDB executable not found: {arguments.duckdb}")
    output = (
        arguments.output.expanduser().resolve()
        if arguments.output is not None
        else package / "jaspar_genome_scan.rebuilt.duckdb"
    )
    reporter = OperationProgress("duckdb_catalog")
    reporter.install()
    build_catalog_database(arguments.duckdb, package, output, reporter)
    expected_files = 2 * sum(task["motif_count"] for task in plan["tasks"])
    query = (
        "SELECT COUNT(*)::BIGINT AS files, "
        "COALESCE(SUM(emitted_hits), 0)::BIGINT AS hits "
        "FROM scan_file_inventory;"
    )
    process = subprocess.run(
        [arguments.duckdb, "-readonly", "-json", str(output), "-c", query],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ScanError(process.stderr.strip() or "rebuilt catalog validation failed")
    rows = json.loads(process.stdout or "[]")
    if len(rows) != 1 or rows[0].get("files") != expected_files:
        raise ScanError("rebuilt catalog inventory does not match the immutable plan")
    print(
        f"I: Built metadata-only DuckDB catalog without opening hit payloads: {output}",
        file=sys.stderr,
    )


def promoted_inventory(
    run_root: Path,
    plan: dict[str, Any],
    plan_sha256: str,
) -> list[tuple[dict[str, Any], dict[str, Any], Path]]:
    package = run_root / "package"
    entries: list[tuple[dict[str, Any], dict[str, Any], Path]] = []
    for task in plan["tasks"]:
        task_root = package / "task_data" / f"task_id={task['task_id']}"
        result_path = task_root / "task_result.json"
        if not result_path.is_file():
            raise ScanError(f"promoted task result is missing: {result_path}")
        result = json_read(result_path)
        validate_promoted_task(task_root, result, task, plan, plan_sha256)
        for row in result["inventory"]:
            entries.append((task, row, task_root / row["output_relative_path"]))
    return entries


def inventory_fingerprint(
    entries: Iterable[tuple[dict[str, Any], dict[str, Any], Path]],
    *,
    include_current_stat: bool = False,
) -> str:
    digest = hashlib.sha256()
    for task, row, path in entries:
        fields: list[object] = [
            task["task_id"], row["output_relative_path"], row["bytes"], row["sha256"]
        ]
        if include_current_stat:
            metadata = path.stat()
            fields.extend([metadata.st_size, metadata.st_mtime_ns])
        digest.update(("\t".join(map(str, fields)) + "\n").encode("utf-8"))
    return digest.hexdigest()


def read_verification_checkpoints(
    checkpoint_directory: Path,
    plan_sha256: str,
    inventory_sha256: str,
) -> dict[tuple[str, str, str, int, int], dict[str, Any]]:
    verified: dict[tuple[str, str, str, int, int], dict[str, Any]] = {}
    if not checkpoint_directory.is_dir():
        return verified
    for path in sorted(checkpoint_directory.iterdir()):
        if not path.is_file() or not path.name.startswith("checksum-") \
                or path.suffix != ".jsonl":
            continue
        with path.open(encoding="utf-8") as stream:
            for line_number, line in enumerate(stream, start=1):
                if not line.strip():
                    continue
                try:
                    row = json.loads(line)
                except json.JSONDecodeError:
                    print(
                        f"W: Ignoring incomplete checkpoint line {line_number} in {path}",
                        file=sys.stderr,
                    )
                    break
                if row.get("plan_sha256") != plan_sha256 \
                        or row.get("inventory_sha256") != inventory_sha256 \
                        or row.get("state") != "verified":
                    continue
                key = (
                    row["task_id"], row["output_relative_path"], row["sha256"],
                    row["bytes"], row["file_mtime_ns"],
                )
                verified[key] = row
    return verified


def verification_status_path(package: Path) -> Path:
    return package / "verification" / "progress.json"


def write_verification_status(
    package: Path,
    reporter: OperationProgress,
    state: str,
    inventory_sha256: str,
) -> None:
    json_write_atomic(
        verification_status_path(package),
        {
            "state": state,
            "updated_at_utc": utc_now(),
            "phase": reporter.phase,
            "files_complete": reporter.files_complete,
            "files_total": reporter.files_total,
            "bytes_complete": reporter.bytes_complete,
            "bytes_total": reporter.bytes_total,
            "current_task": reporter.current_task,
            "current_path": reporter.current_path,
            "inventory_sha256": inventory_sha256,
        },
    )


def verify_checksums(arguments: argparse.Namespace) -> None:
    if not arguments.checksums:
        raise ScanError("verify currently requires --checksums")
    reporter = OperationProgress("checksum_preflight")
    reporter.install()
    run_root = arguments.run_root.expanduser().resolve()
    plan, _, plan_sha256 = load_plan(run_root)
    package = run_root / "package"
    entries = promoted_inventory(run_root, plan, plan_sha256)
    inventory_sha256 = inventory_fingerprint(entries)
    stat_sha256 = inventory_fingerprint(entries, include_current_stat=True)
    audit_path = package / "verification" / "checksum_audit.json"
    if audit_path.is_file():
        audit = json_read(audit_path)
        if audit.get("plan_sha256") == plan_sha256 \
                and audit.get("inventory_sha256") == inventory_sha256 \
                and audit.get("current_stat_sha256") == stat_sha256 \
                and audit.get("state") == "complete":
            print(f"I: Reusing complete checksum audit: {audit_path}", file=sys.stderr)
            return
        raise ScanError(
            f"existing checksum audit no longer matches package state: {audit_path}"
        )

    checkpoint_directory = package / "verification" / "checkpoints"
    checkpoint_directory.mkdir(parents=True, exist_ok=True)
    checkpoints = read_verification_checkpoints(
        checkpoint_directory, plan_sha256, inventory_sha256
    )
    reporter.set_phase("checksum_audit")
    reporter.files_total = len(entries)
    reporter.bytes_total = sum(row["bytes"] for _, row, _ in entries)
    for task, row, path in entries:
        metadata = path.stat()
        key = (
            task["task_id"], row["output_relative_path"], row["sha256"],
            row["bytes"], metadata.st_mtime_ns,
        )
        if key in checkpoints:
            reporter.files_complete += 1
            reporter.bytes_complete += row["bytes"]

    segment = checkpoint_directory / (
        f"checksum-{utc_now().replace(':', '')}-pid-{os.getpid()}.jsonl"
    )
    newly_verified = 0
    print(
        f"I: Checksum audit started with {reporter.files_complete}/"
        f"{reporter.files_total} files already checkpointed.",
        file=sys.stderr,
        flush=True,
    )
    with segment.open("x", encoding="utf-8") as stream:
        for task, row, path in entries:
            metadata = path.stat()
            key = (
                task["task_id"], row["output_relative_path"], row["sha256"],
                row["bytes"], metadata.st_mtime_ns,
            )
            if key in checkpoints:
                continue
            if arguments.max_files is not None \
                    and newly_verified >= arguments.max_files:
                break
            reporter.current_task = task["task_id"]
            reporter.current_path = str(path)
            reporter.current_file_bytes = 0
            observed = sha256_file(
                path,
                reporter,
                arguments.max_read_mib_per_second,
            )
            if observed != row["sha256"]:
                write_verification_status(
                    package, reporter, "failed", inventory_sha256
                )
                raise ScanError(f"checksum mismatch: {path}")
            checkpoint = {
                "state": "verified",
                "verified_at_utc": utc_now(),
                "plan_sha256": plan_sha256,
                "inventory_sha256": inventory_sha256,
                "task_id": task["task_id"],
                "output_relative_path": row["output_relative_path"],
                "bytes": row["bytes"],
                "file_mtime_ns": metadata.st_mtime_ns,
                "sha256": observed,
            }
            stream.write(json.dumps(checkpoint, sort_keys=True, separators=(",", ":")))
            stream.write("\n")
            newly_verified += 1
            reporter.files_complete += 1
            reporter.bytes_complete += row["bytes"]
            reporter.current_file_bytes = 0
            if newly_verified % arguments.checkpoint_files == 0:
                stream.flush()
                os.fsync(stream.fileno())
                write_verification_status(
                    package, reporter, "running", inventory_sha256
                )
        stream.flush()
        os.fsync(stream.fileno())

    if reporter.files_complete != reporter.files_total:
        write_verification_status(package, reporter, "partial", inventory_sha256)
        print(
            f"I: Checksum audit checkpointed {reporter.files_complete}/"
            f"{reporter.files_total} files; rerun the same command to resume.",
            file=sys.stderr,
        )
        return

    reporter.set_phase("checksum_audit_complete", current_path="none", current_task="none")
    audit = {
        "schema_version": 1,
        "state": "complete",
        "verification_mode": "sha256_all_inventory_files",
        "completed_at_utc": utc_now(),
        "plan_sha256": plan_sha256,
        "inventory_sha256": inventory_sha256,
        "current_stat_sha256": stat_sha256,
        "file_count": reporter.files_total,
        "bytes": reporter.bytes_total,
    }
    json_write_new(audit_path, audit)
    write_verification_status(package, reporter, "complete", inventory_sha256)
    print(f"I: Completed checksum audit: {audit_path}", file=sys.stderr)


def status(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan, _, _ = load_plan(run_root)
    package = run_root / "package"
    complete = 0
    emitted = 0
    bytes_written = 0
    for task in plan["tasks"]:
        result_path = (
            package / "task_data" / f"task_id={task['task_id']}" / "task_result.json"
        )
        if result_path.is_file():
            result = json_read(result_path)
            if result.get("state") == "complete":
                complete += 1
                emitted += sum(row["emitted_hits"] for row in result["inventory"])
                bytes_written += sum(row["bytes"] for row in result["inventory"])
    failed_attempts = len(list((run_root / "staging").glob("*/job-*/failure.json")))
    print(f"run_id\t{plan['run_id']}")
    print(f"tasks_complete\t{complete}")
    print(f"tasks_total\t{len(plan['tasks'])}")
    print(f"failed_attempts_retained\t{failed_attempts}")
    print(f"emitted_hits\t{emitted}")
    print(f"parquet_bytes\t{bytes_written}")
    print(f"finalized\t{str((package / 'manifest.json').is_file()).lower()}")
    verification_path = verification_status_path(package)
    if verification_path.is_file():
        verification = json_read(verification_path)
        print(f"checksum_verification_state\t{verification.get('state', 'unknown')}")
        print(
            "checksum_files_complete\t"
            f"{verification.get('files_complete', 0)}"
        )
        print(
            "checksum_files_total\t"
            f"{verification.get('files_total', 0)}"
        )
        print(
            "checksum_bytes_complete\t"
            f"{verification.get('bytes_complete', 0)}"
        )
        print(
            "checksum_bytes_total\t"
            f"{verification.get('bytes_total', 0)}"
        )
    else:
        print("checksum_verification_state\tnot_started")


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Manage requeue-safe direct sparse-Parquet JASPAR genome scans."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser(
        "prepare", help="hash inputs and write immutable motif/chromosome task plan"
    )
    prepare_parser.add_argument("--run-root", required=True, type=Path)
    prepare_parser.add_argument("--run-id", required=True, type=identifier)
    prepare_parser.add_argument("--source-commit", required=True)
    prepare_parser.add_argument("--genome", required=True)
    prepare_parser.add_argument("--fasta-index", required=True)
    prepare_parser.add_argument("--jaspar", required=True)
    prepare_parser.add_argument("--genome-id", required=True, type=identifier)
    prepare_parser.add_argument("--motif-set-id", required=True, type=identifier)
    prepare_parser.add_argument("--taxon-id", required=True, type=positive_integer)
    prepare_parser.add_argument("--species", required=True)
    prepare_parser.add_argument("--assembly-name", required=True)
    prepare_parser.add_argument("--assembly-accession", required=True)
    prepare_parser.add_argument("--ensembl-release", required=True, type=positive_integer)
    prepare_parser.add_argument("--fasta-url", required=True)
    prepare_parser.add_argument("--jaspar-url", required=True)
    prepare_parser.add_argument("--jaspar-version", type=positive_integer, default=2026)
    prepare_parser.add_argument("--chrom", action="append", default=[])
    prepare_parser.add_argument("--motif-batch-size", type=positive_integer, default=128)
    prepare_parser.add_argument(
        "--motif-thresholds", type=Path,
        help="density-capped per-motif threshold TSV; requires its metadata JSON",
    )
    prepare_parser.add_argument(
        "--motif-threshold-metadata", type=Path,
        help="provenance JSON emitted beside --motif-thresholds",
    )
    prepare_parser.add_argument("--special-motif", default="MA0861.2")
    prepare_parser.add_argument(
        "--special-minimum-score", type=finite_float, default=-5.0
    )
    prepare_parser.add_argument(
        "--default-minimum-score", type=finite_float, default=-1.0
    )
    prepare_parser.add_argument(
        "--score-mode", choices=("log2_relative_risk", "log_odds"),
        default="log2_relative_risk",
    )
    prepare_parser.add_argument("--pseudocount", type=nonnegative_float, default=1.0)
    prepare_parser.add_argument("--minimum-free-gib", type=nonnegative_float, default=100.0)
    prepare_parser.set_defaults(function=prepare)

    task_parser = subparsers.add_parser(
        "run-task", help="run and atomically promote one planned Slurm array task"
    )
    task_parser.add_argument("--run-root", required=True, type=Path)
    task_selection = task_parser.add_mutually_exclusive_group(required=True)
    task_selection.add_argument("--task-index", type=int)
    task_selection.add_argument("--task-id")
    task_parser.add_argument("--scanner", required=True)
    task_parser.add_argument("--duckdb", default="duckdb")
    task_parser.add_argument(
        "--duckdb-memory-limit", type=duckdb_memory_limit, default="12GB",
        help="memory ceiling for the in-memory validation database (default: 12GB)",
    )
    task_parser.add_argument(
        "--duckdb-temp-directory", type=Path,
        help="optional job-local spill directory; promoted data never uses it",
    )
    task_parser.add_argument(
        "--task-work-directory", type=Path,
        help=(
            "optional scratch root for scanner output and validation; verified "
            "files are copied once into durable staging before promotion"
        ),
    )
    task_parser.add_argument(
        "--genome-override", type=Path,
        help="verified chromosome-only FASTA staged by run-chromosome",
    )
    task_parser.add_argument(
        "--fasta-index-override", type=Path,
        help="index paired with --genome-override",
    )
    task_parser.add_argument(
        "--staged-input-metadata", type=Path,
        help="hash and stat record written by the chromosome staging step",
    )
    task_parser.add_argument(
        "--jaspar-override", type=Path,
        help="verified JASPAR source staged by run-chromosome",
    )
    task_parser.add_argument(
        "--planned-inputs-verified", action="store_true", help=argparse.SUPPRESS,
    )
    task_parser.add_argument("--allow-local", action="store_true", help=argparse.SUPPRESS)
    task_parser.add_argument(
        "--allow-scanner-provenance-mismatch", action="store_true",
        help=(
            "explicitly permit a dirty build or commit mismatch; the mismatch "
            "is retained in task and run provenance"
        ),
    )
    task_parser.set_defaults(function=run_task)

    chromosome_parser = subparsers.add_parser(
        "run-chromosome",
        help="stage one chromosome in scratch and run all of its motif batches",
    )
    chromosome_parser.add_argument("--run-root", required=True, type=Path)
    chromosome_selection = chromosome_parser.add_mutually_exclusive_group(required=True)
    chromosome_selection.add_argument("--chrom")
    chromosome_selection.add_argument("--chrom-index", type=int)
    chromosome_parser.add_argument("--scanner", required=True)
    chromosome_parser.add_argument("--scratch-directory", required=True, type=Path)
    chromosome_parser.add_argument("--duckdb", default="duckdb")
    chromosome_parser.add_argument(
        "--duckdb-memory-limit", type=duckdb_memory_limit, default="12GB"
    )
    chromosome_parser.add_argument(
        "--minimum-scratch-free-gib", type=nonnegative_float, default=5.0,
        help="free-space reserve checked before staging and every batch (default: 5)",
    )
    chromosome_parser.add_argument(
        "--batch-workers", type=positive_integer, default=1,
        help="concurrent motif-batch scanners sharing the staged FASTA (default: 1)",
    )
    chromosome_parser.add_argument(
        "--scratch-task-output", action="store_true",
        help=(
            "also validate task Parquet in scratch before durable publication; "
            "local copies remain until the cluster clears the job directory"
        ),
    )
    chromosome_parser.add_argument(
        "--allow-local", action="store_true", help=argparse.SUPPRESS
    )
    chromosome_parser.add_argument(
        "--allow-scanner-provenance-mismatch", action="store_true",
        help="forward the explicit run-task provenance override",
    )
    chromosome_parser.set_defaults(function=run_chromosome)

    finalize_parser = subparsers.add_parser(
        "finalize",
        help="check exact task inventories and publish metadata plus DuckDB index",
    )
    finalize_parser.add_argument("--run-root", required=True, type=Path)
    finalize_parser.add_argument("--duckdb", default="duckdb")
    finalize_parser.set_defaults(function=finalize)

    verify_parser = subparsers.add_parser(
        "verify", help="run a resumable payload-integrity audit"
    )
    verify_parser.add_argument("--run-root", required=True, type=Path)
    verify_parser.add_argument(
        "--checksums", action="store_true",
        help="recompute SHA-256 for every exact inventory path",
    )
    verify_parser.add_argument(
        "--max-read-mib-per-second", type=positive_float,
        help="optional process-wide read-rate ceiling for the shared filesystem",
    )
    verify_parser.add_argument(
        "--checkpoint-files", type=positive_integer, default=100,
        help="flush durable resume state after this many new files (default: 100)",
    )
    verify_parser.add_argument(
        "--max-files", type=positive_integer,
        help="verify at most this many new files, then exit with resumable state",
    )
    verify_parser.set_defaults(function=verify_checksums)

    catalog_parser = subparsers.add_parser(
        "build-catalog",
        help="build a new metadata-only DuckDB catalog without opening hit Parquet",
    )
    catalog_parser.add_argument("--run-root", required=True, type=Path)
    catalog_parser.add_argument("--output", type=Path)
    catalog_parser.add_argument("--duckdb", default="duckdb")
    catalog_parser.set_defaults(function=rebuild_catalog)

    status_parser = subparsers.add_parser(
        "status", help="report completed tasks, retained failures, rows, and bytes"
    )
    status_parser.add_argument("--run-root", required=True, type=Path)
    status_parser.set_defaults(function=status)
    return parser


def main() -> int:
    parser = argument_parser()
    arguments = parser.parse_args()
    try:
        arguments.function(arguments)
        return 0
    except (ScanError, StageError, OSError, subprocess.SubprocessError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

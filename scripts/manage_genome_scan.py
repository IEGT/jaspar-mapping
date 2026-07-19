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
from pathlib import Path
from typing import Any, BinaryIO, Iterable


class ScanError(RuntimeError):
    pass


IDENTIFIER = re.compile(r"^[A-Za-z0-9._-]+$")
COMMIT = re.compile(r"^[0-9a-f]{40}$")
PLAN_SCHEMA_VERSION = 1


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
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


def json_read(path: Path) -> Any:
    try:
        with path.open(encoding="utf-8") as stream:
            return json.load(stream)
    except (OSError, json.JSONDecodeError) as error:
        raise ScanError(f"cannot read JSON {path}: {error}") from error


def check_free_space(path: Path, minimum_gib: float) -> None:
    free = shutil.disk_usage(path).free
    required = int(minimum_gib * 1024**3)
    print(
        f"I: Output filesystem has {free / 1024**3:.1f} GiB free; "
        f"required preflight reserve is {minimum_gib:.1f} GiB.",
        file=sys.stderr,
    )
    if free < required:
        raise ScanError(
            f"insufficient free space: {free / 1024**3:.1f} GiB available, "
            f"{minimum_gib:.1f} GiB required"
        )


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


def load_plan(run_root: Path) -> tuple[dict[str, Any], Path, str]:
    path = run_root / "plan" / "scan_plan.json"
    plan = json_read(path)
    if plan.get("schema_version") != PLAN_SCHEMA_VERSION:
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

    def add_batch(batch_id: str, policy_id: str, motif_ids: list[str]) -> None:
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
                "motif_list": str(relative),
                "motif_list_sha256": sha256_text(content),
            }
        )

    add_batch("special-0000", "tp73_calibrated", [arguments.special_motif])
    other_ids = sorted(set(motif_by_id) - {arguments.special_motif})
    for offset in range(0, len(other_ids), arguments.motif_batch_size):
        add_batch(
            f"default-{offset // arguments.motif_batch_size:04d}",
            "default_uncalibrated",
            other_ids[offset : offset + arguments.motif_batch_size],
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
    policy_by_id = {policy["policy_id"]: policy for policy in policies}
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
                    "minimum_score": policy_by_id[batch["policy_id"]]["minimum_score"],
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
            "source_sha256": sha256_file(jaspar),
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

    glob_path = package / "tables" / "jaspar2026" / "motif_hit" / "**" / "*.parquet"
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
    FROM read_parquet({sql_string(glob_path)}, filename=true,
                      hive_partitioning=true, union_by_name=true) h
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
        row["bytes"] = output.stat().st_size
        row["sha256"] = sha256_file(output)
        row["state"] = "complete"
    return records


def validate_promoted_task(
    task_root: Path,
    result: dict[str, Any],
    task: dict[str, Any],
    plan: dict[str, Any],
    plan_sha256: str,
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
        if output.stat().st_size != row.get("bytes") \
                or sha256_file(output) != row.get("sha256"):
            raise ScanError(f"promoted output checksum mismatch: {output}")
        if row.get("run_id") != plan["run_id"] \
                or row.get("task_id") != task["task_id"] \
                or row.get("source_commit") != plan["source_commit"] \
                or row.get("state") != "complete":
            raise ScanError(f"promoted output provenance mismatch: {output}")


def run_task(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan, plan_path, plan_sha256 = load_plan(run_root)
    check_free_space(run_root, plan["minimum_free_gib"])
    task = task_for_arguments(plan, arguments)
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

    job_id = os.environ.get("SLURM_JOB_ID", "local")
    array_id = os.environ.get("SLURM_ARRAY_TASK_ID", str(task["task_index"]))
    restart_count = os.environ.get("SLURM_RESTART_COUNT", "0")
    attempt = f"job-{job_id}-array-{array_id}-restart-{restart_count}-pid-{os.getpid()}"
    staging = run_root / "staging" / task["task_id"] / attempt
    staging.mkdir(parents=True, exist_ok=False)
    task_package = staging / "task_package"
    stats_path = staging / "scan_file_stats.jsonl"
    command = [
        str(scanner),
        "--sparse-parquet",
        "--outdir", str(task_package),
        "--genome", plan["genome"]["fasta_file"],
        "--fasta-index", plan["genome"]["fasta_index_file"],
        "--pssm", plan["motif_set"]["source_file"],
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

    def forward_progress(_signal_number: int, _frame: Any) -> None:
        if child is not None and child.poll() is None:
            child.send_signal(signal.SIGUSR1)

    if hasattr(signal, "SIGUSR1"):
        signal.signal(signal.SIGUSR1, forward_progress)
    try:
        print(f"I: Starting task {task['task_id']}", file=sys.stderr)
        child = subprocess.Popen(command, text=True)
        return_code = child.wait()
        if return_code != 0:
            raise ScanError(f"scanner exited with status {return_code}")
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
                    "scanner_sha256": sha256_file(scanner),
                    "slurm_job_id": job_id,
                    "slurm_array_task_id": array_id,
                    "slurm_restart_count": int(restart_count),
                }
            )
        result = {
            "state": "complete",
            "completed_at_utc": utc_now(),
            "plan_sha256": plan_sha256,
            "plan_file": str(plan_path),
            "source_commit": plan["source_commit"],
            "scanner_file": str(scanner),
            "scanner_sha256": sha256_file(scanner),
            "slurm_job_id": job_id,
            "slurm_array_task_id": array_id,
            "slurm_restart_count": int(restart_count),
            "task": task,
            "inventory": inventory,
        }
        json_write_new(task_package / "task_result.json", result)
        final_task.parent.mkdir(parents=True, exist_ok=True)
        if final_task.exists():
            raise ScanError(f"task appeared concurrently: {final_task}")
        os.replace(task_package, final_task)
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
        existing = existing_root / name / "**" / "*.parquet"
        staged = staged_root / name / "**" / "*.parquet"
        query = f"""
SELECT COUNT(*)::BIGINT AS differences
FROM (
    (SELECT * FROM read_parquet({sql_string(existing)}, union_by_name=true)
     EXCEPT ALL
     SELECT * FROM read_parquet({sql_string(staged)}, union_by_name=true))
    UNION ALL
    (SELECT * FROM read_parquet({sql_string(staged)}, union_by_name=true)
     EXCEPT ALL
     SELECT * FROM read_parquet({sql_string(existing)}, union_by_name=true))
);
"""
        rows = run_duckdb_json(duckdb, query)
        if len(rows) != 1 or rows[0].get("differences") != 0:
            raise ScanError(
                f"existing final metadata differs from staged {name}: {existing_root}"
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


def finalize(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    plan, _, plan_sha256 = load_plan(run_root)
    package = run_root / "package"
    manifest_path = package / "manifest.json"
    if manifest_path.exists():
        manifest = json_read(manifest_path)
        if manifest.get("plan_sha256") != plan_sha256:
            raise ScanError("existing manifest belongs to a different scan plan")
        print(f"I: Reusing finalized package: {package}", file=sys.stderr)
        return
    if shutil.which(arguments.duckdb) is None:
        raise ScanError(f"DuckDB executable not found: {arguments.duckdb}")

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
    if missing:
        raise ScanError(
            f"cannot finalize: {len(missing)} of {len(plan['tasks'])} tasks are missing"
        )

    inventory: list[dict[str, Any]] = []
    seen_files: set[tuple[str, str, str]] = set()
    scanner_sha256s = {result.get("scanner_sha256") for result in task_results}
    if None in scanner_sha256s or len(scanner_sha256s) != 1:
        raise ScanError("promoted tasks were produced by different scanner binaries")
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
    }
    datasets = {
        "motif_set": [motif_set],
        "genome": [genome],
        "sequence_region": sequence_regions,
        "motif_metadata": motif_metadata,
        "scan_run": [scan_run],
        "scan_threshold_policy": threshold_policies,
        "scan_task": scan_tasks,
        "scan_file_inventory": inventory,
    }
    for name, rows in datasets.items():
        source = staging / f"{name}.jsonl"
        destination = metadata_root / name / "part-000000.parquet"
        write_json_lines(source, rows)
        parquet_from_json(arguments.duckdb, source, destination)

    task_data_link = staging_package / "task_data"
    task_data_link.symlink_to((package / "task_data").resolve(), target_is_directory=True)
    schema = Path(__file__).resolve().parent.parent / "sql" / "genome_scan_schema.sql"
    database = staging_package / "jaspar_genome_scan.duckdb"
    process = subprocess.run(
        [arguments.duckdb, str(database), "-bail", "-f", str(schema)],
        cwd=staging_package,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ScanError(process.stderr.strip() or "could not build DuckDB query index")

    manifest = {
        "schema_version": 1,
        "run_id": plan["run_id"],
        "state": "complete",
        "completed_at_utc": completed_at,
        "plan_sha256": plan_sha256,
        "source_commit": plan["source_commit"],
        "genome_id": plan["genome"]["genome_id"],
        "motif_set_id": plan["motif_set"]["motif_set_id"],
        "task_count": len(scan_tasks),
        "file_count": len(inventory),
        "emitted_hit_count": scan_run["emitted_hit_count"],
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
    task_parser.add_argument("--allow-local", action="store_true", help=argparse.SUPPRESS)
    task_parser.set_defaults(function=run_task)

    finalize_parser = subparsers.add_parser(
        "finalize", help="verify every task and publish metadata plus DuckDB index"
    )
    finalize_parser.add_argument("--run-root", required=True, type=Path)
    finalize_parser.add_argument("--duckdb", default="duckdb")
    finalize_parser.set_defaults(function=finalize)

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
    except (ScanError, OSError, subprocess.SubprocessError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

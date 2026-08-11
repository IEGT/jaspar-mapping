#!/usr/bin/env python3

"""Build and publish one chromosome's TP73/control CUT&RUN evidence."""

from __future__ import annotations

import argparse
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


class ProductionError(RuntimeError):
    pass


STARTED = time.monotonic()
CURRENT_PHASE = "startup"
CURRENT_ATTEMPT: Path | None = None
CURRENT_SCRATCH: Path | None = None


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_bytes(path: Path | None) -> int:
    if path is None or not path.exists():
        return 0
    return sum(item.stat().st_size for item in path.rglob("*") if item.is_file())


def progress(_signal_number: int | None = None,
             _frame: object | None = None) -> None:
    print(
        "I: progress "
        f"phase={CURRENT_PHASE} elapsed_seconds={time.monotonic() - STARTED:.1f} "
        f"durable_bytes={tree_bytes(CURRENT_ATTEMPT)} "
        f"scratch_bytes={tree_bytes(CURRENT_SCRATCH)}",
        file=sys.stderr, flush=True,
    )


def set_phase(name: str) -> None:
    global CURRENT_PHASE
    CURRENT_PHASE = name
    progress()


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def run(command: list[str]) -> None:
    print("I: running: " + " ".join(command), file=sys.stderr, flush=True)
    process = subprocess.run(command, check=False)
    if process.returncode != 0:
        raise ProductionError(
            f"command failed with exit code {process.returncode}: {command[0]}"
        )


def query_json(duckdb: Path, database: Path | str,
               query: str) -> list[dict[str, object]]:
    database_arguments = [] if str(database) == ":memory:" else ["-readonly"]
    process = subprocess.run(
        [str(duckdb), "-light-mode", *database_arguments,
         "-json", str(database), "-c", query],
        text=True, capture_output=True, check=False,
    )
    if process.returncode != 0:
        raise ProductionError(process.stderr.strip() or "DuckDB query failed")
    value = json.loads(process.stdout or "[]")
    if not isinstance(value, list):
        raise ProductionError("DuckDB result is not a row array")
    return value


def absolute(path: Path, label: str, *, directory: bool = False) -> Path:
    resolved = path.expanduser().resolve()
    valid = resolved.is_dir() if directory else resolved.is_file()
    if not valid:
        kind = "directory" if directory else "file"
        raise ProductionError(f"{label} {kind} is missing: {resolved}")
    return resolved


def resolve_anchor(duckdb: Path, annotation_run: Path, chrom: str) -> Path:
    database = absolute(
        annotation_run / "final" / "context.duckdb", "annotation catalog"
    )
    rows = query_json(duckdb, database, f"""
SELECT absolute_path
FROM context_file_inventory
WHERE dataset = 'tp73_context_anchor'
  AND CAST(chrom AS VARCHAR) = {sql_string(chrom)}
  AND is_parquet
ORDER BY absolute_path;
""")
    if len(rows) != 1:
        raise ProductionError(
            f"expected one chromosome-{chrom} TP73 anchor file, found {len(rows)}"
        )
    return absolute(Path(str(rows[0]["absolute_path"])), "TP73 anchor source")


def check_free_space(path: Path, minimum_gb: float, label: str) -> None:
    free = shutil.disk_usage(path).free
    required = int(minimum_gb * 1024**3)
    print(
        f"I: {label} free space: {free / 1024**3:.1f} GiB "
        f"(required {minimum_gb:.1f} GiB)",
        file=sys.stderr, flush=True,
    )
    if free < required:
        raise ProductionError(f"insufficient {label} space")


def contract_sha256(contract: dict[str, object]) -> str:
    encoded = json.dumps(contract, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def build_contract(arguments: argparse.Namespace, anchor: Path) -> dict[str, object]:
    annotation_manifest = absolute(
        arguments.annotation_run / "final" / "manifest.json",
        "annotation manifest",
    )
    return {
        "format_version": 1,
        "analysis": "tp73_cutandrun_chromosome_evidence",
        "source": str(arguments.source),
        "source_commit": arguments.source_commit,
        "annotation_run": str(arguments.annotation_run),
        "annotation_manifest_sha256": sha256(annotation_manifest),
        "anchor_source": str(anchor),
        "anchor_source_sha256": sha256(anchor),
        "track_manifest": str(arguments.track_manifest),
        "track_manifest_sha256": sha256(arguments.track_manifest),
        "track_root": str(arguments.track_root),
        "chrom": arguments.chrom,
        "annotation_chrom": arguments.anchor_chrom,
        "chrom_length": arguments.chrom_length,
        "analysis_partition": arguments.analysis_partition,
        "primary_inference": arguments.primary_inference,
        "minimum_anchor_score": arguments.minimum_anchor_score,
        "included_series": ["saos2", "skmel29_2"],
        "excluded_series": ["skmel29_1"],
        "strict_immersion_rule": (
            "merged_component_start < anchor_start AND "
            "merged_component_end > anchor_end"
        ),
        "empty_chromosome_coverage_policy": "zero_support",
        "mitochondrial_alias_policy": (
            "M_MT_25_with_optional_chr_prefix"
            if arguments.analysis_partition == "mitochondrial_bystander_control"
            else "disabled"
        ),
        "mitochondrial_interpretation": (
            "bystander_control_due_to_high_mitochondrial_copy_number"
            if arguments.analysis_partition == "mitochondrial_bystander_control"
            else None
        ),
    }


def canonicalize_provenance(attempt: Path, final: Path,
                            scratch_anchor: Path, durable_anchor: Path) -> None:
    for path in sorted(item for item in attempt.rglob("*")
                       if item.is_file() and item.suffix in {".json", ".tsv"}):
        content = path.read_text(encoding="utf-8")
        updated = content.replace(str(attempt), str(final))
        updated = updated.replace(str(scratch_anchor), str(durable_anchor))
        if updated != content:
            temporary = path.with_name(f".{path.name}.canonicalizing")
            temporary.write_text(updated, encoding="utf-8")
            os.replace(temporary, path)


def validate_output(arguments: argparse.Namespace, anchor: Path,
                    evidence: Path) -> dict[str, int]:
    description = query_json(arguments.duckdb, ":memory:", f"""
DESCRIBE SELECT * FROM read_parquet({sql_string(evidence)});
""")
    columns = {str(row["column_name"]) for row in description}
    support = sorted(name for name in columns if name.startswith("supported_"))
    depth = sorted(name for name in columns if name.startswith("depth_"))
    if len(support) != 12 or len(depth) != 12:
        raise ProductionError(
            f"expected 12 TP73/control support-depth pairs, found "
            f"{len(support)} support and {len(depth)} depth columns"
        )
    if any("skmel29_1" in name for name in support + depth):
        raise ProductionError("excluded skmel29_1 columns entered the evidence")
    predicates = []
    for support_name in support:
        depth_name = support_name.replace("supported_", "depth_", 1)
        if depth_name not in columns:
            raise ProductionError(f"missing depth column for {support_name}")
        predicates.append(
            f"{support_name} IS NULL OR {depth_name} IS NULL OR "
            f"{depth_name} < 0 OR {support_name} <> ({depth_name} > 0)"
        )
    rows = query_json(arguments.duckdb, ":memory:", f"""
WITH source AS (
  SELECT start, "end"
FROM read_parquet({sql_string(anchor)}, hive_partitioning=false)
  WHERE motif_id = 'MA0861.2' AND anchor_selection_class = 'local_peak'
    AND score >= {arguments.minimum_anchor_score:.17g}
), evidence AS (
  SELECT * FROM read_parquet({sql_string(evidence)}, hive_partitioning=false)
)
SELECT
  (SELECT count(DISTINCT (start, "end")) FROM source)::BIGINT AS source_loci,
  (SELECT count(*) FROM evidence)::BIGINT AS evidence_rows,
  (SELECT count(*) - count(DISTINCT (chrom, anchor_start, anchor_end))
     FROM evidence)::BIGINT AS duplicate_keys,
  (SELECT count(*) FROM evidence WHERE anchor_score <
     {arguments.minimum_anchor_score:.17g})::BIGINT AS below_score_floor,
  (SELECT count(*) FROM evidence WHERE {' OR '.join(predicates)})::BIGINT
     AS support_depth_mismatches;
""")
    if len(rows) != 1:
        raise ProductionError("evidence validation returned no summary")
    values = {key: int(value) for key, value in rows[0].items()}
    if (values["source_loci"] != values["evidence_rows"]
            or values["evidence_rows"] == 0
            or values["duplicate_keys"] != 0
            or values["below_score_floor"] != 0
            or values["support_depth_mismatches"] != 0):
        raise ProductionError(f"TP73 evidence validation failed: {values}")
    return values


def output_inventory(attempt: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for path in sorted(item for item in attempt.rglob("*") if item.is_file()):
        if path.name == "manifest.json":
            continue
        rows.append({
            "relative_path": str(path.relative_to(attempt)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    return rows


def parser() -> argparse.ArgumentParser:
    source = Path(__file__).resolve().parent.parent
    result = argparse.ArgumentParser(
        description=(
            "Build a restart-safe one-chromosome TP73/control CUT&RUN evidence "
            "package from schema-7 local-peak anchors and manifest-pinned BigWigs."
        )
    )
    result.add_argument("--run-root", required=True, type=Path)
    result.add_argument("--annotation-run", required=True, type=Path)
    result.add_argument("--track-root", required=True, type=Path)
    result.add_argument("--source", type=Path, default=source)
    result.add_argument("--source-commit", required=True)
    result.add_argument("--track-manifest", type=Path)
    result.add_argument("--duckdb", required=True, type=Path)
    result.add_argument("--bigwig-python", required=True, type=Path)
    result.add_argument("--scratch-root", type=Path)
    result.add_argument("--chrom", required=True)
    result.add_argument(
        "--anchor-chrom",
        help="chromosome label used by the schema-7 annotation catalog",
    )
    result.add_argument("--chrom-length", required=True, type=int)
    result.add_argument("--analysis-partition", required=True)
    result.add_argument("--primary-inference", action="store_true")
    result.add_argument("--minimum-anchor-score", type=float, default=-1.0)
    result.add_argument("--threads", type=int, default=4)
    result.add_argument("--memory-limit", default="28GB")
    result.add_argument("--minimum-free-run-gb", type=float, default=5.0)
    result.add_argument("--minimum-free-scratch-gb", type=float, default=10.0)
    return result


def validate_arguments(arguments: argparse.Namespace) -> None:
    arguments.source = absolute(arguments.source, "source", directory=True)
    arguments.annotation_run = absolute(
        arguments.annotation_run, "annotation run", directory=True
    )
    arguments.track_root = absolute(arguments.track_root, "track root", directory=True)
    if arguments.track_manifest is None:
        arguments.track_manifest = (
            arguments.source / "config" / "h3k4me3_cutandrun_tracks_v1.tsv"
        )
    arguments.track_manifest = absolute(arguments.track_manifest, "track manifest")
    arguments.duckdb = absolute(arguments.duckdb, "DuckDB")
    arguments.bigwig_python = absolute(arguments.bigwig_python, "BigWig Python")
    arguments.run_root = arguments.run_root.expanduser().resolve()
    if not str(arguments.run_root).startswith("/data/sm718/"):
        raise ProductionError("--run-root must be below /data/sm718")
    if re.fullmatch(r"[A-Za-z0-9_.-]+", arguments.chrom) is None:
        raise ProductionError("--chrom contains unsafe characters")
    if arguments.anchor_chrom is None:
        arguments.anchor_chrom = arguments.chrom
    if re.fullmatch(r"[A-Za-z0-9_.-]+", arguments.anchor_chrom) is None:
        raise ProductionError("--anchor-chrom contains unsafe characters")
    if re.fullmatch(r"[A-Za-z0-9_.-]+", arguments.analysis_partition) is None:
        raise ProductionError("--analysis-partition contains unsafe characters")
    if re.fullmatch(r"[0-9a-f]{40}", arguments.source_commit) is None:
        raise ProductionError("--source-commit must be a full lowercase Git hash")
    if arguments.chrom_length <= 0 or arguments.threads <= 0:
        raise ProductionError("chromosome length and threads must be positive")
    if not (-1e100 < arguments.minimum_anchor_score < 1e100):
        raise ProductionError("minimum anchor score must be finite")


def main() -> int:
    global CURRENT_ATTEMPT, CURRENT_SCRATCH
    arguments = parser().parse_args()
    try:
        validate_arguments(arguments)
        signal.signal(signal.SIGUSR1, progress)
        arguments.run_root.mkdir(parents=True, exist_ok=True)
        (arguments.run_root / "attempts").mkdir(exist_ok=True)
        check_free_space(arguments.run_root, arguments.minimum_free_run_gb, "durable")

        set_phase("resolve-schema7-anchor")
        anchor = resolve_anchor(
            arguments.duckdb, arguments.annotation_run, arguments.anchor_chrom
        )
        contract = build_contract(arguments, anchor)
        fingerprint = contract_sha256(contract)
        final = arguments.run_root / "final"
        if final.is_dir():
            manifest = json.loads(
                (final / "manifest.json").read_text(encoding="utf-8")
            )
            if (manifest.get("state") == "complete"
                    and manifest.get("contract_sha256") == fingerprint):
                print(f"I: reusing completed chromosome package: {final}",
                      file=sys.stderr)
                return 0
            raise ProductionError("final package has a different contract")

        job_id = os.environ.get("SLURM_JOB_ID", "manual")
        restart = os.environ.get("SLURM_RESTART_COUNT", "0")
        attempt = arguments.run_root / "attempts" / f"job-{job_id}-restart-{restart}"
        if attempt.exists():
            attempt = attempt.with_name(f"{attempt.name}-pid-{os.getpid()}")
        attempt.mkdir()
        CURRENT_ATTEMPT = attempt

        scratch_root = arguments.scratch_root
        if scratch_root is None:
            scratch_root = Path(os.environ.get(
                "SLURM_TMPDIR", f"/scratch/{os.environ.get('USER', 'sm718')}"
            ))
        scratch_root = scratch_root.expanduser().resolve()
        scratch_root.mkdir(parents=True, exist_ok=True)
        check_free_space(
            scratch_root, arguments.minimum_free_scratch_gb, "scratch"
        )
        scratch = scratch_root / f"tp73-evidence-{job_id}-{restart}-{os.getpid()}"
        scratch.mkdir()
        CURRENT_SCRATCH = scratch
        scratch_anchor = scratch / f"tp73_context_anchor_chr{arguments.chrom}.parquet"

        set_phase("copy-anchor-to-scratch")
        shutil.copy2(anchor, scratch_anchor)
        evidence = attempt / "tp73_anchor_evidence.parquet"
        set_phase("build-strict-tp73-control-evidence")
        evidence_command = [
            str(arguments.bigwig_python),
            str(arguments.source / "scripts" / "build_h3k4me3_anchor_signal.py"),
            "--anchor-source", str(scratch_anchor),
            "--track-manifest", str(arguments.track_manifest),
            "--track-root", str(arguments.track_root),
            "--tp73-output", str(evidence),
            "--evidence-only",
            "--chrom", arguments.chrom,
            "--chrom-length", str(arguments.chrom_length),
            "--minimum-anchor-score", str(arguments.minimum_anchor_score),
            "--threads", str(arguments.threads),
            "--memory-limit", arguments.memory_limit,
            "--duckdb", str(arguments.duckdb),
            "--rscript", str(arguments.bigwig_python),
            "--bigwig-exporter", str(
                arguments.source / "scripts" / "export_bigwig_chrom_bedgraph.py"
            ),
            "--anchor-builder", str(
                arguments.source / "scripts" / "build_tp73_anchor_evidence.py"
            ),
            "--tmpdir", str(scratch),
        ]
        if arguments.analysis_partition == "mitochondrial_bystander_control":
            evidence_command.append("--mitochondrial-chromosome")
        run(evidence_command)

        set_phase("canonicalize-provenance")
        canonicalize_provenance(attempt, final, scratch_anchor, anchor)
        set_phase("validate-and-publish")
        validation = validate_output(arguments, anchor, evidence)
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
        print(f"I: published chromosome evidence: {final}", file=sys.stderr)
        return 0
    except (OSError, ProductionError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        progress()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

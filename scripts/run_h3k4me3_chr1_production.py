#!/usr/bin/env python3

"""Build and publish the chromosome-1 TP73/H3K4me3 production analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
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


def progress(_signal_number: int | None = None, _frame: object | None = None) -> None:
    print(
        "I: progress "
        f"phase={CURRENT_PHASE} elapsed_seconds={time.monotonic() - STARTED:.1f} "
        f"durable_bytes={tree_bytes(CURRENT_ATTEMPT)} "
        f"scratch_bytes={tree_bytes(CURRENT_SCRATCH)}",
        file=sys.stderr,
        flush=True,
    )


def set_phase(name: str) -> None:
    global CURRENT_PHASE
    CURRENT_PHASE = name
    progress()


def run(command: list[str]) -> None:
    print("I: running: " + " ".join(command), file=sys.stderr, flush=True)
    process = subprocess.run(command, check=False)
    if process.returncode != 0:
        raise ProductionError(
            f"command failed with exit code {process.returncode}: {command[0]}"
        )


def query_json(duckdb: Path, database: Path | str, query: str) -> list[dict[str, object]]:
    database_arguments = [] if str(database) == ":memory:" else ["-readonly"]
    process = subprocess.run(
        [str(duckdb), "-light-mode", *database_arguments,
         "-json", str(database), "-c", query],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ProductionError(process.stderr.strip() or "DuckDB query failed")
    try:
        value = json.loads(process.stdout or "[]")
    except json.JSONDecodeError as error:
        raise ProductionError("DuckDB returned invalid JSON") from error
    if not isinstance(value, list):
        raise ProductionError("DuckDB result is not an array")
    return value


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def absolute(path: Path, label: str, *, directory: bool = False) -> Path:
    resolved = path.expanduser().resolve()
    valid = resolved.is_dir() if directory else resolved.is_file()
    if not valid:
        kind = "directory" if directory else "file"
        raise ProductionError(f"{label} {kind} is missing: {resolved}")
    return resolved


def git_value(source: Path, *arguments: str) -> str:
    process = subprocess.run(
        ["git", "-C", str(source), *arguments],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ProductionError(process.stderr.strip() or "git query failed")
    return process.stdout.strip()


def git_success(source: Path, *arguments: str) -> bool:
    return subprocess.run(
        ["git", "-C", str(source), *arguments],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    ).returncode == 0


def read_panel(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"motif_id", "positive_threshold"}
        if not required.issubset(reader.fieldnames or ()):
            raise ProductionError("cofactor panel lacks motif_id/positive_threshold")
        rows = list(reader)
    motifs = [row["motif_id"].strip() for row in rows]
    if not motifs or any(not motif for motif in motifs) or len(set(motifs)) != len(motifs):
        raise ProductionError("cofactor panel motif IDs are empty or duplicated")
    for row in rows:
        try:
            float(row["positive_threshold"])
        except ValueError as error:
            raise ProductionError("cofactor panel contains an invalid threshold") from error
    return motifs


def resolve_anchor(duckdb: Path, annotation_run: Path) -> Path:
    database = absolute(annotation_run / "final" / "context.duckdb", "annotation catalog")
    rows = query_json(
        duckdb,
        database,
        """
SELECT absolute_path
FROM context_file_inventory
WHERE dataset = 'tp73_context_anchor'
  AND CAST(chrom AS VARCHAR) = '1'
  AND is_parquet
ORDER BY absolute_path;
""",
    )
    if len(rows) != 1:
        raise ProductionError(
            f"expected one chromosome-1 tp73_context_anchor file, found {len(rows)}"
        )
    return absolute(Path(str(rows[0]["absolute_path"])), "TP73 anchor source")


def check_free_space(path: Path, minimum_gb: float, label: str) -> None:
    free = shutil.disk_usage(path).free
    required = int(minimum_gb * 1024**3)
    print(
        f"I: {label} free space: {free / 1024**3:.1f} GiB "
        f"(required {minimum_gb:.1f} GiB)",
        file=sys.stderr,
    )
    if free < required:
        raise ProductionError(f"insufficient {label} space")


def stage_scan_inputs(arguments: argparse.Namespace, motifs: list[str]) -> Path:
    output = arguments.run_root / "input" / "chr1_cofactor_panel"
    command = [
        sys.executable,
        str(arguments.source / "scripts" / "stage_motif_context_inputs.py"),
        "--package", str(arguments.scan_package),
        "--output", str(output),
        "--chrom", "1",
        "--duckdb", str(arguments.duckdb),
    ]
    for motif in motifs:
        command.extend(["--motif", motif])
    run(command)
    return output


def copy_scan_inputs(stage: Path, scratch: Path) -> dict[str, dict[str, Path]]:
    manifest_path = stage / "input_manifest.json"
    content = json.loads(manifest_path.read_text(encoding="utf-8"))
    result: dict[str, dict[str, Path]] = {}
    for row in content.get("files", []):
        motif = str(row["motif_id"])
        strand = str(row["strand"])
        if strand not in {"+", "-"}:
            raise ProductionError(f"invalid staged strand for {motif}: {strand}")
        source = absolute(stage / str(row["staged_relative_path"]), "staged motif")
        label = "plus" if strand == "+" else "minus"
        target = scratch / "motifs" / motif / f"{label}.parquet"
        target.parent.mkdir(parents=True, exist_ok=True)
        print(f"I: copying {motif} {label} to scratch", file=sys.stderr, flush=True)
        shutil.copy2(source, target)
        result.setdefault(motif, {})[label] = target
    incomplete = [motif for motif, paths in result.items() if set(paths) != {"plus", "minus"}]
    if incomplete:
        raise ProductionError(f"staged motif inputs lack both strands: {incomplete}")
    return result


def build_contract(arguments: argparse.Namespace, anchor: Path, stage: Path,
                   motifs: list[str], source_commit: str) -> dict[str, object]:
    annotation_manifest = absolute(
        arguments.annotation_run / "final" / "manifest.json", "annotation manifest"
    )
    scan_manifest = absolute(arguments.scan_package / "manifest.json", "scan manifest")
    stage_manifest = absolute(stage / "input_manifest.json", "staging manifest")
    return {
        "format_version": 1,
        "analysis": "tp73_h3k4me3_change_chr1",
        "source": str(arguments.source),
        "source_commit": source_commit,
        "annotation_run": str(arguments.annotation_run),
        "annotation_manifest_sha256": sha256(annotation_manifest),
        "anchor_source": str(anchor),
        "anchor_source_sha256": sha256(anchor),
        "scan_package": str(arguments.scan_package),
        "scan_manifest_sha256": sha256(scan_manifest),
        "staged_input_manifest_sha256": sha256(stage_manifest),
        "track_manifest": str(arguments.track_manifest),
        "track_manifest_sha256": sha256(arguments.track_manifest),
        "track_root": str(arguments.track_root),
        "cofactor_thresholds": str(arguments.thresholds),
        "cofactor_thresholds_sha256": sha256(arguments.thresholds),
        "cofactor_motifs": motifs,
        "chrom": "1",
        "chrom_length": arguments.chrom_length,
        "minimum_anchor_score": arguments.minimum_anchor_score,
        "source_score_floor": arguments.source_score_floor,
        "context_flank_bp": arguments.context_flank,
        "primary_h3k4me3_window": arguments.primary_window,
        "included_series": ["saos2", "skmel29_2"],
        "excluded_series": ["skmel29_1"],
        "negative_reference_thresholds": [-1, 0],
        "integrated_signal_pseudocount": arguments.pseudocount,
        "genomic_block_size_bp": arguments.block_size,
    }


def contract_sha256(contract: dict[str, object]) -> str:
    encoded = json.dumps(contract, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def validate_outputs(arguments: argparse.Namespace, attempt: Path,
                     anchor: Path, motif_count: int) -> dict[str, int]:
    evidence = attempt / "tp73_anchor_evidence.parquet"
    signal_path = attempt / "h3k4me3_anchor_signal.parquet"
    maxima = attempt / "cofactor_maxima.parquet"
    for path in (evidence, signal_path, maxima):
        absolute(path, "production output")
    query = f"""
WITH source AS (
    SELECT CAST(chrom AS VARCHAR) AS chrom, start, \"end\"
    FROM read_parquet({sql_string(anchor)})
), evidence AS (
    SELECT * FROM read_parquet({sql_string(evidence)})
), signal AS (
    SELECT * FROM read_parquet({sql_string(signal_path)})
), maxima AS (
    SELECT * FROM read_parquet({sql_string(maxima)})
)
SELECT
  (SELECT count(DISTINCT (chrom, start, \"end\")) FROM source) AS source_loci,
  (SELECT count(*) FROM evidence) AS evidence_rows,
  (SELECT count(DISTINCT (chrom, anchor_start, anchor_end)) FROM evidence) AS evidence_loci,
  (SELECT count(*) FROM signal) AS signal_rows,
  (SELECT count(DISTINCT window_name) FROM signal) AS windows,
  (SELECT count(DISTINCT series_id) FROM signal) AS series,
  (SELECT count(DISTINCT condition) FROM signal) AS conditions,
  (SELECT count(*) FROM maxima) AS maxima_rows,
  (SELECT count(DISTINCT motif_id) FROM maxima) AS maxima_motifs,
  (SELECT count(*) FROM evidence WHERE anchor_score < {arguments.minimum_anchor_score})
      AS below_anchor_floor,
  (SELECT count(*) FROM signal WHERE series_id = 'skmel29_1') AS excluded_signal_rows;
"""
    rows = query_json(arguments.duckdb, ":memory:", query)
    if len(rows) != 1:
        raise ProductionError("production validation query returned no summary")
    values = {key: int(value) for key, value in rows[0].items()}
    expected_signal = values["evidence_rows"] * values["windows"] * 2 * 3
    expected_maxima = values["evidence_rows"] * motif_count
    if values["source_loci"] != values["evidence_rows"] or \
            values["evidence_rows"] != values["evidence_loci"]:
        raise ProductionError("anchor evidence is not one row per selected physical span")
    if values["series"] != 2 or values["conditions"] != 3 or \
            values["signal_rows"] != expected_signal:
        raise ProductionError("H3K4me3 signal factorial/cardinality is incomplete")
    if values["maxima_motifs"] != motif_count or \
            values["maxima_rows"] != expected_maxima:
        raise ProductionError("cofactor maxima cardinality is incomplete")
    if values["below_anchor_floor"] != 0 or values["excluded_signal_rows"] != 0:
        raise ProductionError("anchor floor or excluded-series invariant failed")
    prefix = attempt / "analysis" / "h3k4me3_chr1"
    for suffix in (
        "_intensity_effect.tsv", "_tp73_interaction.tsv",
        "_series_summary.tsv", "_binding_state_summary.tsv",
        "_occurrence_summary.tsv", "_run_config.tsv",
    ):
        path = Path(str(prefix) + suffix)
        if not path.is_file() or path.stat().st_size == 0:
            raise ProductionError(f"analysis output is absent or empty: {path}")
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


def argument_parser() -> argparse.ArgumentParser:
    source = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(
        description=(
            "Build a restart-safe chromosome-1 TP73/CUT&RUN/H3K4me3 package "
            "from schema-7 local-peak anchors and an exact finalized scan inventory."
        )
    )
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--annotation-run", required=True, type=Path)
    parser.add_argument("--scan-package", required=True, type=Path)
    parser.add_argument("--track-root", required=True, type=Path)
    parser.add_argument("--source", type=Path, default=source)
    parser.add_argument("--track-manifest", type=Path)
    parser.add_argument("--thresholds", type=Path)
    parser.add_argument("--duckdb", required=True, type=Path)
    parser.add_argument("--rscript", required=True, type=Path)
    parser.add_argument("--bigwig-python", required=True, type=Path)
    parser.add_argument("--scratch-root", type=Path)
    parser.add_argument("--chrom-length", type=int, default=248956422)
    parser.add_argument("--minimum-anchor-score", type=float, default=-1.0)
    parser.add_argument("--source-score-floor", type=float, default=-1.0)
    parser.add_argument("--context-flank", type=int, default=150)
    parser.add_argument("--primary-window", default="flank_150_1000")
    parser.add_argument("--pseudocount", type=float, default=1.0)
    parser.add_argument("--block-size", type=int, default=5000000)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--memory-limit", default="28GB")
    parser.add_argument("--max-temp-size", default="100GB")
    parser.add_argument("--minimum-free-run-gb", type=float, default=10.0)
    parser.add_argument("--minimum-free-scratch-gb", type=float, default=30.0)
    return parser


def validate_arguments(arguments: argparse.Namespace) -> None:
    arguments.source = absolute(arguments.source, "source", directory=True)
    if arguments.track_manifest is None:
        arguments.track_manifest = (
            arguments.source / "config" / "h3k4me3_cutandrun_tracks_v1.tsv"
        )
    if arguments.thresholds is None:
        arguments.thresholds = (
            arguments.source / "config" / "h3k4me3_chr1_pilot_cofactors_v1.tsv"
        )
    arguments.annotation_run = absolute(
        arguments.annotation_run, "annotation run", directory=True
    )
    arguments.scan_package = absolute(
        arguments.scan_package, "scan package", directory=True
    )
    arguments.track_root = absolute(arguments.track_root, "track root", directory=True)
    arguments.track_manifest = absolute(arguments.track_manifest, "track manifest")
    arguments.thresholds = absolute(arguments.thresholds, "threshold table")
    arguments.duckdb = absolute(arguments.duckdb, "DuckDB")
    arguments.rscript = absolute(arguments.rscript, "Rscript")
    arguments.bigwig_python = absolute(arguments.bigwig_python, "BigWig Python")
    arguments.run_root = arguments.run_root.expanduser().resolve()
    if not str(arguments.run_root).startswith("/data/sm718/"):
        raise ProductionError("--run-root must be below /data/sm718")
    if arguments.chrom_length <= 0 or arguments.context_flank < 0 or arguments.threads <= 0:
        raise ProductionError("chromosome length/flank/threads are invalid")
    if arguments.minimum_free_run_gb < 0 or arguments.minimum_free_scratch_gb < 0:
        raise ProductionError("minimum free-space limits cannot be negative")


def main() -> int:
    global CURRENT_ATTEMPT, CURRENT_SCRATCH
    arguments = argument_parser().parse_args()
    try:
        validate_arguments(arguments)
        signal.signal(signal.SIGUSR1, progress)
        source_commit = git_value(arguments.source, "rev-parse", "HEAD")
        if not git_success(
            arguments.source, "diff", "--quiet", "--ignore-submodules", "--"
        ) or not git_success(
            arguments.source, "diff", "--cached", "--quiet",
            "--ignore-submodules", "--"
        ):
            raise ProductionError("production execution requires a tracked-clean source tree")
        motifs = read_panel(arguments.thresholds)
        arguments.run_root.mkdir(parents=True, exist_ok=True)
        (arguments.run_root / "input").mkdir(exist_ok=True)
        (arguments.run_root / "attempts").mkdir(exist_ok=True)
        check_free_space(arguments.run_root, arguments.minimum_free_run_gb, "durable")

        set_phase("resolve-schema7-anchor")
        anchor = resolve_anchor(arguments.duckdb, arguments.annotation_run)
        stage = stage_scan_inputs(arguments, motifs)
        contract = build_contract(arguments, anchor, stage, motifs, source_commit)
        fingerprint = contract_sha256(contract)
        final = arguments.run_root / "final"
        if final.is_dir():
            manifest_path = final / "manifest.json"
            if not manifest_path.is_file():
                raise ProductionError("final directory exists without a manifest")
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if manifest.get("state") == "complete" and \
                    manifest.get("contract_sha256") == fingerprint:
                print(f"I: reusing completed production package: {final}", file=sys.stderr)
                return 0
            raise ProductionError("final package exists with a different/incomplete contract")

        job_id = os.environ.get("SLURM_JOB_ID", "manual")
        restart = os.environ.get("SLURM_RESTART_COUNT", "0")
        attempt = arguments.run_root / "attempts" / f"job-{job_id}-restart-{restart}"
        if attempt.exists():
            attempt = arguments.run_root / "attempts" / (
                f"job-{job_id}-restart-{restart}-pid-{os.getpid()}"
            )
        attempt.mkdir()
        CURRENT_ATTEMPT = attempt

        scratch_root = arguments.scratch_root
        if scratch_root is None:
            env_scratch = os.environ.get("SLURM_TMPDIR")
            scratch_root = (
                Path(env_scratch) if env_scratch
                else Path("/scratch") / os.environ.get("USER", "sm718")
            )
        scratch_root = scratch_root.expanduser().resolve()
        scratch_root.mkdir(parents=True, exist_ok=True)
        check_free_space(scratch_root, arguments.minimum_free_scratch_gb, "scratch")
        scratch = scratch_root / f"jaspar-h3k4me3-{job_id}-{restart}-{os.getpid()}"
        scratch.mkdir()
        CURRENT_SCRATCH = scratch

        set_phase("copy-inputs-to-scratch")
        scratch_anchor = scratch / "tp73_context_anchor_chr1.parquet"
        shutil.copy2(anchor, scratch_anchor)
        cofactor_paths = copy_scan_inputs(stage, scratch)

        evidence = attempt / "tp73_anchor_evidence.parquet"
        signal_output = attempt / "h3k4me3_anchor_signal.parquet"
        profile = attempt / "gfp_metaprofile.tsv"
        set_phase("build-cutandrun-and-h3k4me3-evidence")
        run([
            str(arguments.bigwig_python),
            str(arguments.source / "scripts" / "build_h3k4me3_anchor_signal.py"),
            "--anchor-source", str(scratch_anchor),
            "--track-manifest", str(arguments.track_manifest),
            "--track-root", str(arguments.track_root),
            "--tp73-output", str(evidence),
            "--signal-output", str(signal_output),
            "--profile-output", str(profile),
            "--chrom", "1", "--chrom-length", str(arguments.chrom_length),
            "--minimum-anchor-score", str(arguments.minimum_anchor_score),
            "--profile-sample-size", "100000",
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
        ])

        maxima = attempt / "cofactor_maxima.parquet"
        set_phase("build-nine-cofactor-maxima")
        maxima_command = [
            str(arguments.bigwig_python),
            str(arguments.source / "scripts" / "build_sparse_context_maxima.py"),
            "--anchor-parquet", str(evidence),
            "--output", str(maxima),
            "--flank", str(arguments.context_flank),
            "--source-score-floor", str(arguments.source_score_floor),
            "--duckdb", str(arguments.duckdb),
            "--threads", str(arguments.threads),
            "--memory-limit", arguments.memory_limit,
            "--max-temp-size", arguments.max_temp_size,
            "--temp-directory", str(scratch),
        ]
        for motif in motifs:
            maxima_command.extend([
                "--cofactor", motif,
                str(cofactor_paths[motif]["plus"]),
                str(cofactor_paths[motif]["minus"]),
            ])
        run(maxima_command)

        analysis = attempt / "analysis"
        analysis.mkdir()
        prefix = analysis / "h3k4me3_chr1"
        set_phase("fit-gfp-referenced-h3k4me3-change")
        run([
            str(arguments.rscript),
            str(arguments.source / "scripts" / "analyze_h3k4me3_cofactor_change.R"),
            "--signal", str(signal_output),
            "--tp73-evidence", str(evidence),
            "--cofactor-maxima", str(maxima),
            "--thresholds", str(arguments.thresholds),
            "--window", arguments.primary_window,
            "--output-prefix", str(prefix),
            "--series", "saos2", "--series", "skmel29_2",
            "--negative-references", "-1,0",
            "--pseudocount", str(arguments.pseudocount),
            "--block-size", str(arguments.block_size),
            "--duckdb", str(arguments.duckdb),
        ])

        set_phase("validate-and-publish")
        validation = validate_outputs(arguments, attempt, anchor, len(motifs))
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
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(attempt, final)
        CURRENT_ATTEMPT = final
        set_phase("complete")
        print(f"I: published production package: {final}", file=sys.stderr)
        return 0
    except (OSError, ProductionError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        progress()
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

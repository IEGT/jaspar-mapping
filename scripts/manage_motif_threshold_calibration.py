#!/usr/bin/env python3

"""Prepare, inspect, and finalize an all-motif threshold calibration run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


class CalibrationError(RuntimeError):
    pass


JASPAR_HEADER = re.compile(r"^>(\S+)\s+(.+?)\s*$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def run_json(
    command: list[str], *, working_directory: Path | None = None
) -> list[dict[str, object]]:
    process = subprocess.run(
        command,
        cwd=working_directory,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise CalibrationError(process.stderr.strip() or "command failed")
    try:
        value = json.loads(process.stdout or "[]")
    except json.JSONDecodeError as error:
        raise CalibrationError("command returned invalid JSON") from error
    if not isinstance(value, list):
        raise CalibrationError("command did not return a JSON row array")
    return value


def package_database(package: Path) -> tuple[Path, Path]:
    manifest = package / "manifest.json"
    if not manifest.is_file():
        raise CalibrationError(f"scan manifest not found: {manifest}")
    try:
        data = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CalibrationError(f"cannot read scan manifest: {error}") from error
    database_name = data.get("database")
    if not database_name:
        raise CalibrationError("scan manifest does not name its DuckDB catalog")
    database = Path(str(database_name))
    if not database.is_absolute():
        database = package / database
    if not database.is_file():
        raise CalibrationError(f"scan catalog not found: {database}")
    return manifest, database.resolve()


def parse_jaspar(path: Path) -> list[tuple[str, str]]:
    motifs: list[tuple[str, str]] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            match = JASPAR_HEADER.match(line)
            if not match:
                continue
            motif_id, motif_name = match.groups()
            motif_name = " ".join(motif_name.split())
            if motif_id in seen:
                raise CalibrationError(f"duplicate JASPAR motif: {motif_id}")
            seen.add(motif_id)
            motifs.append((motif_id, motif_name))
    if not motifs:
        raise CalibrationError(f"no motifs found in JASPAR source: {path}")
    return motifs


def inventory(duckdb: str, database: Path, chrom: str) -> list[dict[str, object]]:
    query = f"""
SELECT motif_id, strand, task_id, output_relative_path, bytes, sha256,
       emitted_hits
FROM scan_file_inventory
WHERE CAST(chrom AS VARCHAR) = {sql_string(chrom)}
ORDER BY motif_id, strand;
"""
    # Finalized scan catalogs intentionally use package-relative Parquet views.
    # Resolve those views from the catalog directory, independent of the
    # caller's current working directory.
    return run_json(
        [duckdb, "-readonly", "-json", str(database), "-c", query],
        working_directory=database.parent,
    )


def inventory_path(row: dict[str, object]) -> Path:
    relative = Path("task_data") / f"task_id={row['task_id']}" / str(
        row["output_relative_path"]
    )
    if relative.is_absolute() or ".." in relative.parts:
        raise CalibrationError(f"unsafe scan inventory path: {relative}")
    return relative


def immutable_write(path: Path, content: str) -> None:
    encoded = content.encode("utf-8")
    if path.exists():
        if path.read_bytes() != encoded:
            raise CalibrationError(f"existing immutable file differs: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        temporary.write_bytes(encoded)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def git_identity(source: Path) -> tuple[str, bool]:
    commit = subprocess.run(
        ["git", "-C", str(source), "rev-parse", "HEAD"],
        text=True, capture_output=True, check=False,
    )
    if commit.returncode != 0:
        raise CalibrationError(commit.stderr.strip() or "cannot read Git commit")
    dirty = subprocess.run(
        ["git", "-C", str(source), "diff", "--quiet", "--ignore-submodules", "--"],
        check=False,
    ).returncode != 0
    dirty = dirty or subprocess.run(
        ["git", "-C", str(source), "diff", "--cached", "--quiet",
         "--ignore-submodules", "--"],
        check=False,
    ).returncode != 0
    return commit.stdout.strip(), dirty


def prepare(arguments: argparse.Namespace) -> None:
    duckdb = shutil.which(arguments.duckdb)
    if duckdb is None:
        raise CalibrationError(f"DuckDB executable not found: {arguments.duckdb}")
    run_root = arguments.run_root.expanduser().resolve()
    package = arguments.scan_package.expanduser().resolve()
    jaspar = arguments.jaspar.expanduser().resolve()
    source = arguments.source.expanduser().resolve()
    if not package.is_dir() or not jaspar.is_file() or not source.is_dir():
        raise CalibrationError("scan package, JASPAR source, or repository is missing")
    run_root.mkdir(parents=True, exist_ok=True)
    for name in ("plan", "logs", "tasks", "final", "input"):
        (run_root / name).mkdir(exist_ok=True)

    manifest, database = package_database(package)
    motifs = parse_jaspar(jaspar)
    rows = inventory(duckdb, database, arguments.chrom)
    by_motif: dict[str, dict[str, dict[str, object]]] = {}
    for row in rows:
        motif_id = str(row["motif_id"])
        strand = str(row["strand"])
        if strand not in {"+", "-"}:
            raise CalibrationError(f"unexpected strand for {motif_id}: {strand}")
        if strand in by_motif.setdefault(motif_id, {}):
            raise CalibrationError(f"duplicate {motif_id} {strand} inventory row")
        relative = inventory_path(row)
        source_file = package / relative
        if not source_file.is_file() or source_file.stat().st_size != int(row["bytes"]):
            raise CalibrationError(f"missing or size-mismatched scan file: {source_file}")
        row = dict(row)
        row["relative_path"] = str(relative)
        by_motif[motif_id][strand] = row

    jaspar_ids = {motif_id for motif_id, _ in motifs}
    if set(by_motif) != jaspar_ids:
        missing = sorted(jaspar_ids - set(by_motif))
        extra = sorted(set(by_motif) - jaspar_ids)
        raise CalibrationError(
            f"JASPAR/inventory motif mismatch: missing={missing[:5]}, extra={extra[:5]}"
        )
    incomplete = sorted(
        motif_id for motif_id, strands in by_motif.items()
        if set(strands) != {"+", "-"}
    )
    if incomplete:
        raise CalibrationError(f"motifs lack both strands: {incomplete[:5]}")
    if arguments.target_motif not in jaspar_ids:
        raise CalibrationError(f"target motif is absent: {arguments.target_motif}")

    task_motifs = [item for item in motifs if item[0] != arguments.target_motif]
    header = [
        "task_index", "motif_id", "motif_name", "plus_relative_path",
        "minus_relative_path", "plus_bytes", "minus_bytes", "plus_sha256",
        "minus_sha256", "plus_emitted_hits", "minus_emitted_hits",
    ]
    task_lines = ["\t".join(header)]
    for task_index, (motif_id, motif_name) in enumerate(task_motifs):
        plus = by_motif[motif_id]["+"]
        minus = by_motif[motif_id]["-"]
        task_lines.append("\t".join(map(str, [
            task_index, motif_id, motif_name,
            plus["relative_path"], minus["relative_path"],
            plus["bytes"], minus["bytes"], plus["sha256"], minus["sha256"],
            plus["emitted_hits"], minus["emitted_hits"],
        ])))
    immutable_write(run_root / "plan" / "calibration_tasks.tsv",
                    "\n".join(task_lines) + "\n")
    immutable_write(
        run_root / "plan" / "motifs.txt",
        "".join(f"{motif_id}\n" for motif_id, _ in task_motifs),
    )

    target_lines = ["strand\trelative_path\tbytes\tsha256\temitted_hits"]
    for strand in ("+", "-"):
        row = by_motif[arguments.target_motif][strand]
        target_lines.append("\t".join(map(str, [
            strand, row["relative_path"], row["bytes"], row["sha256"],
            row["emitted_hits"],
        ])))
    immutable_write(run_root / "plan" / "target_anchor_files.tsv",
                    "\n".join(target_lines) + "\n")

    commit, dirty = git_identity(source)
    config = {
        "schema_version": 1,
        "run_id": arguments.run_id,
        "threshold_set_id": arguments.threshold_set_id,
        "chrom": arguments.chrom,
        "chrom_size": arguments.chrom_size,
        "scan_package": str(package),
        "scan_manifest_sha256": sha256(manifest),
        "scan_database": str(database),
        "jaspar": str(jaspar),
        "jaspar_sha256": sha256(jaspar),
        "source": str(source),
        "source_commit": commit,
        "source_dirty": dirty,
        "target_motif_id": arguments.target_motif,
        "target_policy": "separate_direct_and_tandem_calibration",
        "cofactor_task_count": len(task_motifs),
        "anchor_evidence": str(arguments.anchor_evidence.expanduser().resolve()),
        "source_minimum_score": -1,
        "context_flank_bp": 150,
        "candidate_grid": "observed_nonnegative_integer_thresholds",
        "minimum_class_fraction": arguments.minimum_class_fraction,
        "folds": 5,
        "fold_definition": "five_equal_width_contiguous_grch38_chr1_spans",
        "selection_metric": "delta_macro_roc_auc",
    }
    immutable_write(
        run_root / "plan" / "run_config.json",
        json.dumps(config, indent=2, sort_keys=True) + "\n",
    )
    print(len(task_motifs))


def load_config(run_root: Path) -> dict[str, object]:
    path = run_root / "plan" / "run_config.json"
    if not path.is_file():
        raise CalibrationError(f"run configuration not found: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise CalibrationError("run configuration is not an object")
    return value


def task_rows(run_root: Path) -> list[dict[str, str]]:
    path = run_root / "plan" / "calibration_tasks.tsv"
    if not path.is_file():
        raise CalibrationError(f"task plan not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise CalibrationError("task plan is empty")
    return rows


def task_directory(run_root: Path, row: dict[str, str]) -> Path:
    return run_root / "tasks" / (
        f"task-{int(row['task_index']):06d}-{row['motif_id']}"
    )


def status(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    rows = task_rows(run_root)
    complete = 0
    failed_contract = 0
    for row in rows:
        marker = task_directory(run_root, row) / "complete.json"
        if not marker.is_file():
            continue
        try:
            value = json.loads(marker.read_text(encoding="utf-8"))
            if value.get("motif_id") != row["motif_id"]:
                failed_contract += 1
            else:
                complete += 1
        except (OSError, json.JSONDecodeError):
            failed_contract += 1
    print("key\tvalue")
    print(f"planned\t{len(rows)}")
    print(f"complete\t{complete}")
    print(f"pending\t{len(rows) - complete - failed_contract}")
    print(f"invalid_complete_markers\t{failed_contract}")


def finalize(arguments: argparse.Namespace) -> None:
    run_root = arguments.run_root.expanduser().resolve()
    config = load_config(run_root)
    source = Path(str(config["source"]))
    finalization_commit, finalization_dirty = git_identity(source)
    rows = task_rows(run_root)
    metric_rows: list[dict[str, str]] = []
    metric_fields: list[str] | None = None
    missing: list[str] = []
    for row in rows:
        directory = task_directory(run_root, row)
        marker = directory / "complete.json"
        metrics = directory / "threshold_metrics.tsv"
        if not marker.is_file() or not metrics.is_file():
            missing.append(row["motif_id"])
            continue
        marker_value = json.loads(marker.read_text(encoding="utf-8"))
        if marker_value.get("motif_id") != row["motif_id"]:
            raise CalibrationError(f"task marker motif mismatch: {marker}")
        with metrics.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise CalibrationError(f"threshold metrics lack a header: {metrics}")
            if metric_fields is None:
                metric_fields = list(reader.fieldnames)
            elif metric_fields != reader.fieldnames:
                raise CalibrationError(f"threshold metric columns differ: {metrics}")
            current = list(reader)
        if not current or {item["motif_id"] for item in current} != {row["motif_id"]}:
            raise CalibrationError(f"threshold metric motif/cardinality error: {metrics}")
        metric_rows.extend(current)
    if missing:
        raise CalibrationError(
            f"cannot finalize; {len(missing)} tasks are incomplete, first={missing[:5]}"
        )
    if metric_fields is None:
        raise CalibrationError("no threshold metrics were collected")

    final = run_root / "final" / "threshold_calibration"
    if final.exists():
        manifest = final / "manifest.json"
        if manifest.is_file():
            print(f"I: Reusing finalized threshold calibration: {final}", file=sys.stderr)
            return
        raise CalibrationError(f"incomplete final output already exists: {final}")
    final.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".threshold-calibration.", dir=final.parent))
    try:
        metrics_tsv = staging / "threshold_metrics.tsv"
        with metrics_tsv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=metric_fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(metric_rows)
        duckdb = shutil.which(arguments.duckdb)
        if duckdb is None:
            raise CalibrationError(f"DuckDB executable not found: {arguments.duckdb}")
        metrics_parquet = staging / "threshold_metrics.parquet"
        process = subprocess.run([
            duckdb, "-batch", ":memory:", "-c",
            "COPY (SELECT * FROM read_csv_auto("
            f"{sql_string(metrics_tsv)}, delim='\\t', header=true, nullstr='NA')) "
            f"TO {sql_string(metrics_parquet)} (FORMAT PARQUET, COMPRESSION ZSTD);",
        ], text=True, capture_output=True, check=False)
        if process.returncode != 0:
            raise CalibrationError(process.stderr.strip() or "metric consolidation failed")

        table = staging / "tables" / "jaspar2026" / "motif_score_threshold"
        table.mkdir(parents=True)
        registry = table / "part-000000.parquet"
        builder = source / "scripts" / "build_motif_score_thresholds.py"
        command = [
            sys.executable, str(builder),
            "--metrics", str(metrics_parquet),
            "--source-metrics-uri", f"threshold-calibration://{config['run_id']}/threshold_metrics.parquet",
            "--jaspar", str(config["jaspar"]),
            "--jaspar-uri", "jaspar://2026/CORE/non-redundant",
            "--motif-list", str(run_root / "plan" / "motifs.txt"),
            "--output", str(registry),
            "--threshold-set-id", str(config["threshold_set_id"]),
            "--calibration-run-id", str(config["run_id"]),
            "--genome-id", "homo_sapiens_grch38_ensembl113_primary",
            "--motif-set-id", "jaspar2026_core_nonredundant",
            "--target-motif-id", str(config["target_motif_id"]),
            "--threshold-role", "tp73_context_binary_feature",
            "--score-mode", "log2_relative_risk",
            "--pseudocount", "1",
            "--background-model-id", "uniform_acgt_v1",
            "--pseudocount-scheme", "additive_per_base",
            "--source-minimum-score", "-1",
            "--context-flank", "150",
            "--context-distance-metric", "signed_interval_edge_distance",
            "--context-max-interval-distance", "150",
            "--context-relation-filter", "any",
            "--calibration-stratum-id", "all_tp73_anchors",
            "--calibration-stratum-json", '{"anchor_motif_id":"MA0861.2","anchor_minimum_score":0}',
            "--calibration-scope", "grch38_chr1",
            "--evidence-dataset-id", "rostock_p73_cutrun_20250602_noDuplicates",
            "--outcome-id", "discordant_anti_p73_only_vs_matched_control_only_strict_immersion",
            "--fold-definition", str(config["fold_definition"]),
            "--folds", str(config["folds"]),
            "--candidate-grid", str(config["candidate_grid"]),
            "--minimum-class-fraction", str(config["minimum_class_fraction"]),
            "--selection-metric", str(config["selection_metric"]),
            "--source-commit", str(config["source_commit"]),
            "--notes", "Exploratory chromosome-1 thresholds; TP73 itself uses separate direct/tandem policies.",
            "--duckdb", duckdb,
        ]
        if bool(config["source_dirty"]):
            command.append("--source-dirty")
        process = subprocess.run(command, text=True, capture_output=True, check=False)
        if process.returncode != 0:
            raise CalibrationError(process.stderr.strip() or "threshold registry build failed")
        summary = run_json([
            duckdb, "-json", ":memory:", "-c",
            "SELECT count(*) AS motifs, count(recommended_threshold) AS recommended, "
            "count(*) FILTER (WHERE calibration_status='pending') AS pending "
            f"FROM read_parquet({sql_string(registry)});",
        ])
        if len(summary) != 1 or int(summary[0]["motifs"]) != len(rows) or \
                int(summary[0]["pending"]) != 0:
            raise CalibrationError(f"final registry validation failed: {summary}")
        final_manifest = {
            "schema_version": 2,
            "run_id": config["run_id"],
            "threshold_set_id": config["threshold_set_id"],
            "motifs": len(rows),
            "recommended_thresholds": summary[0]["recommended"],
            "pending": summary[0]["pending"],
            "threshold_metrics_sha256": sha256(metrics_parquet),
            "threshold_registry_sha256": sha256(registry),
            "source_commit": config["source_commit"],
            "metric_source_commit": config["source_commit"],
            "metric_source_dirty": config["source_dirty"],
            "finalization_source_commit": finalization_commit,
            "finalization_source_dirty": finalization_dirty,
        }
        (staging / "manifest.json").write_text(
            json.dumps(final_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, final)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    print(f"I: Finalized {len(rows)} motif thresholds: {final}", file=sys.stderr)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    subparsers = result.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare", help="write immutable task plan")
    prepare_parser.add_argument("--run-root", type=Path, required=True)
    prepare_parser.add_argument("--scan-package", type=Path, required=True)
    prepare_parser.add_argument("--jaspar", type=Path, required=True)
    prepare_parser.add_argument("--anchor-evidence", type=Path, required=True)
    prepare_parser.add_argument("--source", type=Path, required=True)
    prepare_parser.add_argument("--duckdb", default="duckdb")
    prepare_parser.add_argument("--chrom", default="1")
    prepare_parser.add_argument("--chrom-size", type=int, default=248956422)
    prepare_parser.add_argument("--target-motif", default="MA0861.2")
    prepare_parser.add_argument("--minimum-class-fraction", type=float, default=0.01)
    prepare_parser.add_argument(
        "--run-id", default="jaspar2026_grch38_chr1_tp73_context_thresholds_v1"
    )
    prepare_parser.add_argument(
        "--threshold-set-id",
        default="tp73_chr1_cutrun_context_roc_auc_all_jaspar_v1",
    )
    prepare_parser.set_defaults(function=prepare)

    status_parser = subparsers.add_parser("status", help="report exact task state")
    status_parser.add_argument("--run-root", type=Path, required=True)
    status_parser.set_defaults(function=status)

    finalize_parser = subparsers.add_parser("finalize", help="build final registry")
    finalize_parser.add_argument("--run-root", type=Path, required=True)
    finalize_parser.add_argument("--duckdb", default="duckdb")
    finalize_parser.set_defaults(function=finalize)
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        arguments.function(arguments)
    except (CalibrationError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

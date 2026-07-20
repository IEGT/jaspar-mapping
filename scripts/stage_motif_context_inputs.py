#!/usr/bin/env python3

"""Stage exact scan-inventory files as a metadata-only hard-link tree."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


class StagingError(RuntimeError):
    pass


def sql_string(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def sql_list(values: list[str]) -> str:
    return "[" + ",".join(sql_string(value) for value in values) + "]"


def unique_values(values: list[str]) -> list[str]:
    result: list[str] = []
    for value in values:
        for item in value.split(","):
            item = item.strip()
            if not item:
                raise StagingError("motif/chromosome lists cannot contain empty values")
            if item not in result:
                result.append(item)
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def package_database(package: Path) -> tuple[Path, Path]:
    manifest = package / "manifest.json"
    if not manifest.is_file():
        raise StagingError(f"finalized scan manifest not found: {manifest}")
    try:
        content = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise StagingError(f"cannot read scan manifest: {error}") from error
    database_name = content.get("database")
    if not database_name:
        raise StagingError("scan manifest does not identify its DuckDB catalog")
    database = Path(database_name)
    if not database.is_absolute():
        database = package / database
    database = database.resolve()
    if not database.is_file():
        raise StagingError(f"scan DuckDB catalog not found: {database}")
    return manifest, database


def inventory_rows(duckdb: str, database: Path, motifs: list[str],
                   chromosomes: list[str]) -> list[dict[str, object]]:
    chrom_filter = ""
    if chromosomes:
        chrom_filter = (
            f" AND CAST(chrom AS VARCHAR) IN {sql_list(chromosomes)}"
        )
    query = f"""
SELECT task_id, output_relative_path, motif_id,
       CAST(chrom AS VARCHAR) AS chrom, strand, emitted_hits, bytes, sha256
FROM scan_file_inventory
WHERE motif_id IN {sql_list(motifs)}{chrom_filter}
ORDER BY motif_id, chrom, strand;
"""
    process = subprocess.run(
        [duckdb, "-readonly", "-json", str(database), "-c", query],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise StagingError(process.stderr.strip() or "DuckDB inventory query failed")
    try:
        rows = json.loads(process.stdout or "[]")
    except json.JSONDecodeError as error:
        raise StagingError("DuckDB returned invalid inventory JSON") from error
    if not isinstance(rows, list):
        raise StagingError("DuckDB inventory result is not a row array")
    return rows


def validate_inventory(rows: list[dict[str, object]], motifs: list[str],
                       chromosomes: list[str]) -> None:
    if not rows:
        raise StagingError("no scan inventory files matched the selection")
    observed: dict[tuple[str, str], set[str]] = {}
    for row in rows:
        key = (str(row["motif_id"]), str(row["chrom"]))
        strand = str(row["strand"])
        if strand not in {"+", "-"}:
            raise StagingError(f"unexpected inventory strand for {key}: {strand}")
        if strand in observed.setdefault(key, set()):
            raise StagingError(f"duplicate inventory strand for {key}: {strand}")
        observed[key].add(strand)
    expected_chromosomes = chromosomes or sorted({key[1] for key in observed})
    expected = {(motif, chrom) for motif in motifs for chrom in expected_chromosomes}
    if set(observed) != expected:
        missing = sorted(expected - set(observed))
        extra = sorted(set(observed) - expected)
        raise StagingError(
            f"inventory motif/chromosome mismatch; missing={missing}, extra={extra}"
        )
    incomplete = sorted(key for key, strands in observed.items()
                        if strands != {"+", "-"})
    if incomplete:
        raise StagingError(f"inventory lacks both strands for: {incomplete}")


def staged_relative_path(row: dict[str, object]) -> Path:
    relative = (
        Path("task_data") / f"task_id={row['task_id']}"
        / str(row["output_relative_path"])
    )
    if relative.is_absolute() or ".." in relative.parts:
        raise StagingError(f"unsafe inventory path: {relative}")
    return relative


def verify_existing(output: Path, contract: dict[str, object]) -> bool:
    manifest_path = output / "input_manifest.json"
    if not manifest_path.is_file():
        return False
    try:
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    for key in ("format_version", "source_package", "source_manifest_sha256",
                "motifs", "chromosomes", "files"):
        if existing.get(key) != contract.get(key):
            return False
    files = existing.get("files")
    if not isinstance(files, list):
        return False
    for row in files:
        if not isinstance(row, dict) or "staged_relative_path" not in row:
            return False
        target = output / str(row["staged_relative_path"])
        source = Path(str(row["source_path"]))
        if (not target.is_file() or not source.is_file()
                or target.stat().st_size != int(row["bytes"])
                or not os.path.samefile(source, target)):
            return False
    return True


def stage(arguments: argparse.Namespace) -> None:
    if shutil.which(arguments.duckdb) is None:
        raise StagingError(f"DuckDB executable not found: {arguments.duckdb}")
    package = arguments.package.expanduser().resolve()
    output = arguments.output.expanduser().resolve()
    if not package.is_dir():
        raise StagingError(f"scan package not found: {package}")
    if output == package or package in output.parents:
        raise StagingError("staged input must not be created inside the scan package")

    motifs = unique_values(arguments.motif)
    chromosomes = unique_values(arguments.chrom)
    manifest, database = package_database(package)
    rows = inventory_rows(arguments.duckdb, database, motifs, chromosomes)
    validate_inventory(rows, motifs, chromosomes)

    contract: dict[str, object] = {
        "format_version": 1,
        "source_package": str(package),
        "source_manifest_sha256": sha256(manifest),
        "motifs": motifs,
        "chromosomes": chromosomes,
    }
    file_manifest: list[dict[str, object]] = []
    for row in rows:
        relative = staged_relative_path(row)
        source = (package / relative).resolve()
        if package not in source.parents or not source.is_file():
            raise StagingError(f"inventory payload is missing or unsafe: {source}")
        recorded_bytes = int(row["bytes"])
        if source.stat().st_size != recorded_bytes:
            raise StagingError(f"inventory size mismatch: {source}")
        file_manifest.append({
            **row,
            "source_path": str(source),
            "staged_relative_path": str(relative),
        })
    contract["files"] = file_manifest
    if output.exists():
        if verify_existing(output, contract):
            print(f"I: Reusing verified context input tree: {output}", file=sys.stderr)
            return
        raise StagingError(f"output already exists but does not match: {output}")

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        for row in file_manifest:
            relative = Path(str(row["staged_relative_path"]))
            source = Path(str(row["source_path"]))
            target = staging / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            try:
                os.link(source, target)
            except OSError as error:
                raise StagingError(
                    f"cannot hard-link {source} to {target}; keep staging on "
                    f"the same filesystem as the scan package: {error}"
                ) from error
        (staging / "input_manifest.json").write_text(
            json.dumps(contract, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    print(
        f"I: Staged {len(rows)} exact inventory files as hard links: {output}",
        file=sys.stderr,
    )


def argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create a compact hard-link tree for exact motif/chromosome files "
            "selected from a finalized genome-scan inventory."
        )
    )
    parser.add_argument("--package", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--motif", required=True, action="append",
        help="motif accession; repeat or comma-separate",
    )
    parser.add_argument(
        "--chrom", action="append", default=[],
        help="chromosome; repeat or comma-separate (default: all inventoried)",
    )
    parser.add_argument("--duckdb", default="duckdb")
    return parser


def main() -> int:
    try:
        stage(argument_parser().parse_args())
        return 0
    except (OSError, StagingError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

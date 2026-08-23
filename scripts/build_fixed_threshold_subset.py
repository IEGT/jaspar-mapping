#!/usr/bin/env python3
"""Build a fixed-score sensitivity subset from a motif-threshold registry."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any


class FixedThresholdError(RuntimeError):
    """Raised when the requested threshold subset violates its contract."""


def finite_number(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise argparse.ArgumentTypeError("expected a finite number")
    return result


def identifier(value: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", value):
        raise argparse.ArgumentTypeError(
            "expected an identifier containing letters, digits, dot, underscore, or dash"
        )
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def read_motif_ids(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or "motif_id" not in reader.fieldnames:
            raise FixedThresholdError(f"motif list lacks motif_id: {path}")
        motif_ids = [row["motif_id"].strip() for row in reader]
    if not motif_ids or any(not motif_id for motif_id in motif_ids):
        raise FixedThresholdError("motif list contains no motifs or a blank motif_id")
    if len(motif_ids) != len(set(motif_ids)):
        raise FixedThresholdError("motif list contains duplicate motif_id values")
    return motif_ids


def run_duckdb(executable: str, sql: str) -> None:
    process = subprocess.run(
        [executable, "-batch", ":memory:"], input=sql, text=True,
        capture_output=True, check=False,
    )
    if process.returncode != 0:
        raise FixedThresholdError(
            process.stderr.strip() or "DuckDB threshold-subset build failed"
        )


def write_json(path: Path, value: dict[str, Any]) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--source-registry", type=Path, required=True)
    result.add_argument("--source-threshold-set-id", type=identifier, required=True)
    result.add_argument(
        "--motif-list", type=Path, required=True,
        help="tab-separated file containing one unique motif_id per requested motif",
    )
    result.add_argument("--fixed-threshold", type=finite_number, required=True)
    result.add_argument("--threshold-set-id", type=identifier, required=True)
    result.add_argument(
        "--output-dir", type=Path, required=True,
        help="new package directory receiving motif_score_threshold.parquet and manifest.json",
    )
    result.add_argument("--duckdb", default="duckdb")
    return result


def build(arguments: argparse.Namespace) -> None:
    source = arguments.source_registry.expanduser().resolve()
    motif_list = arguments.motif_list.expanduser().resolve()
    output_dir = arguments.output_dir.expanduser().resolve()
    if not source.is_file() or not motif_list.is_file():
        raise FixedThresholdError("source registry or motif list is missing")
    if output_dir.exists():
        raise FixedThresholdError(f"output directory already exists: {output_dir}")
    if arguments.threshold_set_id == arguments.source_threshold_set_id:
        raise FixedThresholdError("output and source threshold-set IDs must differ")
    if shutil.which(arguments.duckdb) is None:
        raise FixedThresholdError(f"DuckDB executable not found: {arguments.duckdb}")
    motif_ids = read_motif_ids(motif_list)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(
        prefix=f".{output_dir.name}-", dir=output_dir.parent
    ))
    try:
        output = staging / "motif_score_threshold.parquet"
        sql = f"""
CREATE TEMP TABLE requested AS
SELECT (row_number() OVER () - 1)::BIGINT AS motif_order,
       trim(motif_id)::VARCHAR AS motif_id
FROM read_csv_auto({sql_string(motif_list)}, delim='\\t', header=true,
                   all_varchar=true);
SELECT CASE WHEN (SELECT count(*) FROM requested) <> {len(motif_ids)}
  THEN error('motif-list cardinality changed during parsing') END;
CREATE TEMP TABLE source_threshold AS
SELECT motif_id::VARCHAR AS motif_id,
       recommended_threshold::DOUBLE AS recommended_threshold,
       source_minimum_score::DOUBLE AS source_minimum_score,
       calibration_status::VARCHAR AS calibration_status
FROM read_parquet({sql_string(source)})
WHERE threshold_set_id = {sql_string(arguments.source_threshold_set_id)};
SELECT CASE WHEN (SELECT count(*) FROM source_threshold s
                  JOIN requested r USING (motif_id)) <> {len(motif_ids)}
  THEN error('source registry does not contain every requested motif exactly once') END;
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM source_threshold s JOIN requested r USING (motif_id)
  WHERE {arguments.fixed_threshold} < s.source_minimum_score
) THEN error('fixed threshold is below a requested motif source score floor') END;
COPY (
  SELECT {sql_string(arguments.threshold_set_id)}::VARCHAR AS threshold_set_id,
         s.motif_id,
         {arguments.fixed_threshold}::DOUBLE AS recommended_threshold,
         s.source_minimum_score,
         'fixed_threshold_sensitivity'::VARCHAR AS calibration_status,
         {sql_string(arguments.source_threshold_set_id)}::VARCHAR
           AS source_threshold_set_id,
         s.recommended_threshold AS source_recommended_threshold,
         s.calibration_status AS source_calibration_status,
         r.motif_order
  FROM requested r JOIN source_threshold s USING (motif_id)
  ORDER BY r.motif_order
) TO {sql_string(output)} (FORMAT PARQUET, COMPRESSION ZSTD);
"""
        run_duckdb(arguments.duckdb, sql)
        manifest = {
            "analysis": "fixed_motif_threshold_sensitivity_subset",
            "schema_version": 1,
            "builder": "scripts/build_fixed_threshold_subset.py",
            "builder_sha256": sha256(Path(__file__).resolve()),
            "threshold_set_id": arguments.threshold_set_id,
            "fixed_threshold": arguments.fixed_threshold,
            "motif_count": len(motif_ids),
            "motif_list": str(motif_list),
            "motif_list_sha256": sha256(motif_list),
            "source_registry": str(source),
            "source_registry_sha256": sha256(source),
            "source_threshold_set_id": arguments.source_threshold_set_id,
            "output": output.name,
            "output_sha256": sha256(output),
            "created_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        }
        write_json(staging / "manifest.json", manifest)
        os.replace(staging, output_dir)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def main() -> int:
    arguments = parser().parse_args()
    try:
        build(arguments)
    except (FixedThresholdError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

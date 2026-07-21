#!/usr/bin/env python3

"""Build a TP73 anchor table with strict CUT&RUN immersion labels."""

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


class AnchorEvidenceError(RuntimeError):
    pass


SAMPLE_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def parse_coverage(value: str) -> tuple[str, Path]:
    sample_id, separator, filename = value.partition("=")
    if not separator or not SAMPLE_ID.fullmatch(sample_id) or not filename:
        raise argparse.ArgumentTypeError(
            "coverage must be COLUMN_ID=FILE with a SQL-safe column ID"
        )
    return sample_id, Path(filename)


def leading_metadata_lines(path: Path) -> int:
    count = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#") or \
                    stripped.startswith("track") or stripped.startswith("browser"):
                count += 1
                continue
            break
    return count


def coverage_ctes(index: int, sample_id: str, path: Path, chrom: str,
                   skip_lines: int) -> str:
    prefix = f"coverage_{index}"
    return f"""
{prefix}_text AS (
    SELECT trim(column0) AS chrom_text,
           try_cast(column1 AS BIGINT) AS start,
           try_cast(column2 AS BIGINT) AS interval_end,
           try_cast(column3 AS DOUBLE) AS depth
    FROM read_csv(
        {sql_string(path)}, delim='\\t', header=false,
        columns={{'column0':'VARCHAR', 'column1':'VARCHAR',
                 'column2':'VARCHAR', 'column3':'VARCHAR'}},
        skip={skip_lines}, comment='#', auto_detect=false,
        null_padding=true, strict_mode=false,
        quote='', escape=''
    )
    WHERE trim(column0) <> ''
      AND lower(trim(column0)) NOT IN ('track', 'browser')
),
{prefix}_selected AS (
    SELECT start, interval_end, depth
    FROM {prefix}_text
    WHERE regexp_replace(lower(chrom_text), '^chr', '') = {sql_string(chrom)}
),
{prefix}_ordered AS MATERIALIZED (
    SELECT *, lag(interval_end) OVER (ORDER BY start, interval_end) AS previous_end
    FROM {prefix}_selected
),
{prefix}_validation AS MATERIALIZED (
    SELECT CASE
        WHEN count(*) = 0 THEN error('coverage has no positive rows on chromosome {chrom}: {sample_id}')
        WHEN count(*) FILTER (
            WHERE start IS NULL OR interval_end IS NULL OR depth IS NULL
               OR start < 0 OR interval_end <= start
               OR NOT isfinite(depth) OR depth <= 0
        ) > 0 THEN error('coverage contains an invalid positive-depth row: {sample_id}')
        WHEN count(*) FILTER (WHERE previous_end IS NOT NULL AND start < previous_end) > 0
            THEN error('coverage rows overlap or are not coordinate-sorted: {sample_id}')
        ELSE true
    END AS valid
    FROM {prefix}_ordered
),
{prefix}_grouped AS MATERIALIZED (
    SELECT start, interval_end,
           sum(CASE WHEN previous_end IS NULL OR start > previous_end
                    THEN 1 ELSE 0 END)
               OVER (ORDER BY start, interval_end ROWS UNBOUNDED PRECEDING)
               AS component_id
    FROM {prefix}_ordered, {prefix}_validation
    WHERE {prefix}_validation.valid
),
{prefix}_components AS MATERIALIZED (
    SELECT min(start) AS component_start,
           max(interval_end) AS component_end
    FROM {prefix}_grouped
    GROUP BY component_id
)
""".strip()


def build_sql(arguments: argparse.Namespace, staging_output: Path) -> str:
    anchor_statement = f"""
CREATE TABLE anchors AS
    SELECT {sql_string(arguments.chrom)}::VARCHAR AS chrom,
           start::BIGINT AS anchor_start,
           \"end\"::BIGINT AS anchor_end,
           max(score)::FLOAT AS anchor_score
    FROM (
        SELECT start, \"end\", score
        FROM read_parquet({sql_string(arguments.anchor_plus)},
                          hive_partitioning=false)
        UNION ALL
        SELECT start, \"end\", score
        FROM read_parquet({sql_string(arguments.anchor_minus)},
                          hive_partitioning=false)
    )
    GROUP BY start, \"end\"
    HAVING max(score) >= {arguments.minimum_anchor_score:.17g};
""".strip()
    updates: list[str] = []
    for index, (sample_id, path, skip_lines) in enumerate(arguments.coverage):
        updates.append(
            f"ALTER TABLE anchors ADD COLUMN {sample_id} BOOLEAN DEFAULT false;\n"
            f"WITH {coverage_ctes(index, sample_id, path, arguments.chrom, skip_lines)}\n"
            f"UPDATE anchors AS a SET {sample_id} = true\n"
            f"FROM coverage_{index}_components AS c\n"
            "WHERE c.component_start < a.anchor_start\n"
            "  AND c.component_end > a.anchor_end;"
        )
    update_sql = "\n\n".join(updates)
    return f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;

{anchor_statement}

{update_sql}

COPY (
SELECT * FROM anchors
ORDER BY anchor_start, anchor_end
) TO {sql_string(staging_output)}
  (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 131072);
"""


def run_duckdb(duckdb: str, sql: str) -> None:
    process = subprocess.run(
        [duckdb, "-batch", ":memory:"],
        input=sql,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise AnchorEvidenceError(
            process.stderr.strip() or "DuckDB anchor-evidence build failed"
        )


def query_json(duckdb: str, query: str) -> list[dict[str, object]]:
    process = subprocess.run(
        [duckdb, "-json", ":memory:", "-c", query],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise AnchorEvidenceError(
            process.stderr.strip() or "DuckDB anchor-evidence validation failed"
        )
    value = json.loads(process.stdout or "[]")
    if not isinstance(value, list):
        raise AnchorEvidenceError("DuckDB validation did not return a row array")
    return value


def write_tsv_atomic(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def build(arguments: argparse.Namespace) -> None:
    duckdb = shutil.which(arguments.duckdb)
    if duckdb is None:
        raise AnchorEvidenceError(f"DuckDB executable not found: {arguments.duckdb}")

    arguments.anchor_plus = arguments.anchor_plus.expanduser().resolve()
    arguments.anchor_minus = arguments.anchor_minus.expanduser().resolve()
    normalized_coverage: list[tuple[str, Path, int]] = []
    seen_samples: set[str] = set()
    for sample_id, path in arguments.coverage:
        if sample_id in seen_samples:
            raise AnchorEvidenceError(f"duplicate coverage column: {sample_id}")
        seen_samples.add(sample_id)
        resolved = path.expanduser().resolve()
        normalized_coverage.append((
            sample_id, resolved, leading_metadata_lines(resolved)
        ))
    arguments.coverage = normalized_coverage
    for path in [arguments.anchor_plus, arguments.anchor_minus] + [
        path for _, path, _ in arguments.coverage
    ]:
        if not path.is_file():
            raise AnchorEvidenceError(f"input file not found: {path}")

    output = arguments.output.expanduser().resolve()
    run_config = Path(str(output) + ".run_config.tsv")
    coverage_manifest = Path(str(output) + ".coverage_manifest.tsv")
    for path in (output, run_config, coverage_manifest):
        if path.exists():
            raise AnchorEvidenceError(f"refusing to replace existing output: {path}")
    output.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f".{output.name}.", dir=output.parent) as temp:
        staging_output = Path(temp) / output.name
        run_duckdb(duckdb, build_sql(arguments, staging_output))
        summary = query_json(
            duckdb,
            "SELECT count(*) AS anchors, min(anchor_score) AS minimum_score, "
            "max(anchor_score) AS maximum_score FROM read_parquet("
            f"{sql_string(staging_output)});",
        )
        if len(summary) != 1 or int(summary[0].get("anchors", 0)) == 0:
            raise AnchorEvidenceError("anchor-evidence output is empty")
        os.replace(staging_output, output)

    write_tsv_atomic(
        run_config,
        [
            "schema_version", "chrom", "anchor_motif_id", "anchor_score_rule",
            "minimum_anchor_score", "strict_immersion_rule", "anchor_count",
            "observed_minimum_score", "observed_maximum_score", "anchor_plus",
            "anchor_plus_sha256", "anchor_minus", "anchor_minus_sha256",
            "coverage_track_count", "duckdb",
        ],
        [{
            "schema_version": 1,
            "chrom": arguments.chrom,
            "anchor_motif_id": "MA0861.2",
            "anchor_score_rule": "maximum_score_across_orientation_records_at_same_span",
            "minimum_anchor_score": arguments.minimum_anchor_score,
            "strict_immersion_rule": "component_start < anchor_start AND component_end > anchor_end",
            "anchor_count": summary[0]["anchors"],
            "observed_minimum_score": summary[0]["minimum_score"],
            "observed_maximum_score": summary[0]["maximum_score"],
            "anchor_plus": arguments.anchor_plus,
            "anchor_plus_sha256": sha256(arguments.anchor_plus),
            "anchor_minus": arguments.anchor_minus,
            "anchor_minus_sha256": sha256(arguments.anchor_minus),
            "coverage_track_count": len(arguments.coverage),
            "duckdb": duckdb,
        }],
    )
    write_tsv_atomic(
        coverage_manifest,
        ["support_column", "path", "bytes", "mtime_ns", "sha256"],
        [{
            "support_column": sample_id,
            "path": path,
            "bytes": path.stat().st_size,
            "mtime_ns": path.stat().st_mtime_ns,
            "sha256": sha256(path),
        } for sample_id, path, _ in arguments.coverage],
    )
    print(
        f"I: wrote {summary[0]['anchors']} TP73 anchors with "
        f"{len(arguments.coverage)} CUT&RUN support columns: {output}",
        file=sys.stderr,
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Collapse plus/minus sparse TP73 records by alignment span and label "
            "each retained anchor by strict immersion in merged positive-depth "
            "CUT&RUN bedGraph components."
        )
    )
    result.add_argument("--anchor-plus", type=Path, required=True)
    result.add_argument("--anchor-minus", type=Path, required=True)
    result.add_argument(
        "--coverage", action="append", type=parse_coverage, required=True,
        metavar="COLUMN_ID=FILE",
        help="support-column name and sorted positive-depth bedGraph; repeat",
    )
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--chrom", default="1")
    result.add_argument("--minimum-anchor-score", type=float, default=0.0)
    result.add_argument("--duckdb", default="duckdb")
    result.add_argument("--threads", type=int, default=4)
    result.add_argument("--memory-limit", default="12GB")
    return result


def main() -> int:
    arguments = parser().parse_args()
    if arguments.threads <= 0:
        parser().error("--threads must be positive")
    if not (-1e100 < arguments.minimum_anchor_score < 1e100):
        parser().error("--minimum-anchor-score must be finite")
    try:
        build(arguments)
    except (AnchorEvidenceError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

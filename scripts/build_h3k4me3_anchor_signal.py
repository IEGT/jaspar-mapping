#!/usr/bin/env python3

"""Build GFP-blinded H3K4me3 profiles and windowed TP73-anchor signal."""

from __future__ import annotations

import argparse
import csv
import errno
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path


class SignalBuildError(RuntimeError):
    pass


IDENTIFIER = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
CONDITIONS = ("GFP", "TA", "DN")
CHANNELS = ("h3k4me3", "input", "tp73", "negative_control")
WINDOW_DEFAULTS = (
    ("central_150", 0, 150),
    ("central_500", 0, 500),
    ("central_1000", 0, 1000),
    ("flank_150_1000", 150, 1000),
)


@dataclass(frozen=True)
class Track:
    series_id: str
    cell_line: str
    condition: str
    channel: str
    replicate: str
    filename: str
    included: bool
    exclusion_reason: str
    reference_build: str
    signal_scale: str


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_identifier(value: str) -> str:
    if not IDENTIFIER.fullmatch(value):
        raise SignalBuildError(f"unsafe SQL identifier: {value}")
    return '"' + value + '"'


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise SignalBuildError(f"invalid boolean value: {value}")


def parse_window(value: str) -> tuple[str, int, int]:
    fields = value.split(":")
    if len(fields) != 3 or not IDENTIFIER.fullmatch(fields[0]):
        raise argparse.ArgumentTypeError("window must be NAME:INNER_BP:OUTER_BP")
    try:
        inner, outer = int(fields[1]), int(fields[2])
    except ValueError as error:
        raise argparse.ArgumentTypeError("window bounds must be integers") from error
    if inner < 0 or outer <= inner:
        raise argparse.ArgumentTypeError("window requires 0 <= INNER_BP < OUTER_BP")
    return fields[0], inner, outer


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> list[Track]:
    required = {
        "series_id", "cell_line", "condition", "channel", "replicate",
        "filename", "analysis_included", "exclusion_reason",
        "reference_build", "signal_scale",
    }
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise SignalBuildError(
                "track manifest lacks columns: " + ", ".join(sorted(missing))
            )
        tracks = [Track(
            series_id=row["series_id"].strip(),
            cell_line=row["cell_line"].strip(),
            condition=row["condition"].strip(),
            channel=row["channel"].strip(),
            replicate=row["replicate"].strip(),
            filename=row["filename"].strip(),
            included=parse_bool(row["analysis_included"]),
            exclusion_reason=row["exclusion_reason"].strip(),
            reference_build=row["reference_build"].strip(),
            signal_scale=row["signal_scale"].strip(),
        ) for row in reader]
    if not tracks:
        raise SignalBuildError("track manifest is empty")
    seen: set[tuple[str, str, str]] = set()
    for track in tracks:
        key = (track.series_id, track.condition, track.channel)
        if key in seen:
            raise SignalBuildError(f"duplicate track manifest key: {key}")
        seen.add(key)
        if not IDENTIFIER.fullmatch(track.series_id):
            raise SignalBuildError(f"invalid series_id: {track.series_id}")
        if track.condition not in CONDITIONS or track.channel not in CHANNELS:
            raise SignalBuildError(f"invalid condition/channel: {key}")
        if not track.filename or not track.reference_build or not track.signal_scale:
            raise SignalBuildError(f"incomplete track metadata: {key}")
        if track.included and track.exclusion_reason:
            raise SignalBuildError(f"included track has exclusion reason: {key}")
        if not track.included and not track.exclusion_reason:
            raise SignalBuildError(f"excluded track lacks a reason: {key}")
    included_series = sorted({track.series_id for track in tracks if track.included})
    if not included_series:
        raise SignalBuildError("track manifest includes no analysis series")
    for series_id in included_series:
        expected = {(series_id, condition, channel)
                    for condition in CONDITIONS for channel in CHANNELS}
        actual = {key for key in seen if key[0] == series_id and next(
            track.included for track in tracks
            if (track.series_id, track.condition, track.channel) == key
        )}
        if actual != expected:
            missing = expected - actual
            raise SignalBuildError(
                f"included series {series_id} is not factorial; missing {sorted(missing)}"
            )
        series_tracks = [track for track in tracks
                         if track.series_id == series_id and track.included]
        for field in ("cell_line", "replicate", "reference_build", "signal_scale"):
            observed = {getattr(track, field) for track in series_tracks}
            if len(observed) != 1:
                raise SignalBuildError(
                    f"included series {series_id} has inconsistent {field}: "
                    f"{sorted(observed)}"
                )
    return tracks


def leading_metadata_lines(path: Path) -> int:
    count = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "track", "browser")):
                count += 1
            else:
                break
    return count


def run(command: list[str], *, input_text: str | None = None) -> None:
    process = subprocess.run(
        command, input=input_text, text=True, capture_output=True, check=False
    )
    if process.returncode != 0:
        raise SignalBuildError(
            process.stderr.strip() or process.stdout.strip() or
            f"command failed: {' '.join(command)}"
        )


def promote_file(source: Path, target: Path) -> None:
    """Publish a complete file atomically, including across filesystems."""
    try:
        os.replace(source, target)
        return
    except OSError as error:
        if error.errno != errno.EXDEV:
            raise
    temporary = target.with_name(f".{target.name}.tmp-{os.getpid()}")
    if temporary.exists():
        raise SignalBuildError(f"promotion temporary already exists: {temporary}")
    try:
        shutil.copy2(source, temporary)
        os.replace(temporary, target)
    finally:
        if temporary.exists():
            temporary.unlink()


def stage_track(track: Track, root: Path, scratch: Path, chrom: str,
                rscript: str, exporter: Path) -> tuple[Path, Path]:
    source = (root / track.filename).resolve()
    if not source.is_file():
        raise SignalBuildError(f"track not found: {source}")
    lower = source.name.lower()
    if lower.endswith((".bedgraph", ".bdg")):
        return source, source
    if not lower.endswith((".bigwig", ".bw")):
        raise SignalBuildError(f"unsupported track format: {source}")
    output = scratch / (
        f"{track.series_id}_{track.condition}_{track.channel}_chr{chrom}.bedGraph"
    )
    run([
        rscript, str(exporter), "--input", str(source), "--output", str(output),
        "--chrom", chrom,
    ])
    return source, output


def coverage_table_sql(index: int, path: Path, chrom: str) -> str:
    skip = leading_metadata_lines(path)
    name = f"coverage_{index}"
    return f"""
CREATE OR REPLACE TEMP TABLE {name} AS
WITH text_rows AS (
    SELECT trim(column0) AS chrom_text,
           try_cast(column1 AS BIGINT) AS start,
           try_cast(column2 AS BIGINT) AS interval_end,
           try_cast(column3 AS DOUBLE) AS signal
    FROM read_csv(
        {sql_string(path)}, delim='\\t', header=false,
        columns={{'column0':'VARCHAR', 'column1':'VARCHAR',
                 'column2':'VARCHAR', 'column3':'VARCHAR'}},
        skip={skip}, comment='#', auto_detect=false, null_padding=true,
        strict_mode=false, quote='', escape=''
    )
), selected AS (
    SELECT start, interval_end, signal
    FROM text_rows
    WHERE regexp_replace(lower(chrom_text), '^chr', '') = {sql_string(chrom)}
), checked AS (
    SELECT *, lag(interval_end) OVER (ORDER BY start, interval_end) AS previous_end
    FROM selected
)
SELECT start, interval_end, signal
FROM checked
WHERE CASE
    WHEN start IS NULL OR interval_end IS NULL OR signal IS NULL
      OR start < 0 OR interval_end <= start
      OR NOT isfinite(signal) OR signal <= 0
      THEN error('invalid positive signal row in {path.name}')
    WHEN previous_end IS NOT NULL AND start < previous_end
      THEN error('overlapping or unsorted signal rows in {path.name}')
    ELSE true
END;
SELECT CASE WHEN count(*) = 0
    THEN error('no positive chromosome {chrom} rows in {path.name}')
    ELSE true END FROM {name};
""".strip()


def anchor_sql(arguments: argparse.Namespace, source: Path | None = None) -> str:
    if source is not None:
        return f"""
CREATE TABLE anchors AS
SELECT CAST(chrom AS VARCHAR) AS chrom, anchor_start, anchor_end, anchor_score
FROM read_parquet({sql_string(source)}, hive_partitioning=false)
ORDER BY anchor_start, anchor_end;
""".strip()
    if arguments.anchor_source is not None:
        return f"""
CREATE TABLE anchors AS
SELECT {sql_string(arguments.chrom)}::VARCHAR AS chrom,
       start::BIGINT AS anchor_start, "end"::BIGINT AS anchor_end,
       max(score)::FLOAT AS anchor_score
FROM read_parquet(
    {sql_string(arguments.anchor_source)}, hive_partitioning=false)
WHERE motif_id = 'MA0861.2'
  AND anchor_selection_class = 'local_peak'
GROUP BY start, "end"
HAVING max(score) >= {arguments.minimum_anchor_score:.17g}
ORDER BY anchor_start, anchor_end;
""".strip()
    return f"""
CREATE TABLE anchors AS
SELECT {sql_string(arguments.chrom)}::VARCHAR AS chrom,
       start::BIGINT AS anchor_start, "end"::BIGINT AS anchor_end,
       max(score)::FLOAT AS anchor_score
FROM (
    SELECT start, "end", score FROM read_parquet(
        {sql_string(arguments.anchor_plus)}, hive_partitioning=false)
    UNION ALL
    SELECT start, "end", score FROM read_parquet(
        {sql_string(arguments.anchor_minus)}, hive_partitioning=false)
)
GROUP BY start, "end"
HAVING max(score) >= {arguments.minimum_anchor_score:.17g}
ORDER BY anchor_start, anchor_end;
""".strip()


def profile_sql(arguments: argparse.Namespace, tracks: list[tuple[Track, Path]],
                output: Path, anchor_source: Path | None) -> str:
    sample_size = arguments.profile_sample_size
    bins = arguments.profile_flank * 2 // arguments.profile_bin_size
    statements = [
        f"SET threads={arguments.threads};",
        f"SET memory_limit={sql_string(arguments.memory_limit)};",
        "SET preserve_insertion_order=false;",
        anchor_sql(arguments, anchor_source),
        "SELECT CASE WHEN count(*) = 0 THEN error('anchor set is empty') "
        "ELSE true END FROM anchors;",
        f"""
CREATE TABLE profile_anchors AS
WITH ordered AS (
    SELECT *, row_number() OVER (ORDER BY anchor_start, anchor_end) AS rn,
           count(*) OVER () AS total
    FROM anchors
), spaced AS (
    SELECT *
    FROM ordered
    WHERE ((rn - 1) % greatest(1, ceil(total::DOUBLE / {sample_size})::BIGINT)) = 0
    ORDER BY rn
    LIMIT {sample_size}
)
SELECT chrom, anchor_start, anchor_end, anchor_score,
       floor((anchor_start + anchor_end) / 2.0)::BIGINT AS anchor_center
FROM spaced;
CREATE TABLE profile_segments AS
SELECT a.chrom, a.anchor_start, a.anchor_end, a.anchor_score,
       b.offset_start_bp, b.offset_start_bp + {arguments.profile_bin_size}
           AS offset_end_bp,
       greatest(0, a.anchor_center + b.offset_start_bp)::BIGINT AS segment_start,
       least({arguments.chrom_length},
             a.anchor_center + b.offset_start_bp + {arguments.profile_bin_size})::BIGINT
           AS segment_end
FROM profile_anchors AS a
CROSS JOIN range(-{arguments.profile_flank}, {arguments.profile_flank},
                 {arguments.profile_bin_size}) AS b(offset_start_bp)
WHERE greatest(0, a.anchor_center + b.offset_start_bp) <
      least({arguments.chrom_length},
            a.anchor_center + b.offset_start_bp + {arguments.profile_bin_size});
CREATE TABLE profile_denominator AS
SELECT offset_start_bp, offset_end_bp,
       count(*)::BIGINT AS anchors_profiled,
       sum(segment_end - segment_start)::BIGINT AS effective_bp
FROM profile_segments
GROUP BY offset_start_bp, offset_end_bp;
CREATE TABLE profile_long (
    series_id VARCHAR, cell_line VARCHAR, channel VARCHAR,
    offset_start_bp BIGINT, offset_end_bp BIGINT,
    anchors_profiled BIGINT, effective_bp BIGINT,
    signal_area DOUBLE, mean_signal_per_bp DOUBLE,
    positive_signal_bp BIGINT, positive_bp_fraction DOUBLE
);
""".strip(),
    ]
    for index, (track, path) in enumerate(tracks):
        statements.append(coverage_table_sql(index, path, arguments.chrom))
        statements.append(f"""
INSERT INTO profile_long
WITH observed AS (
    SELECT p.offset_start_bp, p.offset_end_bp,
           sum((least(c.interval_end, p.segment_end) -
                greatest(c.start, p.segment_start)) * c.signal)::DOUBLE
               AS signal_area,
           sum(least(c.interval_end, p.segment_end) -
               greatest(c.start, p.segment_start))::BIGINT AS positive_signal_bp
    FROM profile_segments AS p
    JOIN coverage_{index} AS c
      ON c.start < p.segment_end AND c.interval_end > p.segment_start
    GROUP BY p.offset_start_bp, p.offset_end_bp
)
SELECT {sql_string(track.series_id)}, {sql_string(track.cell_line)},
       {sql_string(track.channel)}, d.offset_start_bp, d.offset_end_bp,
       d.anchors_profiled, d.effective_bp,
       coalesce(o.signal_area, 0),
       coalesce(o.signal_area, 0) / d.effective_bp,
       coalesce(o.positive_signal_bp, 0),
       coalesce(o.positive_signal_bp, 0)::DOUBLE / d.effective_bp
FROM profile_denominator AS d
LEFT JOIN observed AS o USING (offset_start_bp, offset_end_bp);
DROP TABLE coverage_{index};
""".strip())
    if output.suffix.lower() == ".parquet":
        copy = (f"COPY (SELECT * FROM profile_long ORDER BY series_id, channel, "
                f"offset_start_bp) TO {sql_string(output)} "
                "(FORMAT PARQUET, COMPRESSION ZSTD);")
    else:
        copy = (f"COPY (SELECT * FROM profile_long ORDER BY series_id, channel, "
                f"offset_start_bp) TO {sql_string(output)} "
                "(FORMAT CSV, DELIMITER '\\t', HEADER true);")
    statements.extend([
        f"SELECT CASE WHEN count(*) <> {len(tracks) * bins} "
        "THEN error('metaprofile row count is incomplete') ELSE true END "
        "FROM profile_long;",
        copy,
    ])
    return "\n\n".join(statements) + "\n"


def window_segments_sql(windows: tuple[tuple[str, int, int], ...],
                        chrom_length: int) -> str:
    values = ",\n".join(
        f"({sql_string(name)}, {inner}, {outer})" for name, inner, outer in windows
    )
    return f"""
CREATE TABLE window_specs(window_name VARCHAR, inner_bp BIGINT, outer_bp BIGINT);
INSERT INTO window_specs VALUES
{values};
CREATE TABLE anchor_centers AS
SELECT *, floor((anchor_start + anchor_end) / 2.0)::BIGINT AS anchor_center
FROM anchors;
CREATE TABLE window_segments AS
SELECT chrom, anchor_start, anchor_end, anchor_score,
       'motif_span'::VARCHAR AS window_name, 0::BIGINT AS inner_bp,
       (anchor_end - anchor_start)::BIGINT AS outer_bp,
       0::INTEGER AS segment_index, anchor_start AS segment_start,
       anchor_end AS segment_end
FROM anchor_centers
UNION ALL
SELECT a.chrom, a.anchor_start, a.anchor_end, a.anchor_score,
       w.window_name, w.inner_bp, w.outer_bp, 0,
       greatest(0, a.anchor_center - w.outer_bp)::BIGINT,
       least({chrom_length}, a.anchor_center + w.outer_bp)::BIGINT
FROM anchor_centers AS a
JOIN window_specs AS w ON w.inner_bp = 0
UNION ALL
SELECT a.chrom, a.anchor_start, a.anchor_end, a.anchor_score,
       w.window_name, w.inner_bp, w.outer_bp, side.segment_index,
       CASE side.segment_index
         WHEN 0 THEN greatest(0, a.anchor_center - w.outer_bp)::BIGINT
         ELSE greatest(0, a.anchor_center + w.inner_bp)::BIGINT END,
       CASE side.segment_index
         WHEN 0 THEN least({chrom_length}, a.anchor_center - w.inner_bp)::BIGINT
         ELSE least({chrom_length}, a.anchor_center + w.outer_bp)::BIGINT END
FROM anchor_centers AS a
JOIN window_specs AS w ON w.inner_bp > 0
CROSS JOIN (VALUES (0), (1)) AS side(segment_index);
DELETE FROM window_segments WHERE segment_end <= segment_start;
CREATE TABLE window_denominator AS
SELECT chrom, anchor_start, anchor_end, anchor_score, window_name,
       min(inner_bp)::BIGINT AS inner_bp, max(outer_bp)::BIGINT AS outer_bp,
       count(*)::INTEGER AS segment_count,
       sum(segment_end - segment_start)::BIGINT AS effective_window_bp
FROM window_segments
GROUP BY chrom, anchor_start, anchor_end, anchor_score, window_name;
""".strip()


def signal_sql(arguments: argparse.Namespace, tracks: list[tuple[Track, Path]],
               anchor_source: Path, output: Path,
               windows: tuple[tuple[str, int, int], ...]) -> str:
    signal_tracks = [pair for pair in tracks
                     if pair[0].channel in {"h3k4me3", "input"}]
    series_count = len({track.series_id for track, _ in signal_tracks})
    window_count = len(windows) + 1
    statements = [
        f"SET threads={arguments.threads};",
        f"SET memory_limit={sql_string(arguments.memory_limit)};",
        f"SET temp_directory={sql_string(arguments.scratch_database.parent / 'duckdb_tmp')};",
        "SET preserve_insertion_order=false;",
        anchor_sql(arguments, anchor_source),
        window_segments_sql(windows, arguments.chrom_length),
        """
CREATE TABLE signal_long (
    chrom VARCHAR, anchor_start BIGINT, anchor_end BIGINT, anchor_score FLOAT,
    series_id VARCHAR, cell_line VARCHAR, condition VARCHAR, replicate VARCHAR,
    window_name VARCHAR, inner_bp BIGINT, outer_bp BIGINT,
    segment_count INTEGER, effective_window_bp BIGINT, channel VARCHAR,
    signal_area DOUBLE, mean_signal_per_bp DOUBLE, maximum_signal DOUBLE,
    positive_signal_bp BIGINT, positive_bp_fraction DOUBLE
);
""".strip(),
    ]
    for index, (track, path) in enumerate(signal_tracks):
        statements.append(coverage_table_sql(index, path, arguments.chrom))
        statements.append(f"""
INSERT INTO signal_long
WITH observed AS (
    SELECT w.chrom, w.anchor_start, w.anchor_end, w.window_name,
           sum((least(c.interval_end, w.segment_end) -
                greatest(c.start, w.segment_start)) * c.signal)::DOUBLE
               AS signal_area,
           max(c.signal)::DOUBLE AS maximum_signal,
           sum(least(c.interval_end, w.segment_end) -
               greatest(c.start, w.segment_start))::BIGINT AS positive_signal_bp
    FROM window_segments AS w
    JOIN coverage_{index} AS c
      ON c.start < w.segment_end AND c.interval_end > w.segment_start
    GROUP BY w.chrom, w.anchor_start, w.anchor_end, w.window_name
)
SELECT d.chrom, d.anchor_start, d.anchor_end, d.anchor_score,
       {sql_string(track.series_id)}, {sql_string(track.cell_line)},
       {sql_string(track.condition)}, {sql_string(track.replicate)},
       d.window_name, d.inner_bp, d.outer_bp, d.segment_count,
       d.effective_window_bp, {sql_string(track.channel)},
       coalesce(o.signal_area, 0),
       coalesce(o.signal_area, 0) / d.effective_window_bp,
       coalesce(o.maximum_signal, 0), coalesce(o.positive_signal_bp, 0),
       coalesce(o.positive_signal_bp, 0)::DOUBLE / d.effective_window_bp
FROM window_denominator AS d
LEFT JOIN observed AS o USING (chrom, anchor_start, anchor_end, window_name);
DROP TABLE coverage_{index};
""".strip())
    expected_factor = series_count * len(CONDITIONS) * window_count * 2
    statements.extend([
        f"""
SELECT CASE WHEN (SELECT count(*) FROM signal_long) <>
    (SELECT count(*) * {expected_factor} FROM anchors)
THEN error('window signal row count is incomplete') ELSE true END;
COPY (
    SELECT chrom, anchor_start, anchor_end, anchor_score,
           series_id, cell_line, condition, replicate, window_name,
           min(inner_bp)::BIGINT AS inner_bp,
           max(outer_bp)::BIGINT AS outer_bp,
           max(segment_count)::INTEGER AS segment_count,
           max(effective_window_bp)::BIGINT AS effective_window_bp,
           max(signal_area) FILTER (channel = 'h3k4me3') AS h3k4me3_area,
           max(mean_signal_per_bp) FILTER (channel = 'h3k4me3') AS h3k4me3_mean,
           max(maximum_signal) FILTER (channel = 'h3k4me3') AS h3k4me3_max,
           max(positive_signal_bp) FILTER (channel = 'h3k4me3')
               AS h3k4me3_positive_bp,
           max(positive_bp_fraction) FILTER (channel = 'h3k4me3')
               AS h3k4me3_positive_fraction,
           max(signal_area) FILTER (channel = 'input') AS input_area,
           max(mean_signal_per_bp) FILTER (channel = 'input') AS input_mean,
           max(maximum_signal) FILTER (channel = 'input') AS input_max,
           max(positive_signal_bp) FILTER (channel = 'input')
               AS input_positive_bp,
           max(positive_bp_fraction) FILTER (channel = 'input')
               AS input_positive_fraction
    FROM signal_long
    GROUP BY chrom, anchor_start, anchor_end, anchor_score,
             series_id, cell_line, condition, replicate, window_name
    HAVING count(*) = 2
    ORDER BY series_id, condition, window_name, anchor_start, anchor_end
) TO {sql_string(output)}
  (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 131072);
""".strip()
    ])
    return "\n\n".join(statements) + "\n"


def write_sidecars(base: Path, arguments: argparse.Namespace, tracks: list[Track],
                   used: list[tuple[Track, Path, Path]],
                   windows: tuple[tuple[str, int, int], ...], mode: str) -> None:
    used_keys = {(track.series_id, track.condition, track.channel)
                 for track, _, _ in used}
    manifest_path = Path(str(base) + ".track_manifest.tsv")
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        fields = [
            "series_id", "cell_line", "condition", "channel", "replicate",
            "analysis_included", "used_in_run", "exclusion_reason", "filename",
            "resolved_source", "bytes", "mtime_ns", "sha256",
            "reference_build", "signal_scale",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        source_by_key = {
            (track.series_id, track.condition, track.channel): source
            for track, source, _ in used
        }
        for track in tracks:
            key = (track.series_id, track.condition, track.channel)
            source = source_by_key.get(key)
            stat = source.stat() if source else None
            writer.writerow({
                "series_id": track.series_id,
                "cell_line": track.cell_line,
                "condition": track.condition,
                "channel": track.channel,
                "replicate": track.replicate,
                "analysis_included": str(track.included).lower(),
                "used_in_run": str(key in used_keys).lower(),
                "exclusion_reason": track.exclusion_reason,
                "filename": track.filename,
                "resolved_source": str(source) if source else "",
                "bytes": stat.st_size if stat else "",
                "mtime_ns": stat.st_mtime_ns if stat else "",
                "sha256": sha256(source) if source else "",
                "reference_build": track.reference_build,
                "signal_scale": track.signal_scale,
            })
    run_config = {
        "schema_version": 1,
        "analysis": "h3k4me3_gfp_referenced_tp73_anchor_signal",
        "mode": mode,
        "chrom": arguments.chrom,
        "chrom_length": arguments.chrom_length,
        "minimum_anchor_score": arguments.minimum_anchor_score,
        "anchor_source_mode": (
            "schema7_local_peak_context_anchor"
            if arguments.anchor_source is not None
            else "orientation_sparse_pair"
        ),
        "anchor_source": (
            str(arguments.anchor_source) if arguments.anchor_source else None
        ),
        "included_series": sorted({track.series_id for track in tracks
                                    if track.included}),
        "excluded_series": sorted({track.series_id for track in tracks
                                    if not track.included}),
        "windows": [
            {"window_name": "motif_span", "geometry": "anchor_alignment_span"},
            *[
                {"window_name": name, "inner_bp": inner, "outer_bp": outer}
                for name, inner, outer in windows
            ],
        ],
        "profile": {
            "condition": "GFP",
            "flank_bp": arguments.profile_flank,
            "bin_size_bp": arguments.profile_bin_size,
            "sample_size": arguments.profile_sample_size,
            "selection": "deterministic_genome_order_spaced_sample",
        },
        "persistent_intermediate_bedgraph": False,
        "track_manifest": str(arguments.track_manifest),
    }
    Path(str(base) + ".run_config.json").write_text(
        json.dumps(run_config, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parser() -> argparse.ArgumentParser:
    repository = Path(__file__).resolve().parent
    result = argparse.ArgumentParser(
        description=(
            "Build a GFP-only H3K4me3/input metaprofile and, in full mode, "
            "strict TP73 support plus windowed H3K4me3/input Parquet. BigWig "
            "is converted only in temporary scratch; no persistent BED layer "
            "is created."
        )
    )
    result.add_argument(
        "--anchor-source", type=Path,
        help=(
            "one-chromosome schema-7 tp73_context_anchor Parquet already "
            "selected as local peaks; chromosome is supplied by --chrom "
            "because partition columns need not be stored in the file; "
            "mutually exclusive with --anchor-plus/--anchor-minus"
        ),
    )
    result.add_argument("--anchor-plus", type=Path)
    result.add_argument("--anchor-minus", type=Path)
    result.add_argument("--track-manifest", type=Path, required=True)
    result.add_argument("--track-root", type=Path, required=True)
    result.add_argument("--tp73-output", type=Path)
    result.add_argument("--signal-output", type=Path)
    result.add_argument("--profile-output", type=Path)
    result.add_argument("--profile-only", action="store_true")
    result.add_argument("--chrom", default="1")
    result.add_argument("--chrom-length", type=int, required=True)
    result.add_argument("--minimum-anchor-score", type=float, default=0.0)
    result.add_argument(
        "--window", type=parse_window, action="append", metavar="NAME:INNER:OUTER",
        help="repeatable symmetric central/annular window; defaults include 150/500/1000 bp",
    )
    result.add_argument("--profile-flank", type=int, default=2000)
    result.add_argument("--profile-bin-size", type=int, default=50)
    result.add_argument("--profile-sample-size", type=int, default=20000)
    result.add_argument("--threads", type=int, default=4)
    result.add_argument("--memory-limit", default="8GB")
    result.add_argument("--duckdb", default="duckdb")
    result.add_argument("--rscript", default="Rscript")
    result.add_argument(
        "--bigwig-exporter", type=Path,
        default=repository / "export_bigwig_chrom_bedgraph.R",
    )
    result.add_argument(
        "--anchor-builder", type=Path,
        default=repository / "build_tp73_anchor_evidence.py",
    )
    result.add_argument("--tmpdir", type=Path)
    return result


def validate_arguments(arguments: argparse.Namespace) -> None:
    source_mode = arguments.anchor_source is not None
    pair_mode = arguments.anchor_plus is not None or arguments.anchor_minus is not None
    if source_mode == pair_mode:
        raise SignalBuildError(
            "provide exactly --anchor-source or both --anchor-plus and --anchor-minus"
        )
    if pair_mode and (arguments.anchor_plus is None or arguments.anchor_minus is None):
        raise SignalBuildError("--anchor-plus and --anchor-minus must be provided together")
    anchor_inputs = (
        (arguments.anchor_source,)
        if source_mode else (arguments.anchor_plus, arguments.anchor_minus)
    )
    for path in (*anchor_inputs, arguments.track_manifest):
        if not path.is_file():
            raise SignalBuildError(f"input not found: {path}")
    if not arguments.track_root.is_dir():
        raise SignalBuildError(f"track root not found: {arguments.track_root}")
    if arguments.chrom_length <= 0 or not math.isfinite(arguments.minimum_anchor_score):
        raise SignalBuildError("chromosome length and anchor score must be finite/positive")
    if arguments.profile_flank <= 0 or arguments.profile_bin_size <= 0 or \
            2 * arguments.profile_flank % arguments.profile_bin_size != 0:
        raise SignalBuildError("profile bin size must divide twice the profile flank")
    if arguments.profile_sample_size <= 0 or arguments.threads <= 0:
        raise SignalBuildError("profile sample size and threads must be positive")
    if arguments.profile_only:
        if arguments.profile_output is None:
            raise SignalBuildError("--profile-only requires --profile-output")
    elif arguments.tp73_output is None or arguments.signal_output is None:
        raise SignalBuildError("full mode requires --tp73-output and --signal-output")
    outputs = [path for path in (
        arguments.tp73_output, arguments.signal_output, arguments.profile_output
    ) if path is not None]
    if len(set(outputs)) != len(outputs):
        raise SignalBuildError("output paths must be distinct")
    for output in outputs:
        if output.exists():
            raise SignalBuildError(f"refusing to replace existing output: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
    if shutil.which(arguments.duckdb) is None:
        raise SignalBuildError(f"DuckDB CLI not found: {arguments.duckdb}")


def main() -> int:
    arguments = parser().parse_args()
    try:
        validate_arguments(arguments)
        tracks = read_manifest(arguments.track_manifest)
        windows = tuple(arguments.window or WINDOW_DEFAULTS)
        if len({name for name, _, _ in windows}) != len(windows):
            raise SignalBuildError("window names must be unique")
        included = [track for track in tracks if track.included]
        if arguments.profile_only:
            required = [track for track in included
                        if track.condition == "GFP" and
                        track.channel in {"h3k4me3", "input"}]
        else:
            required = included
        temp_parent = arguments.tmpdir.resolve() if arguments.tmpdir else None
        with tempfile.TemporaryDirectory(
            prefix="jaspar-h3k4me3-", dir=temp_parent
        ) as temporary:
            scratch = Path(temporary)
            used: list[tuple[Track, Path, Path]] = []
            for index, track in enumerate(required, 1):
                print(
                    f"I: staging track {index}/{len(required)}: "
                    f"{track.series_id} {track.condition} {track.channel}",
                    file=sys.stderr,
                )
                source, bedgraph = stage_track(
                    track, arguments.track_root, scratch, arguments.chrom,
                    arguments.rscript, arguments.bigwig_exporter,
                )
                used.append((track, source, bedgraph))

            anchor_source: Path | None = None
            staged_tp73 = scratch / "tp73_anchor_evidence.parquet"
            if not arguments.profile_only:
                coverage: list[str] = []
                for track, _, bedgraph in used:
                    if track.channel not in {"tp73", "negative_control"}:
                        continue
                    column = f"supported_{track.channel}_{track.series_id}_{track.condition}"
                    coverage.extend(["--coverage", f"{column}={bedgraph}"])
                phase_start = time.monotonic()
                print(
                    "I: phase strict-evidence started: "
                    f"{len(coverage) // 2} TP73/control tracks",
                    file=sys.stderr,
                )
                anchor_arguments = (
                    ["--anchor-source", str(arguments.anchor_source)]
                    if arguments.anchor_source is not None
                    else [
                        "--anchor-plus", str(arguments.anchor_plus),
                        "--anchor-minus", str(arguments.anchor_minus),
                    ]
                )
                run([
                    str(arguments.anchor_builder),
                    *anchor_arguments,
                    *coverage,
                    "--output", str(staged_tp73),
                    "--chrom", arguments.chrom,
                    "--minimum-anchor-score", f"{arguments.minimum_anchor_score:.17g}",
                    "--threads", str(arguments.threads),
                    "--memory-limit", arguments.memory_limit,
                    "--duckdb", arguments.duckdb,
                ])
                print(
                    "I: phase strict-evidence finished: "
                    f"{time.monotonic() - phase_start:.1f} seconds",
                    file=sys.stderr,
                )
                anchor_source = staged_tp73

            if arguments.profile_output is not None:
                staged_profile = scratch / arguments.profile_output.name
                profile_tracks = [
                    (track, bedgraph) for track, _, bedgraph in used
                    if track.condition == "GFP" and
                    track.channel in {"h3k4me3", "input"}
                ]
                sql = profile_sql(
                    arguments, profile_tracks, staged_profile, anchor_source
                )
                run([arguments.duckdb, "-batch", ":memory:"], input_text=sql)
                promote_file(staged_profile, arguments.profile_output)
                write_sidecars(
                    arguments.profile_output, arguments, tracks, used, windows,
                    "profile_only" if arguments.profile_only else "full",
                )

            if not arguments.profile_only:
                staged_signal = scratch / "h3k4me3_anchor_signal.parquet"
                arguments.scratch_database = scratch / "signal.duckdb"
                (scratch / "duckdb_tmp").mkdir()
                phase_start = time.monotonic()
                signal_track_count = sum(
                    track.channel in {"h3k4me3", "input"}
                    for track, _, _ in used
                )
                print(
                    "I: phase window-signal started: "
                    f"{signal_track_count} H3K4me3/input tracks, "
                    f"{len(windows) + 1} windows",
                    file=sys.stderr,
                )
                sql = signal_sql(
                    arguments,
                    [(track, bedgraph) for track, _, bedgraph in used],
                    staged_tp73, staged_signal, windows,
                )
                run([
                    arguments.duckdb, "-batch", str(arguments.scratch_database)
                ], input_text=sql)
                print(
                    "I: phase window-signal finished: "
                    f"{time.monotonic() - phase_start:.1f} seconds",
                    file=sys.stderr,
                )
                promote_file(staged_signal, arguments.signal_output)
                promote_file(staged_tp73, arguments.tp73_output)
                for suffix in (".run_config.tsv", ".coverage_manifest.tsv"):
                    source = Path(str(staged_tp73) + suffix)
                    if source.exists():
                        promote_file(
                            source, Path(str(arguments.tp73_output) + suffix)
                        )
                write_sidecars(
                    arguments.signal_output, arguments, tracks, used, windows, "full"
                )
        print("I: H3K4me3 anchor signal build completed.", file=sys.stderr)
        return 0
    except (OSError, SignalBuildError, ValueError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3

"""Evaluate strict CUT&RUN coverage containment across motif-score thresholds.

Fragment BED rows each contribute unit depth.  bedGraph rows instead contribute
their fourth-column signal, preserving run-length encoded coverage depth.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from bisect import bisect_left
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


@dataclass(frozen=True)
class Motif:
    index: int
    chrom: str
    start: int
    end: int
    name: str
    score: float
    strand: str


@dataclass(frozen=True)
class CoverageInterval:
    index: int
    chrom: str
    start: int
    end: int
    name: str
    depth: float


@dataclass(frozen=True)
class CoverageComponent:
    index: int
    chrom: str
    start: int
    end: int
    member_names: tuple[str, ...]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Measure how motif scores predict strict immersion in merged "
            "CUT&RUN coverage."
        ),
        epilog=(
            "Overlapping and directly adjacent positive-coverage intervals "
            "are merged. A component supports [motif_start,motif_end) only "
            "when component_start < motif_start and component_end > motif_end."
        ),
    )
    parser.add_argument(
        "--motifs",
        required=True,
        type=Path,
        help="BED6 motif intervals; column 5 is a floating-point motif score",
    )
    parser.add_argument(
        "--coverage-bed",
        "--fragments",
        dest="coverage_bed",
        required=True,
        type=Path,
        help=(
            "BED3+ positive CUT&RUN coverage intervals; adjacent rows are "
            "merged (optionally gzip-compressed)"
        ),
    )
    parser.add_argument(
        "--coverage-format",
        choices=("auto", "fragments", "bedgraph"),
        default="auto",
        help=(
            "depth interpretation (default: auto; .bedGraph[.gz] uses column "
            "4, other BED files contribute one per row)"
        ),
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for audit tables and summary.json",
    )
    parser.add_argument(
        "--chrom",
        help="Analyze one chromosome (chr1 and 1 are treated as equivalent)",
    )
    parser.add_argument(
        "--sample-id",
        default="synthetic_tp73",
        help="Sample identifier recorded in summary.json",
    )
    parser.add_argument(
        "--score-mode",
        default="synthetic_score",
        help="Score provenance label recorded in summary.json",
    )
    return parser.parse_args()


def canonical_chromosome(chrom: str) -> str:
    chrom = chrom.strip()
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def resolve_coverage_format(requested: str, path: Path) -> str:
    if requested != "auto":
        return requested
    filename = path.name.lower()
    if filename.endswith(".gz"):
        filename = filename[:-3]
    return "bedgraph" if filename.endswith((".bedgraph", ".bdg")) else "fragments"


def coverage_depth_semantics(coverage_format: str) -> str:
    return (
        "column_4_signal"
        if coverage_format == "bedgraph"
        else "unit_per_fragment_interval"
    )


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def data_rows(path: Path) -> Iterable[tuple[int, list[str]]]:
    if not path.is_file():
        raise ValueError(f"input file not found: {path}")
    with open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            if line.startswith("track ") or line.startswith("browser "):
                continue
            fields = line.split("\t")
            if len(fields) == 1:
                raise ValueError(
                    f"{path}:{line_number}: BED fields must be tab-separated"
                )
            yield line_number, fields


def parse_coordinate(value: str, path: Path, line_number: int, field: str) -> int:
    try:
        coordinate = int(value)
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number}: invalid integer {field}: {value!r}"
        ) from error
    return coordinate


def validate_interval(
    path: Path, line_number: int, start: int, end: int
) -> None:
    if start < 0 or end <= start:
        raise ValueError(
            f"{path}:{line_number}: invalid BED interval [{start},{end})"
        )


def read_motifs(path: Path, selected_chrom: str | None) -> list[Motif]:
    motifs: list[Motif] = []
    for line_number, fields in data_rows(path):
        if len(fields) < 6:
            raise ValueError(f"{path}:{line_number}: motif BED requires 6 columns")
        chrom = canonical_chromosome(fields[0])
        if selected_chrom is not None and chrom != selected_chrom:
            continue
        start = parse_coordinate(fields[1], path, line_number, "start")
        end = parse_coordinate(fields[2], path, line_number, "end")
        validate_interval(path, line_number, start, end)
        try:
            score = float(fields[4])
        except ValueError as error:
            raise ValueError(
                f"{path}:{line_number}: invalid motif score: {fields[4]!r}"
            ) from error
        if not math.isfinite(score):
            raise ValueError(f"{path}:{line_number}: motif score must be finite")
        strand = fields[5]
        if strand not in ("+", "-"):
            raise ValueError(
                f"{path}:{line_number}: motif strand must be '+' or '-'"
            )
        motifs.append(
            Motif(
                index=len(motifs) + 1,
                chrom=chrom,
                start=start,
                end=end,
                name=fields[3],
                score=score,
                strand=strand,
            )
        )
    if not motifs:
        raise ValueError("no motif intervals remain after input filtering")
    lengths_by_chrom: dict[str, set[int]] = defaultdict(set)
    for motif in motifs:
        lengths_by_chrom[motif.chrom].add(motif.end - motif.start)
    inconsistent = [
        chrom for chrom, lengths in lengths_by_chrom.items() if len(lengths) != 1
    ]
    if inconsistent:
        raise ValueError(
            "motif BED must use one motif length per chromosome; inconsistent: "
            + ", ".join(sorted(inconsistent))
        )
    return motifs


def read_coverage_intervals(
    path: Path, selected_chrom: str | None, coverage_format: str
) -> list[CoverageInterval]:
    intervals: list[CoverageInterval] = []
    for line_number, fields in data_rows(path):
        if len(fields) < 3:
            raise ValueError(
                f"{path}:{line_number}: CUT&RUN fragment BED requires 3 columns"
            )
        chrom = canonical_chromosome(fields[0])
        if selected_chrom is not None and chrom != selected_chrom:
            continue
        start = parse_coordinate(fields[1], path, line_number, "start")
        end = parse_coordinate(fields[2], path, line_number, "end")
        validate_interval(path, line_number, start, end)
        if coverage_format == "bedgraph":
            if len(fields) < 4:
                raise ValueError(
                    f"{path}:{line_number}: bedGraph requires depth in column 4"
                )
            try:
                depth = float(fields[3])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: invalid bedGraph depth: {fields[3]!r}"
                ) from error
            if not math.isfinite(depth) or depth <= 0:
                raise ValueError(
                    f"{path}:{line_number}: bedGraph depth must be finite and positive"
                )
            name = f"run_{line_number}"
        else:
            depth = 1.0
            name = (
                fields[3]
                if len(fields) >= 4 and fields[3]
                else f"interval_{line_number}"
            )
        intervals.append(
            CoverageInterval(
                index=len(intervals) + 1,
                chrom=chrom,
                start=start,
                end=end,
                name=name,
                depth=depth,
            )
        )
    if not intervals:
        raise ValueError("no CUT&RUN coverage intervals remain after input filtering")
    return intervals


def merge_coverage_intervals(
    intervals: list[CoverageInterval],
) -> list[CoverageComponent]:
    by_chrom: dict[str, list[CoverageInterval]] = defaultdict(list)
    for interval in intervals:
        by_chrom[interval.chrom].append(interval)

    components: list[CoverageComponent] = []
    for chrom in sorted(by_chrom):
        chromosome_intervals = sorted(
            by_chrom[chrom],
            key=lambda interval: (interval.start, interval.end, interval.index),
        )
        component_start = chromosome_intervals[0].start
        component_end = chromosome_intervals[0].end
        member_names = [chromosome_intervals[0].name]
        for interval in chromosome_intervals[1:]:
            if interval.start <= component_end:
                component_end = max(component_end, interval.end)
                member_names.append(interval.name)
                continue
            components.append(
                CoverageComponent(
                    index=len(components) + 1,
                    chrom=chrom,
                    start=component_start,
                    end=component_end,
                    member_names=tuple(member_names),
                )
            )
            component_start = interval.start
            component_end = interval.end
            member_names = [interval.name]
        components.append(
            CoverageComponent(
                index=len(components) + 1,
                chrom=chrom,
                start=component_start,
                end=component_end,
                member_names=tuple(member_names),
            )
        )
    return components


def ordinary_max_depth(
    motif: Motif, intervals: Iterable[CoverageInterval]
) -> float:
    # Fragment BED rows carry depth 1; bedGraph runs carry their encoded depth.
    deltas: dict[int, float] = defaultdict(float)
    for interval in intervals:
        overlap_start = max(motif.start, interval.start)
        overlap_end = min(motif.end, interval.end)
        if overlap_start < overlap_end:
            deltas[overlap_start] += interval.depth
            deltas[overlap_end] -= interval.depth
    depth = 0.0
    maximum = 0.0
    for position in sorted(deltas):
        depth += deltas[position]
        maximum = max(maximum, depth)
    return maximum


def evaluate_containment(
    motifs: list[Motif],
    intervals: list[CoverageInterval],
    components: list[CoverageComponent],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    motifs_by_chrom: dict[str, list[Motif]] = defaultdict(list)
    intervals_by_chrom: dict[str, list[CoverageInterval]] = defaultdict(list)
    components_by_chrom: dict[str, list[CoverageComponent]] = defaultdict(list)
    for motif in motifs:
        motifs_by_chrom[motif.chrom].append(motif)
    for interval in intervals:
        intervals_by_chrom[interval.chrom].append(interval)
    for component in components:
        components_by_chrom[component.chrom].append(component)

    motif_component: dict[int, CoverageComponent | None] = {}
    component_scores: dict[int, list[float]] = defaultdict(list)

    for chrom, chromosome_motifs in motifs_by_chrom.items():
        chromosome_components = sorted(
            components_by_chrom.get(chrom, []),
            key=lambda component: (component.start, component.end),
        )
        component_starts = [component.start for component in chromosome_components]
        for motif in chromosome_motifs:
            component_index = bisect_left(component_starts, motif.start) - 1
            component = (
                chromosome_components[component_index]
                if component_index >= 0
                else None
            )
            if component is None or component.end <= motif.end:
                component = None
            motif_component[motif.index] = component
            if component is not None:
                component_scores[component.index].append(motif.score)

    motif_rows: list[dict[str, object]] = []
    for motif in sorted(motifs, key=lambda item: item.index):
        component = motif_component[motif.index]
        maximum_depth = ordinary_max_depth(
            motif, intervals_by_chrom.get(motif.chrom, [])
        )
        immersed = component is not None
        motif_rows.append(
            {
                "motif_index": motif.index,
                "chrom": motif.chrom,
                "start": motif.start,
                "end": motif.end,
                "motif_name": motif.name,
                "score": motif.score,
                "strand": motif.strand,
                "ordinary_max_depth": maximum_depth,
                "effective_max_depth": maximum_depth if immersed else 0.0,
                "supported": immersed,
                "coverage_component_start": component.start if component else None,
                "coverage_component_end": component.end if component else None,
                "coverage_component_members": (
                    ",".join(component.member_names) if component else None
                ),
            }
        )

    component_rows: list[dict[str, object]] = []
    for component in sorted(components, key=lambda item: item.index):
        scores = component_scores.get(component.index, [])
        component_rows.append(
            {
                "component_index": component.index,
                "chrom": component.chrom,
                "start": component.start,
                "end": component.end,
                "source_interval_count": len(component.member_names),
                "source_intervals": ",".join(component.member_names),
                "immersed_motif_count": len(scores),
                "best_immersed_motif_score": max(scores) if scores else None,
            }
        )
    return motif_rows, component_rows


def divide(numerator: int | float, denominator: int | float) -> float | None:
    return numerator / denominator if denominator else None


def threshold_curve(
    motif_rows: list[dict[str, object]],
    component_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    thresholds = sorted({float(row["score"]) for row in motif_rows}, reverse=True)
    total_positive = sum(bool(row["supported"]) for row in motif_rows)
    total_negative = len(motif_rows) - total_positive
    prevalence = divide(total_positive, len(motif_rows))
    result: list[dict[str, object]] = []

    for threshold in thresholds:
        selected = [row for row in motif_rows if float(row["score"]) >= threshold]
        true_positive = sum(bool(row["supported"]) for row in selected)
        false_positive = len(selected) - true_positive
        selected_depth = sum(float(row["effective_max_depth"]) for row in selected)
        precision = divide(true_positive, len(selected))
        recall = divide(true_positive, total_positive)
        false_positive_rate = divide(false_positive, total_negative)
        component_hits = sum(
            row["best_immersed_motif_score"] is not None
            and float(row["best_immersed_motif_score"]) >= threshold
            for row in component_rows
        )
        f1 = (
            2.0 * precision * recall / (precision + recall)
            if precision is not None
            and recall is not None
            and precision + recall > 0
            else 0.0
        )
        result.append(
            {
                "threshold": threshold,
                "selected_motifs": len(selected),
                "supported_selected_motifs": true_positive,
                "unsupported_selected_motifs": false_positive,
                "motif_precision": precision,
                "motif_recall": recall,
                "motif_false_positive_rate": false_positive_rate,
                "coverage_component_recall": divide(
                    component_hits, len(component_rows)
                ),
                "mean_effective_depth": divide(selected_depth, len(selected)),
                "precision_enrichment": (
                    precision / prevalence
                    if precision is not None and prevalence not in (None, 0.0)
                    else None
                ),
                "f1": f1,
                "youden_j": (
                    recall - false_positive_rate
                    if recall is not None and false_positive_rate is not None
                    else None
                ),
            }
        )
    return result


def roc_auc(motif_rows: list[dict[str, object]]) -> float | None:
    ranked = sorted(
        (float(row["score"]), bool(row["supported"])) for row in motif_rows
    )
    positives = sum(label for _, label in ranked)
    negatives = len(ranked) - positives
    if positives == 0 or negatives == 0:
        return None

    positive_rank_sum = 0.0
    index = 0
    while index < len(ranked):
        group_end = index + 1
        while group_end < len(ranked) and ranked[group_end][0] == ranked[index][0]:
            group_end += 1
        average_rank = ((index + 1) + group_end) / 2.0
        positive_rank_sum += average_rank * sum(
            label for _, label in ranked[index:group_end]
        )
        index = group_end

    return (
        positive_rank_sum - positives * (positives + 1) / 2.0
    ) / (positives * negatives)


def average_precision(curve: list[dict[str, object]]) -> float | None:
    if not curve or curve[0]["motif_recall"] is None:
        return None
    area = 0.0
    previous_recall = 0.0
    for row in curve:
        recall = float(row["motif_recall"])
        precision = float(row["motif_precision"])
        area += max(0.0, recall - previous_recall) * precision
        previous_recall = recall
    return area


def best_threshold(
    curve: list[dict[str, object]], metric: str
) -> dict[str, object] | None:
    eligible = [row for row in curve if row[metric] is not None]
    if not eligible:
        return None
    return max(eligible, key=lambda row: (float(row[metric]), float(row["threshold"])))


def format_tsv_value(value: object) -> str:
    if value is None:
        # An explicit marker keeps optional terminal columns visible and avoids
        # ambiguous trailing delimiters in tables inspected outside Python.
        return "NA"
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot write an empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: format_tsv_value(value) for key, value in row.items()})


def build_summary(
    arguments: argparse.Namespace,
    motifs: list[Motif],
    intervals: list[CoverageInterval],
    components: list[CoverageComponent],
    motif_rows: list[dict[str, object]],
    component_rows: list[dict[str, object]],
    curve: list[dict[str, object]],
) -> dict[str, object]:
    supported = sum(bool(row["supported"]) for row in motif_rows)
    youden = best_threshold(curve, "youden_j")
    f1 = best_threshold(curve, "f1")
    return {
        "sample_id": arguments.sample_id,
        "score_mode": arguments.score_mode,
        "coverage_format": arguments.coverage_format,
        "coverage_depth_semantics": coverage_depth_semantics(
            arguments.coverage_format
        ),
        "containment_rule": "component_start < motif_start and component_end > motif_end",
        "coverage_merge_rule": "overlapping or directly adjacent intervals",
        "coordinate_mode": "BED 0-based half-open",
        "n_motifs": len(motifs),
        "n_supported_motifs": supported,
        "n_coverage_intervals": len(intervals),
        "n_coverage_components": len(components),
        "n_components_immersing_any_motif": sum(
            int(row["immersed_motif_count"]) > 0 for row in component_rows
        ),
        "motif_support_prevalence": divide(supported, len(motifs)),
        "roc_auc": roc_auc(motif_rows),
        "average_precision": average_precision(curve),
        "youden_threshold": youden["threshold"] if youden else None,
        "youden_j": youden["youden_j"] if youden else None,
        "f1_threshold": f1["threshold"] if f1 else None,
        "f1": f1["f1"] if f1 else None,
    }


def main() -> int:
    arguments = parse_arguments()
    arguments.coverage_format = resolve_coverage_format(
        arguments.coverage_format, arguments.coverage_bed
    )
    selected_chrom = (
        canonical_chromosome(arguments.chrom) if arguments.chrom else None
    )
    motifs = read_motifs(arguments.motifs, selected_chrom)
    intervals = read_coverage_intervals(
        arguments.coverage_bed, selected_chrom, arguments.coverage_format
    )
    components = merge_coverage_intervals(intervals)
    motif_rows, component_rows = evaluate_containment(
        motifs, intervals, components
    )
    curve = threshold_curve(motif_rows, component_rows)
    summary = build_summary(
        arguments,
        motifs,
        intervals,
        components,
        motif_rows,
        component_rows,
        curve,
    )

    output_dir = arguments.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    output_files = [
        output_dir / "motif_evidence.tsv",
        output_dir / "coverage_component_evidence.tsv",
        output_dir / "threshold_curve.tsv",
        output_dir / "summary.json",
    ]
    existing = [path for path in output_files if path.exists()]
    if existing:
        raise ValueError(
            "refusing to replace existing output(s): "
            + ", ".join(str(path) for path in existing)
        )

    write_tsv(output_files[0], motif_rows)
    write_tsv(output_files[1], component_rows)
    write_tsv(output_files[2], curve)
    with output_files[3].open("x", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        raise SystemExit(f"E: {error}") from error

#!/usr/bin/env python3
"""Render the preserved informative-threshold distribution as SVG."""

from __future__ import annotations

import argparse
import csv
from html import escape
from pathlib import Path


EXPECTED_COLUMNS = {
    "threshold",
    "original_v1_motifs",
    "provisional_motifs",
    "negative_optimum_higher_auc",
    "auc_gain_gt_0_001",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Render the aggregate JASPAR 2026 informative-threshold counts as "
            "a dependency-free SVG figure."
        )
    )
    parser.add_argument("input_tsv", type=Path, help="aggregate threshold-count TSV")
    parser.add_argument("output_svg", type=Path, help="SVG path to create")
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, int | None]]:
    with path.open(newline="", encoding="ascii") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != EXPECTED_COLUMNS:
            raise ValueError(
                "input columns must be exactly: " + ", ".join(sorted(EXPECTED_COLUMNS))
            )

        rows: list[dict[str, int | None]] = []
        for source in reader:
            threshold = None if source["threshold"] == "NA" else int(source["threshold"])
            rows.append(
                {
                    "threshold": threshold,
                    "original": int(source["original_v1_motifs"]),
                    "provisional": int(source["provisional_motifs"]),
                    "improved": int(source["negative_optimum_higher_auc"]),
                    "material": int(source["auc_gain_gt_0_001"]),
                }
            )

    if sum(int(row["original"]) for row in rows) != 2632:
        raise ValueError("original_v1_motifs must sum to 2,632")
    if sum(int(row["provisional"]) for row in rows) != 2632:
        raise ValueError("provisional_motifs must sum to 2,632")
    if sum(int(row["improved"]) for row in rows) != 295:
        raise ValueError("negative_optimum_higher_auc must sum to 295")
    if sum(int(row["material"]) for row in rows) != 220:
        raise ValueError("auc_gain_gt_0_001 must sum to 220")

    numeric = [int(row["threshold"]) for row in rows if row["threshold"] is not None]
    if numeric != list(range(-20, 15)):
        raise ValueError("numeric threshold rows must cover every integer from -20 to 14")
    return rows


def text(x: float, y: float, value: str, css_class: str, anchor: str = "start",
         extra: str = "") -> str:
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" class="{css_class}" '
        f'text-anchor="{anchor}" {extra}>{escape(value)}</text>'
    )


def line(x1: float, y1: float, x2: float, y2: float, css_class: str) -> str:
    return (
        f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" '
        f'y2="{y2:.2f}" class="{css_class}" />'
    )


def rect(x: float, y: float, width: float, height: float, css_class: str) -> str:
    return (
        f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(0.0, width):.2f}" '
        f'height="{max(0.0, height):.2f}" class="{css_class}" />'
    )


def render(rows: list[dict[str, int | None]]) -> str:
    width = 1200
    height = 1180
    numeric = [row for row in rows if row["threshold"] is not None]
    negative = [row for row in numeric if int(row["threshold"]) < 0]

    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" role="img" aria-labelledby="title description">',
        '<title id="title">JASPAR 2026 informative-threshold distribution</title>',
        '<desc id="description">Provisional distribution for 2,632 non-TP73 motifs. '
        'The lower panel details the 295 motifs whose negative threshold improved held-out '
        'macro ROC AUC over threshold zero.</desc>',
        """<style>
        .background { fill: #ffffff; }
        .title { fill: #111111; font: 500 28px system-ui, sans-serif; }
        .subtitle { fill: #4d4d4d; font: 400 17px system-ui, sans-serif; }
        .panel { fill: #111111; font: 500 19px system-ui, sans-serif; }
        .axis { fill: #4d4d4d; font: 400 15px system-ui, sans-serif; }
        .value { fill: #111111; font: 500 15px system-ui, sans-serif; }
        .grid { stroke: #d9d9d9; stroke-width: 1; shape-rendering: crispEdges; }
        .baseline { stroke: #737373; stroke-width: 1; shape-rendering: crispEdges; }
        .negative { fill: #0072b2; }
        .nonnegative { fill: #009e73; }
        .material { fill: #d55e00; }
        .original { fill: none; stroke: #111111; stroke-width: 2; }
        .legend { fill: #333333; font: 400 15px system-ui, sans-serif; }
        </style>""",
        rect(0, 0, width, height, "background"),
        text(60, 42, "JASPAR 2026 informative-threshold distribution", "title"),
        text(
            60,
            69,
            "Chromosome-1 TP73-context calibration; provisional negative-threshold sensitivity update",
            "subtitle",
        ),
    ]

    legend = [
        (60, "negative", "Selected below zero"),
        (290, "nonnegative", "Retained at zero or above"),
        (570, "material", "ROC AUC gain > 0.001"),
    ]
    for x, css_class, label in legend:
        swatch_height = 5 if css_class == "material" else 12
        out.append(rect(x, 86 + (12 - swatch_height) / 2, 24, swatch_height, css_class))
        out.append(text(x + 33, 98, label, "legend"))
    out.append(rect(860, 87, 24, 14, "original"))
    out.append(text(893, 99, "Original zero bin", "legend"))

    out.append(text(60, 137, "Provisional all-motif distribution", "panel"))
    plot_x = 80
    plot_y = 160
    plot_width = 1080
    plot_height = 300
    max_y = 500
    for tick in range(0, max_y + 1, 100):
        y = plot_y + plot_height - tick / max_y * plot_height
        out.append(line(plot_x, y, plot_x + plot_width, y, "grid"))
        out.append(text(plot_x - 10, y + 5, str(tick), "axis", "end"))

    step = plot_width / len(numeric)
    bar_width = step * 0.72
    threshold_to_x: dict[int, float] = {}
    for index, row in enumerate(numeric):
        threshold = int(row["threshold"])
        count = int(row["provisional"])
        x = plot_x + index * step + (step - bar_width) / 2
        threshold_to_x[threshold] = x
        bar_height = count / max_y * plot_height
        css_class = "negative" if threshold < 0 else "nonnegative"
        out.append(rect(x, plot_y + plot_height - bar_height, bar_width, bar_height, css_class))

    zero_x = threshold_to_x[0]
    original_height = 438 / max_y * plot_height
    out.append(rect(zero_x, plot_y + plot_height - original_height,
                    bar_width, original_height, "original"))

    for threshold in (-20, -15, -10, -5, 0, 5, 10, 14):
        center = threshold_to_x[threshold] + bar_width / 2
        out.append(line(center, plot_y + plot_height, center, plot_y + plot_height + 6,
                        "baseline"))
        out.append(text(center, plot_y + plot_height + 26, str(threshold), "axis", "middle"))
    out.append(line(plot_x, plot_y + plot_height, plot_x + plot_width,
                    plot_y + plot_height, "baseline"))
    out.append(text(plot_x + plot_width / 2, 505, "Informative threshold", "axis", "middle"))
    out.append(text(20, plot_y + plot_height / 2, "Number of motifs", "axis", "middle",
                    f'transform="rotate(-90 20 {plot_y + plot_height / 2:.2f})"'))

    for threshold, count in ((0, 143), (6, 454)):
        center = threshold_to_x[threshold] + bar_width / 2
        y = plot_y + plot_height - count / max_y * plot_height - 8
        out.append(text(center, y, str(count), "value", "middle"))
    out.append(text(zero_x + bar_width / 2 + 8,
                    plot_y + plot_height - original_height - 8,
                    "438 before", "value"))

    out.append(text(60, 555, "Detail: 295 motifs shifted below zero", "panel"))
    detail_x = 80
    detail_y = 580
    detail_width = 900
    row_height = 25
    detail_height = len(negative) * row_height
    max_x = 70
    for tick in (0, 20, 40, 60):
        x = detail_x + tick / max_x * detail_width
        out.append(line(x, detail_y, x, detail_y + detail_height, "grid"))
        out.append(text(x, detail_y + detail_height + 25, str(tick), "axis", "middle"))

    for index, row in enumerate(negative):
        threshold = int(row["threshold"])
        improved = int(row["improved"])
        material = int(row["material"])
        center_y = detail_y + index * row_height + row_height / 2
        out.append(text(detail_x - 12, center_y + 5, str(threshold), "axis", "end"))
        if improved:
            out.append(rect(detail_x, center_y - 7, improved / max_x * detail_width,
                            14, "negative"))
        if material:
            out.append(rect(detail_x, center_y - 2, material / max_x * detail_width,
                            4, "material"))
        out.append(text(detail_x + detail_width + 18, center_y + 5,
                        f"{improved} / {material}", "value"))

    out.append(line(detail_x, detail_y + detail_height,
                    detail_x + detail_width, detail_y + detail_height, "baseline"))
    out.append(text(detail_x + detail_width / 2, 1134,
                    "Number of motifs (any gain / gain > 0.001)", "axis", "middle"))
    out.append(text(60, 1165,
                    "Seventeen motifs without an evaluable recommendation are excluded from the bars.",
                    "subtitle"))
    out.append("</svg>")
    return "\n".join(out) + "\n"


def main() -> None:
    args = parse_arguments()
    rows = read_rows(args.input_tsv)
    args.output_svg.parent.mkdir(parents=True, exist_ok=True)
    args.output_svg.write_text(render(rows), encoding="ascii")
    print(f"Wrote {args.output_svg}")


if __name__ == "__main__":
    main()

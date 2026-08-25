#!/usr/bin/env python3
"""Finalize exact TP73 motif-distance curves, periodicity, and peak matches."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Iterable


class ExactDistanceSummaryError(RuntimeError):
    """Raised when exact-distance sufficient statistics cannot be finalized."""


Z_975 = 1.959963984540054

TABLE_SCHEMAS: dict[str, tuple[tuple[str, str], ...]] = {
    "distance_response.parquet": (
        ("motif_id", "VARCHAR"), ("motif_name", "VARCHAR"),
        ("anchor_motif_id", "VARCHAR"), ("distance_frame", "VARCHAR"),
        ("relative_orientation", "VARCHAR"), ("center_offset_twice", "BIGINT"),
        ("center_offset_bp", "DOUBLE"), ("isoform", "VARCHAR"),
        ("source_score_floor", "DOUBLE"), ("positive_threshold", "DOUBLE"),
        ("numerator", "DOUBLE"), ("denominator", "DOUBLE"),
        ("positive_anti_discordant", "BIGINT"),
        ("positive_control_discordant", "BIGINT"),
        ("negative_anti_discordant", "BIGINT"),
        ("negative_control_discordant", "BIGINT"),
        ("genomic_blocks", "BIGINT"), ("transaction_scope", "VARCHAR"),
        ("conditioning_scope", "VARCHAR"), ("distance_metric", "VARCHAR"),
        ("overlap_policy", "VARCHAR"), ("transaction_count", "BIGINT"),
        ("context_flank_bp", "BIGINT"), ("block_size_bp", "BIGINT"),
        ("maximum_anchor_span_bp", "BIGINT"),
        ("maximum_neighbor_span_bp", "BIGINT"),
        ("anchors_source_present", "BIGINT"), ("itemset_count", "BIGINT"),
        ("anchors_intermediate", "BIGINT"), ("anchors_negative", "BIGINT"),
        ("physical_locus_pair_count", "BIGINT"), ("source_support", "DOUBLE"),
        ("support", "DOUBLE"), ("adjusted_odds_ratio", "DOUBLE"),
        ("jackknife_blocks", "BIGINT"), ("jackknife_se", "DOUBLE"),
        ("saos2_or", "DOUBLE"), ("skmel29_2_or", "DOUBLE"),
        ("adjusted_log_odds", "DOUBLE"), ("z_value", "DOUBLE"),
        ("p_value", "DOUBLE"), ("confidence_interval_95_lower", "DOUBLE"),
        ("confidence_interval_95_upper", "DOUBLE"),
        ("evaluation_status", "VARCHAR"),
        ("q_value_bh_declared_offsets", "DOUBLE"),
        ("series_direction_consistent", "BOOLEAN"),
    ),
    "isoform_contrast.parquet": (
        ("motif_id", "VARCHAR"), ("motif_name", "VARCHAR"),
        ("distance_frame", "VARCHAR"), ("relative_orientation", "VARCHAR"),
        ("center_offset_twice", "BIGINT"), ("center_offset_bp", "DOUBLE"),
        ("ta_adjusted_odds_ratio", "DOUBLE"),
        ("dn_adjusted_odds_ratio", "DOUBLE"),
        ("ta_vs_dn_log_odds_difference", "DOUBLE"),
        ("ta_vs_dn_odds_ratio_ratio", "DOUBLE"),
        ("genomic_blocks", "BIGINT"), ("jackknife_blocks", "BIGINT"),
        ("jackknife_se", "DOUBLE"),
        ("saos2_ta_vs_dn_odds_ratio_ratio", "DOUBLE"),
        ("skmel29_2_ta_vs_dn_odds_ratio_ratio", "DOUBLE"),
        ("z_value", "DOUBLE"), ("p_value", "DOUBLE"),
        ("confidence_interval_95_lower", "DOUBLE"),
        ("confidence_interval_95_upper", "DOUBLE"),
        ("evaluation_status", "VARCHAR"),
        ("q_value_bh_declared_offsets", "DOUBLE"),
        ("series_contrast_direction_consistent", "BOOLEAN"),
    ),
    "periodicity.parquet": (
        ("motif_id", "VARCHAR"), ("motif_name", "VARCHAR"),
        ("distance_frame", "VARCHAR"), ("relative_orientation", "VARCHAR"),
        ("isoform", "VARCHAR"), ("fixed_period_bp", "DOUBLE"),
        ("estimable_offsets", "BIGINT"), ("minimum_offset_bp", "DOUBLE"),
        ("maximum_offset_bp", "DOUBLE"), ("sine_coefficient", "DOUBLE"),
        ("cosine_coefficient", "DOUBLE"),
        ("harmonic_amplitude_log_odds", "DOUBLE"),
        ("harmonic_phase_peak_bp_modulo_period", "DOUBLE"),
        ("baseline_weighted_rss", "DOUBLE"),
        ("harmonic_weighted_rss", "DOUBLE"),
        ("weighted_rss_fraction_reduced", "DOUBLE"),
        ("inference_status", "VARCHAR"),
    ),
    "periodicity_isoform_contrast.parquet": (
        ("motif_id", "VARCHAR"), ("motif_name", "VARCHAR"),
        ("distance_frame", "VARCHAR"), ("relative_orientation", "VARCHAR"),
        ("fixed_period_bp", "DOUBLE"),
        ("ta_amplitude_log_odds", "DOUBLE"),
        ("dn_amplitude_log_odds", "DOUBLE"),
        ("dn_minus_ta_amplitude_log_odds", "DOUBLE"),
        ("ta_phase_peak_bp_modulo_period", "DOUBLE"),
        ("dn_phase_peak_bp_modulo_period", "DOUBLE"),
        ("dn_minus_ta_phase_shift_bp_wrapped", "DOUBLE"),
        ("inference_status", "VARCHAR"),
    ),
    "peak.parquet": (
        ("motif_id", "VARCHAR"), ("motif_name", "VARCHAR"),
        ("distance_frame", "VARCHAR"), ("relative_orientation", "VARCHAR"),
        ("isoform", "VARCHAR"), ("peak_direction", "VARCHAR"),
        ("peak_center_bp", "DOUBLE"), ("smoothed_log_odds", "DOUBLE"),
        ("prominence_log_odds", "DOUBLE"),
        ("raw_adjusted_odds_ratio", "DOUBLE"),
        ("raw_q_value_bh_declared_offsets", "DOUBLE"),
        ("itemset_count", "BIGINT"), ("smoothing_sigma_bp", "DOUBLE"),
        ("peak_status", "VARCHAR"), ("peak_index", "BIGINT"),
    ),
    "peak_match.parquet": (
        ("motif_id", "VARCHAR"), ("motif_name", "VARCHAR"),
        ("distance_frame", "VARCHAR"), ("relative_orientation", "VARCHAR"),
        ("peak_direction", "VARCHAR"), ("ta_peak_index", "BIGINT"),
        ("dn_peak_index", "BIGINT"), ("ta_peak_center_bp", "DOUBLE"),
        ("dn_peak_center_bp", "DOUBLE"),
        ("dn_minus_ta_peak_shift_bp", "DOUBLE"),
        ("group_shift_hypothesis_bp", "DOUBLE"),
        ("match_residual_to_group_shift_bp", "DOUBLE"),
        ("matched_peaks_in_group", "BIGINT"),
        ("group_shift_nearest_helical_turns", "BIGINT"),
        ("group_shift_helical_residual_bp", "DOUBLE"),
        ("same_side_of_anchor", "BOOLEAN"),
        ("dn_peak_is_closer_to_anchor", "BOOLEAN"),
        ("nearest_helical_turns", "BIGINT"),
        ("helical_shift_residual_bp", "DOUBLE"),
        ("fixed_period_bp", "DOUBLE"), ("match_status", "VARCHAR"),
    ),
    "peak_pair_match.parquet": (
        ("motif_id", "VARCHAR"), ("motif_name", "VARCHAR"),
        ("distance_frame", "VARCHAR"), ("relative_orientation", "VARCHAR"),
        ("peak_direction", "VARCHAR"), ("ta_left_peak_bp", "DOUBLE"),
        ("ta_right_peak_bp", "DOUBLE"), ("dn_left_peak_bp", "DOUBLE"),
        ("dn_right_peak_bp", "DOUBLE"), ("ta_peak_spacing_bp", "DOUBLE"),
        ("dn_peak_spacing_bp", "DOUBLE"),
        ("spacing_difference_dn_minus_ta_bp", "DOUBLE"),
        ("ta_nearest_helical_turns", "BIGINT"),
        ("ta_helical_spacing_residual_bp", "DOUBLE"),
        ("dn_nearest_helical_turns", "BIGINT"),
        ("dn_helical_spacing_residual_bp", "DOUBLE"),
        ("ta_peak_midpoint_bp", "DOUBLE"), ("dn_peak_midpoint_bp", "DOUBLE"),
        ("dn_minus_ta_midpoint_shift_bp", "DOUBLE"),
        ("group_shift_hypothesis_bp", "DOUBLE"),
        ("dn_midpoint_is_closer_to_anchor", "BOOLEAN"),
        ("midpoint_shift_nearest_helical_turns", "BIGINT"),
        ("midpoint_shift_helical_residual_bp", "DOUBLE"),
        ("fixed_period_bp", "DOUBLE"), ("pair_status", "VARCHAR"),
    ),
}


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_path_list(paths: Iterable[Path]) -> str:
    return "[" + ",".join(sql_string(path) for path in paths) + "]"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_duckdb(executable: str, database: Path | str, sql: str) -> None:
    process = subprocess.run(
        [executable, "-batch", "-bail", str(database)],
        input=sql,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ExactDistanceSummaryError(
            process.stderr.strip() or "DuckDB exact-distance finalization failed"
        )


def executable_version(executable: str) -> str:
    process = subprocess.run(
        [executable, "--version"],
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise ExactDistanceSummaryError(
            process.stderr.strip() or f"could not identify {executable}"
        )
    return process.stdout.strip()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bh_adjust(
    rows: list[dict[str, Any]], group_fields: tuple[str, ...]
) -> None:
    groups: dict[tuple[Any, ...], list[dict[str, Any]]] = {}
    for row in rows:
        p_value = row.get("p_value")
        if row.get("evaluation_status") == "ok" and p_value is not None:
            groups.setdefault(tuple(row[name] for name in group_fields), []).append(row)
    for group in groups.values():
        ordered = sorted(group, key=lambda item: item["p_value"])
        next_value = 1.0
        family_size = len(ordered)
        for rank in range(family_size, 0, -1):
            row = ordered[rank - 1]
            next_value = min(next_value, row["p_value"] * family_size / rank)
            row["q_value_bh_declared_offsets"] = min(1.0, next_value)


def parse_optional_float(value: str) -> float | None:
    return None if value == "NA" else float(value)


def parse_response_rows(path: Path, minimum_blocks: int) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    float_fields = {
        "center_offset_bp",
        "source_score_floor",
        "positive_threshold",
        "numerator",
        "denominator",
        "adjusted_odds_ratio",
        "jackknife_se",
        "saos2_or",
        "skmel29_2_or",
        "support",
        "source_support",
    }
    int_fields = {
        "center_offset_twice",
        "positive_anti_discordant",
        "positive_control_discordant",
        "negative_anti_discordant",
        "negative_control_discordant",
        "genomic_blocks",
        "jackknife_blocks",
        "transaction_count",
        "context_flank_bp",
        "block_size_bp",
        "maximum_anchor_span_bp",
        "maximum_neighbor_span_bp",
        "anchors_source_present",
        "itemset_count",
        "anchors_intermediate",
        "anchors_negative",
        "physical_locus_pair_count",
    }
    for raw in read_tsv(path):
        row: dict[str, Any] = dict(raw)
        for name in float_fields:
            row[name] = parse_optional_float(raw[name])
        for name in int_fields:
            row[name] = int(raw[name])
        odds_ratio = row["adjusted_odds_ratio"]
        standard_error = row["jackknife_se"]
        estimable = (
            odds_ratio is not None
            and odds_ratio > 0
            and standard_error is not None
            and standard_error > 0
            and row["jackknife_blocks"] == row["genomic_blocks"]
            and row["jackknife_blocks"] >= minimum_blocks
        )
        if estimable:
            log_or = math.log(odds_ratio)
            z_value = log_or / standard_error
            row["adjusted_log_odds"] = log_or
            row["z_value"] = z_value
            row["p_value"] = math.erfc(abs(z_value) / math.sqrt(2.0))
            row["confidence_interval_95_lower"] = math.exp(
                log_or - Z_975 * standard_error
            )
            row["confidence_interval_95_upper"] = math.exp(
                log_or + Z_975 * standard_error
            )
            row["evaluation_status"] = "ok"
        else:
            row["adjusted_log_odds"] = None
            row["z_value"] = None
            row["p_value"] = None
            row["confidence_interval_95_lower"] = None
            row["confidence_interval_95_upper"] = None
            row["evaluation_status"] = "not_estimable"
        row["q_value_bh_declared_offsets"] = None
        series_values = (row["saos2_or"], row["skmel29_2_or"])
        if all(value is not None and value > 0 for value in series_values):
            row["series_direction_consistent"] = (
                all(value > 1 for value in series_values)
                or all(value < 1 for value in series_values)
            )
        else:
            row["series_direction_consistent"] = None
        result.append(row)
    bh_adjust(
        result,
        ("distance_frame", "relative_orientation", "isoform"),
    )
    return result


def parse_contrast_rows(path: Path, minimum_blocks: int) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    float_fields = {
        "center_offset_bp",
        "ta_adjusted_odds_ratio",
        "dn_adjusted_odds_ratio",
        "ta_vs_dn_log_odds_difference",
        "ta_vs_dn_odds_ratio_ratio",
        "jackknife_se",
        "saos2_ta_vs_dn_odds_ratio_ratio",
        "skmel29_2_ta_vs_dn_odds_ratio_ratio",
    }
    int_fields = {"center_offset_twice", "genomic_blocks", "jackknife_blocks"}
    for raw in read_tsv(path):
        row: dict[str, Any] = dict(raw)
        for name in float_fields:
            row[name] = parse_optional_float(raw[name])
        for name in int_fields:
            row[name] = int(raw[name])
        difference = row["ta_vs_dn_log_odds_difference"]
        standard_error = row["jackknife_se"]
        estimable = (
            difference is not None
            and standard_error is not None
            and standard_error > 0
            and row["jackknife_blocks"] == row["genomic_blocks"]
            and row["jackknife_blocks"] >= minimum_blocks
        )
        if estimable:
            z_value = difference / standard_error
            row["z_value"] = z_value
            row["p_value"] = math.erfc(abs(z_value) / math.sqrt(2.0))
            row["confidence_interval_95_lower"] = math.exp(
                difference - Z_975 * standard_error
            )
            row["confidence_interval_95_upper"] = math.exp(
                difference + Z_975 * standard_error
            )
            row["evaluation_status"] = "ok"
        else:
            row["z_value"] = None
            row["p_value"] = None
            row["confidence_interval_95_lower"] = None
            row["confidence_interval_95_upper"] = None
            row["evaluation_status"] = "not_estimable"
        row["q_value_bh_declared_offsets"] = None
        series_values = (
            row["saos2_ta_vs_dn_odds_ratio_ratio"],
            row["skmel29_2_ta_vs_dn_odds_ratio_ratio"],
        )
        if all(value is not None and value > 0 for value in series_values):
            row["series_contrast_direction_consistent"] = (
                all(value > 1 for value in series_values)
                or all(value < 1 for value in series_values)
            )
        else:
            row["series_contrast_direction_consistent"] = None
        result.append(row)
    bh_adjust(result, ("distance_frame", "relative_orientation"))
    return result


def solve_linear(matrix: list[list[float]], vector: list[float]) -> list[float]:
    size = len(vector)
    augmented = [matrix[row][:] + [vector[row]] for row in range(size)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1e-12:
            raise ValueError("singular weighted regression")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                augmented[row][item] - factor * augmented[column][item]
                for item in range(size + 1)
            ]
    return [augmented[row][-1] for row in range(size)]


def weighted_fit(
    x_values: list[float],
    y_values: list[float],
    weights: list[float],
    period: float,
    harmonic: bool,
) -> tuple[list[float], float]:
    scale = max(max(abs(value) for value in x_values), 1.0)
    omega = 2.0 * math.pi / period
    rows: list[list[float]] = []
    for value in x_values:
        normalized = value / scale
        row = [1.0, normalized, normalized * normalized]
        if harmonic:
            row.extend((math.sin(omega * value), math.cos(omega * value)))
        rows.append(row)
    width = len(rows[0])
    normal = [[0.0] * width for _ in range(width)]
    right = [0.0] * width
    for row, outcome, weight in zip(rows, y_values, weights):
        for i in range(width):
            right[i] += weight * row[i] * outcome
            for j in range(width):
                normal[i][j] += weight * row[i] * row[j]
    coefficients = solve_linear(normal, right)
    residual = 0.0
    for row, outcome, weight in zip(rows, y_values, weights):
        predicted = sum(coefficient * value for coefficient, value in zip(coefficients, row))
        residual += weight * (outcome - predicted) ** 2
    return coefficients, residual


def periodicity_rows(
    responses: list[dict[str, Any]], period: float, minimum_offsets: int
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    groups: dict[tuple[str, ...], list[dict[str, Any]]] = {}
    for row in responses:
        if row["evaluation_status"] != "ok":
            continue
        key = (
            row["motif_id"],
            row["motif_name"],
            row["distance_frame"],
            row["relative_orientation"],
            row["isoform"],
        )
        groups.setdefault(key, []).append(row)
    summaries: list[dict[str, Any]] = []
    for key, group in groups.items():
        group.sort(key=lambda row: row["center_offset_bp"])
        if len(group) < minimum_offsets:
            continue
        x_values = [row["center_offset_bp"] for row in group]
        y_values = [row["adjusted_log_odds"] for row in group]
        weights = [1.0 / max(row["jackknife_se"] ** 2, 1e-6) for row in group]
        try:
            _, baseline_residual = weighted_fit(
                x_values, y_values, weights, period, harmonic=False
            )
            coefficients, harmonic_residual = weighted_fit(
                x_values, y_values, weights, period, harmonic=True
            )
        except ValueError:
            continue
        sine, cosine = coefficients[-2:]
        amplitude = math.hypot(sine, cosine)
        phase = (math.atan2(sine, cosine) / (2.0 * math.pi / period)) % period
        reduction = (
            max(0.0, baseline_residual - harmonic_residual) / baseline_residual
            if baseline_residual > 0
            else None
        )
        summaries.append(
            {
                "motif_id": key[0],
                "motif_name": key[1],
                "distance_frame": key[2],
                "relative_orientation": key[3],
                "isoform": key[4],
                "fixed_period_bp": period,
                "estimable_offsets": len(group),
                "minimum_offset_bp": min(x_values),
                "maximum_offset_bp": max(x_values),
                "sine_coefficient": sine,
                "cosine_coefficient": cosine,
                "harmonic_amplitude_log_odds": amplitude,
                "harmonic_phase_peak_bp_modulo_period": phase,
                "baseline_weighted_rss": baseline_residual,
                "harmonic_weighted_rss": harmonic_residual,
                "weighted_rss_fraction_reduced": reduction,
                "inference_status": "descriptive_fixed_period",
            }
        )
    by_key = {
        (
            row["motif_id"],
            row["distance_frame"],
            row["relative_orientation"],
            row["isoform"],
        ): row
        for row in summaries
    }
    contrasts: list[dict[str, Any]] = []
    bases = sorted({key[:3] for key in by_key})
    for motif_id, frame, orientation in bases:
        ta = by_key.get((motif_id, frame, orientation, "TA"))
        dn = by_key.get((motif_id, frame, orientation, "DN"))
        if ta is None or dn is None:
            continue
        raw_shift = (
            dn["harmonic_phase_peak_bp_modulo_period"]
            - ta["harmonic_phase_peak_bp_modulo_period"]
        )
        phase_shift = (raw_shift + period / 2.0) % period - period / 2.0
        contrasts.append(
            {
                "motif_id": motif_id,
                "motif_name": ta["motif_name"],
                "distance_frame": frame,
                "relative_orientation": orientation,
                "fixed_period_bp": period,
                "ta_amplitude_log_odds": ta["harmonic_amplitude_log_odds"],
                "dn_amplitude_log_odds": dn["harmonic_amplitude_log_odds"],
                "dn_minus_ta_amplitude_log_odds": (
                    dn["harmonic_amplitude_log_odds"]
                    - ta["harmonic_amplitude_log_odds"]
                ),
                "ta_phase_peak_bp_modulo_period": ta[
                    "harmonic_phase_peak_bp_modulo_period"
                ],
                "dn_phase_peak_bp_modulo_period": dn[
                    "harmonic_phase_peak_bp_modulo_period"
                ],
                "dn_minus_ta_phase_shift_bp_wrapped": phase_shift,
                "inference_status": "descriptive_fixed_period",
            }
        )
    return summaries, contrasts


def smooth_curve(
    group: list[dict[str, Any]], sigma: float
) -> list[tuple[float, float, dict[str, Any]]]:
    points = sorted(
        (
            row["center_offset_bp"],
            row["adjusted_log_odds"],
            row,
        )
        for row in group
        if row["evaluation_status"] == "ok"
    )
    radius = max(2.0, 3.0 * sigma)
    result: list[tuple[float, float, dict[str, Any]]] = []
    for center, _, source_row in points:
        numerator = 0.0
        denominator = 0.0
        for position, value, _ in points:
            delta = position - center
            if abs(delta) > radius:
                continue
            weight = math.exp(-0.5 * (delta / sigma) ** 2)
            numerator += weight * value
            denominator += weight
        if denominator > 0:
            result.append((center, numerator / denominator, source_row))
    return result


def detect_peaks(
    responses: list[dict[str, Any]],
    sigma: float,
    minimum_prominence: float,
    minimum_separation: float,
    maximum_peaks: int,
) -> list[dict[str, Any]]:
    groups: dict[tuple[str, ...], list[dict[str, Any]]] = {}
    for row in responses:
        key = (
            row["motif_id"],
            row["motif_name"],
            row["distance_frame"],
            row["relative_orientation"],
            row["isoform"],
        )
        groups.setdefault(key, []).append(row)
    result: list[dict[str, Any]] = []
    for key, group in groups.items():
        curve = smooth_curve(group, sigma)
        candidates: list[dict[str, Any]] = []
        for index in range(1, len(curve) - 1):
            center, value, source = curve[index]
            left_value = curve[index - 1][1]
            right_value = curve[index + 1][1]
            for direction, is_peak, prominence in (
                ("enriched", value > left_value and value >= right_value,
                 value - max(left_value, right_value)),
                ("depleted", value < left_value and value <= right_value,
                 min(left_value, right_value) - value),
            ):
                if not is_peak:
                    continue
                local = [
                    item[1]
                    for item in curve
                    if 2.0 <= abs(item[0] - center) <= max(6.0, minimum_separation)
                ]
                if local:
                    baseline = sum(local) / len(local)
                    prominence = (
                        value - baseline if direction == "enriched" else baseline - value
                    )
                if prominence < minimum_prominence:
                    continue
                candidates.append(
                    {
                        "motif_id": key[0],
                        "motif_name": key[1],
                        "distance_frame": key[2],
                        "relative_orientation": key[3],
                        "isoform": key[4],
                        "peak_direction": direction,
                        "peak_center_bp": center,
                        "smoothed_log_odds": value,
                        "prominence_log_odds": prominence,
                        "raw_adjusted_odds_ratio": source["adjusted_odds_ratio"],
                        "raw_q_value_bh_declared_offsets": source[
                            "q_value_bh_declared_offsets"
                        ],
                        "itemset_count": source["itemset_count"],
                        "smoothing_sigma_bp": sigma,
                        "peak_status": "descriptive_candidate",
                    }
                )
        for direction in ("enriched", "depleted"):
            selected: list[dict[str, Any]] = []
            directional = (
                candidate
                for candidate in candidates
                if candidate["peak_direction"] == direction
            )
            for candidate in sorted(
                directional,
                key=lambda item: item["prominence_log_odds"],
                reverse=True,
            ):
                if len(selected) >= maximum_peaks:
                    break
                if any(
                    abs(item["peak_center_bp"] - candidate["peak_center_bp"])
                    < minimum_separation
                    for item in selected
                ):
                    continue
                selected.append(candidate)
            selected.sort(key=lambda item: item["peak_center_bp"])
            for peak_index, row in enumerate(selected, start=1):
                row["peak_index"] = peak_index
                result.append(row)
    return result


def helical_distance(value: float, period: float) -> tuple[int, float]:
    turns = max(0, round(abs(value) / period))
    return turns, abs(abs(value) - turns * period)


def ordered_peak_alignment(
    ta_rows: list[dict[str, Any]],
    dn_rows: list[dict[str, Any]],
    shift_hypothesis: float,
    tolerance: float,
    maximum_shift: float,
) -> tuple[list[tuple[dict[str, Any], dict[str, Any]]], float, float]:
    """Return an order-preserving maximum-cardinality peak alignment."""

    ta_ordered = sorted(ta_rows, key=lambda row: row["peak_center_bp"])
    dn_ordered = sorted(dn_rows, key=lambda row: row["peak_center_bp"])
    # Each state is (matches, negative residual, prominence, pairs).  Peak order
    # is preserved so a shifted doublet cannot be cross-matched to itself.
    states: list[list[tuple[int, float, float, list[tuple[Any, Any]]]]] = [
        [(0, 0.0, 0.0, []) for _ in range(len(dn_ordered) + 1)]
        for _ in range(len(ta_ordered) + 1)
    ]

    def better(
        left: tuple[int, float, float, list[tuple[Any, Any]]],
        right: tuple[int, float, float, list[tuple[Any, Any]]],
    ) -> tuple[int, float, float, list[tuple[Any, Any]]]:
        return left if left[:3] >= right[:3] else right

    for ta_index in range(len(ta_ordered) + 1):
        for dn_index in range(len(dn_ordered) + 1):
            if ta_index:
                states[ta_index][dn_index] = better(
                    states[ta_index][dn_index], states[ta_index - 1][dn_index]
                )
            if dn_index:
                states[ta_index][dn_index] = better(
                    states[ta_index][dn_index], states[ta_index][dn_index - 1]
                )
            if not ta_index or not dn_index:
                continue
            ta = ta_ordered[ta_index - 1]
            dn = dn_ordered[dn_index - 1]
            observed_shift = dn["peak_center_bp"] - ta["peak_center_bp"]
            if abs(observed_shift) > maximum_shift:
                continue
            residual = abs(observed_shift - shift_hypothesis)
            if residual > tolerance:
                continue
            previous = states[ta_index - 1][dn_index - 1]
            candidate = (
                previous[0] + 1,
                previous[1] - residual,
                previous[2]
                + ta["prominence_log_odds"]
                + dn["prominence_log_odds"],
                previous[3] + [(ta, dn)],
            )
            states[ta_index][dn_index] = better(
                states[ta_index][dn_index], candidate
            )
    selected = states[-1][-1]
    return selected[3], -selected[1], selected[2]


def coherent_peak_alignment(
    ta_rows: list[dict[str, Any]],
    dn_rows: list[dict[str, Any]],
    period: float,
    maximum_shift: float,
) -> tuple[float, list[tuple[dict[str, Any], dict[str, Any]]]]:
    """Choose one shift hypothesis before pairing peaks within an isoform group."""

    if not ta_rows or not dn_rows:
        return 0.0, []
    hypotheses = {0.0}
    turns = int(maximum_shift // period)
    for turn in range(1, turns + 1):
        hypotheses.update((turn * period, -turn * period))
    for ta in ta_rows:
        for dn in dn_rows:
            shift = dn["peak_center_bp"] - ta["peak_center_bp"]
            if abs(shift) <= maximum_shift:
                hypotheses.add(float(shift))
    tolerance = min(3.0, max(1.0, period / 4.0))
    ranked: list[
        tuple[int, float, float, float, float, list[tuple[dict[str, Any], dict[str, Any]]]]
    ] = []
    for hypothesis in hypotheses:
        pairs, residual, prominence = ordered_peak_alignment(
            ta_rows, dn_rows, hypothesis, tolerance, maximum_shift
        )
        ranked.append(
            (
                len(pairs),
                -residual,
                prominence,
                -abs(hypothesis),
                hypothesis,
                pairs,
            )
        )
    best = max(ranked, key=lambda item: item[:4])
    return best[4], best[5]


def match_peaks(
    peaks: list[dict[str, Any]], period: float, maximum_shift: float
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    groups: dict[tuple[str, ...], dict[str, list[dict[str, Any]]]] = {}
    for row in peaks:
        key = (
            row["motif_id"],
            row["motif_name"],
            row["distance_frame"],
            row["relative_orientation"],
            row["peak_direction"],
        )
        groups.setdefault(key, {}).setdefault(row["isoform"], []).append(row)
    matches: list[dict[str, Any]] = []
    for key, isoforms in groups.items():
        ta_rows = isoforms.get("TA", [])
        dn_rows = isoforms.get("DN", [])
        group_shift, aligned = coherent_peak_alignment(
            ta_rows, dn_rows, period, maximum_shift
        )
        group_turns, group_residual = helical_distance(group_shift, period)
        for ta, dn in aligned:
            shift = dn["peak_center_bp"] - ta["peak_center_bp"]
            turns, residual = helical_distance(shift, period)
            matches.append(
                {
                    "motif_id": key[0],
                    "motif_name": key[1],
                    "distance_frame": key[2],
                    "relative_orientation": key[3],
                    "peak_direction": key[4],
                    "ta_peak_index": ta["peak_index"],
                    "dn_peak_index": dn["peak_index"],
                    "ta_peak_center_bp": ta["peak_center_bp"],
                    "dn_peak_center_bp": dn["peak_center_bp"],
                    "dn_minus_ta_peak_shift_bp": shift,
                    "group_shift_hypothesis_bp": group_shift,
                    "match_residual_to_group_shift_bp": abs(shift - group_shift),
                    "matched_peaks_in_group": len(aligned),
                    "group_shift_nearest_helical_turns": group_turns,
                    "group_shift_helical_residual_bp": group_residual,
                    "same_side_of_anchor": (
                        ta["peak_center_bp"] * dn["peak_center_bp"] > 0
                    ),
                    "dn_peak_is_closer_to_anchor": (
                        abs(dn["peak_center_bp"]) < abs(ta["peak_center_bp"])
                    ),
                    "nearest_helical_turns": turns,
                    "helical_shift_residual_bp": residual,
                    "fixed_period_bp": period,
                    "match_status": "descriptive_candidate",
                }
            )
    pair_matches: list[dict[str, Any]] = []
    match_groups: dict[tuple[str, ...], list[dict[str, Any]]] = {}
    for row in matches:
        key = (
            row["motif_id"],
            row["motif_name"],
            row["distance_frame"],
            row["relative_orientation"],
            row["peak_direction"],
        )
        match_groups.setdefault(key, []).append(row)
    for key, group in match_groups.items():
        group.sort(key=lambda row: row["ta_peak_center_bp"])
        for left, right in zip(group, group[1:]):
            ta_spacing = right["ta_peak_center_bp"] - left["ta_peak_center_bp"]
            dn_spacing = right["dn_peak_center_bp"] - left["dn_peak_center_bp"]
            if ta_spacing > 3.0 * period or dn_spacing > 3.0 * period:
                continue
            ta_turns, ta_residual = helical_distance(ta_spacing, period)
            dn_turns, dn_residual = helical_distance(dn_spacing, period)
            ta_midpoint = (
                left["ta_peak_center_bp"] + right["ta_peak_center_bp"]
            ) / 2.0
            dn_midpoint = (
                left["dn_peak_center_bp"] + right["dn_peak_center_bp"]
            ) / 2.0
            midpoint_shift = dn_midpoint - ta_midpoint
            midpoint_turns, midpoint_residual = helical_distance(midpoint_shift, period)
            pair_matches.append(
                {
                    "motif_id": key[0],
                    "motif_name": key[1],
                    "distance_frame": key[2],
                    "relative_orientation": key[3],
                    "peak_direction": key[4],
                    "ta_left_peak_bp": left["ta_peak_center_bp"],
                    "ta_right_peak_bp": right["ta_peak_center_bp"],
                    "dn_left_peak_bp": left["dn_peak_center_bp"],
                    "dn_right_peak_bp": right["dn_peak_center_bp"],
                    "ta_peak_spacing_bp": ta_spacing,
                    "dn_peak_spacing_bp": dn_spacing,
                    "spacing_difference_dn_minus_ta_bp": dn_spacing - ta_spacing,
                    "ta_nearest_helical_turns": ta_turns,
                    "ta_helical_spacing_residual_bp": ta_residual,
                    "dn_nearest_helical_turns": dn_turns,
                    "dn_helical_spacing_residual_bp": dn_residual,
                    "ta_peak_midpoint_bp": ta_midpoint,
                    "dn_peak_midpoint_bp": dn_midpoint,
                    "dn_minus_ta_midpoint_shift_bp": midpoint_shift,
                    "group_shift_hypothesis_bp": left[
                        "group_shift_hypothesis_bp"
                    ],
                    "dn_midpoint_is_closer_to_anchor": abs(dn_midpoint) < abs(ta_midpoint),
                    "midpoint_shift_nearest_helical_turns": midpoint_turns,
                    "midpoint_shift_helical_residual_bp": midpoint_residual,
                    "fixed_period_bp": period,
                    "pair_status": "descriptive_candidate",
                }
            )
    return matches, pair_matches


def write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("x", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, allow_nan=False) + "\n")


def parquet_from_rows(
    duckdb: str,
    rows: list[dict[str, Any]],
    path: Path,
    schema: tuple[tuple[str, str], ...],
) -> None:
    projection = ", ".join(
        f'CAST("{name}" AS {data_type}) AS "{name}"'
        for name, data_type in schema
    )
    if rows:
        jsonl = path.with_suffix(".jsonl")
        write_jsonl(jsonl, rows)
        run_duckdb(
            duckdb,
            ":memory:",
            f"COPY (SELECT {projection} FROM "
            f"read_json_auto({sql_string(jsonl)})) "
            f"TO {sql_string(path)} (FORMAT PARQUET, COMPRESSION ZSTD);",
        )
        jsonl.unlink()
    else:
        empty_projection = ", ".join(
            f'NULL::{data_type} AS "{name}"' for name, data_type in schema
        )
        run_duckdb(
            duckdb,
            ":memory:",
            f"COPY (SELECT {empty_projection} WHERE false) "
            f"TO {sql_string(path)} (FORMAT PARQUET, COMPRESSION ZSTD);",
        )


def write_report(
    path: Path,
    occurrence_rows: list[dict[str, str]],
    periodicity: list[dict[str, Any]],
    peak_matches: list[dict[str, Any]],
    pair_matches: list[dict[str, Any]],
) -> None:
    lines = [
        "# TP73 exact-distance isoform report",
        "",
        "This candidate-level report uses TP73-oriented motif-center distances.",
        "Support/confidence follow TF-COMB transaction terminology, but the",
        "transactions are explicitly conditional on TP73 anchors; lift is not",
        "reported because it is degenerate in that conditioned basket.",
        "",
        "## TP73-conditioned occurrence",
        "",
        "Forward confidence equals support because every transaction is one",
        "TP73 anchor. These are occurrence summaries, not CUT&RUN effects.",
        "",
        "| Motif | Frame | Orientation | Offset bp | Positive / total | Support | Source support |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in occurrence_rows:
        lines.append(
            f"| {row['motif_name']} ({row['motif_id']}) | "
            f"{row['distance_frame']} | {row['relative_orientation']} | "
            f"{float(row['center_offset_bp']):.3g} | "
            f"{row['itemset_count']} / {row['transaction_count']} | "
            f"{float(row['support']):.4g} | "
            f"{float(row['source_support']):.4g} |"
        )
    lines.extend(
        [
            "",
            "## Fixed-period summaries",
            "",
            "| Motif | Frame | Orientation | Isoform | Amplitude | Phase (bp) | RSS reduction |",
            "|---|---|---|---|---:|---:|---:|",
        ]
    )
    for row in sorted(
        periodicity,
        key=lambda item: (
            item["motif_id"], item["distance_frame"],
            item["relative_orientation"], item["isoform"],
        ),
    ):
        reduction = row["weighted_rss_fraction_reduced"]
        reduction_text = "NA" if reduction is None else f"{reduction:.3g}"
        lines.append(
            f"| {row['motif_name']} ({row['motif_id']}) | "
            f"{row['distance_frame']} | {row['relative_orientation']} | "
            f"{row['isoform']} | {row['harmonic_amplitude_log_odds']:.4g} | "
            f"{row['harmonic_phase_peak_bp_modulo_period']:.3g} | "
            f"{reduction_text} |"
        )
    lines.extend(
        [
            "",
            "## Matched single peaks",
            "",
            "| Motif | Direction | TA bp | DN bp | Shift bp | DN closer | Helical residual bp |",
            "|---|---|---:|---:|---:|---|---:|",
        ]
    )
    for row in sorted(
        peak_matches,
        key=lambda item: (item["motif_id"], item["peak_direction"], item["ta_peak_center_bp"]),
    ):
        lines.append(
            f"| {row['motif_name']} ({row['motif_id']}) | {row['peak_direction']} | "
            f"{row['ta_peak_center_bp']:.3g} | {row['dn_peak_center_bp']:.3g} | "
            f"{row['dn_minus_ta_peak_shift_bp']:.3g} | "
            f"{str(row['dn_peak_is_closer_to_anchor']).lower()} | "
            f"{row['helical_shift_residual_bp']:.3g} |"
        )
    lines.extend(
        [
            "",
            "## Matched peak doublets",
            "",
            "| Motif | Direction | TA spacing | DN spacing | Midpoint shift | DN closer |",
            "|---|---|---:|---:|---:|---|",
        ]
    )
    for row in sorted(
        pair_matches,
        key=lambda item: (item["motif_id"], item["peak_direction"], item["ta_left_peak_bp"]),
    ):
        lines.append(
            f"| {row['motif_name']} ({row['motif_id']}) | {row['peak_direction']} | "
            f"{row['ta_peak_spacing_bp']:.3g} | {row['dn_peak_spacing_bp']:.3g} | "
            f"{row['dn_minus_ta_midpoint_shift_bp']:.3g} | "
            f"{str(row['dn_midpoint_is_closer_to_anchor']).lower()} |"
        )
    lines.extend(
        [
            "",
            "Peak and periodicity rows are descriptive candidates. Block-jackknife",
            "inference applies to the underlying one-base odds-ratio curves; peak",
            "location uncertainty and a block-shift spatial null remain validation",
            "requirements before biological claims.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_raw_tables(
    arguments: argparse.Namespace,
    database: Path,
    response_tsv: Path,
    contrast_tsv: Path,
) -> None:
    inventory_paths = [path.expanduser().resolve() for path in arguments.inventory]
    block_paths = [path.expanduser().resolve() for path in arguments.block_components]
    for path in inventory_paths + block_paths:
        if not path.is_file():
            raise ExactDistanceSummaryError(f"input file is missing: {path}")
    run_duckdb(arguments.duckdb, database, f"""
SET threads={arguments.threads};
SET memory_limit={sql_string(arguments.memory_limit)};
SET preserve_insertion_order=false;

CREATE TABLE inventory_chrom AS
SELECT * FROM read_parquet({sql_path_list(inventory_paths)}, hive_partitioning=false);
CREATE TABLE block_component AS
SELECT * FROM read_parquet({sql_path_list(block_paths)}, hive_partitioning=false);

SELECT CASE WHEN EXISTS (
  SELECT motif_id FROM inventory_chrom
  GROUP BY motif_id
  HAVING min(source_score_floor) <> max(source_score_floor)
      OR min(positive_threshold) <> max(positive_threshold)
      OR min(context_flank_bp) <> max(context_flank_bp)
      OR min(block_size_bp) <> max(block_size_bp)
) THEN error('a motif has inconsistent score policies') END;
SELECT CASE WHEN EXISTS (
  SELECT motif_id FROM block_component
  GROUP BY motif_id
  HAVING min(source_score_floor) <> max(source_score_floor)
      OR min(positive_threshold) <> max(positive_threshold)
      OR min(context_flank_bp) <> max(context_flank_bp)
      OR min(block_size_bp) <> max(block_size_bp)
) THEN error('a motif has inconsistent block score policies') END;
SELECT CASE WHEN EXISTS (
  SELECT motif_id, chrom, distance_frame, relative_orientation,
         center_offset_twice
  FROM inventory_chrom GROUP BY ALL HAVING count(*) <> 1
) THEN error('distance inventory contains duplicate chromosome cells') END;
SELECT CASE WHEN EXISTS (
  SELECT motif_id, chrom, component_type, distance_frame,
         relative_orientation, center_offset_twice, block_index,
         tp73_score_stratum, series_id, isoform
  FROM block_component GROUP BY ALL HAVING count(*) <> 1
) THEN error('block components contain duplicate sufficient-statistic cells') END;
SELECT CASE WHEN EXISTS (
  (SELECT DISTINCT motif_id, chrom FROM inventory_chrom
   EXCEPT SELECT DISTINCT motif_id, chrom FROM block_component)
  UNION ALL
  (SELECT DISTINCT motif_id, chrom FROM block_component
   EXCEPT SELECT DISTINCT motif_id, chrom FROM inventory_chrom)
) THEN error('inventory and block motif/chromosome keys differ') END;
SELECT CASE WHEN EXISTS (
  SELECT 1
  FROM (
    SELECT motif_id, chrom, min(source_score_floor) source_score_floor,
           min(positive_threshold) positive_threshold,
           min(context_flank_bp) context_flank_bp,
           min(block_size_bp) block_size_bp
    FROM block_component GROUP BY motif_id, chrom
  ) b
  JOIN (
    SELECT motif_id, chrom, min(source_score_floor) source_score_floor,
           min(positive_threshold) positive_threshold,
           min(context_flank_bp) context_flank_bp,
           min(block_size_bp) block_size_bp
    FROM inventory_chrom GROUP BY motif_id, chrom
  ) i USING (motif_id, chrom)
  WHERE b.source_score_floor <> i.source_score_floor
     OR b.positive_threshold <> i.positive_threshold
     OR b.context_flank_bp <> i.context_flank_bp
     OR b.block_size_bp <> i.block_size_bp
) THEN error('inventory and block score policies differ') END;
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM inventory_chrom
  WHERE transaction_count <> itemset_count + anchors_intermediate + anchors_negative
     OR itemset_count > anchors_source_present
     OR anchors_source_present <> itemset_count + anchors_intermediate
     OR abs(center_offset_bp * 2.0 - center_offset_twice) > 1e-9
) THEN error('distance inventory class counts are inconsistent') END;
SELECT CASE WHEN EXISTS (
  SELECT 1 FROM block_component
  WHERE component_type NOT IN ('baseline', 'observed')
     OR (component_type = 'baseline'
         AND (relative_orientation IS NOT NULL
              OR center_offset_twice IS NOT NULL))
     OR (component_type = 'observed'
         AND (relative_orientation IS NULL
              OR center_offset_twice IS NULL))
) THEN error('block component type/geometry contract is inconsistent') END;

CREATE TABLE distance_inventory AS
SELECT motif_id, max(motif_name) AS motif_name, max(anchor_motif_id) AS anchor_motif_id,
       max(transaction_scope) AS transaction_scope,
       max(conditioning_scope) AS conditioning_scope,
       max(locus_identity_rule) AS locus_identity_rule,
       max(distance_metric) AS distance_metric,
       max(overlap_policy) AS overlap_policy,
       max(context_flank_bp)::BIGINT AS context_flank_bp,
       max(block_size_bp)::BIGINT AS block_size_bp,
       max(maximum_anchor_span_bp)::BIGINT AS maximum_anchor_span_bp,
       max(maximum_neighbor_span_bp)::BIGINT AS maximum_neighbor_span_bp,
       distance_frame, relative_orientation, center_offset_twice,
       max(center_offset_bp)::DOUBLE AS center_offset_bp,
       max(source_score_floor)::DOUBLE AS source_score_floor,
       max(positive_threshold)::DOUBLE AS positive_threshold,
       sum(transaction_count)::BIGINT AS transaction_count,
       sum(anchors_source_present)::BIGINT AS anchors_source_present,
       sum(itemset_count)::BIGINT AS itemset_count,
       sum(anchors_intermediate)::BIGINT AS anchors_intermediate,
       sum(anchors_negative)::BIGINT AS anchors_negative,
       sum(physical_locus_pair_count)::BIGINT AS physical_locus_pair_count,
       min(minimum_interval_distance_bp)::BIGINT AS minimum_interval_distance_bp,
       max(maximum_interval_distance_bp)::BIGINT AS maximum_interval_distance_bp,
       sum(anchors_source_present) / nullif(sum(transaction_count), 0)::DOUBLE
         AS source_support,
       sum(itemset_count) / nullif(sum(transaction_count), 0)::DOUBLE AS support,
       sum(itemset_count) / nullif(sum(transaction_count), 0)::DOUBLE
         AS confidence_anchor_to_neighbor
FROM inventory_chrom
GROUP BY motif_id, distance_frame, relative_orientation, center_offset_twice;

CREATE TABLE run_config AS
SELECT motif_id, max(motif_name) AS motif_name,
       max(anchor_motif_id) AS anchor_motif_id,
       max(conditioning_scope) AS conditioning_scope,
       max(distance_metric) AS distance_metric,
       max(overlap_policy) AS overlap_policy,
       max(source_score_floor)::DOUBLE AS source_score_floor,
       max(positive_threshold)::DOUBLE AS positive_threshold,
       max(context_flank_bp)::BIGINT AS context_flank_bp,
       max(block_size_bp)::BIGINT AS block_size_bp,
       max(maximum_anchor_span_bp)::BIGINT AS maximum_anchor_span_bp,
       max(maximum_neighbor_span_bp)::BIGINT AS maximum_neighbor_span_bp,
       count(DISTINCT chrom)::BIGINT AS chromosome_count,
       list_sort(list(DISTINCT chrom)) AS chromosomes
FROM inventory_chrom GROUP BY motif_id;

CREATE TEMP TABLE dimension AS
SELECT motif_id, motif_name, anchor_motif_id, distance_frame,
       relative_orientation, center_offset_twice, center_offset_bp,
       source_score_floor, positive_threshold
FROM distance_inventory
WHERE itemset_count >= {arguments.minimum_itemset_count};

CREATE TEMP TABLE baseline AS
SELECT motif_id, max(motif_name) AS motif_name, chrom, block_index,
       tp73_score_stratum, series_id, isoform, distance_frame,
       sum(total_anti_discordant)::BIGINT AS total_anti,
       sum(total_control_discordant)::BIGINT AS total_control
FROM block_component
WHERE component_type='baseline'
GROUP BY motif_id, chrom, block_index, tp73_score_stratum,
         series_id, isoform, distance_frame;

CREATE TEMP TABLE observed AS
SELECT motif_id, chrom, block_index, tp73_score_stratum, series_id, isoform,
       distance_frame, relative_orientation, center_offset_twice,
       sum(source_anti_discordant)::BIGINT AS source_anti,
       sum(source_control_discordant)::BIGINT AS source_control,
       sum(positive_anti_discordant)::BIGINT AS positive_anti,
       sum(positive_control_discordant)::BIGINT AS positive_control
FROM block_component
WHERE component_type='observed'
GROUP BY motif_id, chrom, block_index, tp73_score_stratum, series_id, isoform,
         distance_frame, relative_orientation, center_offset_twice;

CREATE TEMP TABLE mh_component AS
WITH cell AS (
  SELECT d.*, b.chrom, b.block_index, b.tp73_score_stratum,
         b.series_id, b.isoform, b.total_anti, b.total_control,
         coalesce(o.source_anti, 0)::BIGINT AS source_anti,
         coalesce(o.source_control, 0)::BIGINT AS source_control,
         coalesce(o.positive_anti, 0)::BIGINT AS positive_anti,
         coalesce(o.positive_control, 0)::BIGINT AS positive_control,
         (b.total_anti - coalesce(o.source_anti, 0))::BIGINT AS negative_anti,
         (b.total_control - coalesce(o.source_control, 0))::BIGINT
           AS negative_control
  FROM dimension d
  JOIN baseline b USING (motif_id, distance_frame)
  LEFT JOIN observed o
    ON o.motif_id=d.motif_id AND o.chrom=b.chrom AND o.block_index=b.block_index
   AND o.tp73_score_stratum=b.tp73_score_stratum
   AND o.series_id=b.series_id AND o.isoform=b.isoform
   AND o.distance_frame=d.distance_frame
   AND o.relative_orientation=d.relative_orientation
   AND o.center_offset_twice=d.center_offset_twice
), component AS (
  SELECT *, (positive_anti + positive_control + negative_anti + negative_control)
              ::DOUBLE AS stratum_total
  FROM cell
)
SELECT *,
       CASE WHEN stratum_total > 0
         THEN positive_anti * negative_control / stratum_total ELSE 0 END::DOUBLE
         AS mh_numerator,
       CASE WHEN stratum_total > 0
         THEN positive_control * negative_anti / stratum_total ELSE 0 END::DOUBLE
         AS mh_denominator
FROM component;

COPY (
WITH total AS (
  SELECT motif_id, max(motif_name) AS motif_name, max(anchor_motif_id) anchor_motif_id,
         distance_frame, relative_orientation, center_offset_twice,
         max(center_offset_bp)::DOUBLE AS center_offset_bp, isoform,
         max(source_score_floor)::DOUBLE AS source_score_floor,
         max(positive_threshold)::DOUBLE AS positive_threshold,
         sum(mh_numerator)::DOUBLE AS numerator,
         sum(mh_denominator)::DOUBLE AS denominator,
         sum(positive_anti)::BIGINT AS positive_anti_discordant,
         sum(positive_control)::BIGINT AS positive_control_discordant,
         sum(negative_anti)::BIGINT AS negative_anti_discordant,
         sum(negative_control)::BIGINT AS negative_control_discordant,
         count(DISTINCT (chrom, block_index))::BIGINT AS genomic_blocks
  FROM mh_component
  GROUP BY motif_id, distance_frame, relative_orientation,
           center_offset_twice, isoform
), block AS (
  SELECT motif_id, distance_frame, relative_orientation, center_offset_twice,
         isoform, chrom, block_index,
         sum(mh_numerator)::DOUBLE AS block_numerator,
         sum(mh_denominator)::DOUBLE AS block_denominator
  FROM mh_component
  GROUP BY motif_id, distance_frame, relative_orientation,
           center_offset_twice, isoform, chrom, block_index
), leave_one_out AS (
  SELECT t.motif_id, t.distance_frame, t.relative_orientation,
         t.center_offset_twice, t.isoform,
         CASE WHEN t.numerator-b.block_numerator > 0
                AND t.denominator-b.block_denominator > 0
              THEN ln((t.numerator-b.block_numerator)
                      /(t.denominator-b.block_denominator)) ELSE NULL END
           AS leave_one_out_log_or
  FROM total t JOIN block b USING (
    motif_id, distance_frame, relative_orientation, center_offset_twice, isoform
  )
), jackknife AS (
  SELECT motif_id, distance_frame, relative_orientation, center_offset_twice,
         isoform, count(leave_one_out_log_or)::BIGINT AS jackknife_blocks,
         sqrt((count(leave_one_out_log_or)-1)
              * var_pop(leave_one_out_log_or))::DOUBLE AS jackknife_se
  FROM leave_one_out GROUP BY ALL
), series AS (
  SELECT motif_id, distance_frame, relative_orientation, center_offset_twice,
         isoform, series_id,
         CASE WHEN sum(mh_numerator)>0 AND sum(mh_denominator)>0
              THEN sum(mh_numerator)/sum(mh_denominator) ELSE NULL END
           AS series_or
  FROM mh_component GROUP BY ALL
), series_wide AS (
  SELECT motif_id, distance_frame, relative_orientation, center_offset_twice,
         isoform,
         max(series_or) FILTER (WHERE series_id='saos2') AS saos2_or,
         max(series_or) FILTER (WHERE series_id='skmel29_2') AS skmel29_2_or
  FROM series GROUP BY ALL
)
SELECT t.*, i.transaction_scope, i.conditioning_scope, i.distance_metric,
       i.overlap_policy, i.transaction_count, i.anchors_source_present,
       i.context_flank_bp, i.block_size_bp, i.maximum_anchor_span_bp,
       i.maximum_neighbor_span_bp,
       i.itemset_count, i.anchors_intermediate, i.anchors_negative,
       i.physical_locus_pair_count, i.source_support, i.support,
       CASE WHEN t.numerator>0 AND t.denominator>0
            THEN t.numerator/t.denominator ELSE NULL END AS adjusted_odds_ratio,
       j.jackknife_blocks, j.jackknife_se, s.saos2_or, s.skmel29_2_or
FROM total t
JOIN distance_inventory i USING (
  motif_id, distance_frame, relative_orientation, center_offset_twice
)
JOIN jackknife j USING (
  motif_id, distance_frame, relative_orientation, center_offset_twice, isoform
)
JOIN series_wide s USING (
  motif_id, distance_frame, relative_orientation, center_offset_twice, isoform
)
ORDER BY motif_id, distance_frame, relative_orientation, isoform,
         center_offset_twice
) TO {sql_string(response_tsv)}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');

COPY (
WITH total AS (
  SELECT motif_id, max(motif_name) AS motif_name, distance_frame,
         relative_orientation, center_offset_twice,
         max(center_offset_bp)::DOUBLE AS center_offset_bp, isoform,
         sum(mh_numerator)::DOUBLE AS numerator,
         sum(mh_denominator)::DOUBLE AS denominator
  FROM mh_component GROUP BY ALL
), paired AS (
  SELECT motif_id, max(motif_name) motif_name, distance_frame,
         relative_orientation, center_offset_twice, max(center_offset_bp) center_offset_bp,
         max(numerator) FILTER (WHERE isoform='TA') AS ta_numerator,
         max(denominator) FILTER (WHERE isoform='TA') AS ta_denominator,
         max(numerator) FILTER (WHERE isoform='DN') AS dn_numerator,
         max(denominator) FILTER (WHERE isoform='DN') AS dn_denominator
  FROM total GROUP BY motif_id, distance_frame, relative_orientation,
                      center_offset_twice HAVING count(DISTINCT isoform)=2
), block AS (
  SELECT motif_id, distance_frame, relative_orientation, center_offset_twice,
         chrom, block_index,
         coalesce(sum(mh_numerator) FILTER (WHERE isoform='TA'),0)::DOUBLE ta_bn,
         coalesce(sum(mh_denominator) FILTER (WHERE isoform='TA'),0)::DOUBLE ta_bd,
         coalesce(sum(mh_numerator) FILTER (WHERE isoform='DN'),0)::DOUBLE dn_bn,
         coalesce(sum(mh_denominator) FILTER (WHERE isoform='DN'),0)::DOUBLE dn_bd
  FROM mh_component GROUP BY motif_id, distance_frame, relative_orientation,
                              center_offset_twice, chrom, block_index
), leave_one_out AS (
  SELECT p.motif_id,p.distance_frame,p.relative_orientation,p.center_offset_twice,
         CASE WHEN p.ta_numerator-b.ta_bn>0 AND p.ta_denominator-b.ta_bd>0
                AND p.dn_numerator-b.dn_bn>0 AND p.dn_denominator-b.dn_bd>0
              THEN ln((p.ta_numerator-b.ta_bn)/(p.ta_denominator-b.ta_bd))
                   -ln((p.dn_numerator-b.dn_bn)/(p.dn_denominator-b.dn_bd))
              ELSE NULL END AS loo
  FROM paired p JOIN block b USING (
    motif_id,distance_frame,relative_orientation,center_offset_twice
  )
), jackknife AS (
  SELECT motif_id,distance_frame,relative_orientation,center_offset_twice,
         count(*)::BIGINT AS genomic_blocks,
         count(loo)::BIGINT AS jackknife_blocks,
         sqrt((count(loo)-1)*var_pop(loo))::DOUBLE AS jackknife_se
  FROM leave_one_out GROUP BY ALL
), series AS (
  SELECT motif_id,distance_frame,relative_orientation,center_offset_twice,
         series_id,isoform,
         CASE WHEN sum(mh_numerator)>0 AND sum(mh_denominator)>0
              THEN sum(mh_numerator)/sum(mh_denominator) ELSE NULL END series_or
  FROM mh_component GROUP BY ALL
), series_contrast AS (
  SELECT motif_id,distance_frame,relative_orientation,center_offset_twice,
         max(series_or) FILTER (WHERE series_id='saos2' AND isoform='TA')
          /nullif(max(series_or) FILTER (WHERE series_id='saos2' AND isoform='DN'),0)
            AS saos2_ta_vs_dn_odds_ratio_ratio,
         max(series_or) FILTER (WHERE series_id='skmel29_2' AND isoform='TA')
          /nullif(max(series_or) FILTER (WHERE series_id='skmel29_2' AND isoform='DN'),0)
            AS skmel29_2_ta_vs_dn_odds_ratio_ratio
  FROM series GROUP BY ALL
)
SELECT p.motif_id,p.motif_name,p.distance_frame,p.relative_orientation,
       p.center_offset_twice,p.center_offset_bp,
       CASE WHEN p.ta_numerator>0 AND p.ta_denominator>0
            THEN p.ta_numerator/p.ta_denominator ELSE NULL END
         AS ta_adjusted_odds_ratio,
       CASE WHEN p.dn_numerator>0 AND p.dn_denominator>0
            THEN p.dn_numerator/p.dn_denominator ELSE NULL END
         AS dn_adjusted_odds_ratio,
       CASE WHEN p.ta_numerator>0 AND p.ta_denominator>0
              AND p.dn_numerator>0 AND p.dn_denominator>0
            THEN ln(p.ta_numerator/p.ta_denominator)
                 -ln(p.dn_numerator/p.dn_denominator) ELSE NULL END
         AS ta_vs_dn_log_odds_difference,
       CASE WHEN p.ta_numerator>0 AND p.ta_denominator>0
              AND p.dn_numerator>0 AND p.dn_denominator>0
            THEN (p.ta_numerator/p.ta_denominator)
                 /(p.dn_numerator/p.dn_denominator) ELSE NULL END
         AS ta_vs_dn_odds_ratio_ratio,
       j.genomic_blocks,j.jackknife_blocks,j.jackknife_se,
       s.saos2_ta_vs_dn_odds_ratio_ratio,
       s.skmel29_2_ta_vs_dn_odds_ratio_ratio
FROM paired p
JOIN jackknife j USING (
  motif_id,distance_frame,relative_orientation,center_offset_twice
)
JOIN series_contrast s USING (
  motif_id,distance_frame,relative_orientation,center_offset_twice
)
ORDER BY motif_id,distance_frame,relative_orientation,center_offset_twice
) TO {sql_string(contrast_tsv)}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
""")


def finalize(arguments: argparse.Namespace) -> None:
    duckdb = shutil.which(arguments.duckdb)
    if duckdb is None:
        raise ExactDistanceSummaryError(
            f"DuckDB executable not found: {arguments.duckdb}"
        )
    output = arguments.output_dir.expanduser().resolve()
    if output.exists():
        raise ExactDistanceSummaryError(f"refusing to replace output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}-", dir=output.parent))
    try:
        work_database = staging / "work.duckdb"
        response_tsv = staging / "response.raw.tsv"
        contrast_tsv = staging / "contrast.raw.tsv"
        build_raw_tables(arguments, work_database, response_tsv, contrast_tsv)
        inventory_parquet = staging / "distance_inventory.parquet"
        run_config_parquet = staging / "run_config.parquet"
        occurrence_tsv = staging / "occurrence.report.tsv"
        run_duckdb(
            duckdb,
            work_database,
            f"COPY distance_inventory TO {sql_string(inventory_parquet)} "
            "(FORMAT PARQUET, COMPRESSION ZSTD);",
        )
        run_duckdb(
            duckdb,
            work_database,
            f"COPY run_config TO {sql_string(run_config_parquet)} "
            "(FORMAT PARQUET, COMPRESSION ZSTD);",
        )
        run_duckdb(
            duckdb,
            work_database,
            f"""
COPY (
  SELECT motif_id, motif_name, distance_frame, relative_orientation,
         center_offset_bp, itemset_count, transaction_count, support,
         source_support
  FROM distance_inventory
  WHERE itemset_count > 0
  QUALIFY row_number() OVER (
    PARTITION BY motif_id, distance_frame, relative_orientation
    ORDER BY support DESC, itemset_count DESC, abs(center_offset_bp),
             center_offset_bp
  ) <= 5
  ORDER BY motif_id, distance_frame, relative_orientation, support DESC,
           itemset_count DESC, abs(center_offset_bp), center_offset_bp
) TO {sql_string(occurrence_tsv)}
  (FORMAT CSV, DELIMITER '\t', HEADER, NULL 'NA');
""",
        )
        responses = parse_response_rows(response_tsv, arguments.minimum_blocks)
        contrasts = parse_contrast_rows(contrast_tsv, arguments.minimum_blocks)
        periodicity, periodicity_contrast = periodicity_rows(
            responses, arguments.period, arguments.minimum_periodicity_offsets
        )
        peaks = detect_peaks(
            responses,
            arguments.smoothing_sigma,
            arguments.minimum_peak_prominence,
            arguments.minimum_peak_separation,
            arguments.maximum_peaks,
        )
        peak_matches, peak_pairs = match_peaks(
            peaks, arguments.period, arguments.maximum_peak_match_shift
        )

        outputs: dict[str, list[dict[str, Any]]] = {
            "distance_response.parquet": responses,
            "isoform_contrast.parquet": contrasts,
            "periodicity.parquet": periodicity,
            "periodicity_isoform_contrast.parquet": periodicity_contrast,
            "peak.parquet": peaks,
            "peak_match.parquet": peak_matches,
            "peak_pair_match.parquet": peak_pairs,
        }
        for name, rows in outputs.items():
            parquet_from_rows(duckdb, rows, staging / name, TABLE_SCHEMAS[name])

        report = staging / "report.md"
        write_report(
            report,
            read_tsv(occurrence_tsv),
            periodicity,
            peak_matches,
            peak_pairs,
        )
        final_database = staging / "tp73_exact_distance.duckdb"
        table_names = {
            "distance_inventory.parquet": "tp73_exact_distance_inventory",
            "run_config.parquet": "tp73_exact_distance_run_config",
            "distance_response.parquet": "tp73_exact_distance_response",
            "isoform_contrast.parquet": "tp73_exact_distance_isoform_contrast",
            "periodicity.parquet": "tp73_exact_distance_periodicity",
            "periodicity_isoform_contrast.parquet": (
                "tp73_exact_distance_periodicity_isoform_contrast"
            ),
            "peak.parquet": "tp73_exact_distance_peak",
            "peak_match.parquet": "tp73_exact_distance_peak_match",
            "peak_pair_match.parquet": "tp73_exact_distance_peak_pair_match",
        }
        statements = []
        for name, table in table_names.items():
            statements.append(
                f"CREATE TABLE {table} AS SELECT * FROM read_parquet("
                f"{sql_string(staging / name)});"
            )
        run_duckdb(duckdb, final_database, "\n".join(statements))

        response_tsv.unlink()
        contrast_tsv.unlink()
        occurrence_tsv.unlink()
        work_database.unlink()
        wal = work_database.with_suffix(work_database.suffix + ".wal")
        if wal.exists():
            wal.unlink()

        files = {}
        for path in sorted(staging.iterdir()):
            if path.name == "manifest.json":
                continue
            files[path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
        manifest = {
            "schema_version": 1,
            "analysis": "tp73_exact_motif_distance_isoform_response",
            "conditioning_scope": "conditional_on_tp73_anchor",
            "transaction_semantics": "one_physical_tp73_anchor",
            "tf_comb_adaptations": [
                "support",
                "confidence_anchor_to_neighbor",
                "exact_distance",
                "orientation",
                "overlaps_retained",
            ],
            "tf_comb_metrics_omitted": {
                "lift": "degenerate_when_every_transaction_contains_tp73_anchor",
                "reverse_confidence": "identically_one_in_tp73_anchor_transactions",
                "cosine": "redundant_square_root_of_support_in_this_transaction_scope",
                "genome_wide_pair_z": "requires_a_separate_unconditioned_spatial_null",
            },
            "fixed_period_bp": arguments.period,
            "periodicity_inference_status": "descriptive_fixed_period",
            "peak_inference_status": "descriptive_candidate",
            "analysis_parameters": {
                "minimum_itemset_count": arguments.minimum_itemset_count,
                "minimum_blocks": arguments.minimum_blocks,
                "minimum_periodicity_offsets": (
                    arguments.minimum_periodicity_offsets
                ),
                "smoothing_sigma_bp": arguments.smoothing_sigma,
                "minimum_peak_prominence_log_odds": (
                    arguments.minimum_peak_prominence
                ),
                "minimum_peak_separation_bp": (
                    arguments.minimum_peak_separation
                ),
                "maximum_peak_match_shift_bp": (
                    arguments.maximum_peak_match_shift
                ),
                "maximum_peaks_per_direction": arguments.maximum_peaks,
            },
            "software": {
                "python": platform.python_version(),
                "duckdb_cli": executable_version(duckdb),
                "finalizer_script_sha256": sha256(Path(__file__).resolve()),
                "chromosome_builder_script_sha256": sha256(
                    Path(__file__).with_name(
                        "build_tp73_exact_distance_counts.py"
                    ).resolve()
                ),
            },
            "input_inventory": [
                {"path": str(path.resolve()), "sha256": sha256(path.resolve())}
                for path in arguments.inventory
            ],
            "input_block_components": [
                {"path": str(path.resolve()), "sha256": sha256(path.resolve())}
                for path in arguments.block_components
            ],
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "files": files,
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--inventory", type=Path, action="append", required=True)
    result.add_argument(
        "--block-components", type=Path, action="append", required=True
    )
    result.add_argument("--output-dir", type=Path, required=True)
    result.add_argument("--period", type=float, default=10.5)
    result.add_argument("--minimum-itemset-count", type=int, default=20)
    result.add_argument("--minimum-blocks", type=int, default=20)
    result.add_argument("--minimum-periodicity-offsets", type=int, default=20)
    result.add_argument("--smoothing-sigma", type=float, default=1.5)
    result.add_argument("--minimum-peak-prominence", type=float, default=0.05)
    result.add_argument("--minimum-peak-separation", type=float, default=3.0)
    result.add_argument("--maximum-peak-match-shift", type=float, default=16.0)
    result.add_argument("--maximum-peaks", type=int, default=20)
    result.add_argument("--threads", type=int, default=2)
    result.add_argument("--memory-limit", default="16GB")
    result.add_argument("--duckdb", default="duckdb")
    return result


def main() -> int:
    arguments = parser().parse_args()
    positive_values = (
        arguments.period,
        arguments.smoothing_sigma,
        arguments.minimum_peak_separation,
        arguments.maximum_peak_match_shift,
    )
    if any(value <= 0 for value in positive_values):
        print("E: period, smoothing, and distance arguments must be positive", file=sys.stderr)
        return 2
    integer_values = (
        arguments.minimum_itemset_count,
        arguments.minimum_blocks,
        arguments.minimum_periodicity_offsets,
        arguments.maximum_peaks,
        arguments.threads,
    )
    if any(value <= 0 for value in integer_values):
        print("E: count and thread arguments must be positive", file=sys.stderr)
        return 2
    if arguments.minimum_peak_prominence < 0:
        print("E: --minimum-peak-prominence cannot be negative", file=sys.stderr)
        return 2
    try:
        finalize(arguments)
    except (ExactDistanceSummaryError, OSError, ValueError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

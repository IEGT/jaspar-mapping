#!/usr/bin/env python3

"""Build a per-motif scan threshold registry with a physical-locus density cap."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
from pathlib import Path
from typing import Any


class ThresholdError(RuntimeError):
    pass


def sha256_text(value: str) -> str:
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise argparse.ArgumentTypeError("expected a lowercase SHA-256 digest")
    return value


def commit_text(value: str) -> str:
    if len(value) != 40 or any(character not in "0123456789abcdef" for character in value):
        raise argparse.ArgumentTypeError("expected a full lowercase Git commit")
    return value


def finite_number(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise argparse.ArgumentTypeError("expected a finite number")
    return result


def positive_number(value: str) -> float:
    result = finite_number(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("expected a number greater than zero")
    return result


def nonnegative_number(value: str) -> float:
    result = finite_number(value)
    if result < 0:
        raise argparse.ArgumentTypeError("expected a non-negative number")
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def integer_score(value: str, label: str) -> int:
    parsed = finite_number(value)
    rounded = round(parsed)
    if abs(parsed - rounded) > 1e-9:
        raise ThresholdError(f"{label} must be an integer-grid score, got {value!r}")
    return int(rounded)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    try:
        with path.open(encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream, delimiter="\t")
            if reader.fieldnames is None:
                raise ThresholdError(f"TSV has no header: {path}")
            return reader.fieldnames, list(reader)
    except OSError as error:
        raise ThresholdError(f"cannot read {path}: {error}") from error


def read_informative_thresholds(path: Path) -> dict[str, dict[str, Any]]:
    fields, rows = read_tsv(path)
    required = {"motif_id", "recommended_threshold"}
    if not required.issubset(fields):
        raise ThresholdError(f"informative threshold list lacks {sorted(required)}: {path}")
    result: dict[str, dict[str, Any]] = {}
    for row in rows:
        motif_id = row["motif_id"].strip()
        if not motif_id or motif_id in result:
            raise ThresholdError(f"blank or duplicate motif_id in {path}: {motif_id!r}")
        raw = row["recommended_threshold"].strip()
        result[motif_id] = {
            "threshold": None if raw in {"", "NA"} else integer_score(
                raw, f"informative threshold for {motif_id}"
            ),
            "source": "tp73_context_v1" if raw not in {"", "NA"}
            else "default_for_unevaluable",
        }
    return result


def apply_negative_sensitivity(
    thresholds: dict[str, dict[str, Any]], path: Path | None
) -> None:
    if path is None:
        return
    fields, rows = read_tsv(path)
    required = {"motif_id", "best_negative_threshold", "sensitivity_conclusion"}
    if not required.issubset(fields):
        raise ThresholdError(f"negative-sensitivity table lacks {sorted(required)}: {path}")
    expected = {
        motif_id for motif_id, row in thresholds.items() if row["threshold"] == 0
    }
    seen: set[str] = set()
    allowed_conclusions = {
        "negative_threshold_higher_auc", "zero_not_worse",
        "zero_not_evaluable", "no_negative_candidate_evaluable",
    }
    for row in rows:
        motif_id = row["motif_id"].strip()
        if motif_id in seen or motif_id not in thresholds:
            raise ThresholdError(f"unexpected or duplicate sensitivity motif: {motif_id}")
        seen.add(motif_id)
        if thresholds[motif_id]["threshold"] != 0:
            raise ThresholdError(
                f"negative sensitivity was supplied for non-zero v1 motif {motif_id}"
            )
        conclusion = row["sensitivity_conclusion"]
        if conclusion not in allowed_conclusions:
            raise ThresholdError(
                f"unknown sensitivity conclusion for {motif_id}: {conclusion!r}"
            )
        if conclusion == "negative_threshold_higher_auc":
            best = integer_score(
                row["best_negative_threshold"],
                f"negative sensitivity threshold for {motif_id}",
            )
            if best >= 0:
                raise ThresholdError(
                    f"best negative threshold is not negative for {motif_id}"
                )
            thresholds[motif_id] = {
                "threshold": best,
                "source": "negative_threshold_sensitivity_v1",
            }
    if seen != expected:
        raise ThresholdError(
            "negative-sensitivity table does not exactly cover v1 zero motifs; "
            f"missing={sorted(expected - seen)[:10]}, extra={sorted(seen - expected)[:10]}"
        )


def apply_overrides(
    thresholds: dict[str, dict[str, Any]], path: Path | None
) -> None:
    if path is None:
        return
    fields, rows = read_tsv(path)
    required = {"motif_id", "informative_threshold", "informative_source"}
    if not required.issubset(fields):
        raise ThresholdError(f"override table lacks {sorted(required)}: {path}")
    seen: set[str] = set()
    for row in rows:
        motif_id = row["motif_id"].strip()
        if not motif_id or motif_id in seen:
            raise ThresholdError(f"blank or duplicate override motif: {motif_id!r}")
        seen.add(motif_id)
        thresholds[motif_id] = {
            "threshold": integer_score(
                row["informative_threshold"], f"override threshold for {motif_id}"
            ),
            "source": row["informative_source"].strip(),
        }
        if not thresholds[motif_id]["source"]:
            raise ThresholdError(f"override source is blank for {motif_id}")


def read_distribution_manifest(path: Path) -> dict[str, dict[str, str]]:
    fields, rows = read_tsv(path)
    required = {"motif_id", "path", "sha256"}
    if not required.issubset(fields):
        raise ThresholdError(f"distribution manifest lacks {sorted(required)}: {path}")
    result: dict[str, dict[str, str]] = {}
    for row in rows:
        motif_id = row["motif_id"].strip()
        if not motif_id or motif_id in result:
            raise ThresholdError(f"blank or duplicate distribution motif: {motif_id!r}")
        resolved = dict(row)
        distribution_path = Path(row["path"]).expanduser()
        if not distribution_path.is_absolute():
            distribution_path = path.parent / distribution_path
        resolved["path"] = str(distribution_path.resolve())
        result[motif_id] = resolved
    return result


def read_distribution(
    motif_id: str,
    manifest_row: dict[str, str],
    expected_chrom: str,
    expected_score_mode: str,
    expected_pseudocount: float,
) -> dict[str, Any]:
    path = Path(manifest_row["path"]).expanduser().resolve()
    if not path.is_file():
        raise ThresholdError(f"distribution file is absent for {motif_id}: {path}")
    observed_sha = sha256(path)
    if observed_sha != manifest_row["sha256"]:
        raise ThresholdError(f"distribution checksum mismatch for {motif_id}: {path}")
    fields, rows = read_tsv(path)
    required = {
        "MotifID", "Chromosome", "Strand", "ScoreMode", "Pseudocount",
        "BinScheme", "BinWidth", "ValidWindows", "SkippedWindows",
        "ScoreBinStart", "ScoreBinEnd", "BinCount", "OrientationAggregation",
    }
    if not required.issubset(fields) or not rows:
        raise ThresholdError(f"distribution has an incomplete schema for {motif_id}: {path}")

    valid_windows: int | None = None
    skipped_windows: int | None = None
    bins: dict[int, int] = {}
    for row in rows:
        if row["MotifID"] != motif_id:
            raise ThresholdError(f"distribution motif mismatch in {path}")
        if row["Chromosome"] != expected_chrom or row["Strand"] != "both":
            raise ThresholdError(f"distribution chromosome/strand mismatch for {motif_id}")
        if row["ScoreMode"] != expected_score_mode:
            raise ThresholdError(f"distribution score mode mismatch for {motif_id}")
        if abs(float(row["Pseudocount"]) - expected_pseudocount) > 1e-9:
            raise ThresholdError(f"distribution pseudocount mismatch for {motif_id}")
        if row["OrientationAggregation"] != "max_score_per_alignment_span":
            raise ThresholdError(f"distribution does not collapse physical loci for {motif_id}")
        row_valid = int(row["ValidWindows"])
        row_skipped = int(row["SkippedWindows"])
        if valid_windows is None:
            valid_windows = row_valid
            skipped_windows = row_skipped
        elif valid_windows != row_valid or skipped_windows != row_skipped:
            raise ThresholdError(f"distribution totals vary between bins for {motif_id}")
        if row["BinScheme"] == "sentinel":
            continue
        if row["BinScheme"] != "fixed" or abs(float(row["BinWidth"]) - 1.0) > 1e-9:
            raise ThresholdError(f"density calibration requires fixed one-unit bins for {motif_id}")
        start = integer_score(row["ScoreBinStart"], f"score bin for {motif_id}")
        end = integer_score(row["ScoreBinEnd"], f"score bin end for {motif_id}")
        if end != start + 1 or start in bins:
            raise ThresholdError(f"invalid or duplicate score bin for {motif_id}: {start}")
        bins[start] = int(row["BinCount"])

    if valid_windows is None or skipped_windows is None:
        raise ThresholdError(f"distribution totals are absent for {motif_id}")
    if sum(bins.values()) != valid_windows:
        raise ThresholdError(
            f"finite bin counts do not equal valid physical loci for {motif_id}"
        )
    return {
        "path": path,
        "sha256": observed_sha,
        "valid_windows": valid_windows,
        "skipped_windows": skipped_windows,
        "bins": bins,
    }


def retained_at_threshold(bins: dict[int, int], threshold: int) -> int:
    return sum(count for start, count in bins.items() if start >= threshold)


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    columns = [
        "motif_id", "informative_threshold", "informative_source",
        "default_minimum_score", "candidate_minimum_score",
        "density_minimum_spacing_bp", "density_maximum_loci",
        "density_threshold", "final_minimum_score", "density_limited",
        "density_chrom", "valid_locus_starts", "skipped_locus_starts",
        "loci_at_candidate_threshold", "loci_at_final_threshold",
        "mean_spacing_bp_at_final_threshold", "distribution_sha256",
    ]
    if path.exists():
        raise ThresholdError(f"refusing to replace existing threshold TSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("x", encoding="ascii", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row[column] for column in columns})
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def write_json(path: Path, value: dict[str, Any]) -> None:
    if path.exists():
        raise ThresholdError(f"refusing to replace existing threshold JSON: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    with temporary.open("x", encoding="ascii") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def build(arguments: argparse.Namespace) -> None:
    informative_path = arguments.informative_thresholds.expanduser().resolve()
    sensitivity_path = (
        arguments.negative_sensitivity.expanduser().resolve()
        if arguments.negative_sensitivity else None
    )
    override_path = (
        arguments.override_thresholds.expanduser().resolve()
        if arguments.override_thresholds else None
    )
    manifest_path = arguments.distribution_manifest.expanduser().resolve()
    thresholds = read_informative_thresholds(informative_path)
    apply_negative_sensitivity(thresholds, sensitivity_path)
    apply_overrides(thresholds, override_path)
    manifest = read_distribution_manifest(manifest_path)
    if set(thresholds) != set(manifest):
        missing_thresholds = sorted(set(manifest) - set(thresholds))
        missing_distributions = sorted(set(thresholds) - set(manifest))
        raise ThresholdError(
            "threshold/distribution motif sets differ; missing thresholds="
            f"{missing_thresholds[:10]}, missing distributions={missing_distributions[:10]}"
        )

    default_floor = integer_score(
        str(arguments.default_minimum_score), "default minimum score"
    )
    result_rows: list[dict[str, Any]] = []
    for motif_id in sorted(thresholds):
        informative = thresholds[motif_id]["threshold"]
        candidate = default_floor if informative is None else min(default_floor, informative)
        distribution = read_distribution(
            motif_id, manifest[motif_id], arguments.chrom,
            arguments.score_mode, arguments.pseudocount,
        )
        valid = int(distribution["valid_windows"])
        maximum_loci = math.floor(valid / arguments.minimum_spacing_bp)
        bins = distribution["bins"]
        selected = candidate
        retained = retained_at_threshold(bins, selected)
        while retained > maximum_loci:
            selected += 1
            retained = retained_at_threshold(bins, selected)
        candidate_loci = retained_at_threshold(bins, candidate)
        mean_spacing = "NA" if retained == 0 else f"{valid / retained:.9f}"
        result_rows.append(
            {
                "motif_id": motif_id,
                "informative_threshold": "NA" if informative is None else informative,
                "informative_source": thresholds[motif_id]["source"],
                "default_minimum_score": default_floor,
                "candidate_minimum_score": candidate,
                "density_minimum_spacing_bp": f"{arguments.minimum_spacing_bp:.9f}",
                "density_maximum_loci": maximum_loci,
                "density_threshold": selected,
                "final_minimum_score": selected,
                "density_limited": str(selected > candidate).lower(),
                "density_chrom": arguments.chrom,
                "valid_locus_starts": valid,
                "skipped_locus_starts": distribution["skipped_windows"],
                "loci_at_candidate_threshold": candidate_loci,
                "loci_at_final_threshold": retained,
                "mean_spacing_bp_at_final_threshold": mean_spacing,
                "distribution_sha256": distribution["sha256"],
            }
        )

    output_tsv = arguments.output_tsv.expanduser().resolve()
    output_json = arguments.output_json.expanduser().resolve()
    write_tsv(output_tsv, result_rows)
    metadata = {
        "schema_version": 1,
        "threshold_set_id": arguments.threshold_set_id,
        "genome_id": arguments.genome_id,
        "motif_set_id": arguments.motif_set_id,
        "genome_fasta_sha256": arguments.genome_fasta_sha256,
        "density_chrom_sequence_sha256": arguments.density_chrom_sequence_sha256,
        "jaspar_sha256": arguments.jaspar_sha256,
        "source_commit": arguments.source_commit,
        "density_chrom": arguments.chrom,
        "score_mode": arguments.score_mode,
        "pseudocount": arguments.pseudocount,
        "orientation_aggregation": "max_score_per_alignment_span",
        "distribution_bin_width": 1,
        "density_minimum_spacing_bp": arguments.minimum_spacing_bp,
        "default_minimum_score": default_floor,
        "density_definition": (
            "valid physical alignment starts divided by distinct retained starts; "
            "plus/minus orientations at one span use their maximum score"
        ),
        "candidate_formula": "min(informative_threshold, default_minimum_score)",
        "final_formula": "max(candidate_minimum_score, density_threshold)",
        "motifs": len(result_rows),
        "density_limited_motifs": sum(row["density_limited"] == "true" for row in result_rows),
        "total_loci_at_candidate_threshold_density_chrom": sum(
            row["loci_at_candidate_threshold"] for row in result_rows
        ),
        "total_loci_at_final_threshold_density_chrom": sum(
            row["loci_at_final_threshold"] for row in result_rows
        ),
        "maximum_orientation_rows_at_final_threshold_density_chrom": 2 * sum(
            row["loci_at_final_threshold"] for row in result_rows
        ),
        "informative_thresholds_sha256": sha256(informative_path),
        "negative_sensitivity_sha256": sha256(sensitivity_path) if sensitivity_path else None,
        "override_thresholds_sha256": sha256(override_path) if override_path else None,
        "distribution_manifest_sha256": sha256(manifest_path),
        "threshold_tsv_sha256": sha256(output_tsv),
    }
    write_json(output_json, metadata)
    print(f"I: Wrote {len(result_rows)} density-capped motif thresholds: {output_tsv}",
          file=sys.stderr)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--informative-thresholds", required=True, type=Path)
    result.add_argument("--negative-sensitivity", type=Path)
    result.add_argument("--override-thresholds", type=Path)
    result.add_argument("--distribution-manifest", required=True, type=Path)
    result.add_argument("--output-tsv", required=True, type=Path)
    result.add_argument("--output-json", required=True, type=Path)
    result.add_argument("--threshold-set-id", required=True)
    result.add_argument("--genome-id", required=True)
    result.add_argument("--motif-set-id", required=True)
    result.add_argument("--genome-fasta-sha256", type=sha256_text)
    result.add_argument("--density-chrom-sequence-sha256", type=sha256_text)
    result.add_argument("--jaspar-sha256", type=sha256_text)
    result.add_argument("--source-commit", type=commit_text)
    result.add_argument("--chrom", default="1")
    result.add_argument(
        "--default-minimum-score", type=finite_number, default=-1.0,
        help="permissive fallback compared with each informative threshold (default: -1)",
    )
    result.add_argument(
        "--minimum-spacing-bp", type=positive_number, default=200.0,
        help="minimum chromosome-average valid bp per retained locus (default: 200)",
    )
    result.add_argument(
        "--score-mode", choices=("log2_relative_risk", "log_odds"),
        default="log2_relative_risk",
    )
    result.add_argument("--pseudocount", type=nonnegative_number, default=1.0)
    return result


def main() -> int:
    try:
        build(parser().parse_args())
    except (OSError, ValueError, ThresholdError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

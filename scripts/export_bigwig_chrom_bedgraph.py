#!/usr/bin/env python3

"""Export one BigWig chromosome as positive, non-overlapping bedGraph rows."""

from __future__ import annotations

import argparse
import math
import os
import sys
from pathlib import Path


def normalize_chrom(value: str) -> str:
    return value[3:] if value.lower().startswith("chr") else value


def export(arguments: argparse.Namespace) -> None:
    try:
        import pyBigWig  # type: ignore[import-not-found]
    except ImportError as error:
        raise RuntimeError(
            "pyBigWig is required; install the conda-forge pybigwig package"
        ) from error

    source = arguments.input.expanduser().resolve()
    output = arguments.output.expanduser().resolve()
    if not source.is_file():
        raise RuntimeError(f"BigWig input not found: {source}")
    if output.exists():
        raise RuntimeError(f"refusing to replace existing output: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    with pyBigWig.open(str(source)) as bigwig:
        matching = [
            chrom for chrom in bigwig.chroms()
            if normalize_chrom(chrom).lower() == normalize_chrom(arguments.chrom).lower()
        ]
        if len(matching) != 1:
            raise RuntimeError(
                f"BigWig has {len(matching)} chromosomes matching {arguments.chrom}"
            )
        chrom = matching[0]
        intervals = bigwig.intervals(chrom) or ()
        temporary = output.with_name(f".{output.name}.tmp-{os.getpid()}")
        written = 0
        try:
            with temporary.open("w", encoding="ascii", newline="") as handle:
                for start, end, value in intervals:
                    if math.isfinite(value) and value > 0:
                        handle.write(f"{chrom}\t{start}\t{end}\t{value:.17g}\n")
                        written += 1
            if written == 0:
                raise RuntimeError(
                    f"BigWig has no positive intervals on chromosome {chrom}"
                )
            os.replace(temporary, output)
        finally:
            if temporary.exists():
                temporary.unlink()
    print(
        f"I: exported {written} positive {chrom} intervals: {output}",
        file=sys.stderr,
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(
        description=(
            "Export one chromosome from BigWig to sorted positive-depth "
            "bedGraph using 0-based half-open coordinates."
        )
    )
    result.add_argument("--input", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--chrom", default="1")
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        export(arguments)
    except (OSError, RuntimeError, ValueError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

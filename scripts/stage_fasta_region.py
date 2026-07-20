#!/usr/bin/env python3

"""Copy one indexed FASTA sequence into a self-contained scratch FASTA."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Callable


class StageError(RuntimeError):
    pass


@dataclass(frozen=True)
class FastaIndexRecord:
    name: str
    length: int
    offset: int
    line_bases: int
    line_bytes: int


def parse_sha256(value: str) -> str:
    if not re.fullmatch(r"[0-9a-f]{64}", value):
        raise argparse.ArgumentTypeError("expected a lowercase SHA-256 digest")
    return value


def read_fasta_index(path: Path) -> dict[str, FastaIndexRecord]:
    records: dict[str, FastaIndexRecord] = {}
    with path.open(encoding="ascii") as stream:
        for line_number, line in enumerate(stream, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise StageError(f"invalid FASTA index line {line_number}: {path}")
            try:
                record = FastaIndexRecord(
                    name=fields[0],
                    length=int(fields[1]),
                    offset=int(fields[2]),
                    line_bases=int(fields[3]),
                    line_bytes=int(fields[4]),
                )
            except ValueError as error:
                raise StageError(
                    f"invalid numeric FASTA index field at line {line_number}: {path}"
                ) from error
            if record.name in records:
                raise StageError(f"duplicate FASTA index record: {record.name}")
            if min(record.length, record.offset, record.line_bases, record.line_bytes) < 0:
                raise StageError(f"negative FASTA index field for {record.name}")
            if record.length == 0 or record.line_bases == 0 \
                    or record.line_bytes < record.line_bases:
                raise StageError(f"invalid FASTA index geometry for {record.name}")
            records[record.name] = record
    if not records:
        raise StageError(f"empty FASTA index: {path}")
    return records


def _write_wrapped(output: BinaryIO, bases: bytes, carry: bytes,
                   width: int) -> bytes:
    data = carry + bases
    complete = len(data) // width * width
    if complete:
        lines = (data[offset:offset + width] for offset in range(0, complete, width))
        output.write(b"\n".join(lines))
        output.write(b"\n")
    return data[complete:]


def stage_fasta_region(
    genome: Path,
    fasta_index: Path,
    sequence: str,
    output: Path,
    *,
    expected_length: int | None = None,
    expected_sha256: str | None = None,
    line_width: int = 60,
    chunk_bytes: int = 8 * 1024 * 1024,
    progress_callback: Callable[[int, int], None] | None = None,
) -> dict[str, object]:
    """Stage one sequence and return its verified identity and local paths."""

    if genome.suffix == ".gz":
        raise StageError(
            "indexed scratch staging currently requires an uncompressed FASTA; "
            "decompress the immutable genome resource before planning the scan"
        )
    if line_width <= 0 or chunk_bytes <= 0:
        raise StageError("line width and chunk size must be positive")
    records = read_fasta_index(fasta_index)
    if sequence not in records:
        raise StageError(f"sequence {sequence!r} is absent from {fasta_index}")
    record = records[sequence]
    if expected_length is not None and record.length != expected_length:
        raise StageError(
            f"sequence length mismatch for {sequence}: index has {record.length}, "
            f"plan expects {expected_length}"
        )

    output = output.expanduser().resolve()
    output_index = Path(f"{output}.fai")
    if output.exists() or output_index.exists():
        raise StageError(f"refusing to replace staged FASTA or index: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.partial-{os.getpid()}")
    temporary_index = Path(f"{temporary}.fai")
    if temporary.exists() or temporary_index.exists():
        raise StageError(f"staging path already exists: {temporary}")

    digest = hashlib.sha256()
    copied = 0
    carry = b""
    header = f">{sequence}\n".encode("ascii")
    try:
        with genome.open("rb") as source, temporary.open("xb") as destination:
            source.seek(record.offset)
            destination.write(header)
            while copied < record.length:
                raw = source.read(chunk_bytes)
                if not raw:
                    raise StageError(
                        f"FASTA ended after {copied} of {record.length} bases for {sequence}"
                    )
                bases = raw.replace(b"\n", b"").replace(b"\r", b"")
                bases = bases[:record.length - copied]
                if not bases:
                    continue
                digest.update(bases.upper())
                carry = _write_wrapped(destination, bases, carry, line_width)
                copied += len(bases)
                if progress_callback is not None:
                    progress_callback(copied, record.length)
            if carry:
                destination.write(carry)
                destination.write(b"\n")
            destination.flush()
            os.fsync(destination.fileno())

        observed_sha256 = digest.hexdigest()
        if expected_sha256 is not None and observed_sha256 != expected_sha256:
            raise StageError(
                f"sequence checksum mismatch for {sequence}: staged {observed_sha256}, "
                f"plan expects {expected_sha256}"
            )
        staged_line_bases = min(line_width, record.length)
        with temporary_index.open("x", encoding="ascii", newline="\n") as index:
            index.write(
                f"{sequence}\t{record.length}\t{len(header)}\t"
                f"{staged_line_bases}\t{staged_line_bases + 1}\n"
            )
            index.flush()
            os.fsync(index.fileno())
        os.replace(temporary, output)
        os.replace(temporary_index, output_index)
    except Exception:
        # Partial files deliberately remain in scratch for diagnosis. The unique
        # job directory is removed later by the cluster's scratch policy.
        raise

    staged_metadata = output.stat()
    index_metadata = output_index.stat()
    index_sha256 = hashlib.sha256(output_index.read_bytes()).hexdigest()
    return {
        "sequence": sequence,
        "length": record.length,
        "sequence_sha256": observed_sha256,
        "source_fasta": str(genome.resolve()),
        "source_fasta_index": str(fasta_index.resolve()),
        "staged_fasta": str(output),
        "staged_fasta_index": str(output_index),
        "staged_bytes": staged_metadata.st_size,
        "staged_mtime_ns": staged_metadata.st_mtime_ns,
        "staged_fasta_index_bytes": index_metadata.st_size,
        "staged_fasta_index_mtime_ns": index_metadata.st_mtime_ns,
        "staged_fasta_index_sha256": index_sha256,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Use a FASTA .fai byte offset to copy one sequence into a small, "
            "self-contained FASTA and create the .fai required by pssm_scan."
        )
    )
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--fasta-index", required=True, type=Path)
    parser.add_argument("--sequence", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-length", type=int)
    parser.add_argument("--expected-sha256", type=parse_sha256)
    parser.add_argument("--metadata", type=Path, help="optional JSON result path")
    arguments = parser.parse_args()

    try:
        result = stage_fasta_region(
            arguments.genome,
            arguments.fasta_index,
            arguments.sequence,
            arguments.output,
            expected_length=arguments.expected_length,
            expected_sha256=arguments.expected_sha256,
        )
        if arguments.metadata is not None:
            arguments.metadata.parent.mkdir(parents=True, exist_ok=True)
            with arguments.metadata.open("x", encoding="utf-8") as stream:
                json.dump(result, stream, indent=2, sort_keys=True)
                stream.write("\n")
        print(json.dumps(result, sort_keys=True))
        return 0
    except (OSError, StageError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

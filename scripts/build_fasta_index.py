#!/usr/bin/env python3

"""Build a samtools-compatible .fai index for an uncompressed FASTA file."""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
from pathlib import Path
from typing import BinaryIO, TextIO

ASCII_WHITESPACE = b" \t\r\n\v\f"


class FastaIndexError(RuntimeError):
    pass


def without_line_ending(raw_line: bytes) -> bytes:
    if raw_line.endswith(b"\n"):
        raw_line = raw_line[:-1]
        if raw_line.endswith(b"\r"):
            raw_line = raw_line[:-1]
    return raw_line


def write_record(output: TextIO, name: str | None, length: int,
                 offset: int, bases_per_line: int, bytes_per_line: int) -> None:
    if name is None:
        return
    if length == 0 or bases_per_line == 0:
        raise FastaIndexError(f"FASTA record {name!r} has no sequence")
    output.write(
        f"{name}\t{length}\t{offset}\t{bases_per_line}\t{bytes_per_line}\n"
    )


def build_index(fasta: BinaryIO, output: TextIO) -> None:
    names: set[str] = set()
    name: str | None = None
    length = 0
    offset = 0
    bases_per_line = 0
    bytes_per_line = 0
    saw_short_line = False

    while True:
        line_offset = fasta.tell()
        raw_line = fasta.readline()
        if not raw_line:
            break

        if raw_line.startswith(b">"):
            write_record(output, name, length, offset, bases_per_line, bytes_per_line)
            header = without_line_ending(raw_line[1:]).split(None, 1)
            if not header:
                raise FastaIndexError(f"empty FASTA header at byte {line_offset}")
            try:
                name = header[0].decode("utf-8")
            except UnicodeDecodeError as error:
                raise FastaIndexError(
                    f"non-UTF-8 FASTA name at byte {line_offset}"
                ) from error
            if name in names:
                raise FastaIndexError(f"duplicate FASTA record name {name!r}")
            names.add(name)
            length = 0
            offset = fasta.tell()
            bases_per_line = 0
            bytes_per_line = 0
            saw_short_line = False
            continue

        if name is None:
            raise FastaIndexError(
                f"sequence data precedes the first FASTA header at byte {line_offset}"
            )

        sequence = without_line_ending(raw_line)
        if not sequence:
            raise FastaIndexError(f"blank sequence line in FASTA record {name!r}")
        if sequence.translate(None, ASCII_WHITESPACE) != sequence:
            raise FastaIndexError(
                f"whitespace inside a sequence line in FASTA record {name!r}"
            )

        line_bases = len(sequence)
        line_bytes = len(raw_line)
        if bases_per_line == 0:
            bases_per_line = line_bases
            bytes_per_line = line_bytes
        else:
            if saw_short_line:
                raise FastaIndexError(
                    f"sequence line follows a short final line in FASTA record {name!r}"
                )
            if line_bases > bases_per_line:
                raise FastaIndexError(
                    f"sequence line exceeds the first line width in FASTA record {name!r}"
                )
            if line_bases == bases_per_line and line_bytes != bytes_per_line:
                if raw_line.endswith(b"\n"):
                    raise FastaIndexError(
                        f"inconsistent line endings in FASTA record {name!r}"
                    )
                saw_short_line = True
            if line_bases < bases_per_line:
                saw_short_line = True

        length += line_bases

    if name is None:
        raise FastaIndexError("FASTA file contains no records")
    write_record(output, name, length, offset, bases_per_line, bytes_per_line)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a samtools-compatible .fai for an uncompressed FASTA"
    )
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--output", "-o", type=Path)
    args = parser.parse_args()

    fasta_path = args.fasta
    output_path = args.output or Path(f"{fasta_path}.fai")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    temporary_name: str | None = None
    try:
        with fasta_path.open("rb") as fasta, tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            prefix=f".{output_path.name}.",
            dir=output_path.parent,
            delete=False,
        ) as temporary:
            temporary_name = temporary.name
            build_index(fasta, temporary)
            temporary.flush()
            os.fsync(temporary.fileno())
        current_umask = os.umask(0)
        os.umask(current_umask)
        os.chmod(temporary_name, 0o666 & ~current_umask)
        os.replace(temporary_name, output_path)
        temporary_name = None
    except (OSError, FastaIndexError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    finally:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass

    print(f"I: Wrote FASTA index to {output_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

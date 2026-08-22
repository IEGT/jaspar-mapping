#!/usr/bin/env python3
"""Build a normalized, queryable catalog from official JASPAR metadata."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any


REQUIRED_COLUMNS = (
    "collection",
    "tax_group",
    "matrix_id",
    "base_id",
    "version",
    "name",
    "class",
    "family",
    "uniprot_ids",
    "validation",
    "comment",
    "source",
    "type",
    "tax_id",
    "species",
)
MATRIX_ID_PATTERN = re.compile(r"^[A-Z]{2}[0-9]+\.[0-9]+$")


class CatalogError(RuntimeError):
    """Raised when the upstream metadata or requested catalog is invalid."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def parse_pssm_headers(path: Path) -> list[dict[str, Any]]:
    motifs: list[dict[str, Any]] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip("\r\n")
            fields = header.split("\t", 1)
            matrix_id = fields[0].strip()
            motif_name = fields[1].strip() if len(fields) == 2 else ""
            if not MATRIX_ID_PATTERN.fullmatch(matrix_id):
                raise CatalogError(
                    f"invalid matrix ID at {path}:{line_number}: {matrix_id!r}"
                )
            if matrix_id in seen:
                raise CatalogError(f"duplicate matrix ID in PSSM: {matrix_id}")
            seen.add(matrix_id)
            motifs.append(
                {
                    "matrix_id": matrix_id,
                    "motif_name_in_pssm": motif_name,
                    "motif_order": len(motifs) + 1,
                }
            )
    if not motifs:
        raise CatalogError(f"no JASPAR matrix headers found in {path}")
    return motifs


def parse_metadata(
    path: Path,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    matrices: list[dict[str, Any]] = []
    species_rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != list(REQUIRED_COLUMNS):
            raise CatalogError(
                "unexpected metadata columns: "
                f"{reader.fieldnames!r}; expected {list(REQUIRED_COLUMNS)!r}"
            )
        for line_number, row in enumerate(reader, 2):
            matrix_id = row["matrix_id"].strip()
            if not MATRIX_ID_PATTERN.fullmatch(matrix_id):
                raise CatalogError(
                    f"invalid matrix ID at {path}:{line_number}: {matrix_id!r}"
                )
            if matrix_id in seen:
                raise CatalogError(f"duplicate metadata matrix ID: {matrix_id}")
            seen.add(matrix_id)

            tax_ids = row["tax_id"].split("::") if row["tax_id"] else []
            species = row["species"].split("::") if row["species"] else []
            if len(tax_ids) != len(species):
                raise CatalogError(
                    f"tax ID/species cardinality differs for {matrix_id}: "
                    f"{len(tax_ids)} != {len(species)}"
                )
            normalized_tax_ids: list[int] = []
            for ordinal, (tax_id_text, species_name) in enumerate(
                zip(tax_ids, species), 1
            ):
                try:
                    tax_id = int(tax_id_text)
                except ValueError as error:
                    raise CatalogError(
                        f"invalid tax ID for {matrix_id}: {tax_id_text!r}"
                    ) from error
                normalized_tax_ids.append(tax_id)
                species_rows.append(
                    {
                        "matrix_id": matrix_id,
                        "species_ordinal": ordinal,
                        "tax_id": tax_id,
                        "species_scientific_name": species_name,
                    }
                )

            matrices.append(
                {
                    "collection": row["collection"],
                    "tax_group": row["tax_group"],
                    "matrix_id": matrix_id,
                    "base_id": row["base_id"],
                    "version": int(row["version"]),
                    "name": row["name"],
                    "class": row["class"] or None,
                    "family": row["family"] or None,
                    "uniprot_ids": row["uniprot_ids"] or None,
                    "validation": row["validation"] or None,
                    "comment": row["comment"] or None,
                    "source": row["source"] or None,
                    "type": row["type"] or None,
                    "species_count": len(species),
                    "includes_homo_sapiens": 9606 in normalized_tax_ids,
                    "includes_mus_musculus": 10090 in normalized_tax_ids,
                    "includes_rattus_norvegicus": 10116 in normalized_tax_ids,
                }
            )
    if not matrices:
        raise CatalogError(f"metadata table is empty: {path}")
    return matrices, species_rows


def write_json_lines(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("x", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, ensure_ascii=True, sort_keys=True) + "\n")


def run_duckdb(duckdb: str, database: Path, sql: str) -> None:
    process = subprocess.run(
        [duckdb, "-batch", str(database)],
        input=sql,
        text=True,
        capture_output=True,
        check=False,
    )
    if process.returncode != 0:
        raise CatalogError(process.stderr.strip() or "DuckDB catalog build failed")


def build(arguments: argparse.Namespace) -> None:
    metadata = arguments.metadata.expanduser().resolve()
    pssm = arguments.pssm.expanduser().resolve()
    output = arguments.output.expanduser().resolve()
    for path in (metadata, pssm):
        if not path.is_file():
            raise CatalogError(f"input file is missing: {path}")
    if output.exists():
        raise CatalogError(f"refusing to replace existing output: {output}")
    if shutil.which(arguments.duckdb) is None:
        raise CatalogError(f"DuckDB executable not found: {arguments.duckdb}")

    metadata_sha256 = sha256(metadata)
    if (
        arguments.expected_metadata_sha256
        and metadata_sha256 != arguments.expected_metadata_sha256
    ):
        raise CatalogError(
            "metadata SHA-256 differs: "
            f"{metadata_sha256} != {arguments.expected_metadata_sha256}"
        )

    matrices, species_rows = parse_metadata(metadata)
    matrix_by_id = {row["matrix_id"]: row for row in matrices}
    membership = parse_pssm_headers(pssm)
    missing = [
        row["matrix_id"] for row in membership if row["matrix_id"] not in matrix_by_id
    ]
    if missing:
        raise CatalogError(
            f"metadata lacks {len(missing)} PSSM matrices; first={missing[:5]}"
        )
    if arguments.expected_motif_count is not None and (
        len(membership) != arguments.expected_motif_count
    ):
        raise CatalogError(
            f"PSSM motif count differs: {len(membership)} != "
            f"{arguments.expected_motif_count}"
        )
    for row in membership:
        matrix = matrix_by_id[row["matrix_id"]]
        row["motif_set_id"] = arguments.motif_set_id
        row["name_matches_metadata"] = row["motif_name_in_pssm"] == matrix["name"]

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".jaspar-catalog-", dir=output.parent))
    try:
        source_dir = staging / "source"
        table_root = staging / "tables" / "jaspar2026"
        source_dir.mkdir()
        table_root.mkdir(parents=True)
        json_paths = {
            "jaspar_matrix": source_dir / "jaspar_matrix.jsonl",
            "jaspar_matrix_species": source_dir / "jaspar_matrix_species.jsonl",
            "jaspar_motif_set_matrix": source_dir / "jaspar_motif_set_matrix.jsonl",
        }
        write_json_lines(json_paths["jaspar_matrix"], matrices)
        write_json_lines(json_paths["jaspar_matrix_species"], species_rows)
        write_json_lines(json_paths["jaspar_motif_set_matrix"], membership)

        table_paths = {
            name: table_root / name / "part-000000.parquet" for name in json_paths
        }
        for path in table_paths.values():
            path.parent.mkdir()
        database = staging / "jaspar_metadata.duckdb"
        sql_parts = ["SET preserve_insertion_order=false;"]
        for name, source in json_paths.items():
            sql_parts.append(
                f"CREATE TABLE {name} AS SELECT * FROM "
                f"read_json_auto({sql_string(source)}, format='newline_delimited');"
            )
        sql_parts.extend(
            [
                "CREATE UNIQUE INDEX jaspar_matrix_id_idx "
                "ON jaspar_matrix(matrix_id);",
                "CREATE UNIQUE INDEX jaspar_matrix_species_idx "
                "ON jaspar_matrix_species(matrix_id, species_ordinal);",
                "CREATE UNIQUE INDEX jaspar_motif_set_matrix_idx "
                "ON jaspar_motif_set_matrix(motif_set_id, matrix_id);",
            ]
        )
        for name, destination in table_paths.items():
            sql_parts.append(
                f"COPY {name} TO {sql_string(destination)} "
                "(FORMAT PARQUET, COMPRESSION ZSTD);"
            )
        sql_parts.extend(
            [
                f"COPY jaspar_matrix TO {sql_string(staging / 'jaspar_matrix.tsv')} "
                "(FORMAT CSV, DELIMITER '\\t', HEADER, NULL 'NA');",
                f"COPY jaspar_matrix_species TO "
                f"{sql_string(staging / 'jaspar_matrix_species.tsv')} "
                "(FORMAT CSV, DELIMITER '\\t', HEADER, NULL 'NA');",
                "CHECKPOINT;",
            ]
        )
        run_duckdb(arguments.duckdb, database, "\n".join(sql_parts) + "\n")
        shutil.rmtree(source_dir)

        selected_tax_groups: dict[str, int] = {}
        selected_human = 0
        for row in membership:
            matrix = matrix_by_id[row["matrix_id"]]
            tax_group = str(matrix["tax_group"])
            selected_tax_groups[tax_group] = selected_tax_groups.get(tax_group, 0) + 1
            selected_human += int(bool(matrix["includes_homo_sapiens"]))
        manifest = {
            "schema_version": 1,
            "jaspar_release": arguments.jaspar_release,
            "motif_set_id": arguments.motif_set_id,
            "metadata_source_url": arguments.source_url,
            "metadata_source_file": metadata.name,
            "metadata_source_sha256": metadata_sha256,
            "metadata_source_bytes": metadata.stat().st_size,
            "pssm_source_file": pssm.name,
            "pssm_source_sha256": sha256(pssm),
            "metadata_matrix_count": len(matrices),
            "metadata_species_row_count": len(species_rows),
            "motif_set_matrix_count": len(membership),
            "motif_set_human_source_count": selected_human,
            "motif_set_tax_group_counts": dict(sorted(selected_tax_groups.items())),
            "pssm_metadata_name_mismatch_count": sum(
                not row["name_matches_metadata"] for row in membership
            ),
            "generated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "database": database.name,
            "database_sha256": sha256(database),
            "tables": {
                name: {
                    "path": str(path.relative_to(staging)),
                    "rows": (
                        len(matrices)
                        if name == "jaspar_matrix"
                        else len(species_rows)
                        if name == "jaspar_matrix_species"
                        else len(membership)
                    ),
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for name, path in table_paths.items()
            },
        }
        (staging / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)

    print(
        f"I: Built JASPAR {arguments.jaspar_release} metadata catalog with "
        f"{len(membership)} motif-set matrices: {output}",
        file=sys.stderr,
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--metadata", type=Path, required=True)
    result.add_argument("--pssm", type=Path, required=True)
    result.add_argument("--output", type=Path, required=True)
    result.add_argument("--duckdb", default="duckdb")
    result.add_argument("--jaspar-release", default="2026")
    result.add_argument(
        "--motif-set-id", default="jaspar2026_core_nonredundant"
    )
    result.add_argument(
        "--source-url",
        default=(
            "https://mencius.uio.no/JASPAR/JASPAR_metadata/2026/"
            "ultimate_metadata_table_CORE.tsv"
        ),
    )
    result.add_argument("--expected-metadata-sha256")
    result.add_argument("--expected-motif-count", type=int)
    return result


def main() -> int:
    arguments = parser().parse_args()
    try:
        build(arguments)
    except (CatalogError, OSError, ValueError, json.JSONDecodeError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

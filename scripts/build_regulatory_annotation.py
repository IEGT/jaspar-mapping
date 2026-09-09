#!/usr/bin/env python3
"""Import pinned Ensembl GFF and annotate physical TP73 anchors without rescanning.

The regulatory package is independent of the per-TSS promoter dimension.
GFF starts are converted once from 1-based inclusive to BED half-open. Missing
extended bounds stay NULL. Native gene ownership is not a GTF-derived link.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
from urllib.parse import unquote


SUBSETS = (
    "all", "promoter_core", "promoter_extended", "outside_promoter_extended",
    "enhancer", "open_chromatin", "no_regulatory_overlap", "tss_window_promoter",
)


def sql_string(value: object) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sha256(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(chunk)
    return result.hexdigest()


def parquet(path: Path) -> str:
    return f"read_parquet({sql_string(path)}, hive_partitioning=false)"


def run_sql(args: argparse.Namespace, sql: str, temporary: Path) -> None:
    settings = (
        f"SET threads={args.threads}; SET memory_limit={sql_string(args.memory_limit)};"
        f"SET temp_directory={sql_string(temporary / 'spill')};"
        "SET preserve_insertion_order=false;"
    )
    result = subprocess.run(
        [args.duckdb, "-batch", "-bail", ":memory:"], input=settings + sql,
        text=True, capture_output=True, check=False,
    )
    if result.returncode:
        raise ValueError(result.stderr.strip())


def copy_table(name: str, directory: Path) -> str:
    return (f"COPY {name} TO {sql_string(directory / (name + '.parquet'))} "
            "(FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 131072);\n")


def publish(staging: Path, output: Path, manifest: dict) -> None:
    files = sorted(staging.glob("*.parquet"))  # Only this small new package, not input discovery.
    manifest["files"] = [{"path": f.name, "bytes": f.stat().st_size,
                          "sha256": sha256(f)} for f in files]
    manifest["builder_sha256"] = sha256(Path(__file__))
    (staging / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (staging / "schema.sql").write_text("\n".join(
        f"CREATE OR REPLACE VIEW {f.stem} AS SELECT * FROM "
        f"read_parquet('{f.name}', hive_partitioning=false);" for f in files
    ) + "\n")
    if output.exists():
        raise ValueError(f"refusing to replace existing package: {output}")
    os.rename(staging, output)


def validated_package(path: Path) -> dict:
    manifest = json.loads((path / "manifest.json").read_text())
    for record in manifest["files"]:
        relative = Path(record["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("invalid inventory path")
        file = path / relative
        if file.stat().st_size != record["bytes"] or sha256(file) != record["sha256"]:
            raise ValueError(f"package checksum changed: {file}")
    return manifest


def gff_rows(path: Path, assembly: str, release: str, digest: str):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as stream:
        for number, line in enumerate(stream, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GFF line {number}: expected nine fields")
            chrom, source, kind, first, last, _, strand, _, attributes = fields
            attrs = {}
            for item in attributes.split(";"):
                if not item:
                    continue
                key, separator, value = item.partition("=")
                if not separator or key in attrs:
                    raise ValueError(f"GFF line {number}: malformed/duplicate attribute")
                attrs[key] = value
            source_id = unquote(attrs.get("ID", ""))
            if not source_id or not re.fullmatch(r"[A-Za-z0-9_.-]+", chrom):
                raise ValueError(f"GFF line {number}: missing ID or unsafe chromosome")
            start, end = int(first) - 1, int(last)
            if start < 0 or end <= start or strand not in ("+", "-", ".", "?"):
                raise ValueError(f"GFF line {number}: invalid interval/strand")
            has_extended = "extended_start" in attrs or "extended_end" in attrs
            if has_extended and not {"extended_start", "extended_end"} <= attrs.keys():
                raise ValueError(f"GFF line {number}: incomplete extended bounds")
            ext_start = int(attrs["extended_start"]) - 1 if has_extended else None
            ext_end = int(attrs["extended_end"]) if has_extended else None
            if has_extended and not (0 <= ext_start <= start < end <= ext_end):
                raise ValueError(f"GFF line {number}: core not contained in extension")
            key = json.dumps(["ensembl_regulation", assembly, release, digest, source_id])
            yield {
                "regulatory_feature_id": hashlib.sha256(key.encode()).hexdigest(),
                "source_feature_id": source_id, "chrom": chrom, "feature_type": kind,
                "start": start, "end": end, "strand": strand,
                "core_start": start if kind == "promoter" else None,
                "core_end": end if kind == "promoter" else None,
                "extended_start": ext_start, "extended_end": ext_end,
                "source": source, "attributes_json": json.dumps(attrs, sort_keys=True),
                # Split before URL decoding: escaped commas belong to a single identifier.
                "native_gene_ids": [unquote(x) for x in attrs.get("gene_id", "").split(",") if x],
            }


def import_gff(args: argparse.Namespace, staging: Path) -> None:
    digest = sha256(args.gff)
    if args.expected_sha256 and digest != args.expected_sha256:
        raise ValueError("GFF checksum does not match --expected-sha256")
    audit = None
    if args.coordinate_audit:
        audit = json.loads(args.coordinate_audit.read_text())
        if (audit.get("kind") != "regulatory_coordinate_audit"
                or audit.get("status") != "passed"
                or audit.get("gff_sha256") != digest
                or audit.get("assembly") != args.assembly):
            raise ValueError("coordinate audit does not validate this assembly/GFF payload")
    # Bounded streaming adapter, not a retained BED/TSV motif-data layer.
    with tempfile.TemporaryDirectory(prefix="regulatory-import-", dir=staging.parent) as temp:
        temporary = Path(temp)
        rows_path = temporary / "features.jsonl"
        rows = 0
        chromosomes = set()
        with rows_path.open("w") as stream:
            for row in gff_rows(args.gff, args.assembly, args.resolved_release, digest):
                stream.write(json.dumps(row) + "\n")
                rows += 1
                chromosomes.add(row["chrom"])
        if not rows:
            raise ValueError("empty GFF")
        columns = """{
          regulatory_feature_id:'VARCHAR', source_feature_id:'VARCHAR',
          chrom:'VARCHAR', feature_type:'VARCHAR', start:'BIGINT', 'end':'BIGINT',
          strand:'VARCHAR', core_start:'BIGINT', core_end:'BIGINT',
          extended_start:'BIGINT', extended_end:'BIGINT', source:'VARCHAR',
          attributes_json:'VARCHAR', native_gene_ids:'VARCHAR[]'}"""
        sql = f"""
CREATE TABLE imported AS SELECT * FROM read_json(
  {sql_string(rows_path)}, format='newline_delimited', columns={columns});
SELECT CASE WHEN EXISTS (SELECT source_feature_id FROM imported
  GROUP BY source_feature_id HAVING count(*) <> 1)
  THEN error('duplicate regulatory source ID') END;
CREATE TABLE regulatory_feature AS SELECT * EXCLUDE(native_gene_ids),
  {sql_string(args.assembly)}::VARCHAR AS assembly,
  {sql_string(args.requested_release)}::VARCHAR AS requested_regulatory_release,
  {sql_string(args.resolved_release)}::VARCHAR AS regulatory_release,
  {sql_string(args.source_uri)}::VARCHAR AS source_uri,
  {sql_string(digest)}::VARCHAR AS source_sha256
FROM imported ORDER BY chrom, start, "end", source_feature_id;
CREATE TABLE regulatory_feature_gene AS
SELECT DISTINCT regulatory_feature_id, unnest(native_gene_ids)::VARCHAR AS gene_id,
  'ensembl_native'::VARCHAR AS link_source,
  {sql_string(args.resolved_release)}::VARCHAR AS annotation_release
FROM imported;
"""
        sql += copy_table("regulatory_feature", staging)
        sql += copy_table("regulatory_feature_gene", staging)
        run_sql(args, sql, temporary)
    publish(staging, args.output, {
        "schema_version": 1, "kind": "ensembl_regulatory_annotation",
        "assembly": args.assembly, "requested_regulatory_release": args.requested_release,
        "regulatory_release": args.resolved_release, "source_uri": args.source_uri,
        "source_sha256": digest, "source_format": "gff_1based_inclusive",
        "coordinate_mode": "bed_0based_half_open", "feature_count": rows,
        "chromosomes": sorted(chromosomes),
        "coordinate_audit_note": args.coordinate_audit_note,
        "coordinate_audit": audit,
        "activity_in_experimental_cells": "not_assessed",
    })


def annotate(args: argparse.Namespace, staging: Path) -> None:
    manifest = validated_package(args.features)
    if manifest["assembly"] != args.assembly:
        raise ValueError("regulatory feature assembly does not match anchor assembly")
    if args.chrom not in manifest["chromosomes"]:
        raise ValueError("chromosome absent from regulatory export; do not infer negative membership")
    if bool(args.tss) != bool(args.transcript_tss):
        raise ValueError("--tss and --transcript-tss must be supplied together")
    feature_file = args.features / "regulatory_feature.parquet"
    native_file = args.features / "regulatory_feature_gene.parquet"
    inputs = {"anchors": args.anchors, "features_manifest": args.features / "manifest.json"}
    # Both inputs are restricted to the same explicit chromosome before the
    # interval join. A redundant chrom equality makes DuckDB choose a huge
    # chromosome hash bucket instead of its inequality-join implementation.
    sql = f"""
CREATE TABLE anchor AS SELECT CAST(chrom AS VARCHAR) AS chrom,
  anchor_start::BIGINT AS anchor_start, anchor_end::BIGINT AS anchor_end
FROM {parquet(args.anchors)} WHERE CAST(chrom AS VARCHAR)={sql_string(args.chrom)};
SELECT CASE WHEN NOT EXISTS (SELECT 1 FROM anchor) OR EXISTS (
  SELECT 1 FROM anchor WHERE chrom IS NULL OR anchor_start IS NULL
    OR anchor_end IS NULL OR anchor_start<0 OR anchor_end<=anchor_start)
  OR (SELECT count(*)-count(DISTINCT (chrom,anchor_start,anchor_end)) FROM anchor)<>0
  THEN error('anchor input must contain unique nonempty physical spans') END;
CREATE TABLE feature AS SELECT * FROM {parquet(feature_file)}
  WHERE chrom={sql_string(args.chrom)};
CREATE TABLE feature_interval AS
  SELECT regulatory_feature_id, chrom, feature_type, 'feature_extent' AS regulatory_definition_id,
    coalesce(extended_start,start) AS start, coalesce(extended_end,"end") AS "end"
  FROM feature
  UNION ALL SELECT regulatory_feature_id, chrom, feature_type, 'promoter_core',
    core_start, core_end FROM feature WHERE feature_type='promoter' AND core_start IS NOT NULL
  UNION ALL SELECT regulatory_feature_id, chrom, feature_type, 'promoter_extended',
    extended_start, extended_end FROM feature
    WHERE feature_type='promoter' AND extended_start IS NOT NULL;
CREATE TABLE tp73_anchor_regulatory_feature AS
SELECT a.*, f.regulatory_feature_id, f.regulatory_definition_id, f.feature_type,
  least(a.anchor_end,f."end")-greatest(a.anchor_start,f.start) AS overlap_bp,
  a.anchor_start>=f.start AND a.anchor_end<=f."end" AS anchor_fully_within
FROM anchor a JOIN feature_interval f
  ON a.anchor_start<f."end" AND a.anchor_end>f.start;
CREATE TABLE regulatory_feature_tss AS
SELECT NULL::VARCHAR AS regulatory_feature_id, NULL::VARCHAR AS regulatory_definition_id,
  NULL::VARCHAR AS tss_id, NULL::VARCHAR AS genome_id, NULL::VARCHAR AS annotation_release
WHERE false;
CREATE TABLE regulatory_feature_gene AS SELECT * FROM {parquet(native_file)}
WHERE regulatory_feature_id IN (SELECT regulatory_feature_id FROM feature);
"""
    if args.tss:
        inputs.update(tss=args.tss, transcript_tss=args.transcript_tss)
        if not args.annotation_release or not args.genome_id:
            raise ValueError("TSS links require --annotation-release and --genome-id")
        sql += f"""
CREATE TABLE tss AS SELECT * FROM {parquet(args.tss)}
WHERE CAST(chrom AS VARCHAR)={sql_string(args.chrom)};
CREATE TABLE ownership AS SELECT * FROM {parquet(args.transcript_tss)};
SELECT CASE WHEN EXISTS (SELECT 1 FROM tss WHERE
    genome_id IS DISTINCT FROM {sql_string(args.genome_id)} OR
    annotation_release IS DISTINCT FROM {sql_string(args.annotation_release)})
  OR (SELECT count(*)-count(DISTINCT tss_id) FROM tss)<>0
  THEN error('TSS annotation provenance/key mismatch') END;
SELECT CASE WHEN EXISTS (SELECT 1 FROM ownership WHERE
    genome_id IS DISTINCT FROM {sql_string(args.genome_id)} OR
    annotation_release IS DISTINCT FROM {sql_string(args.annotation_release)})
  THEN error('TSS ownership annotation provenance mismatch') END;
INSERT INTO regulatory_feature_tss
SELECT DISTINCT f.regulatory_feature_id, f.regulatory_definition_id, t.tss_id,
  t.genome_id, t.annotation_release
FROM feature_interval f JOIN tss t ON t.start<f."end" AND t."end">f.start;
INSERT INTO regulatory_feature_gene
SELECT DISTINCT r.regulatory_feature_id, o.gene_id,
  'tss_overlap_'||r.regulatory_definition_id, r.annotation_release
FROM regulatory_feature_tss r JOIN ownership o
  ON r.tss_id=o.tss_id AND r.genome_id=o.genome_id
  AND r.annotation_release=o.annotation_release;
"""
    sql += """
CREATE TABLE tp73_anchor_regulatory_membership AS
SELECT a.*, true AS "all",
  coalesce(bool_or(b.regulatory_definition_id='promoter_core'),false) AS promoter_core,
  coalesce(bool_or(b.regulatory_definition_id='promoter_extended'),false) AS promoter_extended,
  NOT coalesce(bool_or(b.regulatory_definition_id='promoter_extended'),false)
    AS outside_promoter_extended,
  coalesce(bool_or(b.feature_type='enhancer'),false) AS enhancer,
  coalesce(bool_or(b.feature_type IN ('open_chromatin','open_chromatin_region')),false)
    AS open_chromatin,
  count(b.regulatory_feature_id)=0 AS no_regulatory_overlap,
  NULL::BOOLEAN AS tss_window_promoter
FROM anchor a LEFT JOIN tp73_anchor_regulatory_feature b
  USING(chrom,anchor_start,anchor_end) GROUP BY a.chrom,a.anchor_start,a.anchor_end;
"""
    if args.promoters:
        if not all((args.promoter_definition_id, args.genome_id, args.annotation_release)):
            raise ValueError("--promoters requires definition ID, genome ID and annotation release")
        inputs["promoters"] = args.promoters
        sql += f"""
CREATE TABLE tss_promoter AS SELECT * FROM {parquet(args.promoters)}
WHERE CAST(chrom AS VARCHAR)={sql_string(args.chrom)}
  AND promoter_definition_id={sql_string(args.promoter_definition_id)};
SELECT CASE WHEN NOT EXISTS (SELECT 1 FROM tss_promoter)
  THEN error('requested promoter definition/chromosome is absent') END;
SELECT CASE WHEN EXISTS (SELECT 1 FROM tss_promoter WHERE
  genome_id IS DISTINCT FROM {sql_string(args.genome_id)} OR
  annotation_release IS DISTINCT FROM {sql_string(args.annotation_release)})
  THEN error('promoter annotation provenance mismatch') END;
CREATE TABLE legacy_promoter_member AS SELECT DISTINCT a.*
FROM anchor a JOIN tss_promoter p
  ON a.anchor_start<p.promoter_end AND a.anchor_end>p.promoter_start;
UPDATE tp73_anchor_regulatory_membership a SET tss_window_promoter=EXISTS (
  SELECT 1 FROM legacy_promoter_member p WHERE p.chrom=a.chrom
    AND a.anchor_start=p.anchor_start AND a.anchor_end=p.anchor_end);
"""
    # Absence of an extended definition is unknown, not evidence of being outside it.
    sql += """
SELECT CASE WHEN EXISTS (SELECT 1 FROM feature
    WHERE feature_type='promoter' AND extended_start IS NULL)
  THEN error('extended promoter bounds missing; cannot classify outside promoters') END;
"""
    for name in ("tp73_anchor_regulatory_feature", "tp73_anchor_regulatory_membership",
                 "regulatory_feature_tss", "regulatory_feature_gene"):
        sql += copy_table(name, staging)
    with tempfile.TemporaryDirectory(prefix="regulatory-join-", dir=staging.parent) as temp:
        run_sql(args, sql, Path(temp))
    publish(staging, args.output, {
        "schema_version": 1, "kind": "tp73_regulatory_membership", "chrom": args.chrom,
        "assembly": args.assembly, "regulatory_release": manifest["regulatory_release"],
        "regulatory_source_sha256": manifest["source_sha256"],
        "coordinate_audit_note": manifest["coordinate_audit_note"],
        "coordinate_audit": manifest.get("coordinate_audit"),
        "annotation_release": args.annotation_release, "genome_id": args.genome_id,
        "promoter_definition_id": args.promoter_definition_id,
        "membership_rule": "positive_half_open_overlap_not_abutment",
        "cofactor_neighborhood": "unchanged; subset anchors only",
        "inputs": {k: {"path": str(v.resolve()), "sha256": sha256(v)} for k,v in inputs.items()},
    })


def select_anchors(args: argparse.Namespace, staging: Path) -> None:
    """Export complete evidence rows for a subset, never a feature-expanded join."""
    manifest = validated_package(args.membership)
    if manifest.get("kind") != "tp73_regulatory_membership":
        raise ValueError("expected a TP73 regulatory membership package")
    membership = args.membership / "tp73_anchor_regulatory_membership.parquet"
    # Subset names are a closed set of identifiers, not user-supplied SQL.
    flag = '"' + args.subset + '"'
    sql = f"""
CREATE TABLE m AS SELECT * FROM {parquet(membership)};
CREATE TABLE a AS SELECT * REPLACE(CAST(chrom AS VARCHAR) AS chrom)
FROM {parquet(args.anchors)} WHERE CAST(chrom AS VARCHAR)={sql_string(manifest['chrom'])};
SELECT CASE WHEN (SELECT count(*)-count(DISTINCT (chrom,anchor_start,anchor_end)) FROM a)<>0
  OR (SELECT count(*)-count(DISTINCT (chrom,anchor_start,anchor_end)) FROM m)<>0
  OR EXISTS (SELECT 1 FROM m WHERE {flag} IS NULL)
  OR EXISTS (SELECT chrom,anchor_start,anchor_end FROM a EXCEPT
             SELECT chrom,anchor_start,anchor_end FROM m)
  OR EXISTS (SELECT chrom,anchor_start,anchor_end FROM m EXCEPT
             SELECT chrom,anchor_start,anchor_end FROM a)
  THEN error('membership must cover the exact physical anchor input without duplicates/NULLs') END;
CREATE TABLE selected_anchors AS SELECT a.* FROM a SEMI JOIN m
  ON a.chrom=m.chrom AND a.anchor_start=m.anchor_start AND a.anchor_end=m.anchor_end
  AND m.{flag};
SELECT CASE WHEN NOT EXISTS (SELECT 1 FROM selected_anchors)
  THEN error('regulatory subset is empty; no estimable analysis') END;
"""
    with tempfile.TemporaryDirectory(prefix="regulatory-select-", dir=staging.parent) as temp:
        run_sql(args, sql + copy_table("selected_anchors", staging), Path(temp))
    publish(staging, args.output, {
        "schema_version": 1, "kind": "regulatory_anchor_subset", "subset": args.subset,
        "assembly": manifest["assembly"], "chrom": manifest["chrom"],
        "membership_manifest_sha256": sha256(args.membership / "manifest.json"),
        "anchor_input_sha256": sha256(args.anchors),
        "coordinate_audit_note": manifest["coordinate_audit_note"],
        "selection": "anchor_only; cofactor neighborhood unchanged",
    })


def audit_coordinates(args: argparse.Namespace, staging: Path) -> None:
    """Compare every exported boundary with UCSC's independent BigBed reader."""
    digest = sha256(args.gff)
    expected = {}
    for row in gff_rows(args.gff, args.assembly, "audit", digest):
        identifier = row["source_feature_id"]
        if identifier in expected:
            raise ValueError(f"duplicate GFF ID: {identifier}")
        expected[identifier] = (
            row["chrom"], row["extended_start"] if row["extended_start"] is not None else row["start"],
            row["extended_end"] if row["extended_end"] is not None else row["end"],
            row["strand"], row["start"], row["end"], row["feature_type"],
        )
    if not expected:
        raise ValueError("empty coordinate audit input")
    count = len(expected)
    checked = 0
    # Stream the reader output; there is no intermediate genome-wide BED layer.
    with tempfile.TemporaryFile(mode="w+") as errors:
        with subprocess.Popen(
            [str(args.bigbed_to_bed.resolve()), str(args.bigbed.resolve()), "stdout"],
            stdout=subprocess.PIPE, stderr=errors, text=True,
        ) as process:
            for line in process.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    raise ValueError("BigBed export lacks BED9+ feature type")
                observed = (fields[0], int(fields[1]), int(fields[2]), fields[5],
                            int(fields[6]), int(fields[7]), fields[9].replace(" ", "_"))
                wanted = expected.pop(fields[3], None)
                if wanted != observed:
                    raise ValueError(f"coordinate/identity mismatch {fields[3]}: {wanted} != {observed}")
                checked += 1
        if process.returncode:
            errors.seek(0)
            raise ValueError("BigBed reader failed: " + errors.read(4096))
    if expected or checked != count:
        raise ValueError("BigBed/GFF feature inventories differ")
    publish(staging, args.output, {
        "schema_version": 1, "kind": "regulatory_coordinate_audit", "status": "passed",
        "assembly": args.assembly, "features_compared": count,
        "gff_sha256": digest, "bigbed_sha256": sha256(args.bigbed),
        "bigbed_reader_sha256": sha256(args.bigbed_to_bed),
        "comparison": "all_IDs_chrom_strand_type_core_and_extended_BED_bounds",
        "gff_conversion": "start_minus_one_end_unchanged_including_extended_attributes",
    })


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    imp = commands.add_parser("import", help="import an official GFF, retaining source provenance")
    imp.add_argument("--gff", type=Path, required=True)
    imp.add_argument("--requested-release", required=True)
    imp.add_argument("--resolved-release", required=True,
                     help="actual payload release, not necessarily the requested URL release")
    imp.add_argument("--source-uri", required=True, help="durable public URL, not a signed redirect")
    imp.add_argument("--expected-sha256")
    imp.add_argument("--coordinate-audit", type=Path,
                     help="passed audit JSON for this exact GFF payload and assembly")
    imp.add_argument("--coordinate-audit-note", required=True,
                     help="record independent coordinate validation or unresolved discrepancies")
    imp.set_defaults(function=import_gff)
    ann = commands.add_parser("annotate", help="join one chromosome of physical anchors to features")
    ann.add_argument("--features", type=Path, required=True)
    ann.add_argument("--anchors", type=Path, required=True,
                     help="Parquet with unique chrom, anchor_start, anchor_end")
    ann.add_argument("--chrom", required=True)
    ann.add_argument("--tss", type=Path)
    ann.add_argument("--transcript-tss", type=Path)
    ann.add_argument("--annotation-release")
    ann.add_argument("--genome-id")
    ann.add_argument("--promoters", type=Path)
    ann.add_argument("--promoter-definition-id")
    ann.set_defaults(function=annotate)
    sel = commands.add_parser("select", help="export subset evidence for existing TP73 evaluators")
    sel.add_argument("--anchors", type=Path, required=True)
    sel.add_argument("--membership", type=Path, required=True)
    sel.add_argument("--subset", choices=SUBSETS, required=True)
    sel.set_defaults(function=select_anchors)
    audit = commands.add_parser("audit", help="compare all GFF coordinates with UCSC BigBed output")
    audit.add_argument("--gff", type=Path, required=True)
    audit.add_argument("--bigbed", type=Path, required=True)
    audit.add_argument("--bigbed-to-bed", type=Path, required=True)
    audit.set_defaults(function=audit_coordinates)
    for subparser in (imp, ann, audit):
        subparser.add_argument("--assembly", required=True)
    for subparser in (imp, ann, sel, audit):
        subparser.add_argument("--output", type=Path, required=True,
                               help="new immutable package directory; never overwrites")
        subparser.add_argument("--duckdb", default="duckdb")
        subparser.add_argument("--threads", type=int, default=2)
        subparser.add_argument("--memory-limit", default="1GB")
    return result


def main() -> int:
    args = parser().parse_args()
    try:
        if args.threads < 1:
            raise ValueError("threads must be positive")
        args.output = args.output.resolve()
        if args.output.exists():
            raise ValueError(f"refusing to replace existing package: {args.output}")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=".regulatory-", dir=args.output.parent) as temp:
            staging = Path(temp) / "package"
            staging.mkdir()
            args.function(args, staging)
        print(f"I: Created {args.output}", file=sys.stderr)
        return 0
    except (OSError, ValueError, KeyError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

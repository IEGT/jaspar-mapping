#!/usr/bin/env python3
"""Tag/export stored motif hits using regulatory and TSS intervals, not TP73.

The small completed regulatory-run plan identifies the chromosome-wide GTF
promoter dimension. Its TP73 membership payload is deliberately never read.
Source hit files remain unchanged; overlapping features cannot multiply hits.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

from build_regulatory_annotation import sha256
from query_genome_scan import QueryError, package_paths, sql_list, sql_string


TAGS = {
    "promoter_core": 1, "promoter_extended": 2, "enhancer": 4,
    "open_chromatin": 8, "ctcf": 16, "emar": 32,
    "other_regulatory": 64, "tss_window": 128,
}
SCOPES = {"regulatory_or_tss": 255, "regulatory": 127,
          "promoter_or_tss": 131, **TAGS}


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def record(path):
    path = Path(path).resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def verify(item):
    path = Path(item["path"]).resolve()
    if path.stat().st_size != item["bytes"] or sha256(path) != item["sha256"]:
        raise QueryError(f"annotation input changed: {path}")
    return path


def sql(args, query, output_watch=None):
    settings = (f"SET threads={args.threads}; SET memory_limit={sql_string(args.memory_limit)};"
                "SET max_temp_directory_size='0B';"
                "SET autoinstall_known_extensions=false; SET autoload_known_extensions=false;")
    process = subprocess.Popen([args.duckdb, "-no-init", "-batch", "-bail", "-json", ":memory:"],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True)
    deadline = time.monotonic() + args.timeout_seconds
    payload = settings + query
    try:
        while True:
            try:
                stdout, stderr = process.communicate(input=payload, timeout=0.25)
                break
            except subprocess.TimeoutExpired:
                payload = None
                if time.monotonic() >= deadline:
                    raise QueryError("DuckDB wall-time limit exceeded")
                if output_watch is not None:
                    if output_watch.exists() and output_watch.stat().st_size > args.max_output_bytes:
                        raise QueryError(f"output byte budget exceeded; stopped writer at {output_watch}")
                    if shutil.disk_usage(output_watch.parent).free < args.minimum_free_bytes:
                        raise QueryError("free-space reserve reached; stopped writer before publication")
        if process.returncode:
            raise QueryError(stderr.strip() or "DuckDB failed")
        return json.loads(stdout or "[]")
    finally:
        if process.poll() is None:
            process.kill()
            process.communicate()


def inputs(args):
    package, database = package_paths(args)
    manifest_path = package / "manifest.json"
    scan = read_json(manifest_path)
    if scan.get("state") != "complete" or scan.get("schema_version") not in (1, 2):
        raise QueryError("a completed genome-scan package is required")
    attach = f"ATTACH {sql_string(database)} AS scan (READ_ONLY);"
    genomes = sql(args, attach + "SELECT genome_id,assembly_name FROM scan.genome;")
    if len(genomes) != 1 or genomes[0]["genome_id"] != scan.get("genome_id"):
        raise QueryError("scan genome identity is not unique/consistent")

    root = args.regulatory_run.resolve()
    plan_path, final_path = root / "plan.json", root / "final/manifest.json"
    plan, final = read_json(plan_path), read_json(final_path)
    if (final.get("kind") != "genome_regulatory_membership"
            or final.get("state") != "complete"
            or final.get("production_plan_sha256") != sha256(plan_path)
            or plan.get("genome_id") != scan["genome_id"]
            or final.get("assembly") != genomes[0]["assembly_name"]
            or final.get("annotation_release") != plan.get("annotation_release")
            or final.get("promoter_definition_id") != plan.get("promoter_definition_id")):
        raise QueryError("regulatory plan/completion/genome identity mismatch")
    tasks = [r for r in plan["tasks"] if r["chrom"] == args.chrom]
    if len(tasks) != 1 or args.chrom not in final["chromosomes"]:
        raise QueryError("chromosome lacks complete regulatory AND TSS coverage")
    promoters = verify(tasks[0]["promoters"])
    features_manifest_path = root / "features/manifest.json"
    features = read_json(features_manifest_path)
    audit = features.get("coordinate_audit") or {}
    if (features.get("kind") != "ensembl_regulatory_annotation"
            or features.get("coordinate_mode") != "bed_0based_half_open"
            or features.get("assembly") != final["assembly"]
            or args.chrom not in features.get("chromosomes", [])
            or audit.get("kind") != "regulatory_coordinate_audit"
            or audit.get("status") != "passed"
            or audit.get("assembly") != final["assembly"]
            or audit.get("gff_sha256") != features.get("source_sha256")
            or audit != plan.get("audit") or audit != final.get("coordinate_audit")):
        raise QueryError("complete, coordinate-audited regulatory features are required")
    feature_entries = [f for f in features["files"] if f["path"] == "regulatory_feature.parquet"]
    if len(feature_entries) != 1:
        raise QueryError("regulatory feature inventory is not unique")
    feature_record = dict(feature_entries[0], path=str(root / "features/regulatory_feature.parquet"))
    feature_path = verify(feature_record)

    motif_filter = "" if args.all_motifs else f"AND i.motif_id IN ({','.join(sql_string(m) for m in args.motif)})"
    inventory = sql(args, attach + f"""
SELECT i.*, m.motif_name FROM scan.scan_file_inventory i
JOIN scan.motif_metadata m USING (motif_set_id,motif_id)
WHERE CAST(i.chrom AS VARCHAR)={sql_string(args.chrom)} {motif_filter}
ORDER BY i.motif_id,i.strand;
""")
    if not inventory:
        raise QueryError("no matching motif/chromosome inventory")
    represented = {r["motif_id"] for r in inventory}
    if not args.all_motifs and represented != set(args.motif):
        raise QueryError("one or more requested motifs are absent from the scan")
    if args.all_motifs:
        catalog = sql(args, attach + "SELECT motif_id FROM scan.motif_metadata;")
        if represented != {r["motif_id"] for r in catalog}:
            raise QueryError("chromosome inventory does not cover all catalog motifs")
    seen = set()
    paths = []
    for row in inventory:
        key = (row["motif_id"], row["strand"])
        if key in seen or row["strand"] not in ("+", "-") or row["state"] != "complete":
            raise QueryError("duplicate/incomplete motif-orientation inventory")
        seen.add(key)
        if (row["genome_id"] != scan["genome_id"]
                or row["motif_set_id"] != scan["motif_set_id"]
                or row["coordinate_mode"] != "bed"
                or row["minimum_pwm_relative_score"] is not None
                or row["maximum_pwm_relative_score"] is not None
                or not math.isfinite(float(row["minimum_score"]))
                or float(row["minimum_score"]) > args.minimum_score):
            raise QueryError("scan identity/filter floor cannot cover the requested score threshold")
        path = (package / "task_data" / f"task_id={row['task_id']}"
                / row["output_relative_path"]).resolve()
        if not path.is_relative_to(package) or not path.is_file() or path.stat().st_size != row["bytes"]:
            raise QueryError(f"missing, escaped or size-changed scan payload: {path}")
        paths.append(str(path))
    if len(seen) != 2 * len(represented):
        raise QueryError("both orientations are required for every requested motif")
    provenance = {
        "scan_manifest": record(manifest_path), "scan_run_id": scan["run_id"],
        "genome_id": scan["genome_id"], "assembly": final["assembly"],
        "motif_set_id": scan["motif_set_id"],
        "regulatory_plan": record(plan_path), "regulatory_completion": record(final_path),
        "regulatory_features_manifest": record(features_manifest_path),
        "regulatory_features": feature_record, "tss_promoters": tasks[0]["promoters"],
        "regulatory_release": features["regulatory_release"],
        "annotation_release": plan["annotation_release"],
        "promoter_definition_id": plan["promoter_definition_id"],
        "scan_payload_validation": "exact_inventory_sizes; recorded SHA256, not payload rehash",
        "coordinate_audit": audit,
    }
    return database, feature_path, promoters, paths, inventory, provenance


def query_sql(args, database, features, promoters, paths, provenance):
    chrom = sql_string(args.chrom)
    region = (f'AND h.start < {args.end} AND h."end" > {args.start}'
              if args.start is not None else "")
    scope = SCOPES[args.scope]
    # Union each class before joining, so redundant transcript/promoter owners
    # do not inflate work. A hit spanning several classes receives their bit-OR.
    return f"""
ATTACH {sql_string(database)} AS scan (READ_ONLY);
CREATE TEMP TABLE feature AS SELECT * FROM read_parquet({sql_string(features)}, hive_partitioning=false)
WHERE CAST(chrom AS VARCHAR)={chrom};
CREATE TEMP TABLE promoter AS SELECT * FROM read_parquet({sql_string(promoters)}, hive_partitioning=false);
SELECT CASE WHEN NOT EXISTS(SELECT 1 FROM feature) OR NOT EXISTS(SELECT 1 FROM promoter)
 OR EXISTS(SELECT 1 FROM feature WHERE assembly IS DISTINCT FROM {sql_string(provenance['assembly'])}
   OR start IS NULL OR "end" IS NULL OR feature_type IS NULL OR start<0 OR "end"<=start OR (feature_type='promoter' AND
     (core_start IS NULL OR core_end IS NULL OR extended_start IS NULL OR extended_end IS NULL
      OR extended_start<0 OR extended_start>core_start OR core_start>=core_end OR core_end>extended_end)))
 OR EXISTS(SELECT 1 FROM promoter WHERE CAST(chrom AS VARCHAR) IS DISTINCT FROM {chrom}
   OR genome_id IS DISTINCT FROM {sql_string(provenance['genome_id'])}
   OR annotation_release IS DISTINCT FROM {sql_string(provenance['annotation_release'])}
   OR promoter_definition_id IS DISTINCT FROM {sql_string(provenance['promoter_definition_id'])}
   OR promoter_start IS NULL OR promoter_end IS NULL OR promoter_start<0 OR promoter_end<=promoter_start)
 THEN error('Invalid or incompatible regulatory/TSS intervals') END;
CREATE TEMP TABLE raw_region AS
SELECT start,"end",CASE feature_type WHEN 'enhancer' THEN 4
 WHEN 'open_chromatin_region' THEN 8 WHEN 'CTCF_binding_site' THEN 16
 WHEN 'EMAR' THEN 32 ELSE 64 END AS tag_mask FROM feature WHERE feature_type<>'promoter'
UNION ALL SELECT core_start,core_end,1 FROM feature WHERE feature_type='promoter'
UNION ALL SELECT extended_start,extended_end,2 FROM feature WHERE feature_type='promoter'
UNION ALL SELECT promoter_start,promoter_end,128 FROM promoter;
CREATE TEMP TABLE regions AS
WITH preceding AS (SELECT *,max("end") OVER(PARTITION BY tag_mask ORDER BY start,"end"
 ROWS BETWEEN UNBOUNDED PRECEDING AND 1 PRECEDING) AS previous_end FROM raw_region),
numbered AS (SELECT *,sum(CASE WHEN previous_end IS NULL OR start>previous_end THEN 1 ELSE 0 END)
 OVER(PARTITION BY tag_mask ORDER BY start,"end" ROWS UNBOUNDED PRECEDING) AS island FROM preceding)
SELECT min(start)::BIGINT AS start,max("end")::BIGINT AS "end",tag_mask
FROM numbered GROUP BY tag_mask,island;
CREATE TEMP VIEW tagged_hits AS
WITH hits AS (
 SELECT h.genome_id,h.motif_set_id,CAST(h.chrom AS VARCHAR) AS chrom,h.start,h."end",h.motif_id,
 m.motif_name,CASE h.strand WHEN 'plus' THEN '+' WHEN 'minus' THEN '-' ELSE h.strand END AS strand,
 h.score,h.pwm_relative_score,h.score_mode,CAST(h.pseudocount AS DOUBLE) AS pseudocount,
 h.background_model_id,h.pseudocount_scheme,CAST(h.minimum_score AS DOUBLE) AS source_minimum_score,
 h.n_policy FROM read_parquet({sql_list(paths)}, hive_partitioning=true) h
 JOIN scan.motif_metadata m USING(motif_set_id,motif_id)
 WHERE CAST(h.chrom AS VARCHAR)={chrom} AND h.score>={args.minimum_score:.17g} {region}
), tagged AS (
 SELECT h.*,coalesce((SELECT bit_or(r.tag_mask) FROM regions r
   WHERE h.start<r."end" AND h."end">r.start),0)::USMALLINT AS regulation_tags
 FROM hits h
)
SELECT *,{','.join(f'(regulation_tags & {bit})<>0 AS overlaps_{name}' for name, bit in TAGS.items())},
 (regulation_tags & 127)<>0 AS overlaps_regulatory,
 regulation_tags<>0 AS regulation_candidate
FROM tagged WHERE (regulation_tags & {scope})<>0;
"""


@contextmanager
def output_lock(output):
    # The lock file remains; the kernel releases ownership on exit/preemption.
    with output.with_name(output.name + ".lock").open("a") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise QueryError("another exporter owns this output") from error
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def stage_inputs(args, data):
    database, features, promoters, paths, inventory, provenance = data
    if args.scratch_directory is None:
        return data
    scratch = args.scratch_directory.expanduser().resolve()
    scratch.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(scratch).free < sum(r["bytes"] for r in inventory) + 1073741824:
        raise QueryError("insufficient scratch space for selected files plus 1 GiB reserve")
    attempt = Path(tempfile.mkdtemp(prefix="regulatory-tfbs-", dir=scratch))

    def copy(source, relative, expected):
        target = attempt / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        verify(dict(expected, path=str(target)))
        return target

    staged = []
    package, _ = package_paths(args)
    for source, entry in zip(paths, inventory):
        staged.append(str(copy(source, Path("scan") / Path(source).relative_to(package), entry)))
    features = copy(features, "regulatory_feature.parquet", provenance["regulatory_features"])
    promoters = copy(promoters, "promoter.parquet", provenance["tss_promoters"])
    print(f"I: Verified selected hit and annotation files on scratch: {attempt}", file=sys.stderr)
    return database, features, promoters, staged, inventory, provenance


def reuse(args, data):
    manifest = read_json(args.output / "manifest.json")
    expected = {
        "kind": "regulatory_tfbs_counts" if args.count_only else "regulatory_tfbs_subset",
        "state": "complete", "scope": args.scope, "minimum_score": args.minimum_score,
        "chromosomes": [args.chrom], "region": {"start": args.start, "end": args.end},
        "motif_ids": sorted({r["motif_id"] for r in data[4]}),
        "inputs": data[5], "source_inventory": data[4],
    }
    if (any(manifest.get(k) != v for k, v in expected.items())
            or manifest.get("builder", {}).get("sha256") != sha256(Path(__file__))
            or manifest.get("rows", args.max_rows + 1) > args.max_rows
            or manifest.get("parquet_bytes", args.max_output_bytes + 1) > args.max_output_bytes):
        raise QueryError("existing export does not match requested scope, inputs, builder or caps")
    for entry in manifest["files"]:
        path = (args.output / entry["path"]).resolve()
        if not path.is_relative_to(args.output):
            raise QueryError("unsafe export inventory")
        verify(dict(entry, path=str(path)))
    print(json.dumps({"output": str(args.output), "reused": True,
                      "rows": manifest["rows"], "parquet_bytes": manifest["parquet_bytes"]}))


def execute(args):
    if args.output.exists() and not args.resume:
        raise QueryError("refusing to replace an existing output directory")
    data = inputs(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with output_lock(args.output):
        if args.output.exists():
            if not args.resume:
                raise QueryError("refusing to replace an existing output directory")
            return reuse(args, data)
        return export_new(args, stage_inputs(args, data))


def export_new(args, data):
    database, features, promoters, paths, inventory, provenance = data
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(args.output.parent).free < args.minimum_free_bytes + args.max_output_bytes:
        raise QueryError("insufficient free space for the output budget plus reserve")
    staging = Path(tempfile.mkdtemp(prefix=f".{args.output.name}.attempt-", dir=args.output.parent))
    print(f"I: {len(inventory)} exact input files; output attempt {staging}", file=sys.stderr)
    base = query_sql(args, database, features, promoters, paths, provenance)
    destination = staging / ("counts.parquet" if args.count_only else "motif_hits.parquet")
    selection = ("SELECT motif_id,motif_name,regulation_tags,count(*) AS orientation_records,"
                 'count(DISTINCT (start,"end")) AS physical_loci '
                 "FROM tagged_hits GROUP BY motif_id,motif_name,regulation_tags"
                 if args.count_only else
                 'SELECT * FROM tagged_hits ORDER BY start,"end",motif_id,strand')
    query = base + f"\nCOPY ({selection} LIMIT {args.max_rows + 1}) TO {sql_string(destination)} " \
                   "(FORMAT PARQUET, COMPRESSION ZSTD);"
    # Save the exact query for audit/replay; consumers need not execute it.
    (staging / "query.sql").write_text(query, encoding="utf-8")
    sql(args, query, output_watch=destination)
    stats = sql(args, f"SELECT count(*) AS rows FROM read_parquet({sql_string(destination)});")[0]
    if stats["rows"] > args.max_rows:
        raise QueryError(f"row limit exceeded ({args.max_rows}); narrow selection or explicitly raise limit; {staging}")
    file_info = record(destination)
    if file_info["bytes"] > args.max_output_bytes:
        raise QueryError(f"byte limit exceeded; unpublished attempt preserved at {staging}")
    # Small annotations can be rehashed cheaply; never reread TB-scale payloads.
    for key in ("scan_manifest", "regulatory_plan", "regulatory_completion",
                "regulatory_features_manifest", "regulatory_features", "tss_promoters"):
        verify(provenance[key])
    (staging / "schema.sql").write_text(
        f"-- Open from this package directory; coordinates are BED half-open.\n"
        f"CREATE OR REPLACE VIEW {'regulatory_tfbs_counts' if args.count_only else 'regulatory_tfbs'} AS "
        f"SELECT * FROM read_parquet('{destination.name}', hive_partitioning=false);\n",
        encoding="utf-8")
    manifest = {
        "schema_version": 1, "kind": "regulatory_tfbs_counts" if args.count_only else "regulatory_tfbs_subset",
        "state": "complete", "complete_genome_scan": False,
        "complete_for_requested_scope_at_threshold": True,
        "requires_tp73": False, "coordinate_mode": "bed_0based_half_open",
        "scope": args.scope, "scope_mask": SCOPES[args.scope], "tag_bits": TAGS,
        "chromosomes": [args.chrom], "motif_ids": sorted({r["motif_id"] for r in inventory}),
        "region": {"start": args.start, "end": args.end}, "minimum_score": args.minimum_score,
        "rows": stats["rows"], "parquet_bytes": file_info["bytes"],
        "max_rows": args.max_rows, "max_output_bytes": args.max_output_bytes,
        "scratch_staged_and_verified": args.scratch_directory is not None,
        "inputs": provenance, "source_inventory": inventory,
        "builder": record(Path(__file__)),
        "files": [dict(record(staging / name), path=name)
                  for name in (destination.name, "query.sql", "schema.sql")],
        "non_claims": ["Annotation overlap is not measured regulatory activity or protein occupancy.",
                       "No annotation outside the declared chromosome/window/motif/score scope is inferred.",
                       "Both orientations are separate rows; no strongest-hit or TP73 selection is applied."],
    }
    (staging / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.output.exists():
        raise QueryError("output appeared during export; refusing to replace it")
    os.rename(staging, args.output)
    print(json.dumps({"output": str(args.output), "rows": stats["rows"],
                      "parquet_bytes": file_info["bytes"], "count_only": args.count_only}))


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--package", type=Path, required=True, help="completed permissive genome-scan package")
    p.add_argument("--database", type=Path, help="optional rebuilt scan catalog")
    p.add_argument("--regulatory-run", type=Path, required=True, help="completed regulatory run with plan, features and TSS dimensions")
    selection = p.add_mutually_exclusive_group(required=True)
    selection.add_argument("--motif", action="append", help="exact JASPAR accession; repeat for a panel")
    selection.add_argument("--all-motifs", action="store_true", help="explicitly process every catalog motif")
    p.add_argument("--chrom", required=True, help="exact chromosome; both annotation sources must cover it")
    p.add_argument("--start", type=int, help="BED inclusive start of query interval")
    p.add_argument("--end", type=int, help="BED exclusive end of query interval")
    p.add_argument("--whole-chromosome", action="store_true", help="explicitly allow an unbounded chromosome interval")
    p.add_argument("--scope", choices=SCOPES, default="regulatory_or_tss")
    p.add_argument("--minimum-score", type=float, default=-1, help="inclusive score; source retention must cover it (default -1)")
    p.add_argument("--count-only", action="store_true", help="write compact per-motif/tag counts for sizing; still reads selected hit files")
    p.add_argument("--output", type=Path, required=True, help="new immutable output directory")
    p.add_argument("--resume", action="store_true", help="verify and reuse an identical completed export; failed attempts remain untouched")
    p.add_argument("--scratch-directory", type=Path, help="stage and checksum exact input files in a unique job-local directory")
    p.add_argument("--max-rows", type=int, default=1000000, help="publication ceiling, never silent truncation")
    p.add_argument("--max-output-bytes", type=int, default=1000000000, help="Parquet publication ceiling (default 1 GB)")
    p.add_argument("--minimum-free-bytes", type=int, default=1073741824, help="free-space reserve in addition to output ceiling")
    p.add_argument("--duckdb", default="duckdb")
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--memory-limit", default="1GB", help="DuckDB memory budget; disk spill is disabled")
    p.add_argument("--timeout-seconds", type=int, default=1800, help="per-DuckDB-call wall time limit")
    return p


def main():
    args = parser().parse_args()
    try:
        region = args.start is not None and args.end is not None and 0 <= args.start < args.end <= 2**63-1
        if args.whole_chromosome:
            if args.start is not None or args.end is not None:
                raise QueryError("whole-chromosome and interval selection are exclusive")
        elif not region:
            raise QueryError("supply --start/--end or explicitly --whole-chromosome")
        if (not math.isfinite(args.minimum_score) or args.max_rows < 1 or args.max_output_bytes < 1
                or args.minimum_free_bytes < 0 or not 1 <= args.threads <= 32
                or args.timeout_seconds < 1 or not re.fullmatch(r"[0-9]+(?:MB|GB)", args.memory_limit)):
            raise QueryError("invalid score/resource/output limit")
        if args.motif and len(args.motif) != len(set(args.motif)):
            raise QueryError("duplicate requested motif")
        if shutil.which(args.duckdb) is None:
            raise QueryError("DuckDB executable is unavailable")
        args.output = args.output.expanduser().resolve()
        execute(args)
        return 0
    except (OSError, ValueError, KeyError, QueryError, subprocess.TimeoutExpired) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

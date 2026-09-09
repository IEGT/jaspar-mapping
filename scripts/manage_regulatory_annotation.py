#!/usr/bin/env python3
"""Prepare, stage and finalize restart-safe chromosome regulatory annotation.

Run setup once, then one run-task per planned chromosome, then finalize with
an afterok dependency. Source, evidence, GTF dimensions and coordinate audit
are pinned in an immutable plan. Heavy work runs in the allocation, not on
the login node. Only completed chromosome packages are published on /data.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import urllib.request

import build_regulatory_annotation as annotation


STARTED = time.monotonic()
PHASE = "startup"
SOURCES = ("scripts/build_regulatory_annotation.py", "scripts/manage_regulatory_annotation.py")
DIMENSIONS = {"tss": "transcription_start_site.parquet",
              "transcript_tss": "transcript_tss.parquet", "promoters": "promoter.parquet"}


def status_signal(signum, frame):
    print(f"I: phase={PHASE} elapsed_seconds={time.monotonic()-STARTED:.1f}",
          file=sys.stderr, flush=True)


def sql(args, text):
    result = subprocess.run([str(args.duckdb), "-light-mode", "-batch", "-bail", "-json", ":memory:"],
                            input="SET threads=2; SET memory_limit='2GB';" + text,
                            text=True, capture_output=True)
    if result.returncode:
        raise ValueError(result.stderr)
    return json.loads(result.stdout) if result.stdout.strip() else []


def record(path):
    path = Path(path).resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": annotation.sha256(path)}


def verify_file(item):
    path = Path(item["path"])
    if path.stat().st_size != item["bytes"] or annotation.sha256(path) != item["sha256"]:
        raise ValueError(f"pinned input changed: {path}")
    return path


def write_new_json(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


def prepare(args):
    if args.run_root.exists():
        raise ValueError("run root already exists; reuse the existing plan, never replace it")
    source = args.source.resolve()
    for command in (["diff", "--quiet"], ["diff", "--cached", "--quiet"]):
        subprocess.run(["git", "-C", str(source), *command], check=True)
    commit = subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
    audit = json.loads(args.audit.read_text())
    if (audit.get("status") != "passed" or audit.get("kind") != "regulatory_coordinate_audit"
            or audit.get("assembly") != "GRCh38"):
        raise ValueError("a passed GRCh38 coordinate audit is required")
    chromosomes = args.chromosomes.split(",")
    if len(chromosomes) != len(set(chromosomes)) or not set(chromosomes) <= set(map(str, range(1,23))):
        raise ValueError("chromosomes must be distinct autosome names")
    evidence_manifest = args.evidence_package / "manifest.json"
    catalog_manifest = args.annotation_catalog / "manifest.json"
    if json.loads(evidence_manifest.read_text()).get("state") != "complete":
        raise ValueError("evidence package is incomplete")
    catalog = json.loads(catalog_manifest.read_text())
    if catalog.get("state") != "complete" or catalog.get("context_schema_version") != 9:
        raise ValueError("completed schema-9 annotation is required")
    database = args.annotation_catalog / catalog.get("database", "context.duckdb")
    query = "SELECT dataset,CAST(chrom AS VARCHAR) chrom,absolute_path FROM context_file_inventory "
    query += "WHERE dataset IN ('transcription_start_site.parquet','transcript_tss.parquet','promoter.parquet');"
    process = subprocess.run([str(args.duckdb), "-light-mode", "-readonly", "-json", str(database), "-c", query],
                             text=True, capture_output=True, check=True)
    dimensions = json.loads(process.stdout)
    evidence_inventory = args.evidence_package / "chromosome_file_inventory.tsv"
    with evidence_inventory.open() as stream:
        evidence = list(csv.DictReader(stream, delimiter="\t"))
    tasks = []
    for chrom in chromosomes:
        candidates = [r for r in evidence if r["chrom"] == chrom]
        if len(candidates) != 1:
            raise ValueError(f"evidence inventory is not unique for {chrom}")
        r = candidates[0]
        relative = Path(r["relative_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError("unsafe evidence inventory path")
        entry = record(args.evidence_package / relative)
        if entry["sha256"] != r["sha256"] or entry["bytes"] != int(r["bytes"]):
            raise ValueError("evidence differs from finalized inventory")
        task = {"chrom": chrom, "evidence": entry}
        for key, dataset in DIMENSIONS.items():
            candidates = [r for r in dimensions if r["chrom"] == chrom and r["dataset"] == dataset]
            if len(candidates) != 1:
                raise ValueError(f"expected one {dataset} for chromosome {chrom}")
            task[key] = record(candidates[0]["absolute_path"])
        tasks.append(task)
    # The plan becomes visible atomically, with no half-written restart state.
    args.run_root.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".reg-plan-", dir=args.run_root.parent) as temp:
        staging = Path(temp) / "run"
        staging.mkdir()
        for name in ("logs", "source", "tasks"):
            (staging / name).mkdir()
        write_new_json(staging / "plan.json", {
            "schema_version": 1, "source": str(source), "source_commit": commit,
            "source_files": [record(source / name) for name in SOURCES],
            "audit": audit, "audit_file": record(args.audit), "tasks": tasks,
            "chromosomes": chromosomes, "duckdb": str(args.duckdb.resolve()),
            "genome_id": args.genome_id, "annotation_release": args.annotation_release,
            "promoter_definition_id": args.promoter_definition_id,
            "inputs": [record(evidence_manifest), record(evidence_inventory), record(catalog_manifest)],
        })
        os.rename(staging, args.run_root)
    print(f"I: Prepared {len(tasks)} chromosome tasks at {args.run_root}")


def load_plan(args):
    plan = json.loads((args.run_root / "plan.json").read_text())
    for item in plan["source_files"] + [plan["audit_file"]]:
        verify_file(item)
    args.duckdb = Path(plan["duckdb"])
    return plan


def invoke(plan, *arguments):
    subprocess.run([sys.executable, str(Path(plan["source"]) / SOURCES[0]),
                    *map(str, arguments), "--duckdb", plan["duckdb"],
                    "--memory-limit", "2GB"], check=True)


def setup(args):
    global PHASE
    plan = load_plan(args)
    destination = args.run_root / "features"
    if destination.exists():
        existing = annotation.validated_package(destination)
        if existing.get("coordinate_audit") != plan["audit"]:
            raise ValueError("existing regulatory features have another audit")
        print("I: Reusing validated regulatory features")
        return
    PHASE = "source-download"
    gff = args.run_root / "source" / "regulatory_features.gff.gz"
    if not gff.exists():
        if shutil.disk_usage(args.run_root).free < 2 * 1024**3:
            raise ValueError("less than 2 GiB available for source import")
        with tempfile.TemporaryDirectory(prefix=".reg-source-", dir=gff.parent) as temp:
            candidate = Path(temp) / "source.gff.gz"
            if args.gff:
                shutil.copyfile(args.gff, candidate)
            else:
                with urllib.request.urlopen(plan["audit"]["source_uri"], timeout=90) as response, candidate.open("wb") as out:
                    total = 0
                    while chunk := response.read(1024 * 1024):
                        total += len(chunk)
                        if total > 64 * 1024**2:
                            raise ValueError("regulatory download exceeds 64 MiB limit")
                        out.write(chunk)
            if annotation.sha256(candidate) != plan["audit"]["gff_sha256"]:
                raise ValueError("download differs from coordinate-audited GFF")
            os.rename(candidate, gff)
    if annotation.sha256(gff) != plan["audit"]["gff_sha256"]:
        raise ValueError("cached GFF differs from coordinate audit")
    PHASE = "source-import"
    invoke(plan, "import", "--gff", gff, "--assembly", "GRCh38",
           "--requested-release", plan["audit"]["requested_regulatory_release"],
           "--resolved-release", plan["audit"]["regulatory_release"],
           "--source-uri", plan["audit"]["source_uri"],
           "--expected-sha256", plan["audit"]["gff_sha256"],
           "--coordinate-audit", plan["audit_file"]["path"],
           "--coordinate-audit-note", f"all {plan['audit']['features_compared']} GFF/BigBed features agree after standard GFF conversion",
           "--output", destination)


def verify_task(path, plan_digest):
    manifest = annotation.validated_package(path)
    if manifest.get("production_plan_sha256") != plan_digest:
        raise ValueError("completed chromosome belongs to a different plan")
    return manifest


def run_task(args):
    global PHASE
    plan = load_plan(args)
    index = args.task_index if args.task_index is not None else int(os.environ["SLURM_ARRAY_TASK_ID"])
    if not 0 <= index < len(plan["tasks"]):
        raise ValueError("task index outside plan")
    task = plan["tasks"][index]
    digest = annotation.sha256(args.run_root / "plan.json")
    target = args.run_root / "tasks" / f"chrom-{task['chrom']}"
    if target.exists():
        verify_task(target, digest)
        print(f"I: Reusing chromosome {task['chrom']}")
        return
    PHASE = f"chrom-{task['chrom']}-staging"
    args.scratch_root.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(args.scratch_root).free < 2 * 1024**3:
        raise ValueError("less than 2 GiB free on scratch")
    # Unique attempt directories are never reused after preemption. Haumea
    # owns scratch cleanup; completed durable packages are never removed.
    scratch = Path(tempfile.mkdtemp(prefix=f"reg-{os.environ.get('SLURM_JOB_ID','local')}-", dir=args.scratch_root))
    features = args.run_root / "features"
    feature_manifest = annotation.validated_package(features)
    if feature_manifest.get("coordinate_audit") != plan["audit"]:
        raise ValueError("feature audit differs from immutable plan")
    shutil.copytree(features, scratch / "features")
    paths = {}
    for key in ("evidence", *DIMENSIONS):
        item = task[key]
        path = scratch / f"{key}.parquet"
        shutil.copyfile(item["path"], path)
        verify_file({**item, "path": str(path)})
        paths[key] = path
    PHASE = f"chrom-{task['chrom']}-annotation"
    invoke(plan, "annotate", "--features", scratch / "features",
           "--anchors", paths["evidence"], "--chrom", task["chrom"], "--assembly", "GRCh38",
           "--tss", paths["tss"], "--transcript-tss", paths["transcript_tss"],
           "--promoters", paths["promoters"], "--genome-id", plan["genome_id"],
           "--annotation-release", plan["annotation_release"],
           "--promoter-definition-id", plan["promoter_definition_id"],
           "--output", scratch / "result")
    PHASE = f"chrom-{task['chrom']}-promotion"
    # Copy to the destination filesystem before rename, then verify every file.
    with tempfile.TemporaryDirectory(prefix=f".chrom-{task['chrom']}-", dir=target.parent) as temp:
        staging = Path(temp) / "package"
        shutil.copytree(scratch / "result", staging)
        result = annotation.validated_package(staging)
        result["production_plan_sha256"] = digest
        result["source_commit"] = plan["source_commit"]
        result["durable_inputs"] = {k: task[k] for k in ("evidence", *DIMENSIONS)}
        result["slurm_job_id"] = os.environ.get("SLURM_JOB_ID")
        (staging / "manifest.json").write_text(json.dumps(result, indent=2) + "\n")
        if target.exists():
            verify_task(target, digest)
        else:
            os.rename(staging, target)
    print(f"I: Completed chromosome {task['chrom']}")


def finalize(args):
    global PHASE
    plan = load_plan(args)
    digest = annotation.sha256(args.run_root / "plan.json")
    target = args.run_root / "final"
    if target.exists():
        verify_task(target, digest)
        print("I: Reusing completed regulatory catalog")
        return
    PHASE = "final-validation"
    files = []
    manifests = []
    for task in plan["tasks"]:
        package = args.run_root / "tasks" / f"chrom-{task['chrom']}"
        verify_task(package, digest)
        manifests.append(record(package / "manifest.json"))
        files.append(package / "tp73_anchor_regulatory_membership.parquet")
    paths = "[" + ",".join(annotation.sql_string(p.resolve()) for p in files) + "]"
    with tempfile.TemporaryDirectory(prefix=".reg-final-", dir=args.run_root) as temp:
        staging = Path(temp) / "package"
        staging.mkdir()
        output = staging / "tp73_anchor_regulatory_membership.parquet"
        sql(args, f"""
CREATE TABLE membership AS SELECT * FROM read_parquet({paths},hive_partitioning=false);
SELECT CASE WHEN (SELECT count(*)-count(DISTINCT (chrom,anchor_start,anchor_end)) FROM membership)<>0
  OR EXISTS (SELECT 1 FROM membership WHERE promoter_core AND NOT promoter_extended)
  THEN error('invalid global regulatory membership') END;
COPY (SELECT * FROM membership ORDER BY try_cast(chrom AS INTEGER),anchor_start,anchor_end)
TO {annotation.sql_string(output)} (FORMAT PARQUET,COMPRESSION ZSTD);
""")
        annotation.publish(staging, target, {
            "schema_version": 1, "kind": "genome_regulatory_membership", "state": "complete",
            "production_plan_sha256": digest, "source_commit": plan["source_commit"],
            "assembly": "GRCh38", "chromosomes": plan["chromosomes"],
            "coordinate_audit": plan["audit"], "chromosome_manifests": manifests,
            "annotation_release": plan["annotation_release"],
            "promoter_definition_id": plan["promoter_definition_id"],
        })
    print(f"I: Finalized {len(files)} chromosomes: {target}")


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    p = commands.add_parser("prepare", help="pin exact evidence and GTF inputs; no jobs submitted")
    for name in ("source", "evidence-package", "annotation-catalog", "audit", "duckdb"):
        p.add_argument("--"+name, type=Path, required=True)
    p.add_argument("--chromosomes", default=",".join(map(str, range(1,23))))
    p.add_argument("--genome-id", default="homo_sapiens_grch38_ensembl113_primary")
    p.add_argument("--annotation-release", default="ensembl_113")
    p.add_argument("--promoter-definition-id", default="tss_upstream_2000_downstream_500_v1")
    p.set_defaults(function=prepare)
    s = commands.add_parser("setup", help="download audited source and build shared features inside a job")
    s.add_argument("--gff", type=Path, help="optional already-local source, mainly for offline tests")
    s.set_defaults(function=setup)
    r = commands.add_parser("run-task", help="stage one chromosome on scratch and publish atomically")
    r.add_argument("--task-index", type=int, help="default: SLURM_ARRAY_TASK_ID")
    r.add_argument("--scratch-root", type=Path, required=True)
    r.set_defaults(function=run_task)
    f = commands.add_parser("finalize", help="validate every planned chromosome and publish combined membership")
    f.set_defaults(function=finalize)
    for command in (p,s,r,f):
        command.add_argument("--run-root", type=Path, required=True)
    return result


def main():
    args = parser().parse_args()
    args.run_root = args.run_root.resolve()
    if hasattr(signal, "SIGUSR1"):
        signal.signal(signal.SIGUSR1, status_signal)
    try:
        args.function(args)
        return 0
    except (OSError, ValueError, KeyError, subprocess.CalledProcessError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Build/query a bounded promoter cofactor package; Python stdlib + DuckDB CLI."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


BANDS = ("overlap", "adjacent_0_5", "gap_6_20", "gap_21_50", "gap_51_100", "gap_101_150")
GENES = {
    "Previously associated": ("E2F1", "SP1", "REST", "TP53", "TP63"),
    "Steffen": ("AHCTF1", "LMX1B", "PATZ1", "POU2F2", "TFAP2C"),
    "Glen": ("SUB1", "HMGB3", "IRF9", "YBX3"),
    "Birgit": ("ZBED1", "TGIF2"),
}
TABLES = ("cofactor_distance_isoform_comparison", "cofactor_distance_enrichment",
          "cofactor_distance_frequency_enrichment", "jaspar_matrix", "jaspar_matrix_species")
KIND = "tp73_promoter_collaboration"


def quote(value):
    return "'" + str(value).replace("'", "''") + "'"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    with Path(path).open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def record(path, root):
    return {"path": str(path.relative_to(root)), "bytes": path.stat().st_size,
            "sha256": digest(path)}


def checked(root, row):
    # Package manifests cannot escape their owning directory, even via symlinks.
    path = (root / row["path"]).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("manifest path escapes package")
    if path.stat().st_size != int(row["bytes"]) or digest(path) != row["sha256"]:
        raise ValueError(f"file differs from manifest: {path}")
    return path


def parquet(paths):
    if isinstance(paths, (str, Path)):
        paths = [paths]
    return "read_parquet([" + ",".join(map(quote, paths)) + "], hive_partitioning=false)"


def copy_sql(query, path):
    return f"COPY ({query}) TO {quote(path)} (FORMAT PARQUET, COMPRESSION ZSTD);\n"


def sql(args, query, database=":memory:", cwd=None, rows=False):
    settings = (f"SET threads=2; SET memory_limit={quote(args.memory_limit)};"
                f" SET temp_directory={quote(args.temp)}; SET max_temp_directory_size='32GB';"
                " SET preserve_insertion_order=false;\n")
    command = [args.duckdb, "-batch", "-no-init"]
    if rows:
        command += ["-json"]
        if str(database) != ":memory:":
            command += ["-readonly"]
    result = subprocess.run(command + [str(database)], input=settings + query,
                            text=True, capture_output=True, cwd=cwd)
    if result.returncode:
        raise ValueError(result.stderr.strip() or "DuckDB failed")
    return json.loads(result.stdout or "[]") if rows else None


@contextmanager
def transaction(root, name):
    """Serialize retries, then publish an immutable directory by atomic rename."""
    root.mkdir(parents=True, exist_ok=True)
    with (root / (name + ".lock")).open("a") as lock:
        if os.name == "nt":
            import msvcrt
            lock.write("0")
            lock.flush()
            lock.seek(0)
            msvcrt.locking(lock.fileno(), msvcrt.LK_LOCK, 1)
        else:
            import fcntl
            fcntl.flock(lock, fcntl.LOCK_EX)
        final = root / name
        if final.exists():
            yield final, True
        else:
            # Interrupted durable attempts are preserved, never interpreted as complete.
            stage = Path(tempfile.mkdtemp(prefix=name + ".attempt-", dir=root))
            yield stage, False
            os.rename(stage, final)


def validate_product(path, identity):
    manifest = read_json(path / "complete.json")
    if manifest["identity"] != identity:
        raise ValueError("checkpoint belongs to another input/configuration")
    for row in manifest["files"]:
        checked(path, row)
    return manifest


def finish_product(stage, identity, paths, **extra):
    write_json(stage / "complete.json", {"identity": identity,
               "files": [record(p, stage) for p in paths], **extra})


def select_panel(rows, limit):
    names = {r["motif_id"]: r["motif_name"] for r in rows}
    selected, reasons = [], {}

    def add(motif, reason):
        if motif not in selected:
            if len(selected) >= limit:
                return
            selected.append(motif)
        reasons.setdefault(motif, []).append(reason)

    requested = []
    for group, genes in GENES.items():
        for gene in genes:
            matches = sorted(m for m, name in names.items() if name.upper() == gene)
            requested.append({"group": group, "gene": gene, "motif_ids": matches,
                              "status": "available" if matches else "no_exact_name_in_result"})
            for motif in matches:
                add(motif, "requested:" + gene)
    required = {m for row in requested for m in row["motif_ids"]}
    if not required <= set(selected):
        raise ValueError("panel limit is smaller than the prespecified available motifs")
    queues = []
    for band in BANDS:
        for rank in ("ta_enrichment_rank", "dn_enrichment_rank", "ta_depletion_rank",
                     "dn_depletion_rank", "isoform_difference_absolute_rank"):
            candidates = [r for r in rows if r["distance_band"] == band and r.get(rank) is not None]
            queues.append((band + ":" + rank, sorted(candidates,
                           key=lambda r: (r[rank], r["motif_id"]))))
    # Balanced coverage of both isoforms, both signs, and all exclusive distances.
    for i in range(max((len(q) for _, q in queues), default=0)):
        for reason, queue in queues:
            if i < len(queue):
                add(queue[i]["motif_id"], reason)
        if len(selected) >= limit:
            break
    return [{"motif_id": m, "motif_name": names[m], "selection": reasons[m]}
            for m in selected], requested


def prepare(args):
    run = args.analysis_run.resolve()
    config = read_json(run / "plan/run_config.json")
    final = run / "final/distance_enrichment"
    upstream = read_json(final / "manifest.json")  # No provisional fallback.
    regulatory = config.get("regulatory_selection") or {}
    if (upstream.get("schema_version") != 6 or upstream["run_id"] != config["run_id"]
            or upstream.get("regulatory_selection") != regulatory
            or regulatory.get("subset") != "promoter_extended"
            or config.get("source_score_floor") != -1
            or config.get("tax_group") != "vertebrates"
            or config.get("context_flank_bp") != 150
            or config.get("distance_bands") != list(BANDS)
            or config.get("excluded_series") != ["skmel29_1"]):
        raise ValueError("requires finalized vertebrate, extended-promoter, floor -1, 150bp results")
    if args.panel_size < 1 or args.panel_size > 200 or args.max_bytes < 1:
        raise ValueError("panel size must be 1..200; byte limit must be positive")
    scan_inventory = run / "plan/scan_files.tsv"
    if digest(scan_inventory) != config["scan_files_sha256"]:
        raise ValueError("scan inventory changed")
    split = run / "input/anchors"
    split_manifest = read_json(split / "complete.json")
    if (split_manifest["regulatory_selection"] != regulatory
            or split_manifest["source_sha256"] != config["anchor_evidence_sha256"]
            or digest(split / "anchor_files.tsv") != split_manifest["inventory_sha256"]):
        raise ValueError("anchor split is not the completed promoter cohort")
    identity = {"upstream": digest(final / "manifest.json"), "config": digest(run / "plan/run_config.json"),
                "anchor_split": digest(split / "complete.json"), "panel_size": args.panel_size,
                "max_bytes": args.max_bytes, "exporter_sha256": digest(Path(__file__)),
                "regulatory_features_manifest_sha256": digest(args.regulatory_features / "manifest.json")}
    with transaction(args.run_root, "prepared") as (stage, exists):
        if exists:
            validate_product(stage, identity)
            print("Reused validated preparation")
            return
        copied = []
        for name in TABLES:
            source = checked(final, upstream["tables"][name])
            target = stage / (name + ".parquet")
            shutil.copyfile(source, target)
            copied.append(target)
        comparison = stage / "cofactor_distance_isoform_comparison.parquet"
        rows = sql(args, f"SELECT * FROM {parquet(comparison)} ORDER BY motif_id,distance_band_order;", rows=True)
        if (not rows or any(float(r["source_score_floor"]) != -1 or float(r["positive_threshold"]) != 0
                            or r["tax_group"] != "vertebrates" for r in rows)
                or len(rows) != upstream["task_count"] * len(BANDS)
                or len({(r["motif_id"], r["distance_band"]) for r in rows}) != len(rows)
                or {r["distance_band"] for r in rows} != set(BANDS)):
            raise ValueError("comparison must contain every motif/band once, at common threshold zero")
        panel, requested = select_panel(rows, args.panel_size)
        if not panel:
            raise ValueError("no candidate panel can be selected")
        with (split / "anchor_files.tsv").open() as handle:
            inventory = list(csv.DictReader(handle, delimiter="\t"))
        if sorted(r["chrom"] for r in inventory) != sorted(config["chromosomes"]):
            raise ValueError("anchor inventory chromosome coverage differs")
        paths = [checked(split, {**r, "path": r["relative_path"]}) for r in inventory]
        evidence = args.temp / "source_anchor_evidence.parquet"
        shutil.copyfile(config["anchor_evidence"], evidence)
        if digest(evidence) != config["anchor_evidence_sha256"]:
            raise ValueError("source CUT&RUN evidence changed")
        columns = sql(args, f"DESCRIBE SELECT * FROM {parquet(paths)};", rows=True)
        support = [r["column_name"] for r in columns if r["column_name"].startswith("supported_")]
        if not support or any(not re.fullmatch(r"[A-Za-z0-9_]+", s) for s in support):
            raise ValueError("invalid CUT&RUN support columns")
        depths = [s.replace("supported_", "depth_", 1) for s in support]
        evidence_guard = " OR ".join(f'a."{s}" IS DISTINCT FROM e."{s}" OR e."{d}" IS NULL OR NOT isfinite(e."{d}") OR e."{d}"<0 OR a."{s}"<>(e."{d}">0)'
                                     for s, d in zip(support, depths))
        query = f"""
CREATE TABLE a AS SELECT row_number() OVER (ORDER BY CAST(chrom AS VARCHAR),anchor_start,anchor_end)::BIGINT AS anchor_id,
  * REPLACE(CAST(chrom AS VARCHAR) AS chrom) FROM {parquet(paths)};
SELECT CASE WHEN (SELECT count(*)-count(DISTINCT (chrom,anchor_start,anchor_end)) FROM a)<>0
  OR EXISTS (SELECT 1 FROM a WHERE anchor_start<0 OR anchor_end<=anchor_start OR anchor_score IS NULL)
  OR (SELECT count(*) FROM a)<>{int(split_manifest['anchors'])}
  THEN error('anchor key/size is invalid') END;
CREATE TABLE e AS SELECT e.* FROM {parquet(evidence)} e
 SEMI JOIN a ON a.chrom=CAST(e.chrom AS VARCHAR) AND a.anchor_start=e.anchor_start AND a.anchor_end=e.anchor_end;
SELECT CASE WHEN (SELECT count(*) FROM e)<>(SELECT count(*) FROM a)
 OR (SELECT count(*)-count(DISTINCT (chrom,anchor_start,anchor_end)) FROM e)<>0
 OR EXISTS (SELECT 1 FROM a JOIN e ON a.chrom=CAST(e.chrom AS VARCHAR)
   AND a.anchor_start=e.anchor_start AND a.anchor_end=e.anchor_end
   WHERE a.anchor_score IS DISTINCT FROM e.anchor_score OR {evidence_guard})
 THEN error('source CUT&RUN depth/support does not agree with the analysis cohort') END;
""" + copy_sql('SELECT a.*,' + ','.join(f'e."{d}"' for d in depths) +
               ' FROM a JOIN e ON a.chrom=CAST(e.chrom AS VARCHAR) AND a.anchor_start=e.anchor_start AND a.anchor_end=e.anchor_end ORDER BY anchor_id', stage / "anchors.parquet")
        sql(args, query)
        copied.append(stage / "anchors.parquet")
        # Copy the exact feature bridges already audited upstream, not a new promoter definition.
        regulatory_manifest = Path(regulatory["manifest"])
        if digest(regulatory_manifest) != regulatory["manifest_sha256"]:
            raise ValueError("regulatory manifest changed")
        reg = read_json(regulatory_manifest)
        if reg.get("state") != "complete" or reg.get("coordinate_audit", {}).get("status") != "passed":
            raise ValueError("regulatory annotation lacks a completed coordinate audit")
        links, owners = [], []
        for item in reg["chromosome_manifests"]:
            manifest_path = Path(item["path"])
            checked(manifest_path.parent, {**item, "path": manifest_path.name})
            part = read_json(manifest_path)
            by_name = {r["path"]: r for r in part["files"]}
            links.append(checked(manifest_path.parent, by_name["tp73_anchor_regulatory_feature.parquet"]))
            owners.append(checked(manifest_path.parent, by_name["regulatory_feature_gene.parquet"]))
        features_root = args.regulatory_features.resolve()
        fm = read_json(features_root / "manifest.json")
        if (fm.get("assembly") != "GRCh38"
                or fm.get("source_sha256") != regulatory["coordinate_audit"]["gff_sha256"]):
            raise ValueError("feature source differs from the coordinate-audited source")
        features = checked(features_root, next(r for r in fm["files"] if r["path"] == "regulatory_feature.parquet"))
        query = f"""
CREATE TABLE a AS SELECT * FROM {parquet(stage / 'anchors.parquet')};
CREATE TABLE bridge AS SELECT DISTINCT a.anchor_id, l.regulatory_feature_id,
  l.regulatory_definition_id, l.overlap_bp, l.anchor_fully_within
FROM a JOIN {parquet(links)} l ON a.chrom=CAST(l.chrom AS VARCHAR)
 AND a.anchor_start=l.anchor_start AND a.anchor_end=l.anchor_end
WHERE l.regulatory_definition_id='promoter_extended';
SELECT CASE WHEN EXISTS (SELECT anchor_id FROM a EXCEPT SELECT anchor_id FROM bridge)
  THEN error('selected anchor has no extended promoter membership') END;
CREATE TABLE promoter AS SELECT f.* FROM {parquet(features)} f
WHERE f.regulatory_feature_id IN (SELECT regulatory_feature_id FROM bridge);
SELECT CASE WHEN EXISTS (SELECT regulatory_feature_id FROM bridge EXCEPT SELECT regulatory_feature_id FROM promoter)
 THEN error('promoter dimension is incomplete') END;
"""
        queries = {"anchor_promoter": "SELECT * FROM bridge ORDER BY anchor_id,regulatory_feature_id",
                   "promoter": "SELECT * FROM promoter ORDER BY chrom,start,regulatory_feature_id",
                   "promoter_gene": f"SELECT DISTINCT g.* FROM {parquet(owners)} g WHERE g.regulatory_feature_id IN (SELECT regulatory_feature_id FROM bridge)"}
        for name, select in queries.items():
            target = stage / (name + ".parquet")
            query += copy_sql(select, target)
            copied.append(target)
        sql(args, query)
        with scan_inventory.open() as handle:
            source_rows = list(csv.DictReader(handle, delimiter="\t"))
        selected_ids = {p["motif_id"] for p in panel}
        source_rows = [r for r in source_rows if r["motif_id"] in selected_ids]
        expected = {(m, c, s) for m in selected_ids for c in config["chromosomes"] for s in ("+", "-")}
        if (len(source_rows) != len(expected)
                or {(r["motif_id"], r["chrom"], r["strand"]) for r in source_rows} != expected
                or any(float(r["minimum_score"]) != -1 for r in source_rows)):
            raise ValueError("selected motif scan is missing, duplicated, or censored above floor -1")
        score_configs = set()
        score_keys = ("genome_id", "motif_set_id", "score_mode", "pseudocount",
                      "background_model_id", "pseudocount_scheme", "n_policy")
        for row in source_rows:
            partitions = dict(part.split("=", 1) for part in Path(row["absolute_path"]).parts[:-1] if "=" in part)
            score_configs.add(tuple(partitions[key] for key in score_keys))
        if len(score_configs) != 1:
            raise ValueError("source scan mixes genome/score configurations")
        score_config = dict(zip(score_keys, score_configs.pop()))
        if score_config["genome_id"] != "homo_sapiens_grch38_ensembl113_primary":
            raise ValueError("source scan genome is not the supported GRCh38 reference")
        plan = {"identity": identity, "panel": panel, "requested_genes": requested,
                "chromosomes": config["chromosomes"], "sources": source_rows,
                "score_configuration": score_config,
                "source_run_config": config, "source_result_manifest": upstream,
                "regulatory_features_manifest_sha256": digest(features_root / "manifest.json")}
        write_json(stage / "plan.json", plan)
        copied.append(stage / "plan.json")
        finish_product(stage, identity, copied)
        print(f"Prepared {len(panel)} motifs and {split_manifest['anchors']} anchors")


def feature_sql(anchors, plus, minus, chrom, motif, output):
    """Lossless bin prefilter followed by exact half-open interval geometry."""
    return f"""
CREATE TABLE a AS SELECT * FROM {parquet(anchors)} WHERE chrom={quote(chrom)};
CREATE TABLE raw AS SELECT start,"end",score, '+' AS strand FROM {parquet(plus)} UNION ALL
 SELECT start,"end",score, '-' AS strand FROM {parquet(minus)};
SELECT CASE WHEN EXISTS (SELECT 1 FROM raw WHERE start IS NULL OR "end" IS NULL
 OR start<0 OR "end"<=start OR score IS NULL OR NOT isfinite(score) OR score < -1)
 OR (SELECT count(*)-count(DISTINCT (start,"end",strand)) FROM raw)<>0
 THEN error('invalid or duplicated source locus') END;
CREATE TABLE h AS SELECT start::BIGINT AS hit_start, "end"::BIGINT AS hit_end,
 max(score)::DOUBLE AS best_score,
 max(score) FILTER (WHERE strand='+')::DOUBLE AS plus_score,
 max(score) FILTER (WHERE strand='-')::DOUBLE AS minus_score
FROM raw GROUP BY start,"end";
CREATE TABLE width AS SELECT greatest(1,150+coalesce((SELECT max(anchor_end-anchor_start) FROM a),0)
 +coalesce((SELECT max(hit_end-hit_start) FROM h),0))::BIGINT AS w;
CREATE TABLE ab AS SELECT a.*,floor(((anchor_start+anchor_end)/2.0)/w)::BIGINT AS bin FROM a,width;
CREATE TABLE hb AS SELECT h.*,floor(((hit_start+hit_end)/2.0)/w)::BIGINT+delta AS bin
 FROM h,width,(VALUES (-1),(0),(1)) offsets(delta);
CREATE TABLE candidates AS
WITH geometry AS (
 SELECT a.anchor_id,h.* EXCLUDE(bin),
 greatest(a.anchor_start,h.hit_start)-least(a.anchor_end,h.hit_end) AS interval_distance_bp,
 CASE WHEN h.hit_end<=a.anchor_start THEN 'left'
      WHEN h.hit_start>=a.anchor_end THEN 'right' ELSE 'overlap' END AS genomic_side
 FROM ab a JOIN hb h USING(bin)
 WHERE greatest(a.anchor_start,h.hit_start)-least(a.anchor_end,h.hit_end)<=150
), classified AS (
 SELECT *,CASE WHEN interval_distance_bp<0 THEN 'overlap'
 WHEN interval_distance_bp<=5 THEN 'adjacent_0_5' WHEN interval_distance_bp<=20 THEN 'gap_6_20'
 WHEN interval_distance_bp<=50 THEN 'gap_21_50' WHEN interval_distance_bp<=100 THEN 'gap_51_100'
 ELSE 'gap_101_150' END AS distance_band FROM geometry
)
SELECT *,count(*) OVER band AS n_source_loci,
 count(*) FILTER (WHERE best_score>=0) OVER band AS n_score_zero_loci,
 row_number() OVER (PARTITION BY anchor_id,distance_band
 ORDER BY best_score DESC,abs(interval_distance_bp),hit_start,hit_end) AS priority
FROM classified WINDOW band AS (PARTITION BY anchor_id,distance_band);
""" + copy_sql(f"""SELECT anchor_id,{quote(motif)}::VARCHAR AS motif_id,distance_band,
 hit_start,hit_end,best_score,plus_score,minus_score,
 CASE WHEN plus_score=minus_score THEN '.' WHEN minus_score IS NULL OR plus_score>minus_score
      THEN '+' ELSE '-' END::VARCHAR AS best_strand,
 interval_distance_bp,genomic_side,n_source_loci::BIGINT AS n_source_loci,
 n_score_zero_loci::BIGINT AS n_score_zero_loci
 FROM candidates WHERE priority=1 ORDER BY anchor_id,distance_band""", output)


def motif_task(args):
    prepared = args.run_root / "prepared"
    plan = read_json(prepared / "plan.json")
    if plan["identity"]["exporter_sha256"] != digest(Path(__file__)):
        raise ValueError("exporter changed after preparation; use the pinned source checkout")
    validate_product(prepared, plan["identity"])
    if args.task_index < 0:
        raise ValueError("negative task index")
    if args.task_index >= len(plan["panel"]):
        print("No selected motif assigned to this array slot")
        return
    motif = plan["panel"][args.task_index]["motif_id"]
    if not re.fullmatch(r"MA[0-9]+\.[0-9]+", motif):
        raise ValueError("invalid motif identifier")
    identity = {"plan_sha256": digest(prepared / "plan.json"), "motif_id": motif}
    for chrom in plan["chromosomes"]:
        with transaction(args.run_root / "tasks" / motif, "chrom-" + chrom) as (stage, exists):
            if exists:
                validate_product(stage, identity)
                continue
            print(f"{motif} chromosome {chrom}: staging/scoring-summary", flush=True)
            inputs = [r for r in plan["sources"] if r["motif_id"] == motif and r["chrom"] == chrom]
            with tempfile.TemporaryDirectory(prefix="motif-input-", dir=args.temp) as local:
                local = Path(local)
                paths = {}
                for r in inputs:
                    target = local / ("plus.parquet" if r["strand"] == "+" else "minus.parquet")
                    shutil.copyfile(r["absolute_path"], target)
                    checked(local, {**r, "path": target.name})
                    paths[r["strand"]] = target
                out = local / "feature.parquet"
                sql(args, feature_sql(prepared / "anchors.parquet", paths["+"], paths["-"], chrom, motif, out))
                shutil.copyfile(out, stage / out.name)
                finish_product(stage, identity, [stage / out.name])
    print(f"{motif}: all chromosome checkpoints complete", flush=True)


def package_schema(feature_files):
    query = ""
    for name in (*TABLES, "anchors", "anchor_promoter", "promoter", "promoter_gene"):
        query += f"CREATE VIEW {name} AS SELECT * FROM {parquet(name + '.parquet')};\n"
    query += f"CREATE VIEW feature AS SELECT * FROM {parquet(feature_files)};\n"
    query += "CREATE TABLE distance_band(distance_band VARCHAR, distance_band_order INTEGER);\n"
    query += "INSERT INTO distance_band VALUES " + ",".join(f"({quote(b)},{i})" for i, b in enumerate(BANDS)) + ";\n"
    # Empty sparse feature rows imply zero only inside this completed panel/cohort.
    query += """
CREATE MACRO anchor_features(selected_motif, chromosome, region_start, region_end) AS TABLE
SELECT a.*,b.distance_band,f.* EXCLUDE(anchor_id,motif_id,distance_band,n_source_loci,n_score_zero_loci),
 selected_motif AS motif_id,coalesce(f.n_source_loci,0)::BIGINT AS n_source_loci,
 coalesce(f.n_score_zero_loci,0)::BIGINT AS n_score_zero_loci
FROM anchors a CROSS JOIN distance_band b
LEFT JOIN feature f ON a.anchor_id=f.anchor_id AND f.distance_band=b.distance_band AND f.motif_id=selected_motif
WHERE a.chrom=chromosome AND a.anchor_start<region_end AND a.anchor_end>region_start
 AND CASE WHEN NOT EXISTS (SELECT 1 FROM selected_motif_panel WHERE motif_id=selected_motif)
 THEN error('motif detail is not included; absence is unknown, not negative') ELSE true END;
"""
    return query


def report(rows, requested):
    def cell(value):
        if isinstance(value, float):
            value = f"{value:.4g}"
        return html.escape("" if value is None else str(value))
    cols = ("motif_name", "motif_id", "distance_band", "ta_adjusted_odds_ratio", "dn_adjusted_odds_ratio",
            "ta_vs_dn_odds_ratio_ratio", "positive_anchor_fraction", "negative_anchor_fraction")
    labels = ("Factor", "JASPAR motif", "Distance band", "TA odds ratio", "DN odds ratio",
              "TA / DN odds ratio", "Frequency at score 0", "Negative fraction (below -1)")
    body = "".join("<tr>" + "".join(f"<td>{cell(r.get(c))}</td>" for c in cols) + "</tr>" for r in rows)
    absent = ", ".join(r["gene"] for r in requested if not r["motif_ids"])
    return """<!doctype html><html lang="en"><meta charset="utf-8"><title>TP73 promoter cofactors</title>
<style>body{font:15px system-ui;margin:24px;color:#202124}table{border-collapse:collapse}th,td{padding:6px 12px;text-align:right;border-bottom:1px solid #ddd}th{position:sticky;top:0;background:#edf1f2}td:first-child,th:first-child{text-align:left}</style>
<h1>TP73 promoter cofactors</h1><p>Extended Ensembl promoters; cofactor score &ge; 0; negative reference: no retained locus at score &ge; -1.
TA/DN odds ratios are associations with CUT&amp;RUN support, not causal effects. All eligible promoter anchors remain in the detailed package.
Frequency is the fraction of anchors with at least one score-zero-positive locus in the indicated exclusive band.
This is a selected detail panel, not a new multiple-testing family. Full estimates, confidence intervals, support and unchanged adjusted p-values are in the query tables.</p>
<p>Exact-name motifs unavailable: """ + html.escape(absent or "none") + "</p><table><thead><tr>" + "".join(
        "<th>" + html.escape(label) + "</th>" for label in labels) + "</tr></thead><tbody>" + body + "</tbody></table></html>\n"


def finalize(args):
    prepared = args.run_root / "prepared"
    plan = read_json(prepared / "plan.json")
    if plan["identity"]["exporter_sha256"] != digest(Path(__file__)):
        raise ValueError("exporter changed after preparation; use the pinned source checkout")
    base = validate_product(prepared, plan["identity"])
    identity = {"plan_sha256": digest(prepared / "plan.json")}
    with transaction(args.run_root, "package") as (stage, exists):
        if exists:
            validate_product(stage, identity)
            print("Reused validated complete package")
            return
        paths = []
        for row in base["files"]:
            if row["path"].endswith(".parquet"):
                target = stage / row["path"]
                shutil.copyfile(prepared / row["path"], target)
                paths.append(target)
        feature_files = []
        for item in plan["panel"]:
            motif = item["motif_id"]
            sources = []
            for chrom in plan["chromosomes"]:
                part = args.run_root / "tasks" / motif / ("chrom-" + chrom)
                validate_product(part, {**identity, "motif_id": motif})
                sources.append(part / "feature.parquet")
            name = "feature_" + motif + ".parquet"
            sql(args, copy_sql(f"SELECT * FROM {parquet(sources)} ORDER BY anchor_id,distance_band", stage / name))
            paths.append(stage / name)
            feature_files.append(name)
            if sum(p.stat().st_size for p in paths) > plan["identity"]["max_bytes"]:
                raise ValueError("package exceeds byte cap; preserved unpublished attempt, select fewer motifs in a new run")
        schema = package_schema(feature_files)
        schema = "CREATE TABLE selected_motif_panel(motif_id VARCHAR,motif_name VARCHAR);\nINSERT INTO selected_motif_panel VALUES " + ",".join(
            f"({quote(r['motif_id'])},{quote(r['motif_name'])})" for r in plan["panel"]) + ";\n" + schema
        (stage / "schema.sql").write_text(schema)
        sql(args, schema, database=stage / "cofactors.duckdb", cwd=stage)
        # Compare reconstructed occurrence counts against finalized upstream estimates.
        sql(args, """
SELECT CASE WHEN EXISTS (
 SELECT 1 FROM feature GROUP BY anchor_id,motif_id,distance_band HAVING count(*)<>1)
 OR EXISTS (SELECT anchor_id FROM feature EXCEPT SELECT anchor_id FROM anchors)
 OR EXISTS (
 SELECT 1 FROM cofactor_distance_isoform_comparison c JOIN selected_motif_panel p USING(motif_id)
 LEFT JOIN (SELECT motif_id,distance_band,count(*) AS source_n,
   count(*) FILTER (WHERE best_score>=0) AS positive_n FROM feature GROUP BY ALL) f
 USING(motif_id,distance_band)
 WHERE c.anchors_total<>(SELECT count(*) FROM anchors)
 OR c.anchors_source_present<>coalesce(f.source_n,0)
 OR c.anchors_positive<>coalesce(f.positive_n,0))
 THEN error('reconstructed features disagree with finalized cofactor frequencies') END;
""", database=stage / "cofactors.duckdb", cwd=stage)
        rows = sql(args, "SELECT c.* FROM cofactor_distance_isoform_comparison c JOIN selected_motif_panel p USING(motif_id) ORDER BY motif_name,motif_id,distance_band_order;",
                   database=stage / "cofactors.duckdb", cwd=stage, rows=True)
        (stage / "overview.html").write_text(report(rows, plan["requested_genes"]))
        shutil.copyfile(Path(__file__), stage / "export_promoter_collaboration.py")
        (stage / "README.md").write_text(PACKAGE_README)
        paths += [stage / name for name in ("schema.sql", "cofactors.duckdb", "overview.html",
                                            "export_promoter_collaboration.py", "README.md")]
        manifest = {"kind": KIND, "schema_version": 1, "state": "complete", "identity": identity,
                    "assembly": "GRCh38", "taxon_id": 9606, "coordinate_mode": "bed_0based_half_open",
                    "chromosomes": plan["chromosomes"], "distance_bands": list(BANDS),
                    "retention": "one strongest physical locus per included anchor/motif/exclusive band",
                    "source_score_floor": -1, "positive_threshold": 0,
                    "score_configuration": plan["score_configuration"],
                    "panel": plan["panel"], "requested_genes": plan["requested_genes"],
                    "h3k4me3_model_effects": "not_included", "complete_genome_scan": False,
                    "source_run_config": plan["source_run_config"],
                    "source_result_manifest": plan["source_result_manifest"],
                    "regulatory_features_manifest_sha256": plan["regulatory_features_manifest_sha256"],
                    "files": [record(p, stage) for p in paths]}
        write_json(stage / "manifest.json", manifest)
        paths.append(stage / "manifest.json")
        finish_product(stage, identity, paths)
        total = sum(p.stat().st_size for p in paths) + (stage / "complete.json").stat().st_size
        if total > plan["identity"]["max_bytes"]:
            raise ValueError("complete package exceeds byte cap; no published package")
        print(f"Published complete portable package: {total} bytes")


PACKAGE_README = """# TP73 promoter cofactor package

Open overview.html without a server. All coordinates are GRCh38, BED 0-based
half-open. These are sequence motif matches, not proof of cofactor occupancy.
Ensembl extended promoter membership selects the TP73 anchors, NOT the cofactor
positions: the latter may fall outside the promoter in a 150 bp neighbourhood.
Promoter-gene ownership is many-to-many. Native and TSS-derived links retain
their provenance; do not multiply observations by joining through genes.

Only the listed chromosomes, anchors and selected motifs have detailed local
coverage. All anchors, including negative CUT&RUN anchors, are retained.
No X/Y, mitochondrial, gene-complete or genome-complete coverage is implied.
All vertebrate motif overview estimates and their original uncertainty/BH
correction are retained. Source species is provenance, not an exclusivity claim.
SAOS2 and SK-MEL-29_2 are distinct cell series, not replicates; SK-MEL-29_1 is
excluded. H3K4me3 model effects and raw signal tracks are NOT in this package.

For each anchor/motif/exclusive distance band, feature stores the strongest
physical locus, both retained orientation scores, counts at >= -1 and >= 0,
interval gap (negative for overlap, zero for abutment), and genomic side.
Ties use smallest absolute gap, then start/end. A tied strand is '.'. Missing
orientation scores are censored below -1, not zero. Sparse absent feature rows
mean no locus >= -1 ONLY for included anchors/motifs; best score is unknown
below that floor. Other thresholds' counts and all individual hits cannot be
reconstructed. The per-band maxima DO support any higher-threshold presence
test. Summing exclusive counts gives all-150 counts, but frequencies must not
be summed. All-150 best is the maximum of the six band maxima.

Python 3.9+ and a DuckDB CLI are sufficient, with no downloads/extensions:

    python3 export_promoter_collaboration.py inspect --package .
    python3 export_promoter_collaboration.py region --package . --motif MA0653.1 --chrom 1 --start 1000000 --end 1100000 --output region_example

The region command exports matching TP73 anchors, their strongest cofactor
loci, raw scores and summarized CUT&RUN evidence. Its coordinates are never
clipped to the viewing window. It fails if the result exceeds --max-rows.
It writes a GENtle workflow for importing the BED6 spans with an explicit
GRCh38 reference. BED score is a placeholder 0, NOT the raw motif score;
use scores/evidence in evidence.json. Importing does not load the genome or
create a complete genome-scanning provider. GENtle can navigate these regions
when its matching genome reference is available. The workflow also exports
GENtle's canonical region JSON, with identities computed by GENtle itself.

Direct queries: change into the package directory before opening
`duckdb -readonly -no-init cofactors.duckdb`. Relative paths survive moving the
whole directory. Use `anchor_features('MA0653.1','1',1000000,1100000)` for the
zero-complete six-band rows at overlapping included anchors. An unlisted motif
is an error, not a negative. Always inspect manifest.json for coverage.
schema.sql rebuilds the tiny query index in a NEW DuckDB database. The supplied
index uses views, not a second copy of the data. No unrestricted SQL service is
exposed. A future GENtle adapter can use this bounded region contract directly.
"""


def inspect_or_region(args):
    package = args.package.resolve()
    manifest = read_json(package / "manifest.json")
    if manifest.get("kind") != KIND or manifest.get("schema_version") != 1 or manifest.get("state") != "complete":
        raise ValueError("not a supported completed collaborator package")
    if args.command == "inspect":
        validate_product(package, manifest["identity"])
        print(json.dumps({k: manifest[k] for k in ("assembly", "chromosomes", "retention", "panel", "requested_genes")}, indent=2))
        return
    if (args.chrom not in manifest["chromosomes"] or args.motif not in {r["motif_id"] for r in manifest["panel"]}
            or args.start < 0 or args.end <= args.start or args.max_rows < 1):
        raise ValueError("invalid region or chromosome/motif outside package coverage")
    # Bounded lookups validate only files they read, not the entire multi-motif package.
    complete = read_json(package / "complete.json")
    checked(package, next(r for r in complete["files"] if r["path"] == "manifest.json"))
    for name in ("anchors.parquet", "cofactors.duckdb", "feature_" + args.motif + ".parquet",
                 "anchor_promoter.parquet", "promoter.parquet", "promoter_gene.parquet"):
        checked(package, next(r for r in manifest["files"] if r["path"] == name))
    query = f"SELECT * FROM anchor_features({quote(args.motif)},{quote(args.chrom)},{args.start},{args.end}) ORDER BY anchor_id,distance_band LIMIT {args.max_rows + 1};"
    rows = sql(args, query, package / "cofactors.duckdb", cwd=package, rows=True)
    if len(rows) > args.max_rows:
        raise ValueError("region exceeds --max-rows; narrow the span (nothing exported)")
    if not rows:
        raise ValueError("no included TP73 anchor overlaps this region; not a genome-wide negative")
    ids = ",".join(str(n) for n in sorted({int(r["anchor_id"]) for r in rows}))
    membership = f"SELECT * FROM anchor_promoter WHERE anchor_id IN ({ids})"
    feature_ids = f"SELECT regulatory_feature_id FROM ({membership})"
    extra = {}
    for name, query in (("promoter_memberships", membership),
                        ("promoters", f"SELECT * FROM promoter WHERE regulatory_feature_id IN ({feature_ids})"),
                        ("promoter_gene_links", f"SELECT * FROM promoter_gene WHERE regulatory_feature_id IN ({feature_ids})")):
        extra[name] = sql(args, query + f" LIMIT {args.max_rows+1};", package / "cofactors.duckdb", cwd=package, rows=True)
    if len(rows) + sum(len(v) for v in extra.values()) > args.max_rows:
        raise ValueError("region plus promoter annotation exceeds --max-rows; narrow the span")
    with transaction(args.output.parent.resolve(), args.output.name) as (out, exists):
        if exists:
            raise ValueError("region destination already exists; choose a new name")
        loci = {}
        for r in rows:
            a = (r["anchor_start"], r["anchor_end"], ".")
            loci[("TP73", *a)] = f"TP73_{a[0]}_{a[1]}"
            if r["hit_start"] is not None:
                h = (r["hit_start"], r["hit_end"], r["best_strand"])
                loci[(args.motif, *h)] = f"{args.motif}_{h[0]}_{h[1]}"
        for promoter in extra["promoters"]:
            name = promoter["regulatory_feature_id"]
            loci[(name, promoter["extended_start"], promoter["extended_end"], ".")] = name
        with (out / "regions.bed").open("x") as handle:
            for (_, start, end, strand), name in sorted(loci.items(), key=lambda x: (x[0][1],x[0][2],x[0][0])):
                handle.write(f"{args.chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
        write_json(out / "evidence.json", {"package_manifest_sha256": digest(package / "manifest.json"),
                  "assembly": manifest["assembly"], "coordinate_mode": manifest["coordinate_mode"],
                  "motif_id": args.motif, "rows": rows, **extra})
        set_id = f"tp73_cofactors_{args.chrom}_{args.start}_{args.end}_{args.motif}"
        workflow = {"run_id": set_id, "ops": [{"ImportGenomicRegionSet": {"request": {
            "path": "regions.bed", "format": "bed", "bare_bed_reference": {
                "species_scientific_name": "Homo sapiens", "taxon_id": 9606,
                "assembly_name": "GRCh38", "contig_name": args.chrom},
            "set_id_override": set_id, "collision_policy": "reject",
            "max_bytes": (out / "regions.bed").stat().st_size, "max_rows": len(loci)}}},
            {"ExportGenomicRegionSet": {"request": {"set_id": set_id, "json_path": "regions.gentle.json"}}}]}
        write_json(out / "gentle_workflow.json", workflow)
        print(f"Exported {len(loci)} spans and {len(rows)} anchor-band rows; run GENtle workflow from {args.output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duckdb", default="duckdb", help="DuckDB CLI executable")
    parser.add_argument("--memory-limit", default="2GB", help="DuckDB memory ceiling")
    parser.add_argument("--scratch-root", type=Path, help="temporary/spill parent; node-local scratch for Slurm")
    commands = parser.add_subparsers(dest="command", required=True)
    prep = commands.add_parser("prepare", help="require finalized promoter results, pin a bounded candidate panel")
    prep.add_argument("--analysis-run", required=True, type=Path)
    prep.add_argument("--regulatory-features", required=True, type=Path)
    prep.add_argument("--panel-size", type=int, default=100)
    prep.add_argument("--max-bytes", type=int, default=2_000_000_000, help="hard publication cap including the query index")
    task = commands.add_parser("motif", help="checkpoint one panel motif, one chromosome at a time")
    task.add_argument("--task-index", type=int, required=True)
    end = commands.add_parser("finalize", help="validate all checkpoints, publish a relocatable size-capped package")
    for p in (prep, task, end):
        p.add_argument("--run-root", type=Path, required=True)
    view = commands.add_parser("inspect", help="validate and describe a portable package")
    region = commands.add_parser("region", help="bounded strongest-locus BED/evidence/GENtle workflow export")
    for p in (view, region):
        p.add_argument("--package", type=Path, required=True)
    region.add_argument("--chrom", required=True)
    region.add_argument("--start", required=True, type=int, help="0-based inclusive")
    region.add_argument("--end", required=True, type=int, help="0-based exclusive")
    region.add_argument("--motif", required=True)
    region.add_argument("--max-rows", type=int, default=2000, help="total anchor-band and annotation row cap; exceeding it fails")
    region.add_argument("--output", type=Path, required=True, help="new directory")
    args = parser.parse_args()
    try:
        with tempfile.TemporaryDirectory(prefix="promoter-export-", dir=args.scratch_root) as temp:
            args.temp = Path(temp)
            if hasattr(args, "run_root"):
                args.run_root = args.run_root.resolve()
            {"prepare": prepare, "motif": motif_task, "finalize": finalize,
             "inspect": inspect_or_region, "region": inspect_or_region}[args.command](args)
    except (ValueError, KeyError, OSError, subprocess.SubprocessError) as error:
        print(f"E: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

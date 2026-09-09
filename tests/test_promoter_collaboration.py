#!/usr/bin/env python3
"""Tiny offline end-to-end collaborator package and GENtle handoff checks."""

import argparse
import csv
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/export_promoter_collaboration.py"
spec = importlib.util.spec_from_file_location("collaboration", SCRIPT)
export = importlib.util.module_from_spec(spec)
spec.loader.exec_module(export)


class CollaborationTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="promoter-collaboration-test-")
        self.root = Path(self.temp.name)
        self.args = argparse.Namespace(duckdb="duckdb", memory_limit="512MB", temp=self.root)

    def tearDown(self):
        self.temp.cleanup()

    def sql(self, query, rows=False, database=":memory:", cwd=None):
        return export.sql(self.args, query, rows=rows, database=database, cwd=cwd)

    def pq(self, path, query):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.sql(export.copy_sql(query, path))
        return path

    def cli(self, *args, ok=True):
        result = subprocess.run([sys.executable, str(SCRIPT), "--memory-limit", "512MB", *map(str,args)],
                                text=True, capture_output=True)
        if ok:
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout)
        return result

    def manifest(self, path, files, **extra):
        export.write_json(path, {"files": [export.record(f, path.parent) for f in files], **extra})

    def fixture(self):
        run = self.root / "upstream"
        plan = run / "plan"
        final = run / "final/distance_enrichment"
        split = run / "input/anchors"
        for directory in (plan, final, split):
            directory.mkdir(parents=True)
        anchors = []
        inventory = []
        for chrom, positions in (("1", (100,500)), ("2", (100,))):
            path = self.pq(split / ("chrom-"+chrom+".parquet"), f"""
SELECT '{chrom}'::VARCHAR AS chrom, position::BIGINT AS anchor_start,
 (position+10)::BIGINT AS anchor_end, -2.0::DOUBLE AS anchor_score,
 position=100 AS supported_tp73_saos2_TA, false AS supported_negative_control_saos2_TA
FROM (VALUES {','.join('('+str(n)+')' for n in positions)}) t(position)""")
            anchors.append(path)
            inventory.append({**export.record(path, split), "chrom": chrom,
                              "relative_path": path.name, "rows": len(positions)})
        with (split / "anchor_files.tsv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(inventory[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(inventory)
        evidence = self.pq(self.root / "source_evidence.parquet", f"SELECT *,CASE WHEN supported_tp73_saos2_TA THEN 4 ELSE 0 END::DOUBLE AS depth_tp73_saos2_TA,0::DOUBLE AS depth_negative_control_saos2_TA FROM {export.parquet(anchors)}")
        regulatory = self.root / "regulatory"
        regulatory.mkdir()
        chrom_manifests = []
        feature_rows = []
        for chrom in ("1", "2"):
            part = regulatory / ("chrom-" + chrom)
            part.mkdir()
            # The same anchor belongs to two promoters and multiple genes; no row multiplication.
            links = self.pq(part / "tp73_anchor_regulatory_feature.parquet", f"""
SELECT a.*, 'p'||'{chrom}'||k::VARCHAR AS regulatory_feature_id,
 'promoter_extended'::VARCHAR AS regulatory_definition_id,10::BIGINT AS overlap_bp,
 true AS anchor_fully_within FROM {export.parquet(anchors[int(chrom)-1])} a,range(2) v(k)""")
            owners = self.pq(part / "regulatory_feature_gene.parquet", f"""
SELECT 'p'||'{chrom}'||k::VARCHAR AS regulatory_feature_id,'gene'||k AS gene_id,
 'ensembl_native' AS link_source,'113' AS annotation_release FROM range(2) v(k)""")
            m = part / "manifest.json"
            self.manifest(m, [links, owners])
            chrom_manifests.append({**export.record(m, self.root), "path": str(m)})
            for k in range(2):
                feature_rows.append(f"('p{chrom}{k}','{chrom}',50,600,50,600)")
        features = self.pq(regulatory / "features/regulatory_feature.parquet", "SELECT * FROM (VALUES " + ",".join(feature_rows) + ") t(regulatory_feature_id,chrom,start,\"end\",extended_start,extended_end)")
        self.manifest(features.parent / "manifest.json", [features], assembly="GRCh38", source_sha256="synthetic")
        self.manifest(regulatory / "manifest.json", [], state="complete",
                      coordinate_audit={"status":"passed"}, chromosome_manifests=chrom_manifests)
        selection = {"subset":"promoter_extended", "manifest":str(regulatory / "manifest.json"),
                     "manifest_sha256":export.digest(regulatory / "manifest.json"),
                     "coordinate_audit":{"gff_sha256":"synthetic"}}
        motifs = (("MA0653.1","IRF9"),("MA1961.2","PATZ1"),("MA0001.1","TEST"))
        hit_rows = []
        for task_index, (motif, name) in enumerate(motifs):
            for chrom in ("1","2"):
                for strand in ("+","-"):
                    if chrom == "1":
                        cells = [(95,.2),(90,-.5),(110,2),(115,1),(116,3),(130,-.5),
                                 (131,4),(160,5),(161,6),(210,7),(211,8),(260,9),(261,99)]
                        if strand == "-":
                            cells = [(110,2)]
                    else:
                        cells = [(110,1)] if strand == "+" else []
                    values = ",".join(f"({start},{start+10},{score})" for start,score in cells)
                    query = f"SELECT s::BIGINT AS start,e::BIGINT AS \"end\",v::DOUBLE AS score FROM (VALUES {values}) t(s,e,v)" if cells else "SELECT 0::BIGINT AS start,10::BIGINT AS \"end\",0.0::DOUBLE AS score WHERE false"
                    hit_root = self.root / "hits/genome_id=homo_sapiens_grch38_ensembl113_primary/motif_set_id=jaspar2026_core_nonredundant/score_mode=log2_relative_risk/pseudocount=1/background_model_id=uniform_acgt_v1/pseudocount_scheme=additive_per_base/n_policy=skip"
                    path = self.pq(hit_root / f"{motif}_{chrom}_{strand}.parquet", query)
                    hit_rows.append({"motif_id":motif,"chrom":chrom,"strand":strand,
                                     "absolute_path":str(path),"minimum_score":-1,
                                     "bytes":path.stat().st_size,"sha256":export.digest(path)})
        with (plan / "scan_files.tsv").open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(hit_rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(hit_rows)
        config = {"run_id":"tiny", "regulatory_selection":selection,"source_score_floor":-1,
                  "context_flank_bp":150,"distance_bands":list(export.BANDS),
                  "excluded_series":["skmel29_1"],"tax_group":"vertebrates","chromosomes":["1","2"],
                  "anchor_evidence":str(evidence),"anchor_evidence_sha256":export.digest(evidence),"scan_files_sha256":export.digest(plan / "scan_files.tsv")}
        export.write_json(plan / "run_config.json", config)
        export.write_json(split / "complete.json", {"regulatory_selection":selection,
                          "source_sha256":export.digest(evidence),"inventory_sha256":export.digest(split / "anchor_files.tsv"),"anchors":3})
        values = []
        for index,(motif,name) in enumerate(motifs):
            for order,band in enumerate(export.BANDS):
                present = 2 if order == 1 else 1
                values.append(f"('{motif}','{name}','{band}',{order},-1.0,0.0,'vertebrates',3,{present},{present},{index+1},{index+1},{index+1},{index+1},{index+1},1.2,1.3,0.92,{present}/3.0,(3-{present})/3.0)")
        names = ("motif_id,motif_name,distance_band,distance_band_order,source_score_floor,positive_threshold,tax_group,"
                 "anchors_total,anchors_source_present,anchors_positive,ta_enrichment_rank,dn_enrichment_rank,ta_depletion_rank,dn_depletion_rank,isoform_difference_absolute_rank,"
                 "ta_adjusted_odds_ratio,dn_adjusted_odds_ratio,ta_vs_dn_odds_ratio_ratio,positive_anchor_fraction,negative_anchor_fraction")
        tables = {}
        for name in export.TABLES:
            path = self.pq(final / (name+".parquet"), "SELECT * FROM (VALUES " + ",".join(values) + f") t({names})")
            tables[name] = export.record(path, final)
        export.write_json(final / "manifest.json", {"schema_version":6,"run_id":"tiny","regulatory_selection":selection,"tables":tables,"task_count":3})
        return run, features.parent

    def test_complete_portable_package(self):
        upstream, features = self.fixture()
        run = self.root / "export"
        prep = ("prepare","--run-root",run,"--analysis-run",upstream,"--regulatory-features",features,"--panel-size","3")
        self.cli(*prep)
        self.cli(*prep)  # Identical rerun changes no files.
        plan = export.read_json(run / "prepared/plan.json")
        self.assertEqual({m["motif_id"] for m in plan["panel"]},{"MA0653.1","MA1961.2","MA0001.1"})
        self.assertEqual(next(r["status"] for r in plan["requested_genes"] if r["gene"]=="SUB1"),"no_exact_name_in_result")
        self.cli("finalize","--run-root",run,ok=False)
        self.assertFalse((run / "package").exists())
        for i in range(3):
            self.cli("motif","--run-root",run,"--task-index",i)
        before = export.digest(run / "tasks/MA0653.1/chrom-1/feature.parquet")
        self.cli("motif","--run-root",run,"--task-index",0)
        self.assertEqual(before,export.digest(run / "tasks/MA0653.1/chrom-1/feature.parquet"))
        self.cli("finalize","--run-root",run)
        self.cli("finalize","--run-root",run)
        package = self.root / "moved package"
        shutil.move(run / "package",package)
        self.cli("inspect","--package",package)
        db = package / "cofactors.duckdb"
        rows = self.sql("SELECT * FROM anchor_features('MA0653.1','1',100,501) ORDER BY anchor_start,distance_band",True,db,package)
        self.assertEqual(len(rows),12)
        bands = {r["distance_band"]:r for r in rows if r["anchor_start"]==100}
        self.assertEqual(bands["adjacent_0_5"]["best_strand"],".")
        self.assertEqual(bands["adjacent_0_5"]["n_source_loci"],3)
        self.assertEqual(bands["adjacent_0_5"]["n_score_zero_loci"],2)
        self.assertEqual(bands["gap_101_150"]["best_score"],9)
        self.assertEqual(bands["gap_101_150"]["interval_distance_bp"],150)
        self.assertEqual(bands["gap_101_150"]["hit_start"],260)
        self.assertEqual(bands["overlap"]["interval_distance_bp"],-5)
        self.assertEqual(bands["overlap"]["depth_tp73_saos2_TA"],4)
        self.assertTrue(all(r["n_source_loci"]==0 and r["best_score"] is None for r in rows if r["anchor_start"]==500))
        second = self.sql("SELECT * FROM anchor_features('MA0653.1','2',100,110) WHERE distance_band='adjacent_0_5'",True,db,package)[0]
        self.assertEqual(second["best_score"],1)
        self.assertEqual(self.sql("SELECT count(*) n FROM anchors",True,db,package)[0]["n"],3)
        self.assertEqual(self.sql("SELECT count(*) n FROM anchor_promoter",True,db,package)[0]["n"],6)
        with self.assertRaises(ValueError):
            self.sql("SELECT * FROM anchor_features('MISSING','1',100,110)",True,db,package)
        out = self.root / "region"
        region = ("region","--package",package,"--motif","MA0653.1","--chrom","1","--start","100","--end","110","--output",out)
        self.cli(*region,"--max-rows","5",ok=False)
        self.assertFalse(out.exists())
        self.cli(*region)
        bed = (out / "regions.bed").read_text()
        self.assertIn("1\t110\t120\tMA0653.1_110_120\t0\t.\n",bed)
        self.assertIn("1\t100\t110\tTP73_100_110\t0\t.\n",bed)
        gentle = ROOT.parent / "gentle_rs/target/debug/gentle_cli"
        if gentle.exists():
            result = subprocess.run([str(gentle),"workflow","gentle_workflow.json"],cwd=out,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            canonical = json.loads((out / "regions.gentle.json").read_text())
            self.assertIn("GRCh38", json.dumps(canonical))
        # Source changes invalidate a completed checkpoint instead of silently reusing it.
        altered = dict(plan["identity"], max_bytes=1)
        with self.assertRaises(ValueError):
            export.validate_product(run / "prepared",altered)

    def test_floor_and_byte_cap(self):
        upstream,features = self.fixture()
        run = self.root / "small"
        self.cli("prepare","--run-root",run,"--analysis-run",upstream,"--regulatory-features",features,"--panel-size","2","--max-bytes","1")
        for i in range(2):
            self.cli("motif","--run-root",run,"--task-index",i)
        self.cli("finalize","--run-root",run,ok=False)
        self.assertFalse((run / "package").exists())
        path = upstream / "plan/scan_files.tsv"
        text = path.read_text().replace("\t-1\t","\t0\t")
        path.write_text(text)
        config_path = upstream / "plan/run_config.json"
        config = json.loads(config_path.read_text())
        config["scan_files_sha256"] = export.digest(path)
        config_path.write_text(json.dumps(config))
        failure = self.cli("prepare","--run-root",self.root / "censored","--analysis-run",upstream,"--regulatory-features",features,"--panel-size","2",ok=False)
        self.assertIn("censored",failure.stderr)

    def test_submission_dry_run(self):
        commands = self.root / "commands"
        commands.mkdir()
        fixtures = {
            "git": '#!/bin/sh\ncase "$*" in *rev-parse*) printf "%s\\n" 0123456789012345678901234567890123456789;; esac\n',
            "realpath": '#!/bin/sh\nprintf "%s\\n" "$2"\n',
            "sbatch": '#!/bin/sh\necho "sbatch must not execute in dry run" >&2; exit 99\n',
        }
        for name, content in fixtures.items():
            path = commands / name
            path.write_text(content)
            path.chmod(0o755)
        spooled = self.root / "slurm_script"
        shutil.copyfile(ROOT / "scripts/submit_promoter_collaboration_slurm.sh", spooled)
        result = subprocess.run(["bash",str(spooled),"--source",str(ROOT),
                                 "--run-root","/data/sm718/promoter_collaboration_test_not_created",
                                 "--analysis-run",str(self.root),"--regulatory-features",str(self.root),
                                 "--duckdb",shutil.which("duckdb"),"--afterok","12345","--dry-run"],
                                 text=True,capture_output=True,env={**os.environ,"PATH":str(commands)+os.pathsep+os.environ["PATH"]})
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertEqual(len(result.stdout.splitlines()),3)
        self.assertIn("afterok:12345",result.stdout)
        self.assertIn("--array=0-99%4",result.stdout)
        self.assertEqual(result.stdout.count("--partition=requeue"),3)
        self.assertIn("--mem=32G",result.stdout)
        self.assertEqual(result.stdout.count("--source " + str(ROOT)),3)


if __name__ == "__main__":
    unittest.main()

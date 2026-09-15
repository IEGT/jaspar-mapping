#!/usr/bin/env python3
"""Synthetic regulatory-hit export, independent of TP73 and external downloads.

All inputs are hand-crafted below: two chromosomes, two artificial matrices,
overlapping promoters and plus/minus TSS windows. Run this file to recreate
them in temporary storage and exercise the real DuckDB/Parquet exporter.
"""

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import export_regulatory_tfbs as exporter


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/export_regulatory_tfbs.py"
DUCKDB = os.environ.get("DUCKDB", "duckdb")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def record(path):
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest(path)}


@unittest.skipUnless(shutil.which(DUCKDB), "DuckDB CLI required")
class RegulatoryTFBS(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="regulatory-tfbs-")
        self.root = Path(self.temp.name)
        self.package = self.root / "scan with spaces"
        self.package.mkdir()
        self.reg = self.root / "regulatory"
        (self.reg / "features").mkdir(parents=True)
        (self.reg / "final").mkdir()
        self.database = self.package / "scan.duckdb"
        self.write_json(self.package / "manifest.json", {
            "state": "complete", "schema_version": 1, "database": "scan.duckdb",
            "genome_id": "synthetic", "motif_set_id": "test", "run_id": "synthetic_scan",
        })
        inventory = []
        statements = []
        self.payloads = []
        for chrom in ("1", "2"):
            for motif in ("MA9000.1", "MA9001.1"):
                for strand in ("+", "-"):
                    task = len(inventory)
                    relative = (f"genome_id=synthetic/motif_set_id=test/chrom={chrom}/motif_id={motif}/"
                                "score_mode=log2_relative_risk/pseudocount=1/background_model_id=uniform_acgt_v1/"
                                "pseudocount_scheme=additive_per_base/minimum_score=-1/n_policy=skip/"
                                f"strand={'plus' if strand == '+' else 'minus'}/part.parquet")
                    path = self.package / "task_data" / f"task_id={task}" / relative
                    path.parent.mkdir(parents=True)
                    # Both orientations stay separate, including overlapping
                    # hits and a negative fractional score at a promoter edge.
                    statements.append(f"""COPY (SELECT s::BIGINT AS start,e::BIGINT AS "end",
                        score::FLOAT AS score,0.5::FLOAT AS pwm_relative_score FROM
                        (VALUES (80,90,2.5),(89,91,-0.5),(100,110,4.25),(120,126,1.25),
                                (149,150,0.0),(150,151,2.0),(205,210,2.0),(255,260,2.0),
                                (305,310,2.0),(355,360,2.0),(405,410,2.0),(495,501,3.0),
                                (550,560,2.0),(605,610,2.0)) r(s,e,score))
                        TO '{path}' (FORMAT PARQUET);""")
                    self.payloads.append(path)
                    inventory.append({"task_id": task, "output_relative_path": relative,
                                      "motif_id": motif, "motif_set_id": "test", "genome_id": "synthetic",
                                      "chrom": chrom, "strand": strand, "state": "complete",
                                      "coordinate_mode": "bed", "minimum_score": -1.0,
                                      "minimum_pwm_relative_score": None, "maximum_pwm_relative_score": None})
        self.sql("".join(statements))
        for row, path in zip(inventory, self.payloads):
            row.update(bytes=path.stat().st_size, sha256=digest(path))
        self.write_json(self.root / "inventory.json", inventory)
        self.sql(f"""
CREATE TABLE genome AS SELECT 'synthetic' AS genome_id,'SYN1' AS assembly_name;
CREATE TABLE motif_metadata AS SELECT 'test' AS motif_set_id,* FROM
 (VALUES ('MA9000.1','TEST_A'),('MA9001.1','TEST_B')) t(motif_id,motif_name);
CREATE TABLE scan_file_inventory AS SELECT * FROM read_json_auto('{self.root}/inventory.json');
""", self.database)
        self.features = self.reg / "features/regulatory_feature.parquet"
        self.sql(f"""COPY (SELECT 'SYN1' AS assembly,'1' AS chrom,*,
            CASE WHEN feature_type='promoter' THEN start END AS core_start,
            CASE WHEN feature_type='promoter' THEN "end" END AS core_end
            FROM (VALUES ('promoter',100::BIGINT,120::BIGINT,90::BIGINT,150::BIGINT),
                         ('promoter',110,125,100,150),('enhancer',200,220,NULL,NULL),
                         ('open_chromatin_region',250,270,NULL,NULL),
                         ('CTCF_binding_site',300,320,NULL,NULL),('EMAR',350,370,NULL,NULL),
                         ('future_feature',600,620,NULL,NULL))
            t(feature_type,start,"end",extended_start,extended_end))
            TO '{self.features}' (FORMAT PARQUET);""")
        self.promoters = self.reg / "promoters.parquet"
        self.sql(f"""COPY (SELECT 'synthetic' AS genome_id,'test_gtf' AS annotation_release,
            'tss_u20_d5_v1' AS promoter_definition_id,'1' AS chrom,* FROM
            (VALUES (380::BIGINT,406::BIGINT,'+','shared1'),
                    (380,406,'+','shared2'),(495,521,'-','reverse'))
            t(promoter_start,promoter_end,strand,promoter_id))
            TO '{self.promoters}' (FORMAT PARQUET);""")
        self.audit = {"kind": "regulatory_coordinate_audit", "status": "passed",
                      "gff_sha256": "a" * 64, "assembly": "SYN1",
                      "note": "Synthetic expected endpoints, not a biological audit"}
        self.plan = {"genome_id": "synthetic", "annotation_release": "test_gtf",
                     "promoter_definition_id": "tss_u20_d5_v1", "audit": self.audit,
                     "tasks": [{"chrom": "1", "promoters": record(self.promoters)}]}
        self.write_json(self.reg / "plan.json", self.plan)
        self.write_json(self.reg / "features/manifest.json", {
            "kind": "ensembl_regulatory_annotation", "assembly": "SYN1",
            "coordinate_mode": "bed_0based_half_open", "chromosomes": ["1"],
            "source_sha256": "a" * 64, "regulatory_release": "synthetic",
            "coordinate_audit": self.audit,
            "files": [dict(record(self.features), path="regulatory_feature.parquet")],
        })
        self.write_json(self.reg / "final/manifest.json", {
            "kind": "genome_regulatory_membership", "state": "complete",
            "assembly": "SYN1", "annotation_release": "test_gtf",
            "promoter_definition_id": "tss_u20_d5_v1", "chromosomes": ["1"],
            "production_plan_sha256": digest(self.reg / "plan.json"), "coordinate_audit": self.audit,
        })

    def tearDown(self):
        self.temp.cleanup()

    def write_json(self, path, obj):
        path.write_text(json.dumps(obj))

    def sql(self, query, database=":memory:"):
        result = subprocess.run([DUCKDB, "-no-init", "-batch", "-bail", "-json", str(database)],
                                input=query, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout or "[]")

    def call(self, name, *extra, success=True, all_motifs=False):
        result = subprocess.run([sys.executable, str(SCRIPT), "--package", str(self.package),
                                 "--regulatory-run", str(self.reg), "--chrom", "1",
                                 *( ["--all-motifs"] if all_motifs else ["--motif", "MA9000.1"] ),
                                 "--start", "0", "--end", "1000", "--minimum-free-bytes", "0",
                                 "--max-output-bytes", "1000000", "--output", str(self.root / name),
                                 "--duckdb", DUCKDB, *map(str, extra)], text=True, capture_output=True)
        self.assertEqual(result.returncode == 0, success, result.stderr)
        return result

    def rows(self, name, filename="motif_hits.parquet"):
        return self.sql(f"SELECT * FROM read_parquet('{self.root}/{name}/{filename}') ORDER BY 1,2,3;")

    def test_geometry_tags_no_tp73_and_no_transcript_multiplication(self):
        hashes = [digest(p) for p in self.payloads]
        self.call("hits")
        rows = self.rows("hits")
        self.assertEqual(len(rows), 22)
        self.assertEqual({r["strand"] for r in rows}, {"+", "-"})
        forward = {r["start"]: r for r in rows if r["strand"] == "+"}
        self.assertNotIn(80, forward)  # abutment is not positive overlap
        self.assertNotIn(150, forward)
        self.assertNotIn(550, forward)
        self.assertEqual(forward[89]["score"], -0.5)
        self.assertFalse(forward[89]["overlaps_promoter_core"])
        self.assertTrue(forward[100]["overlaps_promoter_core"])
        self.assertTrue(forward[100]["overlaps_promoter_extended"])
        self.assertTrue(forward[205]["overlaps_enhancer"])
        self.assertTrue(forward[255]["overlaps_open_chromatin"])
        self.assertTrue(forward[305]["overlaps_ctcf"])
        self.assertTrue(forward[355]["overlaps_emar"])
        self.assertTrue(forward[605]["overlaps_other_regulatory"])
        self.assertTrue(forward[405]["overlaps_tss_window"])
        self.assertTrue(forward[495]["overlaps_tss_window"])
        self.assertEqual([digest(p) for p in self.payloads], hashes)
        manifest = json.loads((self.root / "hits/manifest.json").read_text())
        self.assertFalse(manifest["requires_tp73"])
        self.assertFalse(manifest["complete_genome_scan"])
        self.assertTrue(manifest["complete_for_requested_scope_at_threshold"])

    def test_scopes_scores_and_all_motifs(self):
        self.call("promoters", "--scope", "promoter_or_tss", "--minimum-score", "0", all_motifs=True)
        rows = self.rows("promoters")
        self.assertEqual({r["motif_id"] for r in rows}, {"MA9000.1", "MA9001.1"})
        self.assertEqual({r["start"] for r in rows}, {100, 120, 149, 405, 495})
        self.assertEqual(len(rows), 20)
        self.call("counts", "--count-only")
        counts = self.rows("counts", "counts.parquet")
        self.assertEqual(sum(r["orientation_records"] for r in counts), 22)
        self.assertEqual(sum(r["physical_loci"] for r in counts), 11)

    def test_unavailable_and_censored_sources_fail(self):
        self.assertIn("coverage", self.call("badchr", "--chrom", "2", success=False).stderr)
        self.assertIn("absent", self.call("missing", "--motif", "MA9002.1", success=False).stderr)
        self.sql("UPDATE scan_file_inventory SET minimum_score=6;", self.database)
        self.assertIn("floor", self.call("censored", success=False).stderr)

    def test_caps_existing_output_and_empty_valid_subset(self):
        self.call("small", "--max-rows", "1", success=False)
        self.assertFalse((self.root / "small").exists())
        self.call("bytes", "--max-output-bytes", "1", success=False)
        self.assertFalse((self.root / "bytes").exists())
        self.call("empty", "--start", "700", "--end", "710")
        self.assertEqual(self.rows("empty"), [])
        self.assertIn("replace", self.call("empty", success=False).stderr)

    def test_resume_and_scratch_match_direct_export(self):
        self.call("direct")
        self.call("staged", "--scratch-directory", self.root / "scratch")
        self.assertEqual(self.rows("direct"), self.rows("staged"))
        before = digest(self.root / "staged/motif_hits.parquet")
        result = self.call("staged", "--resume")
        self.assertTrue(json.loads(result.stdout)["reused"])
        self.assertEqual(digest(self.root / "staged/motif_hits.parquet"), before)
        self.assertIn("scope", self.call("staged", "--resume", "--minimum-score", "0", success=False).stderr)
        self.call("retry", "--max-rows", "1", success=False)
        self.call("retry", "--resume")
        self.assertEqual(self.rows("retry"), self.rows("direct"))

    def test_corruption_and_unaudited_annotation_fail(self):
        with self.features.open("ab") as stream:
            stream.write(b"changed")
        self.assertIn("changed", self.call("corrupt", success=False).stderr)

    def test_duplicate_orientation_and_assembly_fail(self):
        self.sql("INSERT INTO scan_file_inventory SELECT * FROM scan_file_inventory LIMIT 1;", self.database)
        self.assertIn("duplicate", self.call("duplicate", success=False).stderr)
        final_path = self.reg / "final/manifest.json"
        final = json.loads(final_path.read_text())
        final["assembly"] = "wrong"
        self.write_json(final_path, final)
        self.assertIn("identity", self.call("assembly", success=False).stderr)

    def test_running_writer_is_stopped_at_byte_budget(self):
        watched = self.root / "growing.parquet"
        runtime = self.root / "synthetic-writer"
        runtime.write_text(f"#!{sys.executable}\nimport sys,time\nfrom pathlib import Path\n"
                           f"sys.stdin.read()\nPath({str(watched)!r}).write_bytes(b'x'*16384)\n"
                           "time.sleep(30)\n")
        runtime.chmod(0o700)
        args = SimpleNamespace(duckdb=str(runtime), threads=1, memory_limit="256MB",
                               timeout_seconds=10, max_output_bytes=1024,
                               minimum_free_bytes=0)
        before = time.monotonic()
        with self.assertRaisesRegex(exporter.QueryError, "byte budget"):
            exporter.sql(args, "SELECT 1;", output_watch=watched)
        self.assertLess(time.monotonic() - before, 5)


if __name__ == "__main__":
    unittest.main()

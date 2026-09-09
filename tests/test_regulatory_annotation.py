#!/usr/bin/env python3
"""Small, inspectable regulatory geometry and evidence-subset contracts."""

import gzip
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
SCRIPT = ROOT / "scripts/build_regulatory_annotation.py"
DUCKDB = os.environ.get("DUCKDB", "duckdb")
spec = importlib.util.spec_from_file_location("regulatory", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@unittest.skipUnless(shutil.which(DUCKDB), "DuckDB CLI is required")
class RegulatoryContract(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="jaspar-regulatory-test-")
        self.root = Path(self.temporary.name)
        self.gff = self.root / "features.gff.gz"
        with gzip.open(self.gff, "wt") as stream:
            stream.write(
                "1\tEnsembl\tpromoter\t101\t120\t.\t.\t.\tID=P1;extended_start=91;extended_end=150;gene_id=A,B\n"
                "1\tEnsembl\tpromoter\t111\t125\t.\t.\t.\tID=P2;extended_start=101;extended_end=150;gene_id=A\n"
                "1\tEnsembl\tenhancer\t201\t220\t.\t.\t.\tID=E1\n"
                "1\tEnsembl\tCTCF_binding_site\t251\t270\t.\t+\t.\tID=C1\n"
                "1\tEnsembl\topen_chromatin_region\t401\t420\t.\t.\t.\tID=O1\n"
                "2\tEnsembl\tenhancer\t101\t120\t.\t.\t.\tID=E2\n"
            )
        self.call("import", "--gff", self.gff, "--assembly", "GRCh38",
                  "--requested-release", "2025-12", "--resolved-release", "2025-05",
                  "--source-uri", "https://example.org/fixture.gff.gz",
                  "--coordinate-audit-note", "synthetic known endpoints",
                  "--output", self.root / "features")
        self.sql(f"""
COPY (SELECT '1' AS chrom, s::BIGINT AS anchor_start, e::BIGINT AS anchor_end,
             7.5::DOUBLE AS anchor_score, true AS supported_tp73_saos2_TA
      FROM (VALUES (80,90),(89,91),(99,100),(100,101),(120,121),(149,150),
                   (150,151),(205,206),(260,261),(300,310),(405,406)) r(s,e))
TO '{self.root}/anchors.parquet' (FORMAT PARQUET);
CREATE TABLE tss AS SELECT 'human' AS genome_id, 'gtf-test' AS annotation_release,
  '1' AS chrom, s::BIGINT AS start, (s+1)::BIGINT AS "end", strand, tss_id
FROM (VALUES (105,'+','T1'),(105,'-','T2'),(145,'+','T3')) r(s,strand,tss_id);
COPY tss TO '{self.root}/tss.parquet' (FORMAT PARQUET);
COPY (SELECT 'human' AS genome_id,'gtf-test' AS annotation_release,*
      FROM (VALUES ('T1','A','TA'),('T1','B','TB'),('T2','C','TC'),('T3','D','TD'))
        r(tss_id,gene_id,transcript_id))
TO '{self.root}/owners.parquet' (FORMAT PARQUET);
COPY (SELECT 'human' AS genome_id,'gtf-test' AS annotation_release,'legacy' AS promoter_definition_id,
      '1' AS chrom, 80::BIGINT AS promoter_start, 121::BIGINT AS promoter_end)
TO '{self.root}/promoters.parquet' (FORMAT PARQUET);
""")
        self.annotate("membership")

    def tearDown(self):
        self.temporary.cleanup()

    def call(self, *args, success=True):
        result = subprocess.run([sys.executable, str(SCRIPT), *map(str, args),
                                 "--duckdb", DUCKDB], text=True, capture_output=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def sql(self, sql):
        result = subprocess.run([DUCKDB, "-batch", "-bail", "-json", ":memory:"],
                                input=sql, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        return json.loads(result.stdout) if result.stdout.strip() else []

    def annotate(self, name, chrom="1", success=True):
        return self.call("annotate", "--features", self.root / "features",
                         "--anchors", self.root / "anchors.parquet", "--chrom", chrom,
                         "--assembly", "GRCh38", "--tss", self.root / "tss.parquet",
                         "--transcript-tss", self.root / "owners.parquet",
                         "--annotation-release", "gtf-test", "--genome-id", "human",
                         "--promoters", self.root / "promoters.parquet",
                         "--promoter-definition-id", "legacy", "--output", self.root / name,
                         success=success)

    def test_coordinates_and_identity(self):
        rows = self.sql(f"SELECT * FROM '{self.root}/features/regulatory_feature.parquet' "
                        "WHERE source_feature_id='P1'")
        self.assertEqual(len(rows), 1)
        r = rows[0]
        self.assertEqual((r["core_start"], r["core_end"]), (100, 120))
        self.assertEqual((r["extended_start"], r["extended_end"]), (90, 150))
        self.assertEqual(r["requested_regulatory_release"], "2025-12")
        self.assertEqual(r["regulatory_release"], "2025-05")
        self.assertEqual(r["strand"], ".")
        self.assertEqual(len(r["source_sha256"]), 64)
        rows = self.sql(f"SELECT extended_start,core_start FROM "
                        f"'{self.root}/features/regulatory_feature.parquet' WHERE source_feature_id='E1'")
        self.assertIsNone(rows[0]["extended_start"])
        self.assertIsNone(rows[0]["core_start"])

    def test_membership_without_multiplication(self):
        rows = self.sql(f"SELECT * FROM '{self.root}/membership/tp73_anchor_regulatory_membership.parquet'")
        by_start = {r["anchor_start"]: r for r in rows}
        self.assertEqual(len(rows), 11)
        for start in (80, 150, 300):
            self.assertFalse(by_start[start]["promoter_extended"])
        self.assertTrue(by_start[89]["promoter_extended"])
        self.assertFalse(by_start[99]["promoter_core"])
        self.assertTrue(by_start[100]["promoter_core"])
        self.assertTrue(by_start[120]["promoter_core"])
        self.assertTrue(by_start[149]["promoter_extended"])
        self.assertTrue(by_start[205]["enhancer"])
        self.assertFalse(by_start[260]["no_regulatory_overlap"])
        self.assertTrue(by_start[300]["no_regulatory_overlap"])
        self.assertTrue(by_start[405]["open_chromatin"])
        self.assertTrue(by_start[80]["tss_window_promoter"])
        rows = self.sql(f"SELECT overlap_bp,anchor_fully_within FROM "
                        f"'{self.root}/membership/tp73_anchor_regulatory_feature.parquet' "
                        "WHERE anchor_start=89 AND regulatory_definition_id='promoter_extended'")
        self.assertEqual(rows, [{"overlap_bp": 1, "anchor_fully_within": False}])

    def test_bidirectional_many_gene_links_and_unlinked_feature(self):
        rows = self.sql(f"SELECT g.gene_id, g.link_source FROM "
                        f"'{self.root}/membership/regulatory_feature_gene.parquet' g JOIN "
                        f"'{self.root}/features/regulatory_feature.parquet' f USING(regulatory_feature_id) "
                        "WHERE f.source_feature_id='P1'")
        self.assertEqual({r["gene_id"] for r in rows if r["link_source"] == "ensembl_native"}, {"A", "B"})
        self.assertEqual({r["gene_id"] for r in rows if r["link_source"] == "tss_overlap_promoter_core"}, {"A", "B", "C"})
        self.assertEqual({r["gene_id"] for r in rows if r["link_source"] == "tss_overlap_promoter_extended"}, {"A", "B", "C", "D"})
        rows = self.sql(f"SELECT count(*) AS n FROM '{self.root}/membership/regulatory_feature_gene.parquet' g "
                        f"JOIN '{self.root}/features/regulatory_feature.parquet' f USING(regulatory_feature_id) "
                        "WHERE f.source_feature_id='E1'")
        self.assertEqual(rows[0]["n"], 0)

    def test_export_denominator_and_immutable_output(self):
        arguments = ("select", "--anchors", self.root / "anchors.parquet",
                     "--membership", self.root / "membership", "--subset", "promoter_extended",
                     "--output", self.root / "selected")
        self.call(*arguments)
        rows = self.sql(f"SELECT count(*) AS n,min(anchor_score) AS score FROM '{self.root}/selected/selected_anchors.parquet'")
        self.assertEqual(rows, [{"n": 5, "score": 7.5}])
        self.assertIn("refusing to replace", self.call(*arguments, success=False).stderr)

    def test_reject_wrong_chromosome_and_corrupted_package(self):
        self.assertIn("chromosome absent", self.annotate("bad", "chr1", success=False).stderr)
        file = self.root / "features/regulatory_feature.parquet"
        with file.open("ab") as stream:
            stream.write(b"bad")
        self.assertIn("checksum changed", self.annotate("corrupt", success=False).stderr)

    def test_parser_rejects_invalid_bound_and_preserves_unknown(self):
        file = self.root / "one.gff"
        file.write_text("1\tE\tpromoter\t101\t120\t.\t.\t.\tID=P\n")
        row = list(module.gff_rows(file, "GRCh38", "r", "hash"))[0]
        self.assertIsNone(row["extended_start"])
        file.write_text("1\tE\tpromoter\t101\t120\t.\t.\t.\tID=P;extended_start=102;extended_end=150\n")
        with self.assertRaisesRegex(ValueError, "core not contained"):
            list(module.gff_rows(file, "GRCh38", "r", "hash"))

    def test_missing_extension_cannot_become_negative_membership(self):
        file = self.root / "no-extension.gff"
        file.write_text("1\tE\tpromoter\t101\t120\t.\t.\t.\tID=P\n")
        self.call("import", "--gff", file, "--assembly", "GRCh38",
                  "--requested-release", "r", "--resolved-release", "r",
                  "--source-uri", "synthetic", "--coordinate-audit-note", "fixture",
                  "--output", self.root / "no-extension")
        result = self.call("annotate", "--features", self.root / "no-extension",
                           "--anchors", self.root / "anchors.parquet", "--chrom", "1",
                           "--assembly", "GRCh38", "--output", self.root / "invalid",
                           success=False)
        self.assertIn("extended promoter bounds missing", result.stderr)
        self.assertFalse((self.root / "invalid").exists())

    def test_select_rejects_same_count_but_wrong_anchor_keys(self):
        self.sql(f"COPY (SELECT * REPLACE ((anchor_start+1)::BIGINT AS anchor_start) "
                 f"FROM '{self.root}/anchors.parquet') TO '{self.root}/wrong.parquet' (FORMAT PARQUET)")
        result = self.call("select", "--anchors", self.root / "wrong.parquet",
                           "--membership", self.root / "membership", "--subset", "promoter_core",
                           "--output", self.root / "invalid", success=False)
        self.assertIn("exact physical anchor input", result.stderr)
        self.assertFalse((self.root / "invalid").exists())

    def test_coordinate_audit_checks_all_features(self):
        reader = self.root / "bigBedToBed"
        rows = list(module.gff_rows(self.gff, "GRCh38", "r", "hash"))
        bed = ""
        for r in rows:
            start = r["extended_start"] if r["extended_start"] is not None else r["start"]
            end = r["extended_end"] if r["extended_end"] is not None else r["end"]
            bed += "\t".join(map(str, (r["chrom"], start, end, r["source_feature_id"], 0,
                                          r["strand"], r["start"], r["end"], 0,
                                          r["feature_type"].replace("_", " ")))) + "\n"
        reader.write_text("#!/usr/bin/env python3\nimport sys\nsys.stdout.write(" + repr(bed) + ")\n")
        reader.chmod(0o755)
        self.call("audit", "--gff", self.gff, "--bigbed", self.gff,
                  "--bigbed-to-bed", reader, "--assembly", "GRCh38", "--output", self.root / "audit")
        audit = json.loads((self.root / "audit/manifest.json").read_text())
        self.assertEqual(audit["features_compared"], len(rows))
        self.assertEqual(audit["status"], "passed")
        self.call("import", "--gff", self.gff, "--assembly", "GRCh38",
                  "--requested-release", "r", "--resolved-release", "r", "--source-uri", "synthetic",
                  "--coordinate-audit-note", "fixture", "--coordinate-audit", self.root / "audit/manifest.json",
                  "--output", self.root / "audited-features")
        self.assertEqual(json.loads((self.root / "audited-features/manifest.json").read_text())["coordinate_audit"], audit)
        reader.write_text(reader.read_text().replace("\\t90\\t150", "\\t89\\t150", 1))
        result = self.call("audit", "--gff", self.gff, "--bigbed", self.gff,
                           "--bigbed-to-bed", reader, "--assembly", "GRCh38", "--output", self.root / "bad-audit",
                           success=False)
        self.assertIn("coordinate/identity mismatch", result.stderr)

    def test_production_cohort_validation(self):
        package = self.root / "membership"
        manifest_path = package / "manifest.json"
        manifest = json.loads(manifest_path.read_text())
        manifest.update(kind="genome_regulatory_membership", state="complete", chromosomes=["1"])
        manifest_path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "audited"):
            module.production_selection(package, "promoter_extended", ["1"])
        manifest["coordinate_audit"] = {
            "kind": "regulatory_coordinate_audit", "status": "passed", "assembly": "GRCh38",
        }
        manifest_path.write_text(json.dumps(manifest))
        selected = module.production_selection(package, "promoter_extended", ["1"])
        module.verify_production_selection(selected)
        sql = module.cohort_validation_sql(selected, self.root / "anchors.parquet", ["1"])
        self.sql(sql + "SELECT count(*) n FROM regulatory_membership")
        mutations = {
            "wrong_keys": "SELECT * REPLACE (anchor_start+1 AS anchor_start) FROM m",
            "duplicate": "SELECT * FROM m UNION ALL SELECT * FROM m LIMIT 12",
            "null_flag": "SELECT * REPLACE (NULL::BOOLEAN AS promoter_extended) FROM m",
            "empty": "SELECT * REPLACE (false AS promoter_core, false AS promoter_extended) FROM m",
            "not_nested": "SELECT * REPLACE (true AS promoter_core, false AS promoter_extended) FROM m",
        }
        for name, query in mutations.items():
            with self.subTest(name=name):
                path = self.root / (name + ".parquet")
                self.sql(f"CREATE TABLE m AS SELECT * FROM '{selected['path']}'; "
                         f"COPY ({query}) TO '{path}' (FORMAT PARQUET)")
                changed = dict(selected, path=str(path))
                result = subprocess.run([DUCKDB, "-batch", "-bail", ":memory:"],
                    input=module.cohort_validation_sql(changed, self.root / "anchors.parquet", ["1"]),
                    text=True, capture_output=True)
                self.assertNotEqual(result.returncode, 0, name)
                self.assertIn("regulatory membership must uniquely cover", result.stderr)

    def test_production_restart_and_finalization(self):
        source = self.root / "repo"
        (source / "scripts").mkdir(parents=True)
        for script in ("build_regulatory_annotation.py", "manage_regulatory_annotation.py"):
            shutil.copyfile(ROOT / "scripts" / script, source / "scripts" / script)
        for arguments in (["init", "-q"], ["add", "scripts"],
                          ["-c", "user.name=Test", "-c", "user.email=test@example.invalid",
                           "commit", "-qm", "fixture"]):
            subprocess.run(["git", "-C", str(source), *arguments], check=True, capture_output=True)
        evidence = self.root / "evidence"
        evidence.mkdir()
        shutil.copyfile(self.root / "anchors.parquet", evidence / "anchors.parquet")
        (evidence / "manifest.json").write_text(json.dumps({"state": "complete"}))
        with (evidence / "chromosome_file_inventory.tsv").open("w") as stream:
            writer = csv.DictWriter(stream, delimiter="\t",fieldnames=["chrom","relative_path","bytes","sha256"])
            writer.writeheader()
            writer.writerow({"chrom":"1","relative_path":"anchors.parquet",
                             "bytes":(evidence / "anchors.parquet").stat().st_size,
                             "sha256":module.sha256(evidence / "anchors.parquet")})
        catalog = self.root / "catalog"
        catalog.mkdir()
        (catalog / "manifest.json").write_text(json.dumps({"state":"complete","context_schema_version":9}))
        result = subprocess.run([DUCKDB, "-batch", "-bail", str(catalog / "context.duckdb")],
                                input=f"CREATE TABLE context_file_inventory AS SELECT * FROM (VALUES "
                                f"('transcription_start_site.parquet','1','{self.root}/tss.parquet'),"
                                f"('transcript_tss.parquet','1','{self.root}/owners.parquet'),"
                                f"('promoter.parquet','1','{self.root}/promoters.parquet')) r(dataset,chrom,absolute_path);",
                                text=True, capture_output=True)
        self.assertEqual(result.returncode,0,result.stderr)
        audit = self.root / "production-audit.json"
        audit.write_text(json.dumps({"kind":"regulatory_coordinate_audit","status":"passed",
                                     "assembly":"GRCh38","features_compared":6,
                                     "gff_sha256":module.sha256(self.gff),"source_uri":"synthetic",
                                     "requested_regulatory_release":"r","regulatory_release":"r"}))
        manager = source / "scripts/manage_regulatory_annotation.py"
        run = self.root / "run"

        def invoke(*arguments, success=True):
            result = subprocess.run([sys.executable,str(manager),*map(str,arguments),"--run-root",str(run)],
                                    text=True,capture_output=True)
            self.assertEqual(result.returncode == 0,success,result.stderr)
            return result

        invoke("prepare","--source",source,"--evidence-package",evidence,"--annotation-catalog",catalog,
               "--audit",audit,"--duckdb",Path(shutil.which(DUCKDB)),"--chromosomes","1",
               "--genome-id","human","--annotation-release","gtf-test","--promoter-definition-id","legacy")
        invoke("setup","--gff",self.gff)
        self.assertIn("Reusing",invoke("setup","--gff",self.gff).stdout)
        invoke("finalize",success=False)
        invoke("run-task","--task-index","0","--scratch-root",self.root / "scratch")
        self.assertIn("Reusing",invoke("run-task","--task-index","0","--scratch-root",self.root / "scratch").stdout)
        invoke("finalize")
        self.assertIn("Reusing",invoke("finalize").stdout)
        rows = self.sql(f"SELECT count(*) n FROM '{run}/final/tp73_anchor_regulatory_membership.parquet'")
        self.assertEqual(rows[0]["n"],11)
        with (run / "tasks/chrom-1/tp73_anchor_regulatory_membership.parquet").open("ab") as stream:
            stream.write(b"changed")
        invoke("run-task","--task-index","0","--scratch-root",self.root / "scratch",success=False)


if __name__ == "__main__":
    unittest.main()

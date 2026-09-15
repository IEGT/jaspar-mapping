#!/usr/bin/env python3
"""Build genome-wide regulatory/TSS intersections and restartable motif exports.

prepare runs on a compute node once. submit creates bounded requeue arrays;
run-task checkpoints each motif independently. finalize builds a movable
Parquet bundle and exact-file DuckDB index. Nothing deletes source or attempts.
"""
from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import tempfile

import export_regulatory_tfbs as export

SOURCES = ("manage_regulatory_tfbs.py", "export_regulatory_tfbs.py",
           "build_regulatory_annotation.py", "query_genome_scan.py")
DEFAULT_ROOT = Path("/scratch/sm718")


def save(path, value):
    with Path(path).open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def load(args):
    path = args.run_root.resolve() / "prepared/plan.json"
    plan = export.read_json(path)
    for entry in plan["source_files"] + plan["inputs"]:
        export.verify(entry)
    if (export.sha256(Path(plan['database'])) != plan['catalog_sha256']
            or export.sha256(Path(plan['annotation_package'])/'manifest.json') != plan['annotation_manifest_sha256']):
        raise ValueError('prepared catalog/annotation manifest changed')
    return plan


def space(path, reserve):
    if shutil.disk_usage(path).free < reserve:
        raise ValueError(f"free-space reserve reached: {path}")


def annotation_sql(args, gtf, features, genome, regions):
    q = export.sql_string
    definition = f"tss_upstream_{args.upstream}_downstream_{args.downstream}_v1"
    return f"""
CREATE TABLE sequence_region AS SELECT * FROM read_json_auto({q(regions)});
CREATE TABLE feature AS SELECT * FROM read_parquet({q(features)},hive_partitioning=false);
CREATE TABLE transcript AS
SELECT column0 AS chrom, column3-1 AS start, column4 AS "end", column6 AS strand,
 trim(regexp_extract(column8,'gene_id "?([^";]+)"?',1)) AS gene_id,
 trim(regexp_extract(column8,'transcript_id "?([^";]+)"?',1)) AS transcript_id,
 trim(regexp_extract(column8,'gene_name "?([^";]+)"?',1)) AS gene_name
FROM read_csv({q(gtf)},delim='\\t',header=false,comment='#',quote='',auto_detect=false,
 columns={{'column0':'VARCHAR','column1':'VARCHAR','column2':'VARCHAR',
 'column3':'BIGINT','column4':'BIGINT','column5':'VARCHAR','column6':'VARCHAR',
 'column7':'VARCHAR','column8':'VARCHAR'}})
WHERE lower(column2)='transcript' AND column0 IN (SELECT chrom FROM sequence_region);
CREATE TEMP TABLE gtf_validation AS SELECT CASE WHEN NOT EXISTS(SELECT 1 FROM transcript) OR EXISTS(
 SELECT 1 FROM transcript t JOIN sequence_region s USING(chrom)
 WHERE t.start IS NULL OR t."end" IS NULL OR t.start<0 OR t.start>=t."end" OR t."end">s.length
 OR strand IS NULL OR strand NOT IN ('+','-') OR gene_id='' OR transcript_id='')
 OR EXISTS(SELECT transcript_id FROM transcript GROUP BY transcript_id HAVING count(*)<>1)
 THEN error('Invalid or nonunique GTF transcript records') END;
CREATE TABLE tss_owner AS SELECT {q(genome['genome_id'])}::VARCHAR AS genome_id,
 {q(args.annotation_release)}::VARCHAR AS annotation_release,*,
 CASE WHEN strand='+' THEN start ELSE "end"-1 END::BIGINT AS tss_start
FROM transcript;
CREATE TABLE transcription_start_site AS SELECT DISTINCT
 md5(concat_ws('|',genome_id,annotation_release,chrom,tss_start::VARCHAR,strand)) AS tss_id,
 genome_id,annotation_release,chrom,tss_start AS start,tss_start+1 AS "end",strand FROM tss_owner;
CREATE TABLE transcript_tss AS SELECT o.genome_id,o.annotation_release,o.gene_id,o.gene_name,
 o.transcript_id,s.tss_id FROM tss_owner o JOIN transcription_start_site s
 ON o.chrom=s.chrom AND o.tss_start=s.start AND o.strand=s.strand;
CREATE TABLE promoter AS SELECT md5(concat_ws('|',tss_id,{q(definition)})) AS promoter_id,
 tss_id,t.genome_id,t.annotation_release,{q(definition)}::VARCHAR AS promoter_definition_id,
 t.chrom,t.strand,
 greatest(0,t.start-CASE WHEN t.strand='+' THEN {args.upstream} ELSE {args.downstream} END)::BIGINT AS promoter_start,
 least(r.length,t.start+1+CASE WHEN t.strand='+' THEN {args.downstream} ELSE {args.upstream} END)::BIGINT AS promoter_end,
 t.start AS tss_start,t."end" AS tss_end,{args.upstream}::INTEGER AS upstream_bp,
 {args.downstream}::INTEGER AS downstream_bp
FROM transcription_start_site t JOIN sequence_region r USING(chrom);
"""


def prepare(args):
    root = args.run_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    with export.output_lock(root / "prepare"):
        if (root / "prepared").exists():
            plan = load(args)
            requested = (str(args.package.resolve()), str(args.features.resolve()),
                         str(args.gtf.resolve()), args.upstream, args.downstream, args.batch_size,
                         args.annotation_release)
            recorded = (plan['package'], plan['features_source'], plan['gtf_source'],
                        plan['upstream_bp'], plan['downstream_bp'], plan['batch_size'],
                        plan['annotation_release'])
            if requested != recorded:
                raise ValueError("existing plan differs from requested inputs/window/batching")
            print("I: Reusing verified prepared run")
            return
        space(root, args.minimum_free_bytes)
        stage = Path(tempfile.mkdtemp(prefix=".prepare-", dir=root))
        annotation = stage / "annotation"
        annotation.mkdir()
        args.scratch_root.mkdir(parents=True, exist_ok=True)
        scratch = Path(tempfile.mkdtemp(prefix="regulatory-tss-prepare-", dir=args.scratch_root))
        os.environ['TMPDIR'] = str(scratch)
        args.database = None
        package, database = export.package_paths(args)
        scan = export.read_json(package / "manifest.json")
        if scan.get('state') != 'complete':
            raise ValueError('source scan is not complete')
        features_manifest = args.features / "manifest.json"
        fm = export.read_json(features_manifest)
        audit = fm.get("coordinate_audit") or {}
        if (fm.get('kind') != 'ensembl_regulatory_annotation' or audit.get('status') != 'passed'
                or audit.get('gff_sha256') != fm.get('source_sha256')
                or audit.get('features_compared') != fm.get('feature_count')):
            raise ValueError('full coordinate-audited regulatory annotation required')
        pinned = [export.record(package / 'manifest.json'), export.record(features_manifest),
                  export.record(args.gtf)]
        for name in ('regulatory_feature.parquet', 'regulatory_feature_gene.parquet'):
            entries = [e for e in fm['files'] if e['path'] == name]
            if len(entries) != 1:
                raise ValueError(f'missing/duplicate feature inventory: {name}')
            entry = dict(entries[0], path=str((args.features / name).resolve()))
            export.verify(entry)
            pinned.append(entry)
            shutil.copyfile(entry['path'], annotation / name)
            export.verify(dict(entry, path=str(annotation / name)))
        shutil.copyfile(args.gtf, scratch / args.gtf.name)
        export.verify(dict(pinned[2], path=str(scratch / args.gtf.name)))
        attach = f"ATTACH {export.sql_string(database)} AS scan (READ_ONLY);"
        genome_rows = export.sql(args, attach + 'SELECT * FROM scan.genome;')
        if len(genome_rows) != 1 or genome_rows[0]['genome_id'] != scan['genome_id']:
            raise ValueError('ambiguous genome identity')
        genome = genome_rows[0]
        if genome['assembly_name'] != fm['assembly'] or audit['assembly'] != fm['assembly']:
            raise ValueError('regulatory/scan assembly mismatch')
        catalog = stage / 'catalog.duckdb'
        export.sql(args, attach + f"ATTACH {export.sql_string(catalog)} AS catalog;" + ''.join(
            f'CREATE TABLE catalog.{table} AS SELECT * FROM scan.{table};'
            for table in ('genome', 'motif_metadata', 'sequence_region', 'scan_file_inventory')))
        regions = export.sql(args, attach + 'SELECT * FROM scan.sequence_region WHERE included_in_scan ORDER BY sequence_order;')
        inventory = export.sql(args, attach + 'SELECT chrom,motif_id,strand,state,coordinate_mode,minimum_score,'
                               'minimum_pwm_relative_score,maximum_pwm_relative_score,bytes,emitted_hits FROM scan.scan_file_inventory;')
        motifs = sorted(r['motif_id'] for r in export.sql(args, attach + 'SELECT motif_id FROM scan.motif_metadata;'))
        if any(not re.fullmatch(r'[A-Za-z0-9_.-]+', name) or name in ('.', '..')
               for name in motifs + [r['chrom'] for r in regions]):
            raise ValueError('unsafe motif/chromosome identifier')
        keys = {(r['chrom'], r['motif_id'], r['strand']) for r in inventory}
        expected = {(r['chrom'], m, s) for r in regions for m in motifs for s in ('+', '-')}
        if len(keys) != len(inventory) or keys != expected:
            raise ValueError('scan inventory is not complete and unique across all motifs/chromosomes/strands')
        if any(r['state'] != 'complete' or r['coordinate_mode'] != 'bed'
               or r['minimum_pwm_relative_score'] is not None or r['maximum_pwm_relative_score'] is not None
               or r['minimum_score'] > (-5 if r['motif_id'] == 'MA0861.2' else -1) for r in inventory):
            raise ValueError('production requires the permissive -1 / TP73 -5 scan without relative-score censoring')
        save(stage / 'sequence_regions.json', regions)
        sql = annotation_sql(args, scratch / args.gtf.name, annotation / 'regulatory_feature.parquet',
                             genome, stage / 'sequence_regions.json')
        sql += f"CREATE TEMP TABLE feature_validation AS SELECT CASE WHEN (SELECT count(*) FROM feature)<>{fm['feature_count']} THEN error('incomplete features') END;"
        for table in ('transcription_start_site', 'transcript_tss', 'promoter'):
            sql += f"COPY {table} TO {export.sql_string(annotation / (table+'.parquet'))} (FORMAT PARQUET, COMPRESSION ZSTD);"
        for region in regions:
            chrom = region['chrom']
            directory = annotation / 'chromosomes' / chrom
            directory.mkdir(parents=True)
            for table, name in [('feature', 'features'), ('promoter', 'promoters')]:
                sql += f"COPY (SELECT * FROM {table} WHERE chrom={export.sql_string(chrom)}) TO "
                sql += f"{export.sql_string(directory / (name+'.parquet'))} (FORMAT PARQUET, COMPRESSION ZSTD);"
        sql += '''SELECT r.chrom,r.length,
          (SELECT count(*) FROM feature f WHERE f.chrom=r.chrom) AS feature_count,
          (SELECT count(*) FROM promoter p WHERE p.chrom=r.chrom) AS promoter_count,
          (SELECT count(*) FROM promoter p JOIN feature f ON p.promoter_start<coalesce(f.extended_end,f."end")
           AND p.promoter_end>coalesce(f.extended_start,f.start) WHERE p.chrom=r.chrom AND f.chrom=r.chrom) AS intersections
          FROM sequence_region r ORDER BY r.sequence_order;'''
        # Persist the small summary so the CLI emits only one JSON result array.
        split = sql.rfind('SELECT r.chrom')
        export.sql(args, sql[:split] + f"COPY ({sql[split:].rstrip(';')}) TO {export.sql_string(stage/'coverage.parquet')} (FORMAT PARQUET);")
        coverage = export.sql(args, f"SELECT * FROM read_parquet({export.sql_string(stage/'coverage.parquet')});")
        def relative_record(path):
            return dict(export.record(path), path=str(path.relative_to(annotation)))
        for row in coverage:
            if row['feature_count'] and not row['promoter_count']:
                raise ValueError(f"regulatory chromosome lacks GTF TSS coverage: {row['chrom']}")
            row['coverage'] = 'annotated_intersection' if int(row.pop('intersections')) else 'known_empty_intersection'
            for key in ('features', 'promoters'):
                row[key] = relative_record(annotation / 'chromosomes' / row['chrom'] / (key+'.parquet'))
        am = {'schema_version': 1, 'kind': 'regulatory_tfbs_annotation', 'state': 'complete',
              'genome_id': genome['genome_id'], 'assembly': genome['assembly_name'],
              'annotation_release': args.annotation_release,
              'promoter_definition_id': f'tss_upstream_{args.upstream}_downstream_{args.downstream}_v1',
              'upstream_bp': args.upstream, 'downstream_bp': args.downstream,
              'coordinate_rule': 'BED_half_open_offsets_include_TSS_base_clamped_to_sequence',
              'regulatory_release': fm['regulatory_release'], 'regulatory_gff_sha256': fm['source_sha256'],
              'coordinate_audit': audit, 'chromosomes': coverage, 'inputs': pinned,
              'files': [relative_record(annotation / name) for name in (
                  'regulatory_feature.parquet','regulatory_feature_gene.parquet','transcription_start_site.parquet',
                  'transcript_tss.parquet','promoter.parquet')] + [r[k] for r in coverage for k in ('features','promoters')]}
        save(annotation/'manifest.json', am)
        hit_counts = {}
        for row in inventory:
            key = (row['chrom'], row['motif_id'])
            hit_counts[key] = hit_counts.get(key, 0) + row['emitted_hits']
        tasks = []
        for offset in range(0, len(motifs), args.batch_size):
            for region in coverage:
                tasks.append({'index': len(tasks), 'chrom': region['chrom'], 'coverage': region['coverage'],
                              'motifs': [{'motif_id': m, 'source_rows': hit_counts[region['chrom'], m]}
                                         for m in motifs[offset:offset+args.batch_size]]})
        source = Path(__file__).resolve().parent
        commit = subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
        clean = subprocess.run(['git','-C',str(source),'diff','--quiet','HEAD']).returncode == 0
        prepared = root / 'prepared'
        plan = {'schema_version': 1, 'source_commit': commit, 'source_clean': clean,
                'source_files': [export.record(source/name) for name in SOURCES], 'inputs': pinned,
                'package': str(package), 'features_source': str(args.features.resolve()), 'gtf_source': str(args.gtf.resolve()),
                'database': str(prepared/'catalog.duckdb'), 'annotation_package': str(prepared/'annotation'),
                'annotation_manifest_sha256': export.sha256(annotation/'manifest.json'),
                'catalog_sha256': export.sha256(catalog), 'upstream_bp': args.upstream, 'downstream_bp': args.downstream,
                'annotation_release': args.annotation_release, 'batch_size': args.batch_size,
                'duckdb': str(Path(shutil.which(args.duckdb)).resolve()), 'tasks': tasks,
                'minimum_free_bytes': args.minimum_free_bytes, 'motif_count': len(motifs),
                'chromosomes': [r['chrom'] for r in regions]}
        save(stage/'plan.json', plan)
        for entry in pinned:
            export.verify(entry)
        os.rename(stage, prepared)
        print(json.dumps({'tasks': len(tasks), 'motifs': len(motifs), 'chromosomes': plan['chromosomes']}))


def run_task(args):
    plan = load(args)
    index = args.task_index if args.task_index is not None else int(os.environ['SLURM_ARRAY_TASK_ID']) + args.task_offset
    if not 0 <= index < len(plan['tasks']):
        raise ValueError('task index outside plan')
    task = plan['tasks'][index]
    root = args.run_root.resolve()
    marker = root / 'tasks' / f'{index:05}.json'
    marker.parent.mkdir(exist_ok=True)
    digest = export.sha256(root/'prepared/plan.json')
    with export.output_lock(marker):
        if marker.exists():
            done = export.read_json(marker)
            if done['plan_sha256'] != digest:
                raise ValueError('completed task belongs to another plan')
            for row in done['outputs']:
                export.verify(row['manifest'])
                export.verify(row['payload'])
            print(f'I: Reusing completed task {index}')
            return
        space(root, plan['minimum_free_bytes'])
        args.scratch_root.mkdir(parents=True, exist_ok=True)
        scratch = Path(tempfile.mkdtemp(prefix=f'regulatory-{os.environ.get("SLURM_JOB_ID","local")}-',dir=args.scratch_root))
        os.environ['TMPDIR'] = str(scratch)
        outputs = []
        if task['coverage'] != 'known_empty_intersection':
            for motif in task['motifs']:
                space(root, plan['minimum_free_bytes'])
                target = root/'exports'/task['chrom']/motif['motif_id']
                cmd = [sys.executable, str(Path(__file__).with_name('export_regulatory_tfbs.py')),
                       '--package',plan['package'],'--database',plan['database'],
                       '--annotation-package',plan['annotation_package'],'--motif',motif['motif_id'],
                       '--chrom',task['chrom'],'--whole-chromosome','--scope','regulatory_and_tss','--source-floor',
                       '--output',str(target),'--resume','--scratch-directory',str(scratch),
                       '--max-rows','0','--max-output-bytes','0','--minimum-free-bytes',str(plan['minimum_free_bytes']),
                       '--memory-limit',args.memory_limit,'--threads',str(args.threads),'--duckdb',plan['duckdb'],
                       '--timeout-seconds','6900']
                if motif['source_rows'] > 10000000:
                    cmd += ['--scan-chunk-bp','5000000']
                print(f"I: task={index} chrom={task['chrom']} motif={motif['motif_id']}",flush=True)
                subprocess.run(cmd,check=True)
                manifest = export.read_json(target/'manifest.json')
                outputs.append({'motif_id': motif['motif_id'], 'chrom': task['chrom'],
                                'rows': manifest['rows'], 'bytes': manifest['parquet_bytes'],
                                'manifest': export.record(target/'manifest.json'),
                                'payload': dict(next(f for f in manifest['files'] if f['path']=='motif_hits.parquet'),
                                                path=str(target/'motif_hits.parquet'))})
        pending = marker.with_name(f'.{marker.name}.{os.getpid()}')
        save(pending, {'plan_sha256':digest,'task':task,'outputs':outputs})
        os.rename(pending,marker)
        print(f'I: Task {index} complete',flush=True)


def finalize(args):
    plan = load(args)
    root = args.run_root.resolve()
    digest = export.sha256(root/'prepared/plan.json')
    with export.output_lock(root/'final'):
        if (root/'final').exists():
            m = export.read_json(root/'final/manifest.json')
            if m['production_plan_sha256'] != digest:
                raise ValueError('another final package already exists')
            for entry in m['files']:
                export.verify(dict(entry,path=str(root/'final'/entry['path'])))
            print('I: Final package already complete')
            return
        stage = Path(tempfile.mkdtemp(prefix='.final-',dir=root))
        shutil.copytree(root/'prepared/annotation',stage/'annotation')
        inventory=[]
        for task in plan['tasks']:
            completed=export.read_json(root/'tasks'/f"{task['index']:05}.json")
            if completed['plan_sha256']!=digest or completed['task']!=task:
                raise ValueError('incompatible task completion')
            expected={m['motif_id'] for m in task['motifs']}
            outputs=completed['outputs']
            if task['coverage']=='known_empty_intersection':
                if outputs: raise ValueError('unexpected files for empty intersection')
                inventory.extend({'chrom':task['chrom'],'motif_id':m,'rows':0,'bytes':0,'path':None,
                                  'sha256':None,'state':'known_empty_intersection'} for m in sorted(expected))
                continue
            if len(outputs)!=len(expected) or {r['motif_id'] for r in outputs}!=expected:
                raise ValueError('task motif coverage is incomplete')
            for row in outputs:
                expected_payload=root/'exports'/task['chrom']/row['motif_id']/'motif_hits.parquet'
                if Path(row['payload']['path'])!=expected_payload:
                    raise ValueError('unexpected task payload path')
                export.verify(row['manifest'])
                p=Path(row['payload']['path'])
                if p.stat().st_size!=row['payload']['bytes']:
                    raise ValueError('payload changed size')
                m=export.read_json(row['manifest']['path'])
                if (m['state']!='complete' or m['scope']!='regulatory_and_tss' or m['score_selection']!='source_retention'
                        or m['rows']!=row['rows'] or m['parquet_bytes']!=row['bytes']
                        or next(f['sha256'] for f in m['files'] if f['path']=='motif_hits.parquet')!=row['payload']['sha256']
                        or m['inputs']['annotation_manifest']['sha256']!=plan['annotation_manifest_sha256']):
                    raise ValueError('wrong export selection/provenance')
                relative=Path('hits')/f"chrom={task['chrom']}"/(row['motif_id']+'.parquet')
                (stage/relative).parent.mkdir(parents=True,exist_ok=True)
                os.link(p,stage/relative)
                inventory.append({'chrom':task['chrom'],'motif_id':row['motif_id'],'rows':row['rows'],
                                  'bytes':row['bytes'],'path':str(relative),'sha256':row['payload']['sha256'],
                                  'state':'complete'})
        if len(inventory)!=plan['motif_count']*len(plan['chromosomes']):
            raise ValueError('incomplete genome/motif cross product')
        save(stage/'file_inventory.json',inventory)
        args.package=Path(plan['package'])
        args.duckdb=plan['duckdb']
        outdb=stage/'regulatory_tfbs.duckdb'
        export.sql(args,f"ATTACH {export.sql_string(plan['database'])} AS source (READ_ONLY);"
                   f"ATTACH {export.sql_string(outdb)} AS bundle;"
                   "CREATE TABLE bundle.genome AS SELECT * FROM source.genome;"
                   "CREATE TABLE bundle.motif_metadata AS SELECT * FROM source.motif_metadata;"
                   "CREATE TABLE bundle.sequence_region AS SELECT * FROM source.sequence_region;"
                   "CREATE TABLE bundle.scan_file_inventory AS SELECT * FROM source.scan_file_inventory;"
                   f"CREATE TABLE bundle.file_inventory AS SELECT * FROM read_json_auto({export.sql_string(stage/'file_inventory.json')});"
                   "CREATE MACRO bundle.motif_hits(files) AS TABLE SELECT * FROM read_parquet(files,hive_partitioning=false);")
        (stage/'schema.sql').write_text('''-- Run from the package directory. Metadata index: regulatory_tfbs.duckdb.
CREATE VIEW promoter AS SELECT * FROM read_parquet('annotation/promoter.parquet');
CREATE VIEW transcript_tss AS SELECT * FROM read_parquet('annotation/transcript_tss.parquet');
CREATE VIEW transcription_start_site AS SELECT * FROM read_parquet('annotation/transcription_start_site.parquet');
CREATE VIEW regulatory_feature AS SELECT * FROM read_parquet('annotation/regulatory_feature.parquet');
''')
        am=export.read_json(stage/'annotation/manifest.json')
        save(stage/'manifest.json',{'schema_version':1,'kind':'genome_regulatory_tfbs_subset','state':'complete',
             'genome_id':am['genome_id'],'assembly':am['assembly'],'annotation_release':am['annotation_release'],
             'regulatory_release':am['regulatory_release'],'promoter_definition_id':am['promoter_definition_id'],
             'production_plan_sha256':digest,'source_commit':plan['source_commit'],'database':'regulatory_tfbs.duckdb',
             'coordinate_mode':'bed_0based_half_open','scope':'regulatory_and_tss','score_selection':'source_retention',
             'upstream_bp':plan['upstream_bp'],'downstream_bp':plan['downstream_bp'],
             'chromosomes':plan['chromosomes'],'motif_count':plan['motif_count'],'requires_tp73':False,
             'rows':sum(r['rows'] for r in inventory),'parquet_bytes':sum(r['bytes'] for r in inventory),
             'complete_for_declared_annotation_and_source_floors':True,'complete_genome_scan':False,
             'payload_verification':'SHA256 verified at scratch staging/export; finalization verifies manifests and byte sizes',
             'files':[dict(export.record(stage/name),path=name) for name in
                      ('file_inventory.json','regulatory_tfbs.duckdb','schema.sql','annotation/manifest.json')]})
        os.rename(stage,root/'final')
        print(f'I: Final package: {root / "final"}')


def submit(args):
    plan=load(args)
    if not plan['source_clean']:
        raise ValueError('refusing production submission from uncommitted source')
    root=args.run_root.resolve()
    if (root/'submissions.jsonl').exists():
        raise ValueError('this plan already has submissions; inspect their state before explicit recovery')
    (root/'logs').mkdir(exist_ok=True)
    groups=math.ceil(len(plan['tasks'])/1000)
    group_size=math.ceil(len(plan['tasks'])/groups)
    concurrency=max(1,args.concurrent//groups)
    if groups>args.concurrent: raise ValueError('too few concurrency slots for arrays')
    ids=[]
    script=str(Path(__file__).resolve())
    def dispatch(label, extra, command):
        result=subprocess.check_output(['sbatch','--parsable','--account=cluster','--partition=requeue','--requeue',
              '--cpus-per-task='+str(args.threads),'--mem='+args.slurm_memory,'--time='+args.slurm_time,
              '--chdir='+str(root),'--open-mode=append','--job-name='+label,
              '--output='+str(root/'logs/%x-%A_%a.log'),*extra,'--wrap='+shlex.join(
              ['env','PYTHONDONTWRITEBYTECODE=1','PYTHONNOUSERSITE=1',sys.executable,script,*command])],text=True).strip().split(';')[0]
        with (root/'submissions.jsonl').open('a') as stream:
            stream.write(json.dumps({'job_id':result,'name':label,'command':command,'sbatch_options':extra})+'\n')
        return result
    for offset in range(0,len(plan['tasks']),group_size):
        count=min(group_size,len(plan['tasks'])-offset)
        ids.append(dispatch('regulatory_tfbs', [f'--array=0-{count-1}%{concurrency}'],
                   ['run-task','--run-root',str(root),'--task-offset',str(offset),'--scratch-root',str(args.scratch_root),
                    '--memory-limit',args.memory_limit,'--threads',str(args.threads)]))
    final=dispatch('regulatory_tfbs_finalize',['--dependency=afterok:'+':'.join(ids)],['finalize','--run-root',str(root)])
    print(json.dumps({'arrays':ids,'finalizer':final,'maximum_concurrent':groups*concurrency}))


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    for name, help_text in [('prepare','Build all physical TSS windows and immutable task plan (compute node)'),
                            ('run-task','Export one chromosome/motif batch with per-motif resume'),
                            ('finalize','Validate completeness and publish portable bundle'),
                            ('submit','Submit prepared tasks and finalizer to requeue')]:
        s=sub.add_parser(name,help=help_text)
        s.add_argument('--run-root',type=Path,required=True,help='dedicated durable run directory; never overwritten')
        s.add_argument('--memory-limit',default='16GB',help='DuckDB memory ceiling (default 16GB; disk spill disabled)')
        s.add_argument('--threads',type=int,default=2,help='DuckDB threads and Slurm CPUs (default 2)')
        s.add_argument('--timeout-seconds',type=int,default=14400,help='preparation/finalization SQL deadline per call (default 14400)')
        s.add_argument('--scratch-root',type=Path,default=DEFAULT_ROOT,help='job-local parent, default /scratch/sm718; cleanup belongs to cluster')
        if name=='prepare':
            s.add_argument('--package',type=Path,required=True,help='complete permissive genome scan, all original motif floors')
            s.add_argument('--features',type=Path,required=True,help='full audited Ensembl regulatory feature package')
            s.add_argument('--gtf',type=Path,required=True,help='whole-genome GTF/GTF.gz with transcript records')
            s.add_argument('--annotation-release',default='ensembl_113',help='pinned GTF release label (default ensembl_113)')
            s.add_argument('--upstream',type=int,default=700,help='upstream offset from the TSS base (default 700)')
            s.add_argument('--downstream',type=int,default=300,help='downstream offset from the TSS base (default 300)')
            s.add_argument('--batch-size',type=int,default=128,help='motifs per chromosome job, with per-motif checkpoints (default 128)')
            s.add_argument('--minimum-free-bytes',type=int,default=50*1024**3,help='durable free-space reserve, not an export-size cap (default 50 GiB)')
            s.add_argument('--duckdb',default='duckdb',help='DuckDB CLI path recorded in the immutable plan')
        if name=='run-task':
            s.add_argument('--task-index',type=int,help='exact plan task index; default SLURM_ARRAY_TASK_ID + offset')
            s.add_argument('--task-offset',type=int,default=0,help='offset for subsequent Slurm arrays')
        if name=='submit':
            s.add_argument('--concurrent',type=int,default=20,help='maximum concurrent workers across arrays (default 20)')
            s.add_argument('--slurm-memory',default='24G',help='Slurm allocation per job, above DuckDB ceiling (default 24G)')
            s.add_argument('--slurm-time',default='04:00:00',help='Slurm wall-time per restartable job (default 4 hours)')
    return p


def main():
    args=parser().parse_args()
    try:
        if args.threads<1 or (args.command=='prepare' and (args.upstream<0 or args.downstream<0 or args.batch_size<1 or args.minimum_free_bytes<0)):
            raise ValueError('invalid dimensions/resources')
        globals()[args.command.replace('-','_')](args)
        return 0
    except (OSError,ValueError,KeyError,export.QueryError,subprocess.SubprocessError) as error:
        print(f'E: {error}',file=sys.stderr)
        return 1


if __name__=='__main__':
    raise SystemExit(main())

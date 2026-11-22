"""Isolated APHCP execution acceptance on Pro5000; no external requests.

python -B verify_holo_execution.py CODE_DIR APHCP_DATA_ROOT SCRATCH_EVIDENCE.json
"""
from copy import deepcopy
from datetime import date, timedelta
import io
import json
from pathlib import Path
import runpy
import shutil
import sqlite3
import subprocess
import sys
from tempfile import TemporaryDirectory
import urllib.request
from unittest.mock import patch

import pyarrow.dataset as ds
import pyarrow.parquet as pq

code, root, evidence = map(lambda p: Path(p).resolve(), sys.argv[1:])
sys.path.insert(0, str(code))
from m1_holo_config import date_filters, read_json as read, sha256, write_json as put
from m1_holo_windows import collect, plan_windows

runs = root/'runs'
sept, octo, pilot = [runs/n for n in ('holo-m1-202109-20261121a','holo-m1-202110-20261121a','holo-m2-202109-20261121a')]
release = root.parent/'raw'/'APHCP_upstream_resources'/'PLINDER'/'2026-09'
registry = root/'registry'/'aphcp_registry.sqlite3'
index_tables = ('entry_metadata','entry_sources','annotation_table','entry_biounit_chains','entry_chains','ligand_pocket_residues')


def invoke(name, *args):
    path=code/(name+'.py')
    with patch.object(sys,'argv',[str(path),*map(str,args)]):
        try:
            runpy.run_path(str(path),run_name='__main__')
        except SystemExit as exc:
            return exc.code
    return 0


class Response(io.BytesIO):
    def __init__(self, payload, status=200, url='fixture', headers=None):
        super().__init__(json.dumps(payload).encode()); self.status=status; self.url=url; self.headers=headers or {}
    def geturl(self): return self.url


def sql(path, query):
    with sqlite3.connect(path.as_uri()+'?mode=ro',uri=True) as db:
        return db.execute(query).fetchall()


def stage(argv, check):
    file=Path(argv[2]); args=argv[3:]
    with patch.object(sys,'argv',[str(file),*args]):
        try:
            runpy.run_path(str(file),run_name='__main__')
        except SystemExit as exc:
            if exc.code:
                raise subprocess.CalledProcessError(2,argv) from exc


def m1_fixture(tmp):
    fixture=tmp/'m1'; run=fixture/'runs'/'tiny'; mini=fixture/'index'
    mini.mkdir(parents=True); copy=fixture/'registry'/registry.name; copy.parent.mkdir()
    with sqlite3.connect(registry.as_uri()+'?mode=ro',uri=True) as src, sqlite3.connect(copy) as dst:
        src.backup(dst)
    seed = read(sept/'data'/'entries.json')[:2]
    ids = [row['rcsb_id'] for row in seed]
    day = date.fromisoformat(seed[0]['rcsb_accession_info']['initial_release_date'][:10]); end=day+timedelta(days=1)
    assert seed[1]['rcsb_accession_info']['initial_release_date'][:10]==str(day)
    seed[0]['rcsb_entry_info']['resolution_combined']=None
    for key in index_tables:
        src=release/'index'/f'{key}.parquet'; dst=mini/src.name
        table=ds.dataset(str(src),format='parquet').to_table(filter=ds.field('entry_pdb_id').isin([i.lower() for i in ids]))
        pq.write_table(table,dst)
        rec=read(src.with_suffix('.parquet.source.json')); size=dst.stat().st_size
        put(dst.with_suffix('.parquet.source.json'),{**rec,'bytes':size,'expected_bytes':size,'fixture_only':True})
    cfg=read(sept/'inputs'/'query.json'); conf=fixture/'query.json'
    bounds=date_filters(cfg['m1_1_request']); bounds['greater_or_equal']['value']=str(day); bounds['less']['value']=str(end)
    cfg.update(page_size=1,batch_id='holo-m1-resume-fixture',run_id='holo-m1-resume-fixture-a',cutoff_date=None,
               plinder_index=str(mini),budget_reference_run=str(runs/'holo-pilot-20261120d'))
    put(conf,cfg)
    acc=seed[0]['rcsb_accession_info']; rev=f"{acc['major_revision']}.{acc['minor_revision']}"
    low=ids[0].lower(); url=f'https://files-versioned.rcsb.org/pdb_versioned/data/entries/{low[1:3]}/pdb_0000{low}/pdb_0000{low}_xyz_v{rev.replace(".","-")}.cif.gz'
    put(run/'budget'/'archive_head.json',[{'pdb_id':ids[0],'revision':rev,'url':url,'status':200,'content_length':1000}])
    head_hash=sha256(run/'budget'/'archive_head.json'); mode=0; events=[]
    def transport(request, **options):
        if request.get_method()=='HEAD':
            assert ids[1].lower() in request.full_url
            events.append(['HEAD',ids[1]])
            return Response({},url=request.full_url,headers={'Content-Length':'1000'})
        body=json.loads(request.data)
        if 'variables' in body:
            wanted=body['variables']['ids']; events.append(['GraphQL',wanted]); assert set(wanted)<=set(ids)
            return Response({'data':{'entries':[seed[ids.index(i)] if i==ids[0] or mode>=2 else None for i in wanted]}})
        page=body['request_options'].get('paginate')
        if not page:
            key=date_filters(body)['greater_or_equal']['value']; events.append(['count',key])
            return Response({'total_count':2 if str(day)==key and mode<4 else 0})
        offset=page['start']; events.append(['Search',offset]); assert offset in (0,1)
        if mode>=4: return Response(None,status=204)
        return Response({'total_count':2,'result_set':[{'identifier':ids[offset]}] if offset==0 or mode>=1 else []})
    log=[]; args=('m1_holo_inventory',run,'--config',conf,'--register')
    with patch.object(subprocess,'run',stage),patch.object(urllib.request,'urlopen',transport):
        assert invoke(*args)==2 and read(run/'progress.json')['refresh_offsets']==[1]
        log.append(list(events)); events.clear(); mode=1
        assert invoke(*args)==2 and read(run/'progress.json')['missing_ids']==[ids[1]]
        assert events==[['Search',1],['GraphQL',ids],['GraphQL',[ids[1]]]]
        assert read(run/'data'/'entries.json')[0]['rcsb_entry_info']['resolution_combined'] is None
        assert not (run/'run_receipt.json').exists() and not sql(copy,"SELECT id FROM batch WHERE id='holo-m1-resume-fixture'")
        assert 'Unresolved Search or Entry metadata gaps' in str(invoke('m1_holo_receipt',run,'--no-register'))
        log.append(list(events)); events.clear()
        unavailable=fixture/'runs'/'confirmed-unavailable'
        for rel in ('search/summary.json','data/summary.json','data/entries.json'):
            path=unavailable/rel; path.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(run/rel,path)
        proof={'source':'fixture archive notice','reference':'fixture:notice-1','reason':'simulated documented withdrawal'}
        put(unavailable/'inputs'/'confirmed_unavailable.json',{ids[1]:proof})
        assert invoke('m1_holo_fetch',unavailable,'metadata')==0 and not events
        closed=read(unavailable/'data'/'summary.json')
        assert closed['received']==1 and closed['missing_ids']==[] and closed['confirmed_unavailable']=={ids[1]:proof}
        mode=2
        put(run/'inventory'/'inventory_summary.json',{'entry_metadata_sha256':'pre-fix'})
        put(run/'budget'/'budget_summary.json',{'stale':True})
        assert invoke(*args)==0
        assert events==[['GraphQL',[ids[1]]],['HEAD',ids[1]]],events
        assert sha256(run/'budget'/'archive_head.json')!=head_hash
        log.append(list(events)); events.clear(); mode=3
        assert invoke(*args)==0 and not events
        log.append(list(events))
        cfg.update(batch_id='holo-m1-resume-empty',run_id='holo-m1-resume-empty-a')
        put(conf,cfg); empty=fixture/'runs'/'empty'; mode=4
        assert invoke('m1_holo_inventory',empty,'--config',conf)==0
        assert read(empty/'run_receipt.json')['empty_window'] and not (empty/'data').exists()
        empty_events=list(events); events.clear()
        assert invoke('m1_holo_inventory',empty,'--config',conf)==0 and not events
    final=read(run/'run_receipt.json'); assert final['entry_metadata']['missing_ids']==[]
    states=sql(copy,"SELECT status FROM run WHERE id LIKE 'holo-m1-resume-fixture-a-%'")
    assert final['status']=='INVENTORY_COMPLETE' and states==[('SUCCEEDED',)]*3
    return {'ids':ids,'requests_by_attempt':log,'empty_requests_first_only':empty_events,'confirmed_unavailable_fixture':closed['confirmed_unavailable'],'receipt':final['entry_metadata'],'sqlite_stage_states':states}


def cross_fixture(tmp):
    plan=tmp/'plan.json'; out=tmp/'cross.json'
    def spec(run, cfg=None): return {'run_directory':str(run),'config':str(cfg or run/'inputs'/'query.json')}
    query=deepcopy(read(sept/'inputs'/'query.json')['m1_1_request']); bounds=date_filters(query)
    bounds['less']['value']='2021-11-01'; query['request_options'].pop('paginate'); query['request_options']['return_counts']=True
    baseline={'scope':{'start_inclusive':'2021-09-01','end_exclusive':'2021-11-01'},
              'denominator':{'total_count':2040,'queried_utc':read(octo/'search'/'summary.json')['finished_utc'],
                             'request':query,'source':'fixture derived from the two existing Search receipts'},
              'windows':[spec(sept),spec(octo)]}
    put(plan,baseline)
    with patch.object(subprocess,'run',side_effect=AssertionError('completed window was dispatched')):
        assert invoke('m1_holo_inventory',sept,'--register')==0 and invoke('m1_holo_inventory',octo,'--register')==0
        assert invoke('m1_holo_windows',plan,'--output',out,'--register','--jobs','2')==0
    healthy=read(out); assert healthy['unique_ids']==2040 and not healthy['dispatched']
    missing=tmp/'no-october'; put(plan,{**baseline,'windows':[spec(sept),spec(missing,octo/'inputs'/'query.json')]})
    launches=[]
    def only_missing(args):
        assert Path(args[3])==missing; launches.append(args[3]); return subprocess.CompletedProcess(args,2)
    with patch.object(subprocess,'run',only_missing):
        assert invoke('m1_holo_windows',plan,'--output',out,'--jobs','2')==2
    miss=read(out); assert len(miss['missing_windows'])==1 and miss['unique_ids']==1079 and len(launches)==1
    broken=tmp/'broken-october'
    for file in ('search/summary.json','data/summary.json','inventory/candidate_entries.jsonl','run_receipt.json'):
        dst=broken/file; dst.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(octo/file,dst)
    put(plan,{**baseline,'windows':[spec(sept),spec(broken,octo/'inputs'/'query.json')]})
    g=read(broken/'data'/'summary.json'); search=read(broken/'search'/'summary.json'); old=search['entry_ids'][0]
    put(broken/'data'/'summary.json',{**g,'received':g['received']-1,'missing_ids':[old]})
    assert invoke('m1_holo_windows',plan,'--output',out,'--check-only')==2
    meta=read(out); assert meta['unique_ids']==2040 and meta['windows'][1]['metadata_missing_ids']==[old]
    other=read(sept/'search'/'summary.json')['entry_ids'][0]; search['entry_ids'][0]=other
    put(broken/'data'/'summary.json',g); put(broken/'search'/'summary.json',search)
    c=broken/'inventory'/'candidate_entries.jsonl'; c.write_text(c.read_text().replace('"pdb_id":"'+old+'"','"pdb_id":"'+other+'"',1))
    assert invoke('m1_holo_windows',plan,'--output',out,'--check-only')==2
    dup=read(out); assert dup['unique_ids']==2039 and set(dup['duplicate_ids'])=={other}
    return {'healthy_unique':2040,'healthy_dispatches':healthy['dispatched'],'missing_windows':miss['missing_windows'],
            'missing_window_dispatches':len(launches),'metadata_missing_id':old,'unique_with_metadata_gap':2040,
            'duplicate_id':other,'unique_with_duplicate':2039}


def m2_fixture(tmp):
    inputs=tmp/'m2-inputs'; run=tmp/'m2'; (run/'raw').mkdir(parents=True); (run/'source_proofs').mkdir()
    picked=next(n for n in read(pilot/'bootstrap'/'m2_nominations.json') if n['pdb_id']=='7BOC')
    put(inputs/'m2_nominations.json',[picked]); proof_file=pilot/'source_proofs'/'7BOC.json'; proof=read(proof_file)
    raw=proof['relative_path']; shutil.copyfile(pilot/raw,run/raw); shutil.copyfile(proof_file,run/'source_proofs'/proof_file.name)
    args=('holo_m2_pilot',run,'--inputs',inputs,'--plinder-release',release,'--git-head','isolated-acceptance',
          '--batch-id','holo-m2-resume-fixture','--smoke','7BOC')
    with patch.object(urllib.request,'urlopen',side_effect=AssertionError('coordinate download attempted')) as network:
        assert invoke(*args)==0
        q=read(run/'qc'/'7BOC.json'); digest=sha256(run/'qc'/'7BOC.json')
        assert q['source']['sha256']==proof['sha256'] and q['source']['entry_revision']=='pass'
        assert q['status']=='REVIEW' and q['source_observations'] and set(q['native_roles'])=={'pending'}
        assert invoke(*args)==0 and sha256(run/'qc'/'7BOC.json')==digest and network.call_count==0
    historical=read(pilot/'run_receipt.json')
    assert all(sha256(pilot/name)==digest for name,digest in historical['output_sha256'].items())
    return {'historical_m2_receipt_status':historical['status'],'planned':1,'smoke_ids':['7BOC'],'pinned_revision':q['revision'],'raw_sha256':proof['sha256'],
            'role_rules_present':False,'roles':q['native_roles'],'status':q['status'],'qc_sha256_after_reuse':digest,'download_requests':0}


with TemporaryDirectory(prefix='aphcp-holo-execution-',dir=evidence.parent) as work:
    w=Path(work)
    result={'m1_recovery':m1_fixture(w),'m1_cross_window':cross_fixture(w),'m2':m2_fixture(w)}
    put(evidence,result)
    print(json.dumps(result),flush=True)

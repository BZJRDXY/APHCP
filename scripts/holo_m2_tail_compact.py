"""Seal even unusually large readback ledgers without copying shared Assembly evidence per system."""
import argparse
from collections import Counter
import hashlib
import json
import mmap
import os
from pathlib import Path

from holo_candidate_batch import atomic_json, now, sha256


def rd(p):
    return json.loads(p.read_bytes())


def cut_objects(path, shared, commit):
    with path.open('rb') as f, mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) as buf:
        cursor=0
        while (start:=buf.find(b'\n  {\n',cursor))>=0:
            end=buf.find(b'\n  }',start+4)+4
            assert end>start and buf[end:end+1] in (b',',b'\n')
            cursor=end; hit=buf.find(b'\n    "complete_assembly": {',start,end)
            if hit==-1:
                yield json.loads(buf[start+1:end]); continue
            a=buf.find(b'{',hit,end); b=buf.find(b'\n    }',a,end)+6
            assert b>a
            piece=memoryview(buf)[a:b]; digest=hashlib.sha256(piece).hexdigest()
            if digest not in shared:
                target=commit/'shared_assembly_receipts'/f'{digest}.json'; target.parent.mkdir(exist_ok=True)
                with target.open('wb') as out: out.write(piece)
                item=json.loads(bytes(piece)); proof={'status':item['status'],'shared_evidence':str(target.relative_to(commit)),'shared_evidence_sha256':digest}
                if item['status']=='PASS':
                    art=item['artifact']; proof['artifact']={k:art[k] for k in ('path','bytes','sha256')}; assert proof['artifact']['path'].startswith('assemblies/')
                else:
                    proof.update({k:item[k] for k in ('reason','detail') if k in item})
                shared[digest]=proof
            piece.release()
            prefix=buf[start+1:hit]; suffix=buf[b:end]; assert prefix.endswith(b',')
            if suffix.startswith(b','): compact=prefix+suffix[1:]
            else: compact=prefix[:-1]+suffix
            row=json.loads(compact); row['complete_assembly']=shared[digest]; yield row


def combine(run,p):
    tail=run/'parallel_tail'; shard_dirs=sorted((tail/'shards'/p).iterdir()); final=tail/'compact_work'/p; final.mkdir(parents=True,exist_ok=True)
    records=[]; shared={}; by_scope={}; base=None; primary=None; roster=None; index=None; n=len(shard_dirs); proofs=0
    for d in shard_dirs:
        rr=rd(d/'entry_results'/f'{p}.json'); disc=d/'discovery'/f'{p}.json'; ds=rd(disc); del ds['finished_utc']; digest=hashlib.sha256(json.dumps(ds,separators=(',',':'),sort_keys=True).encode()).hexdigest()
        if base is None:
            base=rr; primary=disc; roster={x['candidate_key']:i for i,x in enumerate(ds['candidate_selection']) if x['eligible_for_structure_review']}; index=digest
        assert digest==index and rr['source_sha256']==base['source_sha256'] and rr['status']!='TECHNICAL_FAILURE'
        del ds
        rows=list(cut_objects(d/'system_results'/f'{p}.json',shared,final)); assert len(rows)==rr['chosen_systems']
        for row in rows:
            key=row['candidate_key']; assert int(key[:16],16)%n==int(d.name)
            art=row.get('artifact'); q=row.get('complete_assembly'); scope=(row['assembly'],row['representative_model'])
            if art:
                for path,check in [(art['path'],art['sha256']),(row['mapping'],None)]:
                    f=d/path; dst=final/path; dst.parent.mkdir(parents=True,exist_ok=True)
                    if check: assert sha256(f)==check; proofs+=1
                    if not dst.exists(): os.link(f,dst)
            if q:
                if scope in by_scope: assert by_scope[scope]==q
                by_scope[scope]=q
                if q['status']=='PASS':
                    f=d/q['artifact']['path']; dst=final/q['artifact']['path']; dst.parent.mkdir(parents=True,exist_ok=True)
                    if not dst.exists():
                        assert sha256(f)==q['artifact']['sha256']; os.link(f,dst)
            records.append(row)
        print('checked',p,d.name,len(records),flush=True)
    assert len(records)==len({s['candidate_key'] for s in records})==len(roster) and {s['candidate_key'] for s in records}==set(roster)
    for s in records:
        if s.pop('shard_assembly_reference',None):
            scope=(s['assembly'],s['representative_model']); s['complete_assembly']=by_scope[scope]
    records.sort(key=lambda x:roster[x['candidate_key']]); result=final/'system_results'/f'{p}.json'; discovery=final/'discovery'/f'{p}.json'
    discovery.parent.mkdir(exist_ok=True)
    if not discovery.exists(): os.link(primary,discovery)
    atomic_json(result,records)
    verified=sum(x['status']=='EXPORTED_READBACK' for x in records); positive=sum(x['candidate_intake']=='PASS' for x in records)
    base.update(status='HAS_QC_PASS_CANDIDATE' if positive else 'REVIEW',chosen_systems=len(records),structure_readback_pass=verified,locally_supported_candidates=positive,
        system_status_counts=dict(Counter(x['status'] for x in records)),discovery_sha256=sha256(discovery),systems_sha256=sha256(result),
        finished_utc=now(),candidate_partition={'slices':n,'shared_assembly_evidence_saved_once':len(shared),'shard_files_preserved':True},
        review_reason=None if positive else 'Candidate roles or topology require review.')
    atomic_json(final/'entry_results'/f'{p}.json',base)
    atomic_json(final/'receipt.json',{'entry':p,'slices':n,'eligible':len(records),'readback':verified,'source_sha256':base['source_sha256'],
        'identical_discovery_across_slices':True,'hash_verified_artifacts':proofs,'shared_assembly_proofs':len(shared),'bytes_compact_system_ledger':result.stat().st_size})
    print(p,'complete with compact evidence',len(records),verified,len(shared),result.stat().st_size,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('run',type=Path); ap.add_argument('entry'); a=ap.parse_args(); combine(a.run.resolve(),a.entry)

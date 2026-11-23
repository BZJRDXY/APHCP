"""Run disjoint candidate slices for large entries using the unchanged first-pass discovery criteria."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import os
from pathlib import Path
import shutil
import time

from holo_candidate_batch import atomic_json, now, sha256
from holo_m2_blind_structure import one, read


def dispatch(run, tail, pin, index, signature, i, total):
    target=tail/'shards'/pin['pdb_id']/str(i); target.mkdir(parents=True,exist_ok=True)
    for label in ('raw','source_proofs'):
        s=target/label
        if not s.exists(): s.symlink_to(run/label,target_is_directory=True)
    return one(target,pin,index,signature,(i,total))


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('run',type=Path); ap.add_argument('entries',nargs='+'); ap.add_argument('--slices',type=int,default=10);a=ap.parse_args()
    r=a.run.resolve(); tail=r/'parallel_tail'; tail.mkdir(exist_ok=True); pins={x['pdb_id']:x for x in read(r/'frozen.json')['files']}; ix=read(r/'index_snapshot.json')['pdb']; sig=read(r/'build_inputs.json')['fingerprint']; start=time.monotonic()
    with ProcessPoolExecutor(max_workers=len(a.entries)*a.slices) as pool:
        fs={pool.submit(dispatch,r,tail,pins[p],ix[p],sig,i,a.slices):(p,i) for p in a.entries for i in range(a.slices)}
        complete=defaultdict(list)
        for f in as_completed(fs):
            p,i=fs[f]; result=f.result(); print(p,i,result,round(time.monotonic()-start,2),flush=True)
            if result[1]=='built' or result[1]=='reused_verified': complete[p].append(i)
            if len(complete[p])==a.slices: assemble(r,tail,p,a.slices,sig)
    print('slice worker pool finished',round(time.monotonic()-start,2),flush=True)


def assemble(r, tail, p, total, sig):
    done=tail/'assembled'/p; done.mkdir(parents=True,exist_ok=True); reviewed=[]; assemblies={}; contexts=set(); s0=None; ledger=None; base=None; digest=None
    for i in range(total):
        shard=tail/'shards'/p/str(i); location=shard/'discovery'/f'{p}.json'; d=read(location); d.pop('finished_utc'); k=read(shard/'entry_results'/f'{p}.json'); rows=read(shard/'system_results'/f'{p}.json')
        if s0 is None: s0=location; ledger=d; base=k; digest=k['source_sha256']
        assert d==ledger and k['source_sha256']==digest and k['signature']==sig and k['status']!='TECHNICAL_FAILURE'
        for row in rows:
            key=row['candidate_key']; assert int(key[:16],16)%total==i
            art=row.get('artifact'); q=row.get('complete_assembly'); context=(row['assembly'],row['representative_model'])
            if art:
                for name,expected in ((art['path'],art['sha256']),(row['mapping'],None)):
                    src=shard/name; dst=done/name; dst.parent.mkdir(parents=True,exist_ok=True)
                    if expected: assert sha256(src)==expected
                    if not dst.exists(): os.link(src,dst)
                if q or row.get('shard_assembly_reference'): contexts.add(context)
            row.pop('shard_assembly_reference',None)
            if q:
                assert context not in assemblies or assemblies[context]==q, 'Assembly evidence differs across slices'
                assemblies[context]=q
                if q['status']=='PASS':
                    for name,expected in ((q['artifact']['path'],q['artifact']['sha256']),
                                          (q['shared_evidence'],q['shared_evidence_sha256'])):
                        src=shard/name; dst=done/name; dst.parent.mkdir(parents=True,exist_ok=True)
                        assert sha256(src)==expected
                        if dst.exists(): assert sha256(dst)==expected
                        else: os.link(src,dst)
            reviewed.append(row)
    assert contexts<=set(assemblies), ('no elected Assembly writer',contexts,set(assemblies))
    for row in reviewed:
        ctx=(row['assembly'],row['representative_model']); a=assemblies.get(ctx)
        if a and row.get('artifact'):
            row['complete_assembly']=a
    ranks={x['candidate_key']:i for i,x in enumerate(ledger['candidate_selection']) if x['eligible_for_structure_review']}
    assert len(reviewed)==len({s['candidate_key'] for s in reviewed})==len(ranks) and {s['candidate_key'] for s in reviewed}==set(ranks)
    reviewed.sort(key=lambda x:ranks[x['candidate_key']]); system=done/'system_results'/f'{p}.json'; discover=done/'discovery'/f'{p}.json'; entry=done/'entry_results'/f'{p}.json'
    discover.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(s0,discover)
    atomic_json(system,reviewed); status=sum(x['candidate_intake']=='PASS' for x in reviewed)
    base.update(status='HAS_QC_PASS_CANDIDATE' if status else 'REVIEW',chosen_systems=len(reviewed),
        structure_readback_pass=sum(x['status']=='EXPORTED_READBACK' for x in reviewed),locally_supported_candidates=status,
        system_status_counts=dict(Counter(x['status'] for x in reviewed)),discovery_sha256=sha256(discover),systems_sha256=sha256(system),
        finished_utc=now(),candidate_partition={'slices':total,'materialized_in':str(done)},
        review_reason=None if status else 'Candidate roles or topology require review.')
    atomic_json(entry,base); atomic_json(done/'receipt.json',{'entry':p,'slices':total,'eligible':len(ranks),'readback':base['structure_readback_pass'],
        'same_discovery_and_input':True,'source_sha256':digest,'assembled_utc':now(),'parent_core_sha256':sig})
    print('assembled',p,len(ranks),base['structure_readback_pass'],flush=True)


def promote(run, tail, p, seconds):
    """Publish verified slice outputs; the Entry receipt is the final write."""
    done=tail/'assembled'/p
    for folder in ('systems','assemblies'):
        for src in (done/folder).rglob('*'):
            if not src.is_file(): continue
            dst=run/src.relative_to(done); dst.parent.mkdir(parents=True,exist_ok=True)
            if dst.exists():
                if sha256(src)!=sha256(dst): raise ValueError('Slice promotion found different existing artifact')
            else: os.link(src,dst)
    for folder in ('discovery','system_results'):
        src=done/folder/f'{p}.json'; dst=run/folder/src.name
        dst.parent.mkdir(parents=True,exist_ok=True)
        # Copy the exact bytes that the merged Entry receipt hashes.
        temp=dst.with_suffix('.promoting'); shutil.copyfile(src,temp); os.replace(temp,dst)
    entry=read(done/'entry_results'/f'{p}.json')
    entry['candidate_partition']['dispatch_to_promotion_s']=seconds
    atomic_json(run/'entry_results'/f'{p}.json',entry)


if __name__=='__main__': main()

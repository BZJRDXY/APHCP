"""Adopt completed disjoint candidate slices only after deterministic overlap verification."""
import argparse
import json
import os
from pathlib import Path

from holo_candidate_batch import atomic_json, now, sha256
from holo_m2_blind_structure import summarize


def load(p):
    return json.loads(p.read_text())


def chosen(r,p):
    paths=[r/'parallel_tail'/v/p for v in ('compact_work','assembled')]
    return next(p for p in paths if (p/'receipt.json').is_file())

def materialized(r,p):
    d=chosen(r,p)
    assert (d/'receipt.json').is_file()
    artifacts=load(d/'system_results'/f'{p}.json'); by_key={a['candidate_key']:a for a in artifacts}; overlap=0
    for f in (r/'systems'/p).glob('*/system.cif'):
        companion=f.with_name('source_mapping.json')
        if companion.is_file():
            prior=load(companion); new=by_key[f.parent.name]['artifact']; assert prior['source_sha256']==load(d/'receipt.json')['source_sha256']
            assert prior['artifact']['sha256']==new['sha256']==sha256(f), (p,f.parent.name,'changed structure bytes')
            overlap+=1
    assert overlap, f'No structures to verify across the original and sliced paths: {p}'
    return overlap


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('run',type=Path); ap.add_argument('--adopt',action='store_true'); ap.add_argument('entries',nargs='+'); args=ap.parse_args()
    r=args.run.resolve(); plan=load(r/'frozen.json')['files']; missing={s['pdb_id'] for s in plan if not (r/'entry_results'/f'{s["pdb_id"]}.json').exists()}
    if not missing:
        print('No rows missing: retained the complete serial outputs.'); return
    assert missing<=set(args.entries), ('unpartitioned active entries',missing)
    counts={p:materialized(r,p) for p in sorted(missing)}
    print('deterministic overlap matched:',counts,flush=True)
    if not args.adopt: return
    interrupt=r/'parallel_tail'/'serial_partial'; interrupt.mkdir(exist_ok=True); moved=[]
    for p in sorted(missing):
        src=chosen(r,p)
        for folder in ('systems','assemblies','discovery','system_results','entry_results','shared_assembly_receipts'):
            for f in list((src/folder).iterdir()) if (src/folder).exists() else ():
                dst=r/folder/f.name; dst.parent.mkdir(exist_ok=True); old=interrupt/folder/f.name
                if dst.exists():
                    old.parent.mkdir(exist_ok=True)
                    assert not old.exists(); dst.rename(old); moved.append(str(old.relative_to(r)))
                if f.is_dir():
                    dst.mkdir(exist_ok=True)
                    for item in f.rglob('*'):
                        to=dst/item.relative_to(f)
                        if item.is_dir(): to.mkdir(exist_ok=True)
                        else: os.link(item,to)
                else: os.link(f,dst)
    fetch=load(r/'last_fetch_invocation.json')
    atomic_json(r/'parallel_tail'/'adoption.json',{'utc':now(),'entries':sorted(missing),'verified_matching_serial_artifacts':counts,
             'preserved_partial_paths':moved,'recovered_entries':len(missing),'sharding_code_sha256':sha256(r/'parallel_tail/bootstrap/holo_m2_blind_structure.py'),
             'method':'Keep every unique candidate, reuse the same source/QC logic, share the elected Assembly export.'})
    summarize(r,plan,{'finished_utc':now(),'fetch':fetch['counts'],'index_snapshot':'first_pass_existing','structure':{'built':len(plan)-len(missing),'built_in_slices':len(missing)},
                      'whole_wall_time_basis':'See first launch, resumed transport and tail adoption timestamps','tail_adoption_receipt':'parallel_tail/adoption.json'})


if __name__=='__main__': main()

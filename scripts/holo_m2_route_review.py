"""Close source-route gaps in the sealed sample and verify selected M3 readbacks."""
import argparse
from collections import Counter
from pathlib import Path

from holo_candidate_batch import atomic_json, now, sha256
from holo_m2_blind_audit import inspect_output, rd
from holo_m2_final_verify import table


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('run',type=Path); a=ap.parse_args(); r=a.run.resolve()
    m3=rd(r/'m3_first_review.json'); audit=rd(r/'independent_audit.json'); snapshot=rd(r/'blind_first_pass.json'); selection=table(r/'system_results.csv')
    known={x['pdb_id'] for x in audit['entries']}; route=set(); cases=[]; results=[]
    for s in sorted(selection,key=lambda s:(s.get('role_qc')!='PASS',s['pdb_id'],s['candidate_key'])):
        if s['status']!='EXPORTED_READBACK' or not s['system_path'] or s['pdb_id'] in known: continue
        p=s['pdb_id']; d=rd(r/'entry_results'/f'{p}.json'); k=d['route']; tags={k}
        if k!='NATIVE_ONLY_NO_ANNOTATION' or not set(tags)-route: continue
        known.add(p); route|=tags; cases.append((p,s['candidate_key'],k,'missing_route_readback_in_first_independent_sample'))
        if 'NATIVE_ONLY_NO_ANNOTATION' in route:break
    for s in m3['review_candidates'][:3]:
        p=s['pdb_id']
        if p not in known: cases.append((p,s['candidate'],rd(r/'entry_results'/f'{p}.json')['route'],'M3_preliminary_roster_direct_readback')); known.add(p)
    for p,key,route,why in cases:
        path=r/'system_results'/f'{p}.json'; entry=rd(r/'source_proofs'/f'{p}.json'); fingerprint=sha256(path)
        assert snapshot['entry_systems_sha256'][p]==fingerprint==sha256(path)
        data=next(x for x in rd(path) if x['candidate_key']==key)
        try:
            trace=inspect_output(r,data,r/entry['relative_path'])
            status='INDEPENDENT_READBACK_PASS'; err=None
        except (AssertionError,ValueError,KeyError) as e:
            status='AUDIT_CHECK_FAILED'; trace=None; err=repr(e)
        results.append({'pdb_id':p,'candidate_key':key,'route':route,'purpose':why,'readback':trace,'status':status,'error':err,'role_decision':'Not part of this readback'})
    atomic_json(r/'route_readback_checks.json',{'updated_utc':now(),'scope':'Supplementary tests outside the sealed 20-case diagnostic sample',
        'counts':dict(Counter(x['status'] for x in results)),'entries':results,'first_pass_sha256':sha256(r/'blind_first_pass.json')})
    print('Supplementary checks',[(x['pdb_id'],x['route'],x['status']) for x in results],flush=True)


if __name__=='__main__':main()

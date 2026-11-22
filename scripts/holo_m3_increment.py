"""Run fixed-source M2, classify a new batch, then review and stage its systems."""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import subprocess
import sys
import time

from holo_m2_chemtype import annotate_batch
from holo_m3_register import frozen_registry_root, require, read, sha, write, file_ref


def once(path, doc):
    if path.is_file():
        require(read(path)==doc, f'{path.name} changed; use a new downstream Run')
    else:
        write(path,doc)


def prior_inputs(m2, batch):
    frozen=read(m2/'frozen.json')
    entries=[]; candidates=[]; nominations=[]; evidence={}
    for spec in frozen['files']:
        p=spec['pdb_id']; revision=spec['revision']; entry_file=m2/'entry_results'/f'{p}.json'
        e=read(entry_file)
        require(e['pdb_id']==p and e['revision']==revision,'M2 Entry identity does not match frozen manifest')
        entries.append(p); evidence[p]={'entry_result_sha256':sha(entry_file),'state':e['status'],
                                        'route':e.get('route'),'stratum':e.get('stratum')}
        if e['status']=='TECHNICAL_FAILURE':
            continue
        discovery=m2/'discovery'/f'{p}.json'; systems=m2/'system_results'/f'{p}.json'
        file_ref(discovery,(m2,),e['discovery_sha256']);file_ref(systems,(m2,),e['systems_sha256'])
        proof_file=m2/'source_proofs'/f'{p}.json'; proof=read(proof_file)
        require(proof['status']==proof['gzip_crc']=='PASS' and proof['identity']==spec and
                proof['observed_revision']==revision and proof['sha256']==e['source_sha256'], 'M2 source proof mismatch')
        raw=file_ref(m2/proof['relative_path'],(m2,),proof['sha256'])
        for s in read(systems):
            key=s['candidate_key']
            require((s['pdb_id'],s['revision'])==(p,revision),'M2 system identity does not match source')
            b={'batch':batch,'pdb_id':p,'revision':revision,'source_path':raw['path'],'source_sha256':raw['sha256'],
               'candidate_key':key,'assembly':s['assembly'],'model':s['representative_model'],
               'member_instances':s['partner_instances'],'discovery_path':str(discovery),'discovery_sha256':e['discovery_sha256'],
               'route':e['route'],'role_after':s['role_qc'],'role_evidence_after':s['role_basis'],
               'prior_status':s['status'],'intake_after':s['candidate_intake'],'evidence_layer':'new_m2_run',
               **{f:s.get(f) for f in ('coordinate_identity','experimental_fidelity','topology','partner_scope_qc')}}
            if s['status']=='EXPORTED_READBACK':
                cf=file_ref(m2/s['artifact']['path'],(m2,),s['artifact']['sha256']); mp=file_ref(m2/s['mapping'],(m2,))
                b.update(system_path=cf['path'],system_sha256=cf['sha256'],mapping_path=mp['path'],mapping_sha256=mp['sha256'])
                if s['candidate_intake']!='EXCLUDED':
                    nominations.append({k:b[k] for k in ('batch','pdb_id','candidate_key','revision','assembly','source_sha256','system_sha256')} |
                                       {'source_model':b['model'],'frozen_source_proof':str(proof_file)})
            candidates.append(b)
    return {'batch':batch,'source_run':str(m2),'source_manifest_sha256':sha(m2/'frozen.json'),
            'entries':entries,'expected_entries':frozen['plan_count'],'expected_candidates':len(candidates),
            'entry_evidence':evidence,'candidates':candidates}, nominations


def run_pipeline(a):
    import fcntl
    start=time.monotonic(); root=a.root.resolve(strict=True); m2=a.m2_run.resolve(strict=True); m3=a.m3_run.resolve(strict=True)
    require(m2.parent==root/'runs' and m2.name.startswith('holo-m2-'),'not_an_isolated_m2_run')
    frozen_registry_root(root,m3)
    here=Path(__file__).resolve().parent
    with (m3/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        subprocess.run([sys.executable,'-B',str(here/'holo_m2_blind.py'),str(m2),'--index',str(a.index),
                        '--jobs',str(a.jobs),'--fetch-jobs',str(a.jobs)],check=True)
        evidence, nominations=prior_inputs(m2,a.batch)
        require(len(evidence['entries'])==a.expected_entries,'M2 entries do not match the approved input count')
        in_f=m3/'candidate_evidence.json'; once(in_f,evidence)
        out_f=m3/'chemical_labels.jsonl'; therm=annotate_batch(in_f,a.reference_run,out_f,a.jobs)
        ref=lambda path:{'path':str(path),'sha256':sha(path)}
        manifest=m3/'registration_input.json'
        once(manifest,{'schema':'aphcp-m3-increment-v1','candidate_evidence':ref(in_f),'chemistry':ref(out_f),
                       'review_count':len(nominations),'review_candidates':nominations})
        if nominations:
            cmd=[sys.executable,'-B',str(here/'holo_m3_register.py'),'register','--root',str(root),'--run',str(m3),
                 '--manifest',str(manifest),'--git-head',a.git_head,'--jobs',str(a.jobs)]
            if a.seed_run:cmd+=['--seed-run',str(a.seed_run)]
            subprocess.run(cmd,check=True)
        invocation={'entry_statuses':dict(Counter(x['state'] for x in evidence['entry_evidence'].values())),
                    'strata':dict(Counter(x['stratum'] for x in evidence['entry_evidence'].values())),
                    'candidates':len(evidence['candidates']),'nominated':len(nominations),'chemistry':therm,
                    'm3_results_sha256':sha(m3/'review_results.json') if nominations else None,
                    'elapsed_s':round(time.monotonic()-start,3)}
        write(m3/'increment_receipt.json',invocation)
        print(invocation,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('root','m2-run','m3-run','reference-run','index'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--seed-run',type=Path)
    parser.add_argument('--batch',required=True)
    parser.add_argument('--expected-entries',type=int,required=True)
    parser.add_argument('--git-head',required=True)
    parser.add_argument('--jobs',type=int,default=8)
    run_pipeline(parser.parse_args())

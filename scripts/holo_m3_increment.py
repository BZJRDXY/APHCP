"""Run fixed-source M2, classify a new batch, then review and stage its systems."""
from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path
import os
import subprocess
import sys
import time

from holo_m2_chemtype import annotate_batch
from holo_m3_register import frozen_registry_root, require, read, sha, write, file_ref
from holo_m3_observations import append_increment
from holo_archive_runtime import sync_path


def once(path, doc):
    if path.is_file():
        require(read(path)==doc, f'{path.name} changed; use a new downstream Run')
    else:
        write(path,doc)


def prior_inputs(m2, batch, selected=None):
    frozen=read(m2/'frozen.json')
    entries=[]; candidates=[]; nominations=[]; evidence={}
    for spec in frozen['files']:
        p=spec['pdb_id']; revision=spec['revision']; entry_file=m2/'entry_results'/f'{p}.json'
        if selected is not None and p not in selected:
            continue
        e=read(entry_file) if entry_file.is_file() else {'pdb_id':p, 'revision':revision, 'status':'TECHNICAL_FAILURE',
                                                           'error':'missing_worker_checkpoint'}
        require(e['pdb_id']==p and e['revision']==revision,'M2 Entry identity does not match frozen manifest')
        entries.append(p); evidence[p]={'entry_result_sha256':sha(entry_file) if entry_file.is_file() else None,'state':e['status'],
                                        'error':e.get('error'),'failure_stage':e.get('failure_stage'),
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
               'prior_status':s['status'],'export_reason':s.get('reason'),'intake_after':s['candidate_intake'],'evidence_layer':'new_m2_run',
               **{f:s.get(f) for f in ('coordinate_identity','experimental_fidelity','topology','partner_scope_qc')}}
            if s['status']=='EXPORTED_READBACK':
                cf=file_ref(m2/s['artifact']['path'],(m2,),s['artifact']['sha256']); mp=file_ref(m2/s['mapping'],(m2,))
                b.update(system_path=cf['path'],system_sha256=cf['sha256'],mapping_path=mp['path'],mapping_sha256=mp['sha256'])
                if s['candidate_intake']!='EXCLUDED':
                    nominations.append({k:b[k] for k in ('batch','pdb_id','candidate_key','revision','assembly','source_sha256','system_sha256')} |
                                       {'source_model':b['model'],'frozen_source_proof':str(proof_file)})
            candidates.append(b)
    return {'batch':batch,'source_run':str(m2),'source_manifest_sha256':sha(m2/'frozen.json'),
            'entries':entries,'expected_entries':len(entries),'expected_candidates':len(candidates),
            'entry_evidence':evidence,'candidates':candidates}, nominations


def owner_lease(a):
    return tuple(getattr(a, "controller_lease_fd", ()))


def submit_m2(a, m2):
    here = Path(__file__).resolve().parent
    cmd=[sys.executable,'-B',str(here/'holo_m2_blind.py'),str(m2),'--index',str(a.index),
         '--jobs',str(a.jobs),'--fetch-jobs',str(a.fetch_jobs)]
    if getattr(a,'archive_slots',None):
        cmd+=['--archive-slots',str(a.archive_slots),'--archive-limit',str(a.archive_limit),
              '--archive-ancillary',str(a.archive_ancillary),'--primary-budget',str(a.primary_budget)]
    if getattr(a,'stream_build',False):
        cmd.append('--stream-build')
    subprocess.run(cmd,check=True,pass_fds=owner_lease(a))


def link_cached_raw(previous, target, pdb):
    from holo_m2_blind import raw_path

    proof = previous/'source_proofs'/f'{pdb}.json'
    if proof.is_file():
        saved = read(proof)
        if saved['status'] != 'PASS':
            return
        source = file_ref(previous/saved['relative_path'], (previous,), saved['sha256'])
        raw = raw_path(target, saved['identity'])
        raw.parent.mkdir(exist_ok=True)
        if not raw.exists():
            os.link(source['path'], raw)


def prepare_retry(m2, base, expected, previous_by_pdb):
    m2.mkdir(exist_ok=True)
    inputs = m2/'inputs'; inputs.mkdir(exist_ok=True)
    for name in ('blind.json', 'archive_head.json'):
        rows = [row for row in read(base/'inputs'/name) if row['pdb_id'] in expected]
        require({row['pdb_id'] for row in rows} == set(expected), 'retry_inputs_not_in_base_freeze')
        once(inputs/name, rows)
    for p in expected:
        link_cached_raw(previous_by_pdb[p], m2, p)


def stage(a, root, run, m2, seed, keep, previous_successes):
    started = time.monotonic(); evidence, nominations = prior_inputs(m2,a.batch,keep)
    require(set(evidence['entries']) == set(keep), 'stage_entries_do_not_match_plan')
    # Keep prior delivered candidates and their asset links. Reconsider the retry delta only.
    retry_catalogued = [(r['pdb_id'],r['candidate_key']) for r in evidence['candidates']
                        if (r['pdb_id'],r['candidate_key']) in previous_successes]
    evidence['candidates'] = [r for r in evidence['candidates'] if (r['pdb_id'],r['candidate_key']) not in previous_successes]
    evidence['expected_candidates'] = len(evidence['candidates'])
    nominations = [r for r in nominations if (r['pdb_id'],r['candidate_key']) not in previous_successes]
    in_f = run/'candidate_evidence.json'; once(in_f,evidence)
    out_f = run/'chemical_labels.jsonl'; chem = annotate_batch(in_f,a.reference_run,out_f,a.jobs)
    manifest = run/'registration_input.json'
    once(manifest,{'schema':'aphcp-m3-increment-v1','candidate_evidence':file_ref(in_f,(run,)),
                   'chemistry':file_ref(out_f,(run,)), 'review_count':len(nominations),'review_candidates':nominations})
    if nominations:
        here = Path(__file__).resolve().parent
        subprocess.run([sys.executable,'-B',str(here/'holo_m3_register.py'),'register','--root',str(root),
                        '--run',str(run),'--manifest',str(manifest),'--git-head',a.git_head,'--jobs',str(a.jobs),
                        '--seed-run',str(seed)],check=True,pass_fds=owner_lease(a))
    catalog = append_increment(root,run,seed,in_f,out_f,a.git_head)
    review_file = run/'review_results.json'
    reviews = read(review_file) if review_file.is_file() else []
    m2_retry = {p for p, row in evidence['entry_evidence'].items() if row['state']=='TECHNICAL_FAILURE'} | {
        r['pdb_id'] for r in evidence['candidates'] if r['prior_status']=='EXPORT_TECHNICAL_FAILURE'}
    m3_retry = {r['pdb_id'] for r in reviews if r['disposition']=='TECHNICAL_FAILURE'}
    failed_keys = {(r['pdb_id'],r['candidate_key']) for r in reviews if r['disposition']=='TECHNICAL_FAILURE'}
    settled = [(r['pdb_id'],r['candidate_key']) for r in evidence['candidates'] if
               r['prior_status']!='EXPORT_TECHNICAL_FAILURE' and (r['pdb_id'],r['candidate_key']) not in failed_keys]
    receipt = {'selected':sorted(keep), 'm2_run':str(m2), 'registry_run':str(run), 'previous_catalog_run':str(seed),
               'entry_states':dict(Counter(x['state'] for x in evidence['entry_evidence'].values())),
               'routes':dict(Counter(x['stratum'] for x in evidence['entry_evidence'].values())),
               'structural_states':dict(Counter(x['prior_status'] for x in evidence['candidates'])),
               'review_dispositions':dict(Counter(r['disposition'] for r in reviews)), 'chemistry':chem,
               'catalog':catalog, 'm2_pending':sorted(m2_retry), 'm3_pending':sorted(m3_retry), 'settled_keys':settled,
               'skipped_prior_successes':len(retry_catalogued), 'elapsed_s':round(time.monotonic()-started,3)}
    # The catalog transaction has returned and closed before advertising its receipt.
    registry=Path(catalog['registry']); sync_path(registry); sync_path(registry.parent)
    sync_path(registry.parent.parent); sync_path(registry.parent.parent.parent); sync_path(run.parent)
    write(run/'phase_receipt.json', receipt)
    return receipt


def run_pipeline(a):
    import fcntl
    from holo_m2_blind import freeze

    start=time.monotonic(); root=a.root.resolve(strict=True); m2=a.m2_run.resolve(strict=True); m3=a.m3_run.resolve(strict=True)
    require(m2.parent==root/'runs' and m2.name.startswith('holo-m2-'),'not_an_isolated_m2_run')
    frozen_registry_root(root,m3)
    with (m3/'execute.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        plan = read(m2/'frozen.json')['files'] if a.adopt_m2 else freeze(m2)
        require(len(plan)==a.expected_entries,'M2 entries do not match the approved input count')
        once(m3/'increment_plan.json', {'m2_run':str(m2),'m2_manifest_sha256':sha(m2/'frozen.json'),
            'adopt_m2':a.adopt_m2,'entries':len(plan),'batch':a.batch,'seed_run':str(a.seed_run.resolve(strict=True)),
            'reference_run':str(a.reference_run.resolve(strict=True)),'index':str(a.index.resolve()),'git_head':a.git_head,
            **({'stream_build':True} if getattr(a,'stream_build',False) else {})})
        seed=a.seed_run.resolve(strict=True); pending={x['pdb_id'] for x in plan}; cumulative=set(); receipts=[]
        before={}; m2_pending=pending.copy(); m3_pending=set(); executed=[]
        while pending:
            index=len(receipts); label = f'{m3.name}-phase-{index+1:03d}'
            target=root/'runs'/label; target.mkdir(exist_ok=True)
            done=target/'phase_receipt.json'
            if done.is_file():
                result=read(done)
                require(set(result['selected'])<=pending and result['previous_catalog_run']==str(seed), 'phase_sequence_changed')
            else:
                if (index>0 and getattr(a,'defer_pending',False)) or len(executed)>=a.passes:
                    break
                if index==0:
                    selected=pending; source=m2; process=not a.adopt_m2
                elif m2_pending:
                    selected=m2_pending; source=m2.with_name(f'{m2.name}-retry-{index:03d}')
                    process=True; prepare_retry(source,m2,selected,before)
                else:
                    # Retry only registration/reading; keep the original M2 evidence frozen.
                    source=before[sorted(m3_pending)[0]]
                    selected={p for p in m3_pending if before[p]==source}; process=False
                phase_plan=target/'phase_plan.json'
                once(phase_plan, {'selected':sorted(selected),'previous_catalog_run':str(seed),
                                 'm2_run':str(source),'process_m2':process})
                if process and not (target/'candidate_evidence.json').exists():
                    submit_m2(a,source)
                result=stage(a,root,target,source,seed,selected,cumulative)
                executed.append(label)
            pending -= set(result['selected']); m2_pending -= set(result['selected']); m3_pending -= set(result['selected'])
            m2_pending.update(result['m2_pending']); m3_pending.update(result['m3_pending'])
            pending.update(m2_pending | m3_pending)
            before.update({p:Path(result['m2_run']) for p in result['selected']})
            cumulative.update(map(tuple,result['settled_keys']))
            receipts.append(result); seed=target
            # Publish the head only after both the SQLite commit and the phase receipt exist.
            write(m3/'catalog_head.json',{'catalog_run':seed.name,'phase':len(receipts)})
        summary={'catalog_run':str(seed),'frozen_entries':len(plan),'phases':[r['registry_run'] for r in receipts],
                 'executed_phases':executed,'reused_phases':len(receipts)-len(executed),'pending_technical':sorted(pending),
                 'status':'STAGED_PENDING_RETRY' if pending else 'STAGED_COMPLETE','elapsed_s':round(time.monotonic()-start,3)}
        write(m3/'increment_receipt.json',summary);print(summary,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('root','m2-run','m3-run','reference-run','index'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--seed-run',type=Path,required=True)
    parser.add_argument('--adopt-m2',action='store_true',help='Use an existing frozen M2 Run without writing to it')
    parser.add_argument('--passes',type=int,default=2,help='Maximum new phases per invocation')
    parser.add_argument('--fetch-jobs',type=int,default=16)
    parser.add_argument('--defer-pending',action='store_true',help='Replay sealed phases but never start a new retry phase')
    parser.add_argument('--controller-lease-fd',type=int,action='append',default=[],help=argparse.SUPPRESS)
    parser.add_argument('--archive-slots',type=Path); parser.add_argument('--archive-limit',type=int,default=24)
    parser.add_argument('--archive-ancillary',type=int,default=6); parser.add_argument('--primary-budget',type=float,default=30)
    parser.add_argument('--stream-build',action='store_true')
    parser.add_argument('--batch',required=True)
    parser.add_argument('--expected-entries',type=int,required=True)
    parser.add_argument('--git-head',required=True)
    parser.add_argument('--jobs',type=int,default=8)
    args=parser.parse_args()
    if args.stream_build and not args.archive_slots:
        parser.error('--stream-build requires --archive-slots')
    if min(args.jobs,args.fetch_jobs,args.passes,args.expected_entries)<1:
        parser.error('positive counts are required')
    run_pipeline(args)

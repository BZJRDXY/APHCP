"""Check the frozen rolling denominator and catalog independently of exit codes."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

from holo_m3_increment import once
from holo_m3_observations import STAGE, lookup_records
from holo_m3_register import connect, database, digest, read, require, sha
from holo_m3_rolling import SPEC, archived


def run_path(runs, name):
    require(isinstance(name,str) and Path(name).name==name,'terminal_catalog_path_invalid')
    path=runs/name
    require(path.resolve()==path and path.is_dir(),'terminal_run_missing_or_redirected')
    return path


def frozen_inputs(run, expected_sha):
    path=run/'rolling_plan.json'; require(sha(path)==expected_sha,'terminal_plan_hash_differs')
    plan=read(path); provenance=[run]; runs=run.parent
    while parent:=read(provenance[-1]/'rolling_config.json').get('continued_from'):
        ancestor=run_path(runs,Path(parent).name)
        require(str(ancestor)==parent and ancestor not in provenance and sha(ancestor/'rolling_plan.json')==expected_sha,
                'terminal_predecessor_plan_differs')
        provenance.append(ancestor)
    parts={}; seen=set()
    for n,part in enumerate(plan['partitions'],1):
        bid=part['id']; require(bid==f'b{n:04d}','terminal_partition_ordinals_differ')
        file=next((d/'partitions'/f'{bid}.json' for d in provenance if (d/'partitions'/f'{bid}.json').is_file()),None)
        require(file is not None and sha(file)==part['sha256'],'terminal_partition_hash_differs')
        pins=[row['source'] for row in read(file)]
        by_pdb={pin['pdb_id']:pin for pin in pins}
        require(len(by_pdb)==len(pins)==part['entries'] and not seen.intersection(by_pdb) and all(map(archived,pins)),
                'terminal_duplicate_or_unpinned_entry')
        parts[bid]=by_pdb; seen.update(by_pdb)
    excluded=next((a/'excluded.json' for a in provenance if (a/'excluded.json').is_file()),None)
    require(excluded is not None and sha(excluded)==plan['excluded_sha256'],'terminal_exclusion_hash_differs')
    exclusions=read(excluded); ids={x['pdb_id'] for x in exclusions}
    require(len(ids)==len(exclusions) and not seen.intersection(ids) and len(ids)+len(seen)==plan['denominator'] and
            Counter(x['status'] for x in exclusions)==Counter(plan['exclusions']), 'terminal_m1_denominator_differs')
    return parts, provenance, plan['exclusions']


def batch_proof(runs, saved, pins, selected, pending):
    old=saved['previous_catalog']; new=saved['catalog_run']
    if not saved.get('m3_run'):
        require(old==new and selected==pending==set(saved['archive_deferred']),'terminal_unbuilt_batch_evidence_differs')
        return
    deferred=set(saved['archive_deferred']); planned=selected-deferred
    require(deferred<=pending<=selected and planned, 'terminal_archive_defer_evidence_differs')
    seq=[]; seen=set()
    while new!=old:
        require(new not in seen,'terminal_phase_chain_cycle'); seen.add(new)
        path=run_path(runs,new); ph=read(path/'phase_receipt.json')
        require(str(path)==ph['registry_run'],'terminal_phase_identity_differs')
        previous=Path(ph['previous_catalog_run']); new=previous.name
        require(previous==runs/new,'terminal_phase_seed_outside_runs')
        seq.append(ph)
    require(seq and Path(seq[-1]['previous_catalog_run'])==runs/old,'terminal_catalog_chain_broken')
    current=planned; delivered=set()
    for ph in reversed(seq):
        ids=set(ph['selected']); unresolved=set(ph['m2_pending'])|set(ph['m3_pending'])
        require(len(ids)==len(ph['selected']) and unresolved<=ids<=current,'terminal_phase_selection_differs')
        src=run_path(runs,Path(ph['m2_run']).name)
        require(str(src)==ph['m2_run'],'terminal_phase_source_outside_runs')
        rows=read(src/'frozen.json')['files']; evidence={x['pdb_id']:x for x in rows}
        require(len(evidence)==len(rows) and all(p in evidence and all(evidence[p][k]==pins[p][k] for k in SPEC) for p in ids),
                'terminal_phase_fixed_version_differs')
        delivered.update(ids-unresolved); current=(current-ids)|unresolved
    require(delivered==selected-pending and current|deferred==pending,'terminal_per_entry_disposition_differs')


def audit_keys(run, expected_sha, stop, max_retries):
    parts, ancestry, excluded=frozen_inputs(run,expected_sha)
    require(1<=stop<=len(parts),'terminal_stop_outside_plan')
    required=set(list(parts)[:stop]); last=f'b{stop:04d}'
    runs=run.parent; config=read(run/'rolling_config.json'); previous=Path(config['seed_run']).name
    require(str(run_path(runs,previous))==config['seed_run'],'terminal_seed_differs')
    journal=sorted(((f.stem,read(f),f) for f in (run/'receipts').glob('*.json')),key=lambda x:x[1]['sequence'])
    prefix=state_snapshot(run/'active_dispatch.json')
    require(not prefix,'terminal_active_dispatch_remaining')
    require(not any(re.match(r'b\d{4}', f.name) and int(f.name[1:5])>stop
                    for f in (run/'preflight').glob('b*')),'terminal_dispatched_beyond_stop')
    child_stem=run.name.removeprefix('holo-m3-')
    require(not any((m:=re.match(r'-b(\d{4})(?:$|-)', f.name[len(stem):])) and int(m[1])>stop
                    for stem in ('holo-m2-'+child_stem,run.name) for f in runs.glob(stem+'-b*')),
            'terminal_child_run_beyond_stop')
    h=hashlib.sha256(); retry_round=Counter(); primary=set(); pending={}; consumed={}
    for n,(label,report,path) in enumerate(journal,1):
        bid,sep,attempt=label.partition('-r')
        require(bid in required and report['sequence']==n and report['previous_catalog']==previous,
                'terminal_receipt_chain_differs')
        pins=parts[bid]; unresolved=set(report['pending'])
        require(len(unresolved)==len(report['pending']),'terminal_duplicate_pending_key')
        if sep:
            retry_round[bid]+=1
            require(bid in primary and attempt==f'{retry_round[bid]:02d}','terminal_retry_order_differs')
            listed=report['retry_selected']; selected=set(listed)
            require(len(listed)==len(selected) and selected<=pending[bid],'terminal_retry_selection_differs')
            if attempts:=report.get('retry_attempts'):
                require(set(attempts)==selected and all(type(v) is int and v>consumed.get((bid,p),0) for p,v in attempts.items()),
                        'terminal_retry_attempt_evidence_differs')
                consumed.update({(bid,p):v for p,v in attempts.items()})
        else:
            require(bid not in primary,'terminal_duplicate_primary_receipt')
            primary.add(bid); selected=set(pins); pending[bid]=set()
        require(unresolved<=selected and len(selected)==report['total'] and len(selected-unresolved)==report['settled'],
                'terminal_receipt_denominator_differs')
        batch_proof(runs,report,pins,selected,unresolved)
        pending[bid]=(pending[bid]-selected)|unresolved
        previous=report['catalog_run']; h.update(f'{label}:{sha(path)}\n'.encode())
    require(primary==required and last in primary,'terminal_unattempted_partitions')
    require(read(run/'catalog_head.json')=={'catalog_run':previous},'terminal_catalog_head_differs')
    retry=read(run/'retry_progress.json')
    require(all(retry[k]==0 for k in ('running','retryable_not_ready','verified_waiting_for_catalog')),
            'terminal_actionable_retries_remaining')
    historical=ancestry[1] if len(ancestry)>1 else None
    old_ids=set()
    if historical:
        baseline=read(historical/'rolling_progress.json')['unresolved_by_partition']
        old_ids={(bid,pdb) for bid,items in baseline.items() for pdb in items}
    exceptions=[]
    for bid in sorted(pending):
        for pdb in sorted(pending[bid]):
            pin=parts[bid][pdb]; key=digest(tuple(pin[k] for k in SPEC)); rel=Path('retry_cache')/key/'attempt.json'
            f=next((p/rel for p in ancestry if (p/rel).is_file()),None)
            require(f is not None,'terminal_exhausted_attempt_missing')
            state=read(f); attempts=state['attempts']; phase=state['state']
            require(attempts>=max_retries and (phase=='FAILED' or (phase=='READY' and attempts==consumed.get((bid,pdb)))) and
                    all(state['pin'][k]==pin[k] for k in SPEC),'terminal_entry_is_not_proven_exhausted')
            exceptions.append({'partition':bid,'pdb_id':pdb,'revision':pin['revision'],'attempts':attempts,
                               'last_state':phase,'last_error':state.get('check',{}).get('error'),
                               'checkpoint_sha256':sha(f),'checkpoint_run':f.parents[2].name,
                               'historical':(bid,pdb) in old_ids})
    require(len(exceptions)==retry['exhausted']==retry['queued'],'terminal_retry_denominator_differs')
    denominator=sum(len(parts[p]) for p in required)
    return {'absolute_stop':stop,'partitions':len(primary),'rolling_entries':denominator,
            'settled':denominator-len(exceptions),'exhausted':exceptions,'unattempted':0,
            'historical_exhausted':sum(x['historical'] for x in exceptions),
            'new_exhausted':sum(not x['historical'] for x in exceptions),
            'catalog_run':previous,'journal_sha256':h.hexdigest(),'receipt_count':len(journal),
            'frozen_plan_sha256':expected_sha,'m1_exclusions':excluded}


def state_snapshot(path):
    return read(path) if path.is_file() else None


def audit_catalog(root, run, result):
    catalog=run_path(run.parent,result['catalog_run']); dbpath=database(root/'batches'/catalog.name)
    catalog_evidence=read(catalog/'catalog_receipt.json'); frozen=read(catalog/'frozen_catalog.json')
    with connect(dbpath,readonly=True) as db:
        db.execute('PRAGMA mmap_size=536870912')
        require(db.execute('PRAGMA integrity_check').fetchone()[0]=='ok' and not db.execute('PRAGMA foreign_key_check').fetchone(),
                'terminal_catalog_integrity_failed')
        require(db.execute('SELECT COUNT(*) FROM release').fetchone()[0]==0,'terminal_catalog_unexpected_release')
        proof=db.execute('SELECT config_hash,error_summary FROM run WHERE id=?',(catalog.name+'-observations',)).fetchone()
        require(proof and proof[0]==digest(frozen) and json.loads(proof[1])==catalog_evidence,'terminal_catalog_commit_differs')
        examples=[]
        for status in ('REGISTERED','REVIEW'):
            for order in ('ASC','DESC'):
                record=db.execute(f'''SELECT object_id FROM run_item WHERE stage=? AND status=? AND output_uri IS NOT NULL
                                    ORDER BY updated_at {order}, object_id LIMIT 1''',(STAGE,status)).fetchone()
                if record and record[0] not in examples: examples.append(record[0])
    require(examples,'terminal_no_retrievable_catalog_examples')
    checked=[]
    for key in examples:
        rows=lookup_records(catalog,system=key,verify_files=True)
        require(len(rows)==1 and {'system_cif','source_mapping','source_mmcif'}<=set(rows[0].get('files',{})),
                'terminal_files_not_retrievable')
        checked.append({'system_id':key,'status':rows[0]['status'],'assets':sorted(rows[0]['files'])})
    return checked


def accept_terminal(root, run, expected_sha, stop, max_retries, *, publish=True):
    result=audit_keys(run,expected_sha,stop,max_retries)
    if publish:
        require(stop==len(read(run/'rolling_plan.json')['partitions']),'terminal_not_absolute_end')
    result.update(accepted=True,files=audit_catalog(root,run,result))
    if publish: once(run/'terminal_receipt.json',result)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True); p.add_argument('--plan-sha256',required=True)
    p.add_argument('--stop',type=int,required=True); p.add_argument('--max-retries',type=int,default=3)
    a=p.parse_args(); run=a.run.resolve(strict=True)
    print(json.dumps(accept_terminal(run.parent.parent,run,a.plan_sha256,a.stop,a.max_retries,publish=False),indent=2))


if __name__=='__main__': main()

"""Check full APHCP chemistry ledgers against reconciliation, including decision isolation."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import time

from holo_m2_chemtype import REF, jsonl, layer_rows, load, sha256, write

DECISION_FIELDS=('role_after','prior_role','prior_status','intake_after','prior_intake','partner_scope_qc',
                 'coordinate_identity','experimental_fidelity','topology','system_path','system_sha256','mapping_path')


def key(r):return r['batch'], r['pdb_id'], r['candidate_key']


def verify(run,before=None):
    clock=time.monotonic(); ref=run.parent/REF; results=run/'results'; summary=load(run/'summary.json')
    baseline=list(jsonl(ref/'system_attribution.jsonl')); overlay,entries=layer_rows(ref,baseline)
    resulting_entries={}; merge={}; output={}; missing=Counter(); decision_status=Counter()
    for file, expected in (('baseline_system_labels.jsonl',baseline),('supplemental_system_labels.jsonl',overlay)):
        index={key(x):x for x in expected}; assert len(index)==len(expected)
        result_ids=set(); status=Counter()
        for row in jsonl(results/file):
            k=key(row); old=index[k]; assert k not in result_ids; result_ids.add(k)
            assert row['classification_changes_existing_decision'] is False
            for fld in DECISION_FIELDS+('revision','source_sha256','source_path','discovery_path','candidate_is_increment','evidence_layer'):
                assert old.get(fld)==row.get(fld), (k,fld)
            c=row['classification']; s=c['annotation_status']; status[s]+=1
            assert c['chemical_subclass']==sorted(set(c['chemical_subclass']))
            if s in ('PARTIAL','UNKNOWN'):assert c['status_reasons'],k
            if c['primary_class']=='UNKNOWN':assert not c['chemical_subclass'],k
            old_status=row.get('prior_status');
            if c['primary_class']=='UNKNOWN':decision_status[(file, str(old_status))]+=1
            if file.startswith('supplemental'):
                assert bool(row['candidate_is_increment'])==(k not in merge)
            merge[k]={'status':s,'primary':c['primary_class']}
        assert result_ids==set(index)
        output[file]={'rows':len(result_ids),'status':dict(status),'matches_reconciliation_keys_and_decisions':True}
    for row in jsonl(results/'entry_dispositions.jsonl'):
        p=row['pdb_id']; assert p in entries and p not in resulting_entries; old=entries[p];resulting_entries[p]=row
        assert all(old.get(f)==row.get(f) for f in ('batch','revision','source_path','source_sha256'))
        if row['status']=='NO_CANDIDATE_IN_BASELINE_OR_RECOVERY':
            assert row['no_candidate_evidence'] and row['establishes_no_target'] is False;missing[row['batch']]+=1
    assert len(merge)==46864 and set(resulting_entries)==set(entries)
    assert sum(not (r['baseline_annotations'] or r['unique_candidate_increment']) for r in resulting_entries.values())==sum(missing.values())
    codes={x['ccd']:x['baseline_observed'] for x in jsonl(results/'ccd_labels.jsonl')}
    assert {p for p,v in codes.items() if v}=={c for r in baseline for c in r['ccds']}
    assert sum(codes.values())==1028
    artifacts={x.name:{'bytes':x.stat().st_size,'sha256':sha256(x)} for x in sorted(results.glob('*.jsonl'))}
    assert artifacts==summary['artifacts']
    prior=load(before) if before else None
    if prior:
        assert prior['artifacts']==artifacts
        assert summary['input']['source_shards_built']==0 and summary['input']['source_shards_reused']==2007
        assert summary['validation']==prior['validation']
    receipt={'verified_entries':len(entries),'baseline_ccds':1028,'supplement_only_codes':sorted(k for k,v in codes.items() if not v),
             'files':output,'effective_candidate_rows':len(merge),'effective_status':dict(Counter(v['status'] for v in merge.values())),
             'no_candidate_by_batch':dict(missing),'unresolved_primary_and_original_status':{'/'.join(k):v for k,v in sorted(decision_status.items())},
             'artifacts':artifacts,'science_or_decision_fields_changed':0,'verified_after_rerun':bool(prior),
             'seconds':round(time.monotonic()-clock,2)}
    write(run/'verification.json',receipt);print(json.dumps(receipt,ensure_ascii=False,sort_keys=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('--compare-with',type=Path)
    a=p.parse_args();verify(a.run,a.compare_with)

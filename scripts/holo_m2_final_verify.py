"""Check the frozen final M2 roster, evidence ledgers and first-pass/continuation counts."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path


def read(p):
    return json.loads(p.read_text(encoding='utf8'))


def sha(p):
    with p.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()


def table(p):
    with p.open(newline='',encoding='utf8') as f: return list(csv.DictReader(f))


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('run',type=Path); a=ap.parse_args(); r=a.run.resolve()
    f=read(r/'frozen.json'); plan={x['pdb_id']:x for x in f['files']}; out=table(r/'entry_results.csv'); sys=table(r/'system_results.csv'); sel=read(r/'selection_evidence.json')
    sm={x['pdb_id']:x for x in sel}; excluded=read(r/'inputs/history_excluded.json'); ex=set(excluded if isinstance(excluded,list) else excluded['ids'])
    assert len(plan)==len(out)==2000==len(sm)==len(f['files']) and set(plan)=={x['pdb_id'] for x in out}==set(sm)
    assert not set(plan)&ex; assert sha(r/'inputs/blind.json')==f['input_sha256'] and sha(r/'inputs/archive_head.json')==f['head_sha256']
    assert len(sys)==len({(x['pdb_id'],x['candidate_key']) for x in sys})
    assert all(x['revision']==plan[x['pdb_id']]['revision'] and x['pdb_id'] in plan for x in out+sys)
    assert Counter(x['group'] for x in sel)==Counter(x['stratum'] for x in out)==Counter({'INDEX_MATCH':1200,'NO_LIGAND_ANNOTATION':600,'INDEX_VERSION_MISMATCH':200})
    by=defaultdict(list)
    for x in sys:
        by[x['pdb_id']].append(x)
        if x['status']=='EXPORTED_READBACK':
            assert x['coordinate_identity']==x['experimental_fidelity']==x['topology']=='PASS' and x['system_path'] and (r/x['system_path']).is_file()
        elif x['candidate_intake']=='PASS': raise AssertionError('Only exported, QC-passing candidates can pass intake')
        if x['candidate_intake']=='PASS': assert x['role_qc']==x['partner_scope_qc']=='PASS'
    proof=Counter(); qc=Counter(); techniques=Counter(); groups={}; holds=Counter(); gaps=Counter()
    for x in out:
        p=x['pdb_id']; s=read(r/'source_proofs'/f'{p}.json'); proof[s['status']]+=1
        if s['status']=='PASS':
            assert s['identity']==plan[p] and s['gzip_crc']=='PASS' and s['observed_revision']==plan[p]['revision']
            assert s['bytes']==plan[p]['expected_bytes']==(r/s['relative_path']).stat().st_size and x['source_sha256']==s['sha256']
        group=groups.setdefault(x['stratum'],Counter()); group['entries']+=1; group[x['status']]+=1
        candidates=by[p]; ec=Counter(s['status'] for s in candidates); group.update(ec); group['candidate_systems']+=len(candidates)
        group['no_structure_outputs']+=not ec['EXPORTED_READBACK']; group['without_candidate_system']+=not candidates
        if x['status']!='TECHNICAL_FAILURE':
            assert len(candidates)==int(x['chosen_systems']) and ec['EXPORTED_READBACK']==int(x['structure_readback_pass'])
        else: techniques[f"{x.get('failure_stage')}: {x.get('error')}"]+=1
        if not candidates: gaps['no_structural_candidates']+=1
        elif not ec['EXPORTED_READBACK']: gaps['candidates_but_no_readback']+=1
    for s in sys:
        for key in ('coordinate_identity','experimental_fidelity','topology','role_qc','candidate_intake'):
            qc[f'{key}:{s[key]}']+=1
        if s.get('reason'): holds[s['reason']]+=1
    fist=read(r/'blind_first_pass.json'); assert all(sha(r/'first_pass'/f)==h for f,h in fist['inputs_and_tables_sha256'].items())
    assert len(fist['entry_discovery_sha256'])==sum((r/'discovery'/f'{p}.json').is_file() for p in plan)
    changed={d:sorted(p for p,h in fist[q].items() if not (r/d/f'{p}.json').is_file() or sha(r/d/f'{p}.json')!=h)
             for d,q in [('discovery','entry_discovery_sha256'),('system_results','entry_systems_sha256')]}
    review=read(r/'independent_audit.json'); assert review['first_pass_sha256']==sha(r/'blind_first_pass.json') and review['n']==20
    handoff=read(r/'m3_first_review.json'); assert handoff['input_sha256']==sha(r/'frozen.json') and handoff['review_count']==len(handoff['review_candidates'])
    for c in handoff['review_candidates']:
        s=next(s for s in by[c['pdb_id']] if s['candidate_key']==c['candidate']); assert s['candidate_intake']=='PASS' and c['system_sha256']==sha(r/c['system_cif'])
        assert c['status']=='REVIEW_BEFORE_PRODUCTION'
    output={'run':str(r),'frozen_ids':len(plan),'history_ids':len(ex),'overlap':0,'final_source_proofs':dict(proof),
            'unique_system_candidates':len(sys),'qualified_by_stratum':groups,'qc_axis_status':dict(qc),'hold_reasons':dict(holds),
            'technical_failures':dict(techniques),'no_output_breakdown':dict(gaps),'first_pass_changed_evidence':changed,
            'audit_status':review['status_counts'],'m3_manual_review_candidates':handoff['review_count']}
    (r/'delivery_validation.json').write_text(json.dumps(output,ensure_ascii=False,indent=2),encoding='utf8')
    print(json.dumps(output,ensure_ascii=False))


if __name__=='__main__':main()

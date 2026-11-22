"""Replay only new structural candidates; independently read back changed QC files."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import time

import gemmi

from holo_candidate_batch import atomic_json, now, revision_from_cif, sha256
from holo_m2_blind_structure import one, explicit_role_names
from holo_role_evidence import assess, role_source
from holo_m2_rootcause import RUNS, EXPECTED, SCAN_VERSION, load, dump, source_metadata
from holo_v1_export import rows, text


def jsonl(path):
    with path.open(encoding='utf8') as f:
        for line in f:
            yield json.loads(line)


def verify_file(row):
    file=Path(row['system_path']); mapping=load(Path(row['mapping_path']))
    assert sha256(file)==row['system_sha256'] and row['source_sha256']==mapping['source_sha256']
    block=gemmi.cif.read_file(str(file)).sole_block(); st=gemmi.make_structure_from_block(block)
    ids={text(a['label_entity_id']) for a in rows(block,'_atom_site.') if float(a['occupancy'])>0}
    assert set(row['entity_ids'])<=ids and len(st)==1 and st[0].num==row['model']
    if row['ccds']:
        assert set(row['ccds'])<={text(a['label_comp_id']) for a in rows(block,'_atom_site.')}
    return file,block,st,ids


def run_structure(out, selections):
    started=time.monotonic(); partition={}; runs={}; result=[]
    for row in selections:
        p=row['pdb_id']; partition.setdefault((row['batch'],p),[]).append(row)
    for (batch,p), plan in partition.items():
        assert all(x['revision']==plan[0]['revision'] for x in plan)
        root=out.parent/RUNS[batch]; old=load(out/f'scan_{SCAN_VERSION}'/f'{p}.json')['entry']; sub=out/'new_structure'/p
        for dirname in ('raw','source_proofs','discovery','entry_results','system_results'):
            (sub/dirname).mkdir(parents=True,exist_ok=True)
        raw=Path(old['source_path']); proof=root/'source_proofs'/f'{p}.json'
        assert sha256(raw)==old['source_sha256']==load(proof)['sha256']
        for ref in (proof,raw):
            dest=sub/ref.parent.name/ref.name
            if dest.is_symlink():
                assert dest.resolve()==ref.resolve()
            else:
                dest.symlink_to(ref)
        ix=load(root/'index_snapshot.json')['pdb'][p]
        observed_role_qc={s['candidate_key']:s['prior_role'] for s in load(out/f'scan_{SCAN_VERSION}'/f'{p}.json')['systems']}
        if observed_role_qc:
            raise ValueError('Selected Entry has old structures: only a candidate-restricted rebuild should be used')
        signature=hashlib.sha256(dump({'inputs':old['source_sha256'],'code':{f.name:sha256(f) for f in Path(__file__).parent.glob('holo*.py') if f.name!=Path(__file__).name},
                                       'plan':plan,'index':ix}).encode()).hexdigest()
        runs[p]=one(sub,{'pdb_id':p,'revision':old['revision']},ix,signature)[1]
        rr=load(sub/'entry_results'/f'{p}.json'); led=load(sub/'discovery'/f'{p}.json'); sats=load(sub/'system_results'/f'{p}.json')
        choices={r['candidate_key']:r for r in sats}
        assert set(choices)=={r['candidate_key'] for r in plan}, (p,set(choices))
        for r in plan:
            k=r['candidate_key']; s=choices[k]; obs=next(x for x in led['native_observations'] if x['candidate_key']==k)
            assert s['status']=='EXPORTED_READBACK' and s['role_qc']=='REVIEW' and s['candidate_intake']=='REVIEW'
            file=sub/s['artifact']['path']; block=gemmi.cif.read_file(str(file)).sole_block(); st=gemmi.make_structure_from_block(block)
            assert sha256(file)==s['artifact']['sha256'] and len(st)==1 and st[0].num==r['model']
            members={m['entity'] for m in obs['members']}
            native={rr['label_seq'] for m in obs['members'] for rr in m['residues'] if rr['positive_heavy_atoms']}
            written={int(text(x['label_seq_id'])) for x in rows(block,'_atom_site.')
                     if text(x['label_entity_id']) in members and float(x['occupancy'])>0 and text(x['type_symbol']) not in ('H','D','T')}
            assert native==written==set(r['evidence']['positive_observed_label_seq'])
            result.append(r | {'system_path':str(file),'system_sha256':sha256(file),'mapping_path':str(sub/s['mapping']),
                                'readback':'PASS','role':'REVIEW','intake':'REVIEW','source_models':rr['models'],
                                'declared_length':r['evidence']['declared_length'],'observed_fragment_residues':len(written),
                                'system_bytes':file.stat().st_size,'source_sha256':old['source_sha256']})
    return {'invocations':runs,'systems':result,'elapsed_s':round(time.monotonic()-started,3)}


def control_changes(out):
    # One expectation per independent Entry/entity and role from direct reading of the source and CIF.
    expected={('7JG4','literal_ccd_with_explicit_title_role','PASS'):('2','GAR',1),
              ('6XS5','entity_named_cyclic_peptide','PASS'):('2',None,1),
              ('6XS9','entity_named_cyclic_peptide','PASS'):('3',None,4),
              ('1GAG','entity_explicit_peptide_inhibitor','PASS'):('2',None,3),
              ('4CSM','entity_chemical_bound_site','PASS'):('2','TYR',6),
              ('5N4E','explicit_crystal_cryoprotectant','NON_TARGET_EVIDENCED'):('3','GOL',1),
              ('8RNP','explicit_crystal_cryoprotectant','NON_TARGET_EVIDENCED'):('4','GOL',5)}
    changed=Counter(); receipt=[]; blocks={}
    for r in jsonl(out/'role_changes.jsonl'):
        p=r['pdb_id']; b=r['evidence_after']; key=p,b['rule'],r['role_after']
        ent,ccd,_=expected[key]; assert r['entity_ids']==[ent] and r['ccds']==([ccd] if ccd else [])
        if p not in blocks:
            blocks[p]=revision_from_cif(Path(r['source_path']),p)[0]
        raw=blocks[p]; file,block,st,ids=verify_file(r); conns=rows(block,'_struct_conn.')
        assert 'PASS'==r['coordinate_identity']==r['experimental_fidelity']==r['topology']==r['partner_scope_qc']
        match=next(x for x in rows(raw,'_entity.') if text(x['id'])==ent)
        if key[1]=='entity_chemical_bound_site':
            assert 'TYROSINE BOUND TO REGULATORY SITE'==text(match['details']) and 'TYROSINE' not in text(raw.find_value('_struct.title'))
        elif key[1]=='entity_named_cyclic_peptide':
            assert 'Cyclic peptide RT-' in text(match['details'])
            assert any(text(c['conn_type_id'])=='covale' for c in conns)
        elif key[1]=='entity_explicit_peptide_inhibitor':
            assert 'BISUBSTRATE PEPTIDE INHIBITOR'==text(match['pdbx_description'])
            link=[c for c in conns if {text(c['ptnr1_label_comp_id']),text(c['ptnr2_label_comp_id'])}=={'112','PHE'}]
            assert any(text(c['conn_type_id'])=='covale' for c in link) and '4' in ids, file
        elif key[1]=='literal_ccd_with_explicit_title_role':
            assert 'GAR substrate' in text(raw.find_value('_struct.title'))
        else:
            assert 'glycerol' in b['text'].lower() and 'cryoprotected' in b['text'].lower()
            assert b['text'] in ' '.join(text(x) for x in raw.find_values(b['field']))
        changed[key]+=1; receipt.append({'pdb_id':p,'key':r['candidate_key'],'role':r['role_after'],
            'rule':b['rule'],'disposition':'SOURCE_EVIDENCE_AND_EXISTING_CIF_REREAD_VERIFIED',
            'system_path':str(file),'system_sha256':r['system_sha256'],'source_sha256':r['source_sha256']})
    assert changed==Counter({k:v[2] for k,v in expected.items()})
    return receipt


def individual_controls(out):
    checks=[]
    for p,code,why in (('6TDU','CDL','Title absent; primary citation specifically names native bound cardiolipin, but does not classify the role of every instance for this resource.'),
                       ('7JG4','V9V','Depositor title names AGF131, but the fixed CIF has no explicit identity bridge from that experimental name to this deposited CCD.'),
                       ('6XS5','GOL','Geometry and chemical identity without explicit target or preparation-role statement are insufficient.'),
                       ('8RNP','9IR','The original R-limonene title/instance pass must remain unchanged beside the glycerol excluded from intake.')):
        frame=load(out/f'scan_{SCAN_VERSION}'/f'{p}.json'); path=Path(frame['entry']['source_path']); raw=revision_from_cif(path,p)[0]
        meta=source_metadata(raw); title=meta['extra']['_struct.title']; expected='PASS' if p=='8RNP' else 'REVIEW'
        selected=[s for s in frame['systems'] if s['ccds']==[code] and s['prior_status']=='EXPORTED_READBACK']
        assert selected
        for x in selected:
            assert x['prior_role']==expected and assess(x['role_input'],x['prior_role'],title,explicit_role_names(raw),meta['polylen'],role_source(raw)) is None
        verify_file(selected[0])
        checks.append({'pdb_id':p,'ccd':code,'matching_readbacks':len(selected),'status':expected,'reason':why,
                       'source_path':str(path),'system_path':selected[0]['system_path']})
    return checks


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('run',type=Path); out=ap.parse_args().run.resolve()
    start=time.monotonic(); selection=load(out/'structure_retest_selection.json')
    new=run_structure(out,[x for x in selection if not x['historically_verified']]); changed=control_changes(out); controls=individual_controls(out)
    repro=load(out/'attribution_summary.json'); hist=load(out/'historical_deltas.json')['deltas']
    old={b:{m:sum(x[k] for x in hist if x['pdb_id'] in ('6TO3','6X8H','7NDX','7C5Y')) for m,k in [('candidates','delta_candidates'),('readbacks','delta_readbacks')]} if b=='200' else {} for b in EXPECTED}
    old['200']['intake']=sum(x.get('delta_intake',0) for x in hist)
    by={b:{'candidates':a['candidates']+old[b].get('candidates',0)+sum(x['batch']==b for x in new['systems']),
           'readbacks':a['readbacks_before']+old[b].get('readbacks',0)+sum(x['batch']==b for x in new['systems']),
           'program_intake':a['intake_after']+old[b].get('intake',0),
           'zero_candidates':a['zero_candidates']-(2 if b=='200' else len({x['pdb_id'] for x in new['systems']}))} for b,a in repro['by_batch'].items()}
    atomic_json(out/'validation.json',{'finished_utc':now(),'changed_structure_roles_reopened':changed,'positive_negative_and_unresolved_controls':controls,
                 'structure_retest':new,'historical_deltas_outside_sealed_baseline':old,'distinct_combined_results_with_deltas':by,
                 'new_coordinate_download_bytes':0,'elapsed_s':round(time.monotonic()-start,3)})
    print(dump({'new':new,'changed_role_assertions_checked':len(changed),'controls':len(controls),'combined':by}),flush=True)


if __name__=='__main__':
    main()

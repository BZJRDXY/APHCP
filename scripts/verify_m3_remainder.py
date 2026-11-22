"""Validate an M3 combined staging registry, real asset lookup and a repeat control."""
import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import subprocess
import sys

import holo_m3_register as m


def key(x):
    return str(x['batch']), x['pdb_id'], x['candidate_key']


def outcome(rows):
    return {'candidates':len(rows),'entries':len({x['pdb_id'] for x in rows}),
            'dispositions':dict(Counter(x['disposition'] for x in rows)),
            'reasons':dict(Counter(z for x in rows for z in x['reasons']))}


def inspect(row):
    files=row['files']; model=m.gemmi.cif.read_file(files['system_cif']['path']).sole_block()
    raw,revision=m.revision_from_cif(Path(files['source_mmcif']['path']),row['pdb_id'])
    mapping=m.read(files['source_mapping']['path']); atoms=len(model.find_values('_atom_site.id'))
    assert revision==row['revision'] and atoms==row['qc']['atoms']>0
    assert mapping['artifact']['sha256']==files['system_cif']['sha256']
    assert sorted(mapping['partner_member_identities'])==row['identity']['target_residue_instances']
    output={'system_id':row['system_id'],'pdb_id':row['pdb_id'],'revision':revision,'atoms':atoms,'disposition':row['disposition'],
            'reasons':row['reasons'],'target_residue_instances':row['identity']['target_residue_instances'],
            'files':{k:files[k] for k in ('source_mmcif','system_cif','source_mapping')}}
    return output,model,mapping


def verify(first, batch, control):
    one=m.read(first/'review_results.json'); two=m.read(batch/'review_results.json'); diag=m.read(control/'review_results.json')
    frozen=m.read(batch/'frozen_manifest.json'); keys_first={key(x) for x in one}; keys_second={key(x) for x in two}
    assert len(keys_first)==len(one) and len(keys_second)==len(two)==frozen['review_count']
    assert not (keys_first&keys_second) and keys_second=={key(x) for x in frozen['review_candidates']}
    assert len(diag)==1 and key(diag[0]) in keys_first and diag[0]['disposition']=='REGISTERED'
    reg=m.read(batch/'registration_receipt.json'); rerun=m.read(batch/'rerun_receipt.json'); ctrl=m.read(control/'registration_receipt.json')
    assert reg['registry_counts']==rerun['registry_counts'] and reg['results_sha256']==rerun['results_sha256']==m.sha(batch/'review_results.json')
    assert reg['readback_files_reused_on_disk']==rerun['prior_review_checks_reused'] and rerun['database_registration_reused']
    for line in (batch/'protected_inputs_before.sha256').read_text().splitlines()+(batch/'rerun_guard.sha256').read_text().splitlines():
        hash_,src=line.split('  ',1); assert m.sha(src)==hash_,src
    totals={}
    with m.connect(reg['registry_path'],readonly=True) as db:
        assert db.execute('PRAGMA integrity_check').fetchone()[0]=='ok' and not db.execute('PRAGMA foreign_key_check').fetchall()
        assert db.execute('SELECT COUNT(*) FROM release').fetchone()[0]==0
        for source, rows in ((first,one),(batch,two)):
            actual=Counter(r[0] for r in db.execute("SELECT status FROM run_item WHERE stage='M3_REVIEW' AND run_id=?",(source.name+'-run',)))
            assert actual==Counter(x['disposition'] for x in rows)
        totals['work_item_kinds']={k:n for k,n in db.execute('SELECT kind,COUNT(*) FROM work_item GROUP BY kind')}
        totals['unique_fixed_revisions']=db.execute('SELECT COUNT(*) FROM entry_version').fetchone()[0]
        totals['registered_systems'],totals['registered_pdb_entries'],totals['registered_ligand_instances']=tuple(db.execute("""
          SELECT count(DISTINCT ri.object_id),count(DISTINCT ri.pdb_id),count(DISTINCT m.member_id)
          FROM run_item ri JOIN system_member m ON ri.object_id=m.system_id WHERE ri.stage='M3_REVIEW' AND ri.status='REGISTERED' AND m.role='ligand'
          """).fetchone())
    def lookup(run, selector, value):
        cmd=[sys.executable,'-B',str(batch/'code/holo_m3_register.py'),'lookup','--run',str(run),selector,value,'--verify-files']
        return json.loads(subprocess.check_output(cmd,text=True))['systems']
    repeat=diag[0]['system_id']; first_record=next(x for x in one if x.get('system_id')==repeat)
    assert first_record['identity']==diag[0]['identity'] and first_record['files']['system_cif']==diag[0]['files']['system_cif']
    with m.connect(ctrl['registry_path'],readonly=True) as db:
        assert db.execute('PRAGMA integrity_check').fetchone()[0]=='ok' and not db.execute('PRAGMA foreign_key_check').fetchall()
        assert db.execute('SELECT count(*) FROM work_item WHERE id=?',(repeat,)).fetchone()[0]==1
        lineage=[r[0] for r in db.execute("SELECT run_id FROM run_item WHERE stage='M3_REVIEW' AND object_id=? ORDER BY run_id",(repeat,))]
        assert set(lineage)=={first.name+'-run',control.name+'-run'}
        assert all(reg['registry_counts'][k]==ctrl['registry_counts'][k] for k in ('pdb_entry','entry_version','asset_file','work_item','system_member','release'))
    cross=lookup(control,'--system',repeat)
    assert len(cross)==1 and {x['run_id'] for x in cross[0]['provenance']}==set(lineage)
    repeat_proof=inspect(cross[0])[0]|{'provenance_runs':lineage,'canonical_rows':1,'new_canonical_or_asset_rows':0}
    same=lookup(batch,'--pdb','7DCW'); same={x['system_id']:x for x in same if x['disposition']=='REGISTERED'}
    assert len(same)==2
    instances=[]
    for sid in sorted(same):
        found=lookup(batch,'--system',sid); assert len(found)==1 and same[sid]==found[0]
        instances.append(inspect(found[0])[0])
    a,b=[v['target_residue_instances'][0] for v in instances]; assert a[-1]==b[-1] and a!=b
    prior_review=next(x for x in one if x['pdb_id']=='1JPL' and x['disposition']=='REVIEW')
    old=lookup(batch,'--system',prior_review['system_id']); assert len(old)==1 and old[0]['disposition']=='REVIEW'
    retained=inspect(old[0])[0]
    # A missing chemical label must not block an otherwise eligible observation (D-0018).
    prior_seal=m.read(first/'frozen_staging.json'); n, = [x for x in m.read_manifest(Path(prior_seal['manifest'])) if x[1]=='1AVD']
    k=n[:3]; runs=batch.parent; m2=runs/'holo-m2-final2000-20261122a'; ref=runs/'holo-m2-full-reconciliation-20261122a'
    base,_,_=m.pool_rows(ref,runs/'holo-m2-chemtype-20261122a',{k})
    unknown=m.assess_candidate(k,n[3],base[k],None,(m2,ref,runs/'holo-m2-blind200-20261122a'),m2,None,{})
    assert unknown['type']['primary_class']=='UNKNOWN' and unknown['disposition']=='REGISTERED' and unknown['system_id']==repeat
    unknown_control={k:unknown[k] for k in ('disposition','system_id')}|{'primary_class':'UNKNOWN'}
    history=[]
    for row in two:
        if row.get('evidence_layer')!='covalent_copy':continue
        found=lookup(batch,'--system',row['system_id']); assert len(found)==1 and found[0]['disposition']=='REGISTERED'
        proof,block,mapping=inspect(found[0]); endpoint_ids=set(block.find_values('_struct_conn.id'))
        edges=[edge for edge in mapping['artifact']['connections'] if edge['output']['kind'] in {'covale','disulf','modres'}]
        partner_labels={row['output_instances'][v['instance']]['output_label_asym'] for v in row['partner_instances']}
        cross_chain=[edge for edge in edges if len({e['label_asym_id'] for e in edge['output']['ends']})>1 and
                     partner_labels & {e['label_asym_id'] for e in edge['output']['ends']}]
        assert cross_chain and all(v['output']['id'] in endpoint_ids for v in cross_chain)
        history.append(proof | {'partner_copies':sorted(partner_labels),'cross_chain_covalent':cross_chain,
                               'historical_verifier_sha256':row['qc']['historical_verifier_sha256']})
    ext={(str(x['batch']),x['pdb_id'],x['candidate_key']):x for x in one+two if 'complete_partner_membership_requires_confirmation' in x['reasons']}
    type_root=batch.parent/'holo-m2-chemtype-20261122a/results'
    chem={}
    for f in (type_root/'baseline_system_labels.jsonl',type_root/'supplemental_system_labels.jsonl'):
        for x in m.jsonl(f):
            if key(x) in ext:chem[key(x)]=x
    covariance=[]
    for k,row in sorted(ext.items()):
        _,block,mapping=inspect(row); disc=chem[k]['partner_identity_evidence']['external_covalent']; assert disc
        present=set(block.find_values('_struct_asym.id')); ids=set(block.find_values('_struct_conn.id'))
        coverage=[]
        for link in disc:
            ins=mapping['instances'].get(link['foreign_instance']); label=m.output_label(ins) if ins else None
            actual=[e for e in mapping['artifact']['connections'] if e['source']['id']==link['id'] or e['source'].get('source_id')==link['id']]
            edges=[e['output']['id'] for e in actual if e['output']['id'] in ids and label in {x['label_asym_id'] for x in e['output']['ends']}]
            coverage.append({'original_connection_id':link['id'],'foreign_instance':link['foreign_instance'],
                             'mapped_output_chain':label,'chain_observed_in_cif':label in present,
                             'mapped_connections_observed_in_cif':edges})
        covariance.append({'batch':row['batch'],'pdb_id':row['pdb_id'],'candidate_key':row['candidate_key'],'system_id':row['system_id'],
                           'coverage':coverage,'current_disposition':row['disposition']})
    splits={}
    for name, rows in (('first',one),('remaining',two),('combined',one+two)):
        yes=[x for x in rows if x['disposition']=='REGISTERED']; has=[x for x in rows if x.get('system_id')]
        splits[name]=outcome(rows)|{'entry_counts_by_disposition':{label:len({x['pdb_id'] for x in rows if x['disposition']==label})
                    for label in outcome(rows)['dispositions']},'fixed_revisions':len({(x['pdb_id'],x['revision']) for x in has}),
                    'retrievable_systems':len({x['system_id'] for x in has}),'registered_chemical_classes':dict(Counter(x['type']['primary_class'] for x in yes)),
                    'unknown_type_by_disposition':dict(Counter(x['disposition'] for x in has if x['type']['primary_class']=='UNKNOWN'))}
    return {'scope':'M3_STAGING_UNPUBLISHED','splits':splits,'combined_sqlite':totals,'cross_batch_repeat':repeat_proof,
            'distinct_same_ccd_instances':instances,'first_review_retrieved_in_combined':retained,'d0018_unknown_type_control':unknown_control,'reused_historical_observations':history,
            'external_covalent_pending_evidence':covariance,'registration_receipt':reg,'rerun_receipt':rerun,'control_receipt':ctrl,
            'protected_source_hashes_unchanged':True,'rerun_hashes_unchanged':True}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('first','batch','control'):p.add_argument('--'+k,type=Path,required=True)
    args=p.parse_args(); report=verify(args.first,args.batch,args.control)
    m.write(args.batch/'verification.json',report)
    print(m.raw_json(report['splits']|{'combined_sqlite':report['combined_sqlite'],'cross_batch_provenance':report['cross_batch_repeat']['provenance_runs']}))

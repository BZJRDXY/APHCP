"""Post-freeze, independent checks of source observations and written V1 artifacts."""
from collections import Counter
import csv
from hashlib import sha256 as hashfunc
import json
from pathlib import Path

import gemmi
import numpy as np
from scipy.spatial.distance import cdist

from holo_candidate_batch import atomic_json, now, revision_from_cif, sha256


def rd(p):
    return json.loads(p.read_text(encoding='utf8'))


def body(path, kind):
    return gemmi.cif.read(str(path)).sole_block().get_mmcif_category(kind)


def inspect_output(run, system, raw):
    art=system['artifact']; path=run/art['path']; assert sha256(path)==art['sha256']
    source=body(raw,'_atom_site.'); output=body(path,'_atom_site.'); mapping=rd(run/system['mapping'])
    source_byid={(str(m),str(i)):j for j,(m,i) in enumerate(zip(source['pdbx_PDB_model_num'],source['id']))}
    src_labels=set(); errors=[]; alternates=Counter(); pose=Counter()
    child=body(path,'_aphcp_atom_source.'); expanded=bool(child)
    if expanded:
        line_byid={str(i):j for j,i in enumerate(child['atom_site_id'])}
        assert len(line_byid)==len(output['id'])
    want_models={str(system['representative_model'])}; assert set(output['pdbx_PDB_model_num'])==want_models
    metrics=('occupancy','label_alt_id','B_iso_or_equiv','label_atom_id','label_comp_id','label_seq_id','auth_seq_id','pdbx_PDB_ins_code')
    for i,oid in enumerate(output['id']):
        model=output['pdbx_PDB_model_num'][i]
        if expanded:
            l=line_byid[str(oid)]; member=mapping['instances'][child['assembly_instance_id'][l]]
            assert child['source_model_id'][l]==model
            source_id=child['source_atom_site_id'][l]
        else:
            source_id=oid; member=None
        original=source_byid[str(model),str(source_id)]
        fields=[name for name in metrics if name in source]
        assert all(output[name][i]==source[name][original] for name in fields)
        coord=('Cartn_x','Cartn_y','Cartn_z')
        expected=np.array([float(source[name][original]) for name in coord])
        if member:
            assert source['label_asym_id'][original]==member['label_asym']
            expected=np.array(member['matrix'])@expected+np.array(member['vector'])
        else:
            assert output['label_asym_id'][i]==source['label_asym_id'][original]
        errors.append(float(np.max(np.abs(np.asarray([float(output[name][i]) for name in coord])-expected))))
        src_labels.add(source['label_asym_id'][original])
        if output['label_alt_id'][i]:
            alternates[str(output['label_alt_id'][i])]+=1
            pose[(str(output['label_alt_id'][i]),str(output['occupancy'][i]))]+=1
    worst=max(errors); assert worst<=(5.01e-7 if expanded else 1e-12)
    conn=body(path,'_struct_conn.'); pairs=body(path,'_ndb_struct_na_base_pair.')
    return {'file':art['path'],'candidate_key':system['candidate_key'],'role_qc':system['role_qc'],'atom_rows_independently_compared':len(errors),
            'max_axis_residual_A':worst,'altloc_atoms':dict(alternates),'source_chain_labels':sorted(src_labels),
            'altloc_and_occupancy':dict((f'{a}:{b}',n) for (a,b),n in pose.items()),'links_written':len(conn.get('id',[])),
            'base_pair_rows_written':len(pairs.get('model_number',[])), 'assembly_copy':expanded}


def source_observations(st, block, models):
    # Gemmi's direct ASU traversal and bounding-box/Euclidean contact check, independent
    # of the run's assembly expansion, candidate grouping and tree-based selection.
    entity=block.get_mmcif_category('_entity_poly.'); ty=dict(zip(entity.get('entity_id',()),entity.get('type',())))
    types=Counter(); pos=Counter(); components={}; links=block.get_mmcif_category('_struct_conn.')
    for model in st:
        if model.num not in models:
            continue
        receptor=[]; obs=[]
        for chain in model:
            for res in chain:
                ats=[(a.pos.x,a.pos.y,a.pos.z, a.altloc if a.altloc not in ('\0',' ','.','?') else '') for a in res if a.occ>0 and not a.element.is_hydrogen]
                if not ats or res.is_water():
                    continue
                t=str(ty.get(res.entity_id,''))
                if t.startswith('polypeptide'):
                    receptor.extend(ats)
                elif not t:
                    obs.append((res.subchain,res.name,str(res.seqid),ats))
                else:
                    types[t]+=1
        receptor_coords=np.asarray([it[:3] for it in receptor]); receptor_alt=np.asarray([it[3] for it in receptor])
        for chain,ccd,seq,ats in obs:
            key=(model.num,chain,ccd,seq)
            pos[ccd]+=1
            x=np.array([it[:3] for it in ats]); lo=x.min(axis=0)-4; hi=x.max(axis=0)+4
            if not receptor:
                continue
            r=receptor_coords; selection=np.flatnonzero(np.all((lo<=r)&(r<=hi),axis=1))
            if len(selection):
                d=cdist(x,r[selection]); al=np.asarray([it[3] for it in ats]); bl=receptor_alt[selection];
                if np.any((d<=4.) & ((al[:,None]==bl[None,:]) | (al[:,None]=='') | (bl[None,:]==''))):
                    components[key]=len(ats)
    return {'asu_nonpolymer_components':dict(pos),'asu_polymer_partner_residue_counts':dict(types),
            'asu_direct_contact_components':[{'model':k[0],'source_label':k[1],'ccd':k[2],'source_auth_seq':k[3],'heavy_atoms':n} for k,n in components.items()],
            'deposited_connection_types':dict(Counter(links.get('conn_type_id',[]))),'source_base_pair_rows':len(block.get_mmcif_category('_ndb_struct_na_base_pair.').get('model_number',[]))}


def audit(run, plan):
    roster=rd(run/'audit_sample.json'); sample=roster['entries']; first=rd(run/'blind_first_pass.json')
    assert 15<=len(sample)<=20 and len({x['pdb_id'] for x in sample})==len(sample)
    current={r['pdb_id']:r for r in plan}; results=[]
    for spec in sample:
        p=spec['pdb_id']; pin=current[p]; src=rd(run/'source_proofs'/f'{p}.json'); rec=rd(run/'entry_results'/f'{p}.json')
        item={'pdb_id':p,'revision':pin['revision'],'sample_reason':spec['reason'],'entry_status':rec['status'],'stratum':rec['stratum']}
        try:
            if src['status']!='PASS':
                item.update(audit_status='SOURCE_UNAVAILABLE',source_fetch_error=src['error']);results.append(item);continue
            raw=run/src['relative_path']; block,rev=revision_from_cif(raw,p); assert rev==pin['revision'] and sha256(raw)==src['sha256']
            structure=gemmi.make_structure_from_block(block); models=[m.num for m in structure]
            item.update(source_models=models,source_assembly_ids=[a.name for a in structure.assemblies])
            if rec['status']=='TECHNICAL_FAILURE':
                item.update(audit_status='PIPELINE_TECHNICAL_FAILURE',error=rec['error']);results.append(item);continue
            ledger=run/'discovery'/f'{p}.json'; results_path=run/'system_results'/f'{p}.json'
            assert sha256(ledger)==first['entry_discovery_sha256'][p] and sha256(results_path)==first['entry_systems_sha256'][p]
            d=rd(ledger); sys=rd(results_path); eligible=sum(x['eligible_for_structure_review'] for x in d['candidate_selection'])
            native_keys={(r['model'],x['label_asym'],u['ccd'],u['auth_seq']) for r in d['native_observations'] for x in r['members'] for u in x['residues']}
            # ASU-only geometry is diagnostic: contacts may be absent from the deposited
            # assembly or not belong to a target. Record the disagreement rather than force PASS.
            independent=source_observations(structure,block,{models[0],d['representative_model']}); checks=independent['asu_direct_contact_components']
            missing=[x for x in checks if (x['model'],x['source_label'],x['ccd'],x['source_auth_seq']) not in native_keys]
            to_check=[x for x in sys if x.get('artifact')]
            to_check.sort(key=lambda x:(x['candidate_intake']!='PASS',not bool(x['connections_from_source']), x['artifact']['bytes']))
            disk=inspect_output(run,to_check[0],raw) if to_check else None
            item.update(eligible=eligible,system_states=dict(Counter(x['status'] for x in sys)),audit_status='MANUAL_REVIEW' if not disk or missing else 'INDEPENDENT_READBACK_PASS',
                        source_observations=independent,asu_direct_components_missing_in_assembly_observations=missing,
                        disk=disk,role_acceptance='REVIEW_PENDING_INDEPENDENT_CONTEXTUAL_ROLE_REVIEW',
                        interpretation='ASU contacts are an omission diagnostic, not a biological-role decision.')
        except Exception as exc:
            item.update(audit_status='AUDIT_CHECK_FAILED', error=f'{type(exc).__name__}: {exc}')
        results.append(item)
    atomic_json(run/'independent_audit.json',{'checked_utc':now(),'first_pass_sha256':sha256(run/'blind_first_pass.json'),'n':len(results),
        'status_counts':dict(Counter(x['audit_status'] for x in results)),'entries':results})
    print('independent audit',len(results),dict(Counter(x['audit_status'] for x in results)),flush=True)


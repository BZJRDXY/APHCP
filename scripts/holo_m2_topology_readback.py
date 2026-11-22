"""Independently compare deposited topology with first-read V1 chain and atom lineage."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from math import dist
from pathlib import Path

import gemmi
from holo_candidate_batch import atomic_json, now, revision_from_cif, sha256
from holo_m2_blind_audit import inspect_output, rd


def rows(block, prefix):
    data=block.get_mmcif_category(prefix)
    if not data:
        return []
    return [dict(zip(data, values)) for values in zip(*data.values())]


def lookup(rows, key):
    out={row[key]:row for row in rows}
    assert len(out)==len(rows), key
    return out


def check(run, p):
    proof=rd(run/'source_proofs'/f'{p}.json'); raw=run/proof['relative_path']
    original, rev = revision_from_cif(raw,p)
    assert proof['sha256']==sha256(raw) and proof['observed_revision']==rev
    pairs=rows(original,'_ndb_struct_na_base_pair.')
    pair_src={(r['model_number'],r['pair_number']):r for r in pairs}
    assert len(pair_src)==len(pairs)
    conn_src=lookup(rows(original,'_struct_conn.'),'id')
    orig_atoms={ (a['pdbx_PDB_model_num'],a['id']):a for a in rows(original,'_atom_site.') }
    artifacts=[]; states=Counter(); blockers=Counter()
    for rec in rd(run/'system_results'/f'{p}.json'):
        states[rec['status']]+=1
        if not rec.get('artifact'):
            blockers[rec.get('reason')]+=1
            continue
        output=gemmi.cif.read(str(run/rec['artifact']['path'])).sole_block()
        atom=lookup(rows(output,'_atom_site.'),'id'); lineage=lookup(rows(output,'_aphcp_atom_source.'),'atom_site_id')
        copies=lookup(rows(output,'_aphcp_assembly_instance.'),'id'); conn=lookup(rows(output,'_struct_conn.'),'id')
        base_pairs=lookup(rows(output,'_ndb_struct_na_base_pair.'),'pair_number')
        out_links=[]; grouped_pairs=[]
        for t in rows(output,'_aphcp_conn_source.'):
            orig=conn_src[t['source_conn_id']]; written=conn[t['id']]; model=t['source_model_id']; ends=[]; points=[]
            assert orig['conn_type_id']==written['conn_type_id']
            for n in range(1,int(t['endpoint_count'])+1):
                stem=f'ptnr{n}_'; ix=t[stem+'output_atom_site_id']; si=t[stem+'source_atom_site_id']; ci=t[stem+'assembly_instance_id']
                a=atom[ix]; o=orig_atoms[(model,si)]; line=lineage[ix]
                assert (line['source_model_id'],line['source_atom_site_id'],line['assembly_instance_id'])==(model,si,ci)
                assert copies[ci]['source_label_asym_id']==o['label_asym_id'] and copies[ci]['output_label_asym_id']==a['label_asym_id']
                for key in ('label_asym_id','label_comp_id','label_seq_id','label_atom_id','auth_seq_id','label_alt_id'):
                    field=('pdbx_' if key=='label_alt_id' or n==3 else '')+stem+key
                    if orig.get(field) not in (False,None):
                        assert orig[field]==o[key], (p,field,orig[field],o[key])
                    if written.get(field) not in (False,None):
                        assert written[field]==a[key], (p,field,written[field],a[key])
                sym=('pdbx_' if n==3 else '')+stem+'symmetry'; assert (orig.get(sym) or '1_555')==t[stem+'source_symmetry']
                assert written[sym]=='1_555'
                assert all(a[k]==o[k] for k in ('label_alt_id','occupancy','B_iso_or_equiv'))
                points.append([float(a['Cartn_'+k]) for k in 'xyz'])
                ends.append({'instance':ci,'source_atom_site_id':si,'output_atom_site_id':ix,'source_symmetry':t[stem+'source_symmetry'],
                             'altloc':o['label_alt_id'] or None,'occupancy':o['occupancy'],'operator':copies[ci]['source_operator_id']})
            assert len({e['output_atom_site_id'] for e in ends})==len(ends)
            gap=dist(*points[:2]); target=orig.get('pdbx_dist_value')
            out_links.append({'source':orig['id'],'type':orig['conn_type_id'],'endpoints':ends,'distance_A':round(gap,6),
                              'source_distance_A':target or None})
        for t in rows(output,'_aphcp_base_pair_source.'):
            m=t['source_model_id']; orig=pair_src[(m,t['source_pair_number'])]; new=base_pairs[t['output_pair_number']]; ends=[]
            assert new['model_number']==m
            for side in ('i','j'):
                ci=t[side+'_assembly_instance_id']; cp=copies[ci]
                assert new[side+'_label_asym_id']==cp['output_label_asym_id'] and orig[side+'_label_asym_id']==cp['source_label_asym_id']
                assert orig[side+'_symmetry']==t[side+'_source_symmetry'] and new[side+'_symmetry']=='1_555'
                assert all(orig[f'{side}_{k}']==new[f'{side}_{k}'] for k in ('label_comp_id','label_seq_id','auth_seq_id','PDB_ins_code'))
                assert any(a['pdbx_PDB_model_num']==m and all(a[k]==new[f'{side}_{k}'] for k in ('label_asym_id','label_comp_id','label_seq_id')) for a in atom.values())
                ends.append(ci)
            assert all(orig[k]==new[k] for k in ('shear','stretch','stagger','buckle','propeller','opening') if k in orig)
            grouped_pairs.append({'source':orig['pair_number'],'copies':ends,'components':[orig['i_label_comp_id'],orig['j_label_comp_id']]})
        if not (out_links or grouped_pairs):
            continue
        independent=inspect_output(run,rec,raw)
        artifacts.append({'file':rec['artifact']['path'],'hash':rec['artifact']['sha256'],'role':rec['role_qc'],
                          'links':out_links,'pairs':grouped_pairs,'max_axis_residual_A':independent['max_axis_residual_A']})
    return {'pdb_id':p,'revision':rev,'raw_sha256':proof['sha256'],'deposited_connections':len(conn_src),'deposited_pairs':len(pair_src),
            'polymer_types': [x['type'] for x in rows(original,'_entity_poly.')], 'artifact_states':dict(states),'held_reasons':dict(blockers),'artifacts':artifacts}


def scan_raw(run, jobs):
    plan=rd(run/'frozen.json')['files']
    def scan(i):
        p=i['pdb_id']; f=run/'raw'/f'{p.lower()}-v{i["revision"]}.cif.gz'
        assert f.stat().st_size==i['expected_bytes']
        b=gemmi.cif.read(str(f)).sole_block()
        c=rows(b,'_struct_conn.'); a=rows(b,'_pdbx_struct_conn_angle.'); pairs=rows(b,'_ndb_struct_na_base_pair.')
        cov=[t['id'] for t in c if (t['conn_type_id'] in ('disulf','modres') or t['conn_type_id'].startswith('covale'))
              and (x:=t.get('ptnr1_symmetry')) and (y:=t.get('ptnr2_symmetry')) and x!=y]
        third=[t['id'] for t in c if any(t.get(f'pdbx_ptnr3_{k}') for k in ('label_asym_id','label_comp_id','label_atom_id'))]
        differs=sum(bool((x:=t.get('i_symmetry')) and (y:=t.get('j_symmetry')) and x!=y) for t in pairs)
        return {'pdb_id':p,'nonlocal_covalent_rows':cov,'nonlocal_base_pair_rows':differs,'third_endpoint_bond_rows':third,
                'three_atom_angle_rows':len(a),'first_angle':a[0] if a else None}
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        evidence=list(pool.map(scan,plan))
    keys=('nonlocal_covalent_rows','nonlocal_base_pair_rows','third_endpoint_bond_rows','three_atom_angle_rows')
    counts={k:{'entries':sum(bool(x[k]) for x in evidence),'rows':sum(len(x[k]) if isinstance(x[k],list) else x[k] for x in evidence)} for k in keys}
    atomic_json(run/'raw_topology_scan_full.json',{'scanned':len(evidence),'counts':counts,'entries':[x for x in evidence if any(x[k] for k in keys)]})
    print(counts,flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('run',type=Path); p.add_argument('entries',nargs='*'); p.add_argument('--scan-raw',type=int); args=p.parse_args()
    run=args.run.resolve(); checks=[]
    if args.scan_raw:
        scan_raw(run,args.scan_raw); return
    if not args.entries:
        p.error('supply Entry IDs or --scan-raw JOBS')
    for pdb in args.entries:
        one=check(run,pdb); checks.append(one)
        nonlocal_cov=[l for a in one['artifacts'] for l in a['links'] if l['type']=='disulf' or l['type'].startswith('covale')
                      if len({e['source_symmetry'] for e in l['endpoints']})>1]
        print(pdb,one['artifact_states'],'max written pairs',max((len(x['pairs']) for x in one['artifacts']),default=0),
              'nonlocal covalent bonds across outputs',len(nonlocal_cov),flush=True)
    atomic_json(run/'independent_topology_readback.json',{'generated_utc':now(),'entries':checks})


if __name__=='__main__': main()

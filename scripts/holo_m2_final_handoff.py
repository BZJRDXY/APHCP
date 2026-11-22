"""Combine the sealed final M2 results with the narrowly verified topology delta."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import statistics

from holo_m2_final_verify import read, sha, table


def write_csv(path, rows):
    fields=list(dict.fromkeys(k for row in rows for k in row))
    tmp=path.with_suffix('.tmp')
    with tmp.open('w',newline='',encoding='utf8') as f:
        writer=csv.DictWriter(f,fields); writer.writeheader(); writer.writerows(rows)
    tmp.replace(path)


def metrics(entries, systems):
    output=Counter(s['pdb_id'] for s in systems if s['status']=='EXPORTED_READBACK')
    return {'entries':dict(Counter(x['status'] for x in entries)),'system_rows':dict(Counter(x['status'] for x in systems)),
            'candidate_intake':dict(Counter(x['candidate_intake'] for x in systems)),
            'role_axis':dict(Counter(x['role_qc'] for x in systems)),
            'entries_with_readback':len(output),'entries_without_readback':len(entries)-len(output),
            'entries_without_admitted_candidate':len(entries)-len({s['pdb_id'] for s in systems}),
            'readback_reasons':dict(Counter(x['reason'] for x in systems if x['reason'])),
            'assembly_status':dict(Counter(x['complete_assembly_status'] for x in systems if x['complete_assembly_status']))}


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('run',type=Path); a=ap.parse_args(); r=a.run.resolve(); delta=r/'point_symmetry_retest'; result=r/'handoff'; result.mkdir(exist_ok=True)
    first_entries=table(r/'first_pass'/'entry_results.csv'); first_systems=table(r/'first_pass'/'system_results.csv')
    entries=table(r/'entry_results.csv'); systems=table(r/'system_results.csv'); scope={x['pdb_id'] for x in read(delta/'scope.json')['target_entries']}
    es=table(delta/'entry_results.csv'); ys=table(delta/'system_results.csv')
    replaced={s['pdb_id']:s for s in es if s['pdb_id'] in scope}
    prior={p:{z['candidate_key']:z for z in systems if z['pdb_id']==p} for p in scope}
    after={p:{z['candidate_key']:z for z in ys if z['pdb_id']==p} for p in scope}
    proof={z['pdb_id']:z for z in read(delta/'independent_topology_readback.json')['entries']}
    assert set(replaced)==scope
    changes=[]
    for p in sorted(scope):
        assert set(prior[p])==set(after[p]) and proof[p]['artifacts']
        for c,row in after[p].items():
            orig=prior[p][c]; assert orig['role_qc']==row['role_qc'] and orig['revision']==row['revision']
            if row['status']!=orig['status']:
                assert any(x['file']==row['system_path'] and x['hash']==row['system_sha256'] for x in proof[p]['artifacts'])
                changes.append({'pdb_id':p,'candidate_key':c,'previous_status':orig['status'],'new_status':row['status'],
                                'role':row['role_qc'],'previous_candidate_intake':orig['candidate_intake'],'new_candidate_intake':row['candidate_intake'],
                                'original_reason':orig['reason'],'readback':f'point_symmetry_retest/{row["system_path"]}'})
    def tag(row, version):
        return row|{'evidence_layer':version,'result_root':str(delta if version=='point_symmetry_retest' else r)}
    final_e=[tag(replaced[z['pdb_id']],'point_symmetry_retest') if z['pdb_id'] in scope else tag(z,'main') for z in entries]
    final_s=[tag(after[z['pdb_id']][z['candidate_key']],'point_symmetry_retest') if z['pdb_id'] in scope else tag(z,'main') for z in systems]
    assert len({x['pdb_id'] for x in final_e})==len(final_e)==2000 and len(final_s)==len({(x['pdb_id'],x['candidate_key']) for x in final_s})
    write_csv(result/'entry_results.csv',final_e); write_csv(result/'system_results.csv',final_s)
    m3=read(r/'m3_first_review.json'); roster=m3['review_candidates']; paths=set(); raw_size=0; times=[]
    for c in roster:
        p=c['pdb_id']; assert p not in scope; d=read(r/'source_proofs'/f'{p}.json'); assert d['sha256']==c['source_sha256']; raw_size+=d['bytes']
        paths.update((c['system_cif'],c['mapping'])); times.append(float(read(r/'entry_results'/f'{p}.json')['elapsed_seconds']))
    out={'run':str(r),'final_entries_csv':'handoff/entry_results.csv','final_systems_csv':'handoff/system_results.csv',
         'frozen_input_sha256':sha(r/'frozen.json'),'first_pass_seal_sha256':sha(r/'blind_first_pass.json'),
         'baseline_code':read(r/'first_pass'/'build_inputs.json'),'targeted_version':read(delta/'build_inputs.json'),
         'first_pass':metrics(first_entries,first_systems),'final':metrics(final_e,final_s),
         'delta_only_affected_entries':sorted(scope),'delta_changed_candidates':changes,'m3_first_manual_review':'m3_first_review.json',
         'm3_first_manual_review_scope':{'entry_count':len(roster),'planned_additional_download_bytes':0,
            'verified_pinned_compressed_original_bytes':raw_size,'actual_existing_system_and_mapping_bytes':sum((r/p).stat().st_size for p in paths),
            'm2_single_worker_equivalent_seconds':round(sum(times),3),'m2_entry_median_seconds':round(statistics.median(times),3) if times else None,
            'meaning':'Historical per-entry timings cover M2 discovery and all systems for each selected Entry. M3 human decisions are not included.'},
         'scope_note':'Point symmetry refinements only where deposited endpoint identifiers and distances can be independently confirmed. New topology evidence is not a role sign-off.'}
    (result/'summary.json').write_text(json.dumps(out,indent=2),encoding='utf8')
    print('handoff reconciled:',len(final_e),'Entries',len(final_s),'systems',len(changes),'changed systems',len(roster),'M3 manual checks',flush=True)


if __name__=='__main__':main()

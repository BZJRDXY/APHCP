"""Summarize full Holo B pocket checks and reconcile the preselected sample."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3

from holo_b_science_sampling import eligible_candidate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--database',type=Path,required=True)
    p.add_argument('--sites',type=Path,required=True)
    p.add_argument('--sample-sites',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    db = sqlite3.connect(args.database.resolve().as_uri()+'?mode=ro',uri=True)
    db.row_factory = sqlite3.Row
    db.execute('attach database ? as site',(args.sites.resolve().as_uri()+'?mode=ro',))
    expected = db.execute("select count(*) from result where match_level='SAME_REVISION_INSTANCE_CHAIN_MATCH'").fetchone()[0]
    checked = db.execute('select count(*) from site.site_check').fetchone()[0]
    if checked != expected:
        raise ValueError(f'Incomplete full site check: {checked}/{expected}')
    sample = json.loads(args.sample_sites.read_text(encoding='utf-8'))
    compared = 0
    for row in sample:
        if row['state'] != 'COMPARED':
            continue
        record = db.execute('select details from site.site_check where system_id=?',(row['system_id'],)).fetchone()
        if not record:
            raise ValueError('A sampled system is missing from full site output')
        actual = json.loads(record[0])
        for key in ('state','mapped_near_residues','plinder_neighbor_residues','intersection'):
            if row[key] != actual[key]:
                raise ValueError(f'Sampled site differs: {row["system_id"]} {key}')
        compared += 1
    groups = {k:Counter() for k in ('all','role_only','candidate')}
    other_models = {k:Counter() for k in groups}
    strict_candidate = Counter()
    reasons = '["target_role_evidence_insufficient_or_conflicting"]'
    query = '''select result.*,site_check.model_evidence,site_check.state as site_state,
               site_check.pocket_agreement from result join site.site_check using(system_id)'''
    for row in db.execute(query):
        key = (row['model_evidence'],row['site_state'],row['pocket_agreement'])
        groups['all'][key] += 1
        scopes = ['all']
        if row['disposition'] == 'REVIEW' and row['reasons'] == reasons:
            groups['role_only'][key] += 1
            scopes.append('role_only')
            if eligible_candidate(row):
                groups['candidate'][key] += 1
                scopes.append('candidate')
                if row['model_evidence'] == 'SOURCE_SINGLE_MODEL_1' and row['pocket_agreement'] == 'EXACT_NONEMPTY':
                    strict_candidate['systems'] += 1
                    annotations = json.loads(row['annotations'])
                    strict_candidate['covalent'] += any(a.get('ligand_is_covalent') is True for a in annotations)
                    strict_candidate['NAG'] += any(a.get('ligand_ccd_code') == 'NAG' for a in annotations)
                    strict_candidate['competing_interpretation'] += 'competing' in ((json.loads(row['source_evidence']) or {}).get('reason') or '')
        if row['model_evidence'] != 'SOURCE_SINGLE_MODEL_1':
            qc = json.loads(row['new_qc']) if row['new_qc'] else None
            models = qc.get('all_deposited_models') if isinstance(qc,dict) else None
            classification = 'NOT_RECORDED' if models is None else 'MULTIPLE' if len(models)>1 else 'OTHER'
            for scope in scopes:
                other_models[scope][classification] += 1
    result = {'expected_same_revision_links':expected,'checked':checked,'sample_sites_reconciled':compared,
              'groups':{name:[{'model':key[0],'state':key[1],'pocket':key[2],'systems':count}
                               for key,count in sorted(counter.items())] for name,counter in groups.items()},
              'other_model_breakdown':{k:dict(v) for k,v in other_models.items()},
              'single_model_exact_candidates':dict(strict_candidate)}
    with args.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2)
    print(json.dumps(result),flush=True)


if __name__ == '__main__':
    main()

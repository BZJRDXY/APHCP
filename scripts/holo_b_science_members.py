"""Count distinct chain members and residue instances in the frozen Holo B audit."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sqlite3

from holo_b_science_sampling import eligible_candidate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--database', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--chain-check', action='store_true', help='check source-chain uniqueness for all strong links')
    args = p.parse_args()
    db = sqlite3.connect(args.database.resolve().as_uri()+'?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    if args.chain_check:
        keys = defaultdict(set)
        for row in db.execute("select system_id,document from observation join result using(system_id) where match_level='SAME_REVISION_INSTANCE_CHAIN_MATCH'"):
            i = json.loads(row['document'])['identity']; t = i['target_residue_instances'][0]
            keys[(i['pdb_id'],i['revision'],i['assembly'],tuple(i['model_scope']),t[0])].add(row['system_id'])
        identities = defaultdict(set)
        for row in db.execute('select document from observation'):
            i = json.loads(row['document'])['identity']
            base = (i['pdb_id'],i['revision'],i['assembly'],tuple(i['model_scope']))
            for t in i.get('target_residue_instances') or []:
                k = base+(t[0],)
                if k in keys:
                    identities[k].add((str(t[5]),str(t[6])))
        ambiguous = [k for k,v in identities.items() if len(v)>1]
        summary = {'strong_target_chain_keys':len(keys),'multi_residue_identity_keys':len(ambiguous),
                   'affected_systems':sum(len(keys[k]) for k in ambiguous),
                   'examples':[{'key':k,'identities':sorted(identities[k])} for k in ambiguous[:10]]}
        with args.output.open('x',encoding='utf-8') as f:
            json.dump(summary,f,indent=2)
        print(json.dumps(summary))
        return
    groups = {k:{'chains':set(),'residues':set()} for k in ('all','role_only','candidate')}
    reasons = '["target_role_evidence_insufficient_or_conflicting"]'
    query = 'select result.*,observation.document as observation_document from result join observation using(system_id)'
    for row in db.execute(query):
        obs = json.loads(row['observation_document'])
        i = obs['identity']
        base = (i['pdb_id'], i['revision'], str(i['assembly']), tuple(i.get('model_scope') or []))
        scopes = ['all']
        if row['disposition'] == 'REVIEW' and row['reasons'] == reasons:
            scopes.append('role_only')
            if eligible_candidate(row):
                scopes.append('candidate')
        for t in i.get('target_residue_instances') or []:
            if len(t) < 7:
                raise ValueError(f'Incomplete target identity: {row["system_id"]}')
            chain = base + tuple(t[:5])
            residue = base + tuple(t[:7])
            for scope in scopes:
                groups[scope]['chains'].add(chain)
                groups[scope]['residues'].add(residue)
    summary = {k:{'assembly_member_chains':len(v['chains']), 'residue_instances':len(v['residues'])}
               for k,v in groups.items()}
    with args.output.open('x',encoding='utf-8') as f:
        json.dump(summary,f,indent=2)
    print(json.dumps(summary))


if __name__ == '__main__':
    main()

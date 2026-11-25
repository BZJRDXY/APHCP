"""Fix a comparable one-candidate-per-source inspection budget and execute it."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import shutil
import subprocess
import sys

from apo_c_pilot import cache_get, now, read_json, sha, write_json


def latest_search(run, case):
    stem = case['id'].replace(':', '_')
    paths = sorted((run / 'search').glob(stem + '*.json'))
    if not paths:
        raise ValueError(f'no search result for {case["id"]}')
    return paths[-1], read_json(paths[-1])


def rcsb_chain(run, identifier):
    pdb, entity = identifier.split('_', 1)
    body, proof = cache_get(run, f'rcsb_entity/{pdb.lower()}_{entity}',
                            f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{entity}')
    container = json.loads(body)['rcsb_polymer_entity_container_identifiers']
    authors = sorted(set(container.get('auth_asym_ids') or []))
    if not authors:
        raise ValueError(f'RCSB entity {identifier} has no author chain IDs')
    return authors[0], {'all_author_chains': authors, 'data_api': proof}


def plan(run):
    path = run / 'budget' / 'selection.json'
    if path.exists():
        raise FileExistsError(f'budget selection already fixed: {path}')
    selected = []
    for case in read_json(run / 'samples.json'):
        search_path, data = latest_search(run, case)
        rows = []
        pl = data['plinder']
        links = sorted(pl.get('links', []), key=lambda x: (x['rank'], x['linked_structure_id']))
        link = next((x for x in links if x['source_entry_id'].upper() != case['pdb']), None)
        if link and link.get('source_chain_auth_id') and link.get('source_biounit_id'):
            rows.append({'route': 'PLINDER', 'status': 'selected', 'pdb': link['source_entry_id'].upper(),
                         'author_chain': link['source_chain_auth_id'], 'assembly': link['source_biounit_id'],
                         'source_rank': link['rank'], 'source_row': link})
        else:
            rows.append({'route': 'PLINDER', 'status': 'candidate_chain_unavailable' if link else pl['status']})
        ah = data['ahoj']
        hit = next((a for m in ah.get('matches', []) for a in m['apo'] if a['pdb'] != case['pdb']), None)
        if hit:
            rows.append({'route': 'AHoJ', 'status': 'selected', 'pdb': hit['pdb'],
                         'author_chain': hit['auth_chain'], 'assembly': '1', 'source_rank': hit['source_order'],
                         'assembly_note': 'Assembly 1 is an independent inspection choice; the online hit did not specify it.'})
        else:
            rows.append({'route': 'AHoJ', 'status': 'no_apo_sites' if ah.get('matches') else ah['status']})
        rs = data['sifts_rcsb']
        hit = next(((s, h) for s in rs.get('searches', []) for h in s['hits']
                    if h['identifier'].split('_')[0].upper() != case['pdb']), None)
        if hit:
            search, item = hit
            identifier = item['identifier']
            pdb = identifier.split('_')[0].upper()
            try:
                chain, proof = rcsb_chain(run, identifier)
                rows.append({'route': 'RCSB', 'status': 'selected', 'pdb': pdb, 'author_chain': chain,
                             'assembly': '1', 'source_entity': identifier, 'chain_selection': proof,
                             'query_id': search.get('query_id') or search.get('accession'),
                             'assembly_note': 'Assembly 1 is an independent inspection choice.'})
            except Exception as exc:
                rows.append({'route': 'RCSB', 'status': 'candidate_chain_technical_failure',
                             'source_entity': identifier, 'error': f'{type(exc).__name__}: {exc}'})
        else:
            rows.append({'route': 'RCSB', 'status': 'no_other_entry' if rs['status'].startswith('searched') else rs['status']})
        selected.append({'case': case['id'], 'split': case['split'], 'holo_pdb': case['pdb'],
                         'holo_binding_instances': case['identity']['binding_instances'],
                         'search_evidence': str(search_path), 'routes': rows})
    plan_data = {'created_utc': now(), 'budget': 'At most one ranked non-self candidate per source and Holo system.',
                 'selection_rules': ['PLINDER: lowest source rank.', 'AHoJ: first rendered Apo-site row.',
                     'RCSB: first non-self polymer entity in saved response order, first sorted author chain.',
                     'RCSB and AHoJ use Assembly 1 only as a disclosed inspection choice.',
                     'Identical Holo/entry/author-chain/Assembly requests share one inspection.',
                     'Candidate lists remain preserved in search results; budget does not delete them.'],
                 'cases': selected}
    write_json(path, plan_data)
    print(json.dumps({'selection': str(path), 'cases': len(selected),
                      'selected_source_rows': sum(r['status'] == 'selected' for x in selected for r in x['routes'])}))


def execute(run):
    selection = read_json(run / 'budget' / 'selection.json')
    if (run / 'budget' / 'complete.json').exists():
        raise FileExistsError('budget execution already has a completion receipt')
    log = run / 'budget' / 'execution.jsonl'
    prior = [json.loads(x) for x in log.read_text(encoding='utf8').splitlines()] if log.exists() else []
    done = {x['key'] for x in prior}
    tasks = {}
    for case in selection['cases']:
        for row in case['routes']:
            if row['status'] != 'selected':
                continue
            key = '|'.join([case['case'], row['pdb'], row['author_chain'], row['assembly']])
            task = tasks.setdefault(key, {'key': key, 'case': case['case'], 'pdb': row['pdb'],
                                         'author_chain': row['author_chain'], 'assembly': row['assembly'],
                                         'holo_chain_count': len(case['holo_binding_instances']), 'routes': []})
            task['routes'].append(row['route'])
    current_dir = Path(__file__).resolve().parent
    code_dir = run / 'code'
    if not (code_dir / f'{Path(__file__).stem}-{sha(Path(__file__))[:12]}.py').exists():
        shutil.copy2(Path(__file__), code_dir / f'{Path(__file__).stem}-{sha(Path(__file__))[:12]}.py')
    for key, task in tasks.items():
        if key in done:
            print('BUDGET_REUSE', key, flush=True)
            continue
        row = dict(task, started_utc=now())
        if task['holo_chain_count'] > 3:
            row.update({'execution_status': 'NOT_PROCESSED_MULTICHAIN_SCOPE',
                        'reason': 'More than three Holo binding instances; separate bounded mapping required.'})
        else:
            multi = task['holo_chain_count'] > 1
            script = current_dir / ('apo_c_multichain.py' if multi else 'apo_c_pair.py')
            if multi:
                prior_dir = run / 'multichain' / (task['case'].replace(':', '_') + '__' + task['pdb'])
            else:
                prior_dir = run / 'pairs' / (task['case'].replace(':', '_') + '__' + task['pdb'] + '_' + task['author_chain'])
            reusable = None
            for f in sorted(prior_dir.glob('attempt*/pair.json'), reverse=True) if prior_dir.exists() else []:
                evidence = read_json(f)
                req = evidence.get('request', {})
                if (evidence.get('execution_status') == 'INSPECTED'
                        and evidence.get('dependency_sha256', {}).get(script.name) == sha(script)
                        and req.get('assembly') == task['assembly']
                        and (not multi or (req.get('anchor_auth_chain') == task['author_chain']
                                           and not req.get('simulate_missing_chain')))):
                    reusable = (f, evidence)
                    break
            if reusable:
                f, evidence = reusable
                row.update({'execution_status': evidence['execution_status'], 'evidence': str(f),
                            'outcome': evidence.get('geometric_screen') or evidence.get('outcome'), 'reused': True})
            else:
                cmd = [sys.executable, '-B', str(script), '--run', str(run), '--case', task['case'],
                       '--candidate', task['pdb'], '--assembly', task['assembly']]
                if multi:
                    cmd += ['--anchor-auth-chain', task['author_chain']]
                else:
                    cmd += ['--author-chain', task['author_chain'], '--source', '+'.join(task['routes'])]
                try:
                    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
                    text = proc.stdout.strip().splitlines()
                    summary = json.loads(text[-1]) if text else {}
                    evidence_path = str(Path(summary['path']) / 'pair.json') if summary.get('path') else None
                    row.update({'execution_status': summary.get('status', 'TECHNICAL_FAILURE'),
                                'outcome': summary.get('screen') or summary.get('outcome'),
                                'evidence': evidence_path, 'return_code': proc.returncode,
                                'error': summary.get('error') or proc.stderr[-1000:] or None, 'reused': False})
                except Exception as exc:
                    row.update({'execution_status': 'TECHNICAL_FAILURE', 'error': f'{type(exc).__name__}: {exc}'})
        row['finished_utc'] = now()
        with log.open('a', encoding='utf8') as f:
            f.write(json.dumps(row, ensure_ascii=False) + '\n')
        prior.append(row)
        print('BUDGET', key, row['execution_status'], row.get('outcome'), flush=True)
    counts = Counter(x['execution_status'] for x in prior)
    write_json(run / 'budget' / 'complete.json', {'finished_utc': now(), 'unique_tasks': len(tasks),
               'recorded_tasks': len(prior), 'execution_status': dict(counts),
               'reused': sum(bool(x.get('reused')) for x in prior)})


def main():
    p = argparse.ArgumentParser()
    p.add_argument('stage', choices=['plan', 'execute'])
    p.add_argument('--run', required=True)
    args = p.parse_args()
    if args.stage == 'plan':
        plan(Path(args.run))
    else:
        execute(Path(args.run))


if __name__ == '__main__':
    main()

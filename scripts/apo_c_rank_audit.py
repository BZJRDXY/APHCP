"""Freeze source rankings and replay known, explicitly non-blind reference cases."""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
import json
import shutil

from apo_c_pilot import now, read_json, sha, write_json


SELECTED = {
    '3ABX': 'holo:3ABX:913b8553b26ddb99c6fbc12136c1',
    '2SRC': 'holo:2SRC:415139e6df7bbbc2f01d4c41dde7',
    '1AKE': 'holo:1AKE:9086db90a158337553fbbd5ec0d8',
}


def unique_entry_ranking(rows, holo):
    grouped = {}
    for row in rows:
        pdb = row['pdb'].upper()
        if pdb == holo:
            continue
        if pdb not in grouped:
            grouped[pdb] = {'pdb': pdb, 'entry_rank': len(grouped) + 1, 'source_rows': []}
        grouped[pdb]['source_rows'].append(row)
    return list(grouped.values())


def freeze(run):
    output = run / 'rank_comparison' / 'frozen_candidates.json'
    if output.exists():
        raise FileExistsError('candidate comparison already frozen')
    origin = Path(read_json(run / 'run.json')['baseline_run'])
    cases = []
    for holo, system in SELECTED.items():
        files = sorted((origin / 'search').glob(system.replace(':', '_') + '*.json'))
        if not files:
            raise ValueError(f'no first-run search evidence for {system}')
        path = files[-1]
        d = read_json(path)
        pl = [{'pdb': x['source_entry_id'], 'author_chain': x['source_chain_auth_id'],
               'assembly': x['source_biounit_id'], 'source_rank': x['rank'],
               'chain_instance': x['source_chain_instance']}
              for x in sorted(d['plinder'].get('links', []), key=lambda z: (z['rank'], z['linked_structure_id']))]
        ah = [{'pdb': x['pdb'], 'author_chain': x['auth_chain'], 'source_rank': x['source_order'],
               'source_target': match['entry']['entry_key'], 'assembly': None}
              for match in d['ahoj'].get('matches', []) for x in match['apo']]
        rc = [{'pdb': h['identifier'].split('_')[0], 'entity': h['identifier'],
               'query_id': s.get('query_id') or s.get('accession'), 'source_order': i + 1}
              for s in d['sifts_rcsb'].get('searches', []) for i, h in enumerate(s['hits'])]
        rankings = {'PLINDER': unique_entry_ranking(pl, holo), 'AHoJ': unique_entry_ranking(ah, holo),
                    'RCSB': unique_entry_ranking(rc, holo)}
        union = defaultdict(dict)
        for source, items in rankings.items():
            for item in items[:10]:
                union[item['pdb']][source] = item
        cases.append({'holo': holo, 'system_id': system, 'search_evidence': str(path),
                      'search_sha256': sha(path), 'source_status': {
                          'PLINDER': d['plinder']['status'], 'AHoJ': d['ahoj']['status'],
                          'RCSB': d['sifts_rcsb']['status']},
                      'source_unique_entry_counts': {k: len(v) for k, v in rankings.items()},
                      'top10_rankings': {k: v[:10] for k, v in rankings.items()},
                      'top10_union': dict(sorted(union.items()))})
    result = {'created_utc': now(), 'baseline_run': str(origin), 'cases': cases,
              'unit': 'Unique PDB Entry for preliminary ranking only; chain/Assembly observations remain distinct.',
              'rules': ['Exclude the query Entry before ranking.', 'Retain all source rows for each top-ten Entry.',
                        'All original full candidate lists remain in search_evidence.',
                        'AHoJ did not specify the candidate Assembly; no chain-equivalence assumption is made.',
                        'This fixed reference set was known before freezing and is not a blinded evaluation.']}
    write_json(output, result)
    code = run / 'code' / f'{Path(__file__).stem}-{sha(Path(__file__))[:12]}.py'
    if not code.exists():
        shutil.copy2(Path(__file__), code)
    print(json.dumps({'path': str(output), 'cases': len(cases),
                      'union_entries': sum(len(x['top10_union']) for x in cases)}), flush=True)


def replay(run, reference):
    output = run / 'rank_comparison' / 'known_reference_replay.json'
    if output.exists():
        raise FileExistsError('reference ranking replay already exists')
    frozen_path = run / 'rank_comparison' / 'frozen_candidates.json'
    frozen, refs = read_json(frozen_path), read_json(reference)
    ranks = {}
    for case in frozen['cases']:
        ranks[case['holo']] = {source: {v['pdb']: v['entry_rank'] for v in values}
                               for source, values in case['top10_rankings'].items()}
    groups = ['experimental_unbound_supported_pair_scope_pending', 'unbound_evidence_pending',
              'other_ligand_occupied', 'wrong_domain', 'different_species_homolog']
    summaries = []
    for k in [1, 3, 10]:
        by_source = {}
        for source in ['RCSB', 'PLINDER', 'AHoJ']:
            found = [r for r in refs['references'] if ranks[r['holo']][source].get(r['candidate'], 999) <= k]
            hit_sets = {g: sorted(r['holo'] + '>' + r['candidate'] for r in found if r['status'] == g) for g in groups}
            unknown = 0
            for case in frozen['cases']:
                reference_ids = {r['candidate'] for r in refs['references'] if r['holo'] == case['holo']}
                unknown += sum(x['pdb'] not in reference_ids for x in case['top10_rankings'][source][:k])
            by_source[source] = {'known_reference_hits': hit_sets, 'unadjudicated_entries': unknown}
        r = set(by_source['RCSB']['known_reference_hits'][groups[0]])
        pl = set(by_source['PLINDER']['known_reference_hits'][groups[0]])
        ah = set(by_source['AHoJ']['known_reference_hits'][groups[0]])
        summaries.append({'top_k': k, 'sources': by_source,
                          'PLINDER_unbound_refs_beyond_RCSB': sorted(pl - r),
                          'AHoJ_unbound_refs_beyond_RCSB_and_PLINDER': sorted(ah - r - pl)})
    out = {'created_utc': now(), 'frozen_candidates': str(frozen_path), 'frozen_sha256': sha(frozen_path),
           'reference_file': str(reference), 'reference_sha256': sha(reference),
           'reference_count': len(refs['references']), 'reference_method': refs['method'],
           'denominator_unbound_supported': sum(r['status'] == groups[0] for r in refs['references']),
           'summaries': summaries,
           'limits': ['Reference cases were known before the comparison and were not independently blinded.',
                      'Entry-level ranking is a screening proxy; specific chain/Assembly observations require separate checks.',
                      'Most candidate Entries remain unadjudicated; no population accuracy or complete recall is estimated.',
                      'Experimental unbound support does not by itself establish a formally accepted APHCP Pair.']}
    write_json(output, out)
    print(json.dumps({'path': str(output), 'references': len(refs['references']),
                      'unbound_denominator': out['denominator_unbound_supported']}), flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('stage', choices=['freeze', 'replay'])
    p.add_argument('--run', required=True)
    p.add_argument('--reference-file')
    args = p.parse_args()
    run = Path(args.run)
    if args.stage == 'freeze':
        freeze(run)
    else:
        if not args.reference_file:
            p.error('replay requires --reference-file')
        ref = Path(args.reference_file)
        if not ref.resolve().is_relative_to(run.resolve()):
            raise ValueError('reference file must be in the research Run')
        replay(run, ref)


if __name__ == '__main__':
    main()

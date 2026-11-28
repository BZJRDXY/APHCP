"""Compare sealed Apo evaluation runs and independently read delivered structures."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from functools import cache
import json
import re
import subprocess
from pathlib import Path

import gemmi

from verify_apo_c_autopilot import audit_delivery_assembly, check_observation, digest, instances


OBS_KEYS = ('model', 'instance', 'label_asym', 'author_chain', 'assembly')
ASSEMBLY_CHECKS = ('source_identity_revision_and_hash_ok', 'output_model_matches',
                   'atom_ids_and_altloc_ok', 'coordinates_ok', 'occupancy_ok',
                   'applicable_connection_endpoints_ok', 'hash_ok')
GEOMETRY_CHECKS = ('holo_site_matches_source', 'candidate_site_observed_count_matches',
                   'candidate_contact_set_and_distances_match')
RRF_K = 60
RANKING_CONTROLS = (('3Q71', '4ABK', 'SCIENTIFIC_EXCLUSION', 'OCCUPIED', 'OCCUPIED_MAPPED_SUBSITE'),
                    ('2E3N', '2E3S', 'REVIEW', 'INCOMPLETE', 'MAPPED_SITE_NOT_FULLY_OBSERVED'))
RANKING_REFERENCE_REQUESTS = (('1AQ1', '3PXR'), ('1AQ1', '1HCL'), ('1AQ1', '1BUH'), ('1AQ1', '5IF1'))
source_digest = cache(digest)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_batch(root):
    shards = sorted(p for p in root.glob('shard-*') if p.is_dir())
    cases, literature, shard_seconds = {}, Counter(), {}
    for shard in shards:
        receipt = read_json(max((shard / 'receipts').glob('automation_run*.json')))
        literature.update(receipt['literature'])
        shard_seconds[shard.name] = receipt['elapsed_s']
        samples = read_json(shard / 'samples.json')
        require(receipt['holo_dispatched'] == len(samples), f'{shard}: sample count differs from receipt')
        for sample in samples:
            target = sample['id']
            require(target not in cases, f'duplicate Holo target {target}')
            directory = shard / 'automated' / target.replace(':', '_')
            result_path = max([directory / 'result.json', *directory.glob('result.retry*.json')])
            result = read_json(result_path)
            require(result['system_id'] == target and result['holo'] == sample['pdb'],
                    f'{result_path}: result target identity differs')
            require(source_digest(result['search_evidence']) == result['search_sha256'],
                    f'{result_path}: search checksum differs')
            require(source_digest(sample['source_uri']) == sample['identity']['source_sha256'],
                    f'{result_path}: Holo checksum differs')
            holo_system = result['holo_system']
            require(source_digest(holo_system['path']) == holo_system['sha256'] and
                    source_digest(sample['system_uri']) == holo_system['sha256'],
                    f'{result_path}: Holo system checksum differs')
            selected = [row['pdb'] for row in result['selected_candidates']]
            observed = [row['pdb'] for row in result['candidates']]
            require(selected == observed and len(selected) == len(set(selected)) <= 6,
                    f'{result_path}: candidate order, unique count or quota differs')
            candidates = {row['pdb']: read_json(row['evidence_path']) for row in result['candidates']}
            cases[target] = (sample, result, candidates, shard.name)
    return cases, len(shards), dict(literature), shard_seconds


def source_identity(evidence):
    source = evidence.get('candidate_source')
    if source is None:
        return None
    require(source.get('revision') and source_digest(source['path']) == source['sha256'],
            f'{source.get("path")}: fixed candidate revision or checksum invalid')
    return source['revision'], source['sha256']


def compare_cases(baseline, replay):
    require(list(baseline) == list(replay), 'Holo order or identity differs across shards')
    checked_sources = 0
    for target, (old_sample, old, old_candidates, old_shard) in baseline.items():
        sample, result, candidates, shard = replay[target]
        require(old_shard == shard and old_sample['identity'] == sample['identity'],
                f'{target}: frozen Holo identity or shard differs')
        require(old['holo_system']['sha256'] == result['holo_system']['sha256'] and
                old['holo_evidence']['site_residues'] == result['holo_evidence']['site_residues'],
                f'{target}: Holo system or observed site differs')
        require(old['search_sha256'] == result['search_sha256'] and
                old['selected_candidates'] == result['selected_candidates'],
                f'{target}: frozen search or ordered candidates differ')
        for code, old_evidence in old_candidates.items():
            source = source_identity(old_evidence)
            require(source == source_identity(candidates[code]), f'{target}/{code}: frozen candidate source differs')
            checked_sources += source is not None
    return checked_sources


def frozen_source_ranks(search, holo):
    plinder = sorted(search['plinder'].get('links', []), key=lambda row: (row['rank'], row['linked_structure_id']))
    routes = {
        'PLINDER': (row['source_entry_id'] for row in plinder),
        'AHoJ': (row['pdb'] for match in search['ahoj'].get('matches', []) for row in match['apo']),
        'RCSB': (hit['identifier'].split('_', 1)[0] for query in search['sifts_rcsb'].get('searches', [])
                 for hit in query['hits']),
    }
    ranks = defaultdict(dict)
    for route, entries in routes.items():
        unique = dict.fromkeys(entry.upper() for entry in entries if entry.upper() != holo)
        for rank, entry in enumerate(unique, 1):
            ranks[entry][route] = rank
    return ranks


def rrf_score(source_ranks):
    return sum(Fraction(1, RRF_K + rank) for rank in source_ranks.values())


def independent_fusion(ranks):
    return sorted(ranks, key=lambda entry: (-rrf_score(ranks[entry]), min(ranks[entry].values()), entry))


@cache
def check_ranked_source(path, code, revision):
    block = gemmi.cif.read(path).sole_block()
    rows = block.find('_pdbx_audit_revision_history.', ['data_content_type', 'major_revision', 'minor_revision'])
    versions = [(int(row[1]), int(row[2])) for row in rows
                if gemmi.cif.as_string(row[0]).casefold() == 'structure model']
    require(gemmi.cif.as_string(block.find_value('_entry.id') or '').upper() == code and versions and
            '.'.join(map(str, max(versions))) == revision, f'{path}: source Entry or revision differs')


def ranking_configs(baseline_root, replay_root, shards):
    for shard in sorted(shards):
        before, after = baseline_root / shard, replay_root / shard
        old, new = (read_json(p / 'frozen_config.json') for p in (before, after))
        require(new['rerank_from_run'] == str(before) and new['prior_literature_cache'] == old['literature_cache'],
                f'{shard}: ranking input or prior literature cache differs')
        for key in ('catalog', 'index', 'sample_spec_sha256', 'scientific_rule', 'formal_pair_acceptance'):
            require(old[key] == new[key], f'{shard}: frozen configuration changed: {key}')
        old_code, new_code = (dict(config['code_sha256']) for config in (old, new))
        for code in (old_code, new_code):
            del code['apo_c_autopilot.py']
        require(old_code == new_code, f'{shard}: frozen scientific or search dependency changed')


def compare_ranking(baseline, replay):
    require(list(baseline) == list(replay), 'Holo order or identity differs across shards')
    fixed_entries = defaultdict(set)
    for _, _, candidates, _ in baseline.values():
        for code, evidence in candidates.items():
            if source := source_identity(evidence):
                fixed_entries[code].add(source)
    counts, statuses, cases, ranks_by_target = Counter(), Counter(), [], {}
    for target, (old_sample, old, old_candidates, old_shard) in baseline.items():
        sample, result, candidates, shard = replay[target]
        require(old_sample == sample and old_shard == shard, f'{target}: frozen Holo input or shard differs')
        require(old['holo_system']['sha256'] == result['holo_system']['sha256'] and
                old['holo_evidence'] == result['holo_evidence'], f'{target}: frozen Holo evidence differs')
        require(old['search_sha256'] == result['search_sha256'], f'{target}: frozen search differs')
        ranks = frozen_source_ranks(read_json(old['search_evidence']), sample['pdb'])
        fused = independent_fusion(ranks)
        ranks_by_target[target] = {entry: rank for rank, entry in enumerate(fused, 1)}
        quota = len(old_candidates)
        priority = {'method': 'reciprocal_rank_fusion', 'k': RRF_K, 'candidate_quota': quota}
        require(list(candidates) == fused[:quota] and len(candidates) == quota and
                result['candidate_priority'] == priority,
                f'{target}: order or metadata disagrees with independent RRF and frozen quota')
        for selection in result['selected_candidates']:
            actual = {route: row['entry_rank'] for route, row in selection['sources'].items()}
            require(actual == ranks[selection['pdb']],
                    f'{target}/{selection["pdb"]}: selected source ranks differ from frozen search')
        added = [code for code in candidates if code not in old_candidates]
        dropped = [code for code in old_candidates if code not in candidates]
        counts['replaced_slots'] += len(added)
        counts['targets_with_replacements'] += bool(added)
        old_deliveries = {delivery['candidate'] for delivery in old['deliveries']}
        new_deliveries = {delivery['candidate'] for delivery in result['deliveries']}
        rows = []
        for rank, (code, evidence) in enumerate(candidates.items(), 1):
            require(evidence['pdb'] == code, f'{target}/{code}: candidate evidence Entry differs')
            source = source_identity(evidence)
            outcome = relation_outcome(evidence)
            if code in old_candidates:
                previous = old_candidates[code]
                require(source == source_identity(previous), f'{target}/{code}: frozen candidate source differs')
                ignored = {'replayed_from', 'materialized_assemblies'}
                current_fields = {key: evidence[key] for key in evidence.keys() - ignored}
                previous_fields = {key: previous[key] for key in previous.keys() - ignored}
                require(current_fields == previous_fields, f'{target}/{code}: frozen observation evidence differs')
                source_kind = 'overlap_fixed' if source else 'overlap_gap'
            elif source:
                check_ranked_source(evidence['candidate_source']['path'], code, source[0])
                source_kind = 'added_reuses_baseline_entry' if source in fixed_entries[code] else 'added_original'
            else:
                source_kind = 'added_gap'
            counts[source_kind] += 1
            statuses[outcome['scientific_status']] += 1
            rows.append({'pdb': code, 'rank': rank, 'source_ranks': ranks[code],
                         'rrf_score': str(rrf_score(ranks[code])), 'source_kind': source_kind,
                         'delivered': code in new_deliveries, **outcome})
        cases.append({'target': target, 'holo': sample['pdb'], 'quota': quota, 'search_sha256': old['search_sha256'],
                      'baseline_candidates': list(old_candidates), 'selected': rows, 'added': added, 'dropped': dropped,
                      'retained_pairs': sorted(old_deliveries & new_deliveries),
                      'new_pairs': sorted(new_deliveries - old_deliveries), 'lost_pairs': sorted(old_deliveries - new_deliveries)})
    return {'counts': dict(counts), 'scientific_status_counts': dict(statuses), 'cases': cases}, ranks_by_target


def atom_names(residue):
    return {atom.name for atom in residue if atom.occ > 0 and not atom.element.is_hydrogen}


def check_mapped_atoms(sample, candidate_source, observation):
    holo = gemmi.read_structure(sample['source_uri'])
    candidate = gemmi.read_structure(candidate_source['path'])
    identity = sample['identity']
    _, holo_rows = instances(holo, identity['assembly'])[identity['binding_instances'][0][0]]
    _, candidate_rows = instances(candidate, observation['assembly'])[observation['instance']]
    require(candidate[0].num == observation['model'], 'delivered model differs from verified source model')
    holo_residues = {res.label_seq: res for _, res in holo_rows if res.label_seq is not None}
    candidate_residues = {res.label_seq: (author, res) for author, res in candidate_rows if res.label_seq is not None}
    mapped = observation['mapped_site']
    require(mapped and len({row['candidate_label_seq'] for row in mapped}) == len(mapped),
            'site mapping is empty or repeats candidate residues')
    for row in mapped:
        holo_residue = holo_residues[row['holo_label_seq']]
        author, candidate_residue = candidate_residues[row['candidate_label_seq']]
        require(candidate_residue.subchain == observation['label_asym'] and author == observation['author_chain'] and
                candidate_residue.name == holo_residue.name and str(candidate_residue.seqid) == row['author_seq'],
                f'wrong chain copy or residue identity at Holo label {holo_residue.label_seq}')
        names = atom_names(holo_residue)
        require(names and names <= atom_names(candidate_residue),
                f'missing heavy atoms at Holo label {holo_residue.label_seq}')
    return len(mapped)


def check_experimental_citation(fact, delivery):
    """Independent rejection of an absent result and verification of source-row direction."""
    basis = fact.get('basis')
    text = ' '.join(str(fact.get('excerpt') or '').split())
    entry = delivery['candidate']
    if basis == 'MATCHED_APO_DEPOSITION_AND_PRIMARY_ABSTRACT':
        absence = re.compile(r'\b(?:no|neither|never|lack(?:ed|ing)?)\b.{0,65}\b(?:structures?|available|reported)\b|'
                             r'\b(?:structures?|data)\b.{0,60}\b(?:unavailable|unknown|not (?:yet )?(?:known|available|reported))\b', re.I)
        require(bool(text) and not absence.search(text), 'cited abstract describes missing or unavailable evidence')
    if basis == 'FIXED_HOLO_DEPOSITION_NAMES_CANDIDATE_STATE':
        source = Path(delivery['holo_original'])
        require(fact['entry'] == entry and source_digest(source) == fact['source_sha256'] == delivery['holo_identity']['source_sha256']
                and fact['source_revision'] == delivery['holo_identity']['revision'], 'wrong fixed deposit for the experimental state')
        block = gemmi.cif.read(str(source)).sole_block()
        rows = block.get_mmcif_category('_pdbx_database_related.')
        reported = [(i, str(d)) for i, (db, k, d) in enumerate(zip(rows['db_name'], rows['db_id'], rows['details']), 1)
                    if str(db).upper() == 'PDB' and str(k).upper() == entry]
        require(str(block.find_value('_entry.id')).upper() == fact['source_entry'] and reported and
                [i for i, _ in reported] == fact['row_numbers'] and text == ' '.join(reported[0][1].split()) and
                all(re.search(r'\b(?:apo|unbound|unliganded|uninhibited|(?:ligand|inhibitor)[- ]free)\b', d, re.I)
                    and not re.search(r'\b(?:not|negated|template|complexed|bound to|inhibited|occupied)\b', d, re.I) for _, d in reported),
                'deposition row does not identify this candidate as experimentally unbound')
    if basis == 'AUTHOR_SUPPLEMENT_ENUMERATES_UNBOUND_STRUCTURES':
        loc = fact['location']; src = fact['source_receipt']; pdf = Path(src['path'])
        require(src['doi'] == fact['doi'] and source_digest(pdf) == src['sha256'], 'supplement hash or paper identity mismatch')
        content = subprocess.run(['pdftotext', '-layout', '-f', str(loc['page']), '-l', str(loc['page']), str(pdf), '-'],
                                 check=True, capture_output=True, timeout=30).stdout.decode('utf8')
        normalized = ' '.join(content.split())
        require(text in normalized and text.startswith(loc['figure'] + ' ') and re.search(r'(?:unligated|unbound|apo|unliganded) (?:form|state)', text, re.I)
                and 'PDB' in text and entry in re.findall(r'\b[0-9][A-Z0-9]{3}\b', text.split('PDB')[-1]),
                'figure and candidate group cannot be reproduced from the cited page')


def check_delivery(sample, result, evidence, delivery):
    code = delivery['candidate']
    source = evidence['candidate_source']
    identity = sample['identity']
    copy = delivery['apo_observation']
    require(delivery['system_id'] == sample['id'] and delivery['holo_identity'] == identity and
            delivery['apo_original'] == source and delivery['holo_original'] == sample['source_uri'],
            'delivery target or source differs from frozen evidence')
    experiment = delivery['candidate_experiment']
    require(experiment == evidence['scientific_decision']['candidate_experiment'] and
            experiment['status'] == 'EXPLICIT_UNBOUND' and experiment['source_entry'].upper() == code and
            experiment['fixed_revision'] == source['revision'] and experiment['source_sha256'] == source['sha256'] and
            experiment['evidence'] and not experiment.get('literature', {}).get('conflict'),
            'unbound evidence is missing, conflicting, or belongs to another candidate')
    for fact in experiment['evidence']:
        check_experimental_citation(fact, delivery)
    matches = [(entity, obs) for entity in evidence['entities'] for obs in entity['assembly_observations']
               if all(obs.get(key) == copy[key] for key in OBS_KEYS)]
    require(len(matches) == 1, 'delivery does not identify exactly one inspected chain copy')
    entity, observation = matches[0]
    decision = observation['decision']
    require(decision['delivery_eligible'] and decision['status'] == 'SUPPORTED_RESEARCH_PAIR' and
            decision['site_observation']['status'] == 'OBSERVED_NO_COMPETING_OCCUPANCY' and
            not decision['unknown_holo_positions'] and decision == delivery['assessment'] and
            observation['mapped_site'] == delivery['mapped_site'], 'ineligible or altered observation')
    holo = result['holo_evidence']
    require(sorted(row['holo_label_seq'] for row in observation['mapped_site']) ==
            sorted(row['label_seq'] for row in holo['site_residues']), 'mapping does not cover the complete Holo site')
    geometry = check_observation(sample, holo, evidence, entity, observation)
    require(all(geometry[key] for key in GEOMETRY_CHECKS), 'source geometry differs from the delivered observation')
    mapped_count = check_mapped_atoms(sample, source, observation)
    holo_check = audit_delivery_assembly(
        delivery['holo_original'], {'pdb': identity['pdb_id'], 'revision': identity['revision'],
                                    'sha256': identity['source_sha256']}, identity['assembly'], holo['model'],
        delivery['holo_assembly'], {edge['id'] for edge in holo['target_connections']})
    apo_check = audit_delivery_assembly(
        source['path'], {'pdb': code, 'revision': source['revision'], 'sha256': source['sha256']},
        copy['assembly'], copy['model'], delivery['apo_assembly'])
    require(all(item['atoms'] and all(item[key] for key in ASSEMBLY_CHECKS) for item in (holo_check, apo_check)),
            'assembly readback differs from fixed originals')
    return {'target': sample['id'], 'holo': sample['pdb'], 'candidate': code, 'observation': copy,
            'mapped_residues_verified': mapped_count, 'geometry': geometry,
            'holo_assembly': holo_check, 'apo_assembly': apo_check, 'experiment': experiment['evidence']}


def delivered_pairs(cases):
    return {(target, delivery['candidate']) for target, (_, result, _, _) in cases.items()
            for delivery in result['deliveries']}


def delivered_copy_count(cases):
    return sum(len(result['deliveries']) for _, result, _, _ in cases.values())


def target_outcomes(cases, pairs, requests):
    outcomes = []
    for holo, code in requests:
        hits = []
        for target, (sample, result, candidates, _) in cases.items():
            if sample['pdb'] != holo:
                continue
            evidence = candidates.get(code, {})
            decisions = [row['decision'] for entity in evidence.get('entities', [])
                         for row in entity['assembly_observations'] if 'decision' in row]
            hits.append({'target': target, 'selected': code in candidates, 'delivered': (target, code) in pairs,
                         'status': evidence.get('scientific_decision', {}).get('status'),
                         'site_statuses': sorted({row['site_observation']['status'] for row in decisions}),
                         'reasons': sorted({reason for row in decisions for reason in row['reasons']}),
                         'delivery_copies': sum(d['candidate'] == code for d in result['deliveries'])})
        outcomes.append({'holo': holo, 'candidate': code, 'matches': hits})
    return outcomes


def check_ranking_controls(baseline, replay, pairs, ranks_by_target):
    controls = []
    for holo, code, status, site_status, reason in RANKING_CONTROLS:
        matching = [target for target, (sample, _, candidates, _) in baseline.items()
                    if sample['pdb'] == holo and code in candidates]
        require(len(matching) == 1, f'{holo}/{code}: frozen control is not unique')
        target = matching[0]
        selected = code in replay[target][2]
        sample, result, candidates, _ = (replay if selected else baseline)[target]
        evidence = candidates[code]
        require((target, code) not in pairs and evidence['scientific_decision']['status'] == status,
                f'{target}/{code}: control was delivered or changed scientific status')
        observations = []
        for entity in evidence['entities']:
            for observation in entity['assembly_observations']:
                decision = observation.get('decision', {})
                if decision.get('site_observation', {}).get('status') != site_status or reason not in decision['reasons']:
                    continue
                geometry = check_observation(sample, result['holo_evidence'], evidence, entity, observation)
                require(all(geometry[key] for key in GEOMETRY_CHECKS), f'{target}/{code}: control geometry differs')
                recorded_missing = decision['unknown_holo_positions']
                if site_status == 'OCCUPIED':
                    require(geometry['actual_nonwater_within_4A'] > 0, f'{target}/{code}: no actual occupied site')
                else:
                    candidate = gemmi.read_structure(evidence['candidate_source']['path'])
                    _, rows = instances(candidate, observation['assembly'])[observation['instance']]
                    residues = {res.label_seq: res for _, res in rows if res.label_seq is not None}
                    missing = [row['holo_label_seq'] for row in observation['mapped_site']
                               if row['candidate_label_seq'] not in residues or
                               not atom_names(residues[row['candidate_label_seq']])]
                    require(132 in recorded_missing and 132 in missing, f'{target}/{code}: position 132 is not absent')
                observations.append({'geometry': geometry, 'unknown_holo_positions': recorded_missing,
                                     'reasons': decision['reasons']})
        require(observations, f'{target}/{code}: no observation established control status {site_status}')
        controls.append({'target': target, 'holo': holo, 'candidate': code, 'selected': selected,
                         'rrf_rank': ranks_by_target[target].get(code), 'delivered': False, 'status': status,
                         'site_status': site_status, 'geometry_input': 'replay' if selected else 'baseline_only',
                         'observations': observations})
    return controls


def relation_outcome(evidence):
    scientific = evidence['scientific_decision']
    decisions = [row['decision'] for entity in evidence.get('entities', [])
                 for row in entity['assembly_observations'] if 'decision' in row]
    return {'entry_status': evidence['status'], 'scientific_status': scientific['status'],
            'experiment': scientific.get('candidate_experiment', {}).get('status'),
            'chain_copy_reasons': dict(Counter(reason for decision in decisions for reason in decision['reasons']))}


def relation_changes(baseline, replay):
    transitions = Counter()
    changes = []
    for target, (sample, _, old_candidates, _) in baseline.items():
        for code, evidence in replay[target][2].items():
            old, new = relation_outcome(old_candidates[code]), relation_outcome(evidence)
            transitions[(old['scientific_status'], new['scientific_status'])] += 1
            if old != new:
                changes.append({'target': target, 'holo': sample['pdb'], 'candidate': code, 'baseline': old, 'replay': new})
    return {'status_transitions': [{'baseline': before, 'replay': after, 'count': count}
                                   for (before, after), count in sorted(transitions.items())],
            'changed': len(changes), 'relations': changes}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True, type=Path)
    parser.add_argument('--replay', required=True, type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--ranking-replay', action='store_true', help='audit an equal-quota RRF selection replay')
    parser.add_argument('--target', action='append', nargs=2, metavar=('HOLO', 'CANDIDATE'), default=[])
    parser.add_argument('--expected-shards', type=int, default=10)
    parser.add_argument('--expected-holo', type=int, default=100)
    parser.add_argument('--expected-candidates', type=int, default=444)
    parser.add_argument('--expected-sources', type=int, default=415)
    args = parser.parse_args()
    baseline_root, replay_root = args.baseline.resolve(), args.replay.resolve()
    default_output = 'ranking_replay_audit.json' if args.ranking_replay else 'evidence_replay_audit.json'
    output = (args.output or replay_root / 'quality' / default_output).resolve()
    if baseline_root == replay_root or output.is_relative_to(baseline_root) or output.exists():
        parser.error('choose a distinct replay Run and an unused audit path outside the baseline')
    receipt = {'created_utc': datetime.now(timezone.utc).isoformat(), 'code_sha256': digest(__file__),
               'baseline': str(baseline_root), 'replay': str(replay_root), 'errors': [], 'delivery_checks': []}
    errors = receipt['errors']
    checked_errors = (KeyError, ValueError, RuntimeError, OSError, StopIteration)
    try:
        baseline, old_shards, baseline_literature, _ = read_batch(baseline_root)
        replay, shards, literature, seconds = read_batch(replay_root)
        require(old_shards == shards == args.expected_shards, 'unexpected shard count')
        selected_count = sum(len(candidates) for _, _, candidates, _ in replay.values())
        require(len(baseline) == len(replay) == args.expected_holo and selected_count == args.expected_candidates,
                'unexpected target or inspected candidate count')
        if args.ranking_replay:
            ranking_configs(baseline_root, replay_root, seconds)
            receipt['ranking'], ranks_by_target = compare_ranking(baseline, replay)
            source_count = receipt['ranking']['counts'].get('overlap_fixed', 0)
        else:
            source_count = compare_cases(baseline, replay)
            require(source_count == args.expected_sources, f'unexpected frozen candidate source count: {source_count}')
        old_pairs, pairs = delivered_pairs(baseline), delivered_pairs(replay)
        old_holo, now_holo = {h for h, _ in old_pairs}, {h for h, _ in pairs}
        targets = [(holo.upper(), candidate.upper()) for holo, candidate in args.target]
        if args.ranking_replay:
            targets = list(dict.fromkeys([*targets, *RANKING_REFERENCE_REQUESTS]))
        receipt.update(shards=shards, holo=len(replay), candidates=selected_count, fixed_sources=source_count,
                       prior_holo=len(old_holo), replay_holo=len(now_holo), prior_pairs=len(old_pairs),
                       replay_pairs=len(pairs), prior_delivery_copies=delivered_copy_count(baseline),
                       replay_delivery_copies=delivered_copy_count(replay),
                       new_holo=sorted(now_holo - old_holo), lost_holo=sorted(old_holo - now_holo),
                       new_pairs=sorted(pairs - old_pairs), lost_pairs=sorted(old_pairs - pairs),
                       targets=target_outcomes(replay, pairs, targets), literature=literature,
                       baseline_literature=baseline_literature, shard_seconds=seconds)
        if args.ranking_replay:
            receipt.update(retained_holo=sorted(now_holo & old_holo), net_holo=len(now_holo) - len(old_holo),
                           retained_pairs=sorted(pairs & old_pairs),
                           controls=check_ranking_controls(baseline, replay, pairs, ranks_by_target))
            for row in receipt['targets']:
                for match in row['matches']:
                    match['rrf_rank'] = ranks_by_target[match['target']].get(row['candidate'])
        else:
            receipt['relation_changes'] = relation_changes(baseline, replay)
        for sample, result, candidates, _ in replay.values():
            seen = set()
            for delivery in result['deliveries']:
                code = delivery['candidate']
                try:
                    copy = (code, *(delivery['apo_observation'][key] for key in OBS_KEYS))
                    require(copy not in seen, 'duplicate delivered chain copy')
                    seen.add(copy)
                    require(all(Path(delivery[key]['path']).resolve().is_relative_to(replay_root)
                                for key in ('holo_assembly', 'apo_assembly')),
                            'delivered assembly is outside the replay Run')
                    receipt['delivery_checks'].append(check_delivery(sample, result, candidates[code], delivery))
                except checked_errors as exc:
                    errors.append(f'{sample["id"]}/{code}: {exc}')
    except checked_errors as exc:
        errors.append(str(exc))
    receipt['all_pass'] = not errors
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8') as stream:
        json.dump(receipt, stream, ensure_ascii=False, indent=2)
        stream.write('\n')
    print(json.dumps({'output': str(output), 'all_pass': receipt['all_pass'], 'errors': errors,
                      'new_holo': receipt.get('new_holo'), 'lost_holo': receipt.get('lost_holo')}, ensure_ascii=False))
    if errors:
        raise SystemExit(1)


if __name__ == '__main__':
    main()

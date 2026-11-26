"""Bounded Holo-only candidate search and evidence pilot; no formal acceptance."""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import shutil
import time

import gemmi
import numpy as np
import pyarrow.parquet as pq
import requests
from scipy.spatial import cKDTree

import apo_c_pair
import apo_c_pilot
import apo_c_rank_audit
import apo_c_science
import apo_c_triplet_regression
import apo_c_union_triage
import holo_source_mapping
from apo_c_pilot import now, read_json, sha, write_json


def rankings(search, holo):
    pl = [{'pdb': x['source_entry_id'], 'author_chain': x['source_chain_auth_id'],
           'assembly': x['source_biounit_id'], 'source_rank': x['rank'],
           'chain_instance': x['source_chain_instance'],
           'target_match_granularity': x.get('target_match_granularity', 'LEGACY_UNVERIFIED')}
          for x in sorted(search['plinder'].get('links', []), key=lambda z: (z['rank'], z['linked_structure_id']))]
    ah = [{'pdb': x['pdb'], 'author_chain': x['auth_chain'], 'source_rank': x['source_order'],
           'source_target': match['entry']['entry_key'], 'assembly': None}
          for match in search['ahoj'].get('matches', []) for x in match['apo']]
    rc = [{'pdb': h['identifier'].split('_')[0], 'entity': h['identifier'],
           'query_id': s.get('query_id') or s.get('accession'), 'source_order': i + 1}
          for s in search['sifts_rcsb'].get('searches', []) for i, h in enumerate(s['hits'])]
    return {k: apo_c_rank_audit.unique_entry_ranking(v, holo)
            for k, v in [('PLINDER', pl), ('AHoJ', ah), ('RCSB', rc)]}


def select_budget(ranked, per_source, max_unique):
    # Equal rank-depth first, then source order; provenance from all ranked lists is retained.
    selected = []
    for i in range(per_source):
        for name in ('PLINDER', 'AHoJ', 'RCSB'):
            items = ranked[name]
            if i < len(items) and items[i]['pdb'] not in selected and len(selected) < max_unique:
                selected.append(items[i]['pdb'])
    return [{'pdb': pdb, 'sources': {name: item for name, items in ranked.items()
                                    for item in items[:10] if item['pdb'] == pdb}}
            for pdb in selected]


def prepare_holo(case):
    identity = case['identity']
    if len(identity['binding_instances']) != 1:
        return {'status': 'MULTICHAIN_SCOPE_REQUIRES_REVIEW',
                'reason': 'All required chain copies must be matched together.'}, None
    path = Path(case['source_uri'])
    if sha(path) != identity['source_sha256']:
        raise ValueError('Holo source hash changed')
    block = gemmi.cif.read(str(path)).sole_block()
    revision, first_date = apo_c_pair.revision(block, case['pdb'])
    if revision != identity['revision']:
        raise ValueError('Holo source revision changed')
    structure = gemmi.make_structure_from_block(block)
    if not structure or structure[0].num not in identity['model_scope']:
        raise ValueError('first deposited model is not within the requested model scope')
    instances = holo_source_mapping.assembly_instances(
        structure, str(identity['assembly']), structure[0], include_water=True)
    key = identity['binding_instances'][0][0]
    operator, label, rows = instances[key]
    residues = apo_c_triplet_regression.residue_index(rows)
    members, points = apo_c_triplet_regression.target_members(block, instances, identity)
    tree = cKDTree(points)
    site = {}
    for number, (author, residue) in residues.items():
        atoms = apo_c_triplet_regression.positive_heavy(residue, operator)
        if atoms and float(np.min(tree.query(np.array([a['xyz'] for a in atoms]))[0])) <= 4:
            site[number] = {'author_seq': str(residue.seqid), 'author_chain': author,
                            'component': residue.name, 'atom_names': {a['name'] for a in atoms}}
    if not site:
        raise ValueError('no observed Holo site residues')
    entity_ids = {str(r.entity_id) for _, r in rows if r.entity_type == gemmi.EntityType.Polymer}
    if len(entity_ids) != 1:
        raise ValueError('Holo binding chain has no unique polymer entity')
    acc = apo_c_union_triage.accessions(block, next(iter(entity_ids)))
    edges = apo_c_triplet_regression.connections_for_target(block, identity)
    external = [x for x in edges if x['kind'] == 'covale' and sum(x['target_endpoint_membership']) == 1]
    public = {'status': 'HOLO_OBSERVED_SITE_READY', 'system_id': case['id'], 'holo_status': case['holo_status'],
              'qc': case.get('qc') or {}, 'role': case.get('role'), 'role_evidence': case.get('role_evidence'),
              'identity': identity, 'source': str(path), 'source_sha256': sha(path),
              'source_revision': revision, 'first_public_model_date': first_date, 'model': structure[0].num,
              'binding_instance': key, 'binding_label': label, 'deposited_accessions': acc,
              'accession_basis': 'fixed mmCIF entity-level reference, not residue-level SIFTS',
              'target_members': members, 'target_connections': edges, 'external_covalent_connections': external,
              'site_definition': 'positive-occupancy target and protein heavy atoms within 4 A',
              'site_residues': [{'label_seq': n, **{k: sorted(v) if isinstance(v, set) else v
                                for k, v in item.items()}} for n, item in sorted(site.items())]}
    private = (block, structure, apo_c_pair.get_sequence(structure, rows), site, acc)
    return public, private


def observation_status(entity, observation):
    if not entity.get('common_uniprot'):
        return 'IDENTITY_REQUIRES_REVIEW'
    n = entity.get('holo_site_residues', 0)
    if not n or observation.get('site_observed_residues') != n:
        return 'SITE_NOT_FULLY_OBSERVED'
    if observation.get('site_complete_relative_holo_residues') != n:
        return 'SITE_ATOMS_OR_IDENTITY_DIFFER'
    contacts = observation.get('actual_component_contacts', {})
    if contacts.get('status') != 'MEASURED':
        return 'CONTACT_MEASUREMENT_UNAVAILABLE'
    if any(x['atoms_within_4A_of_mapped_site'] > 0 for x in contacts.get('nonwater', [])):
        return 'NONWATER_CONTACT_REQUIRES_REVIEW'
    return 'OBSERVED_SCOPE_NO_NONWATER_CONTACT_UNBOUND_UNPROVEN'


def observation_flags(entity, observation):
    contacts = observation.get('actual_component_contacts', {})
    return {'identity_has_common_accession': bool(entity.get('common_uniprot')),
            'holo_site_residues': entity.get('holo_site_residues', 0),
            'candidate_observed_residues': observation.get('site_observed_residues'),
            'candidate_complete_relative_holo_residues': observation.get('site_complete_relative_holo_residues'),
            'nonwater_within_4A': [x for x in contacts.get('nonwater', [])
                                   if x['atoms_within_4A_of_mapped_site'] > 0],
            'contact_measurement_status': contacts.get('status'),
            'unbound_experiment_evidence': 'NOT_ESTABLISHED_BY_GEOMETRY'}


class DownloadBudgetExceeded(RuntimeError):
    def __init__(self, reason, evidence):
        super().__init__(reason)
        self.evidence = evidence


def download_bounded(url, max_bytes, deadline_seconds):
    """Limit declared and observed size and elapsed streaming time for a raw original."""
    started = time.monotonic()
    evidence = {'url': url, 'started_utc': now(), 'bytes_received': 0}
    try:
        response = requests.head(url, allow_redirects=True, timeout=(5, 8))
        evidence['head_status'] = response.status_code
        declared = response.headers.get('Content-Length')
        evidence['head_declared_bytes'] = int(declared) if declared and declared.isdigit() else None
        response.close()
        if evidence['head_declared_bytes'] is not None and evidence['head_declared_bytes'] > max_bytes:
            raise DownloadBudgetExceeded('declared compressed size exceeds budget', evidence)
    except requests.RequestException as exc:
        evidence['head_error'] = f'{type(exc).__name__}: {exc}'
    with requests.get(url, stream=True, timeout=(5, 8),
                      headers={'User-Agent': 'APHCP-apo-research-pilot/0.1'}) as response:
        response.raise_for_status()
        evidence['status'] = response.status_code
        declared = response.headers.get('Content-Length')
        evidence['declared_bytes'] = int(declared) if declared and declared.isdigit() else None
        if evidence['declared_bytes'] is not None and evidence['declared_bytes'] > max_bytes:
            raise DownloadBudgetExceeded('declared compressed size exceeds budget', evidence)
        chunks = []
        for chunk in response.iter_content(chunk_size=65536):
            evidence['bytes_received'] += len(chunk)
            if evidence['bytes_received'] > max_bytes:
                raise DownloadBudgetExceeded('streamed compressed size exceeds budget', evidence)
            if time.monotonic() - started > deadline_seconds:
                raise DownloadBudgetExceeded('download elapsed-time budget exceeded', evidence)
            chunks.append(chunk)
    evidence['elapsed_s'] = round(time.monotonic() - started, 3)
    return b''.join(chunks), evidence


def fetch_candidate(run, pdb, origins, config, revision=None):
    key = pdb if revision is None else f"{pdb}.v{revision.replace('.', '-')}"
    base = run / 'fetch' / f'{key}.json'
    previous = apo_c_pilot.attempts(base)
    if previous:
        result = read_json(previous[-1])
        if result['status'] != 'FETCHED':
            if not config.get('retry_failed') or result.get('permanent'):
                return None
        elif revision is not None and result['revision'] != revision:
            raise ValueError('previously fetched revision differs from pinned revision')
        else:
            if sha(Path(result['path'])) != result['sha256']:
                raise ValueError('previously fetched original has changed')
            return result
    budget = config.setdefault('download_session', set())
    if key in budget or len(budget) >= config['max_candidate_download_attempts']:
        return None
    budget.add(key)
    receipt = apo_c_pilot.next_attempt(base)
    try:
        if revision:
            major, minor = map(int, revision.split('.'))
            basis = 'PLINDER 2026-09 entry_sources' if 'PLINDER' in origins else 'pinned comparison revision'
        else:
            metadata, info = apo_c_pilot.cache_get(run, f'rcsb_entry/{pdb.lower()}',
                                                   f'https://data.rcsb.org/rest/v1/core/entry/{pdb}')
            accession = json.loads(metadata)['rcsb_accession_info']
            major, minor = int(accession['major_revision']), int(accession['minor_revision'])
            basis = 'RCSB current fixed revision; AHoJ revision not established'
        rev, ident = f'{major}.{minor}', pdb.lower()
        url = (f'https://files-versioned.rcsb.org/pdb_versioned/data/entries/{ident[1:3]}/'
               f'pdb_0000{ident}/pdb_0000{ident}_xyz_v{major}-{minor}.cif.gz')
        body, request = download_bounded(url, config['max_compressed_bytes'], config['download_deadline_seconds'])
        block = gemmi.cif.read_string(gzip.decompress(body).decode('utf8')).sole_block()
        observed, first_date = apo_c_pair.revision(block, pdb)
        if observed != rev:
            raise ValueError('fixed candidate revision differs from selected revision')
        path = run / 'raw' / f'{ident}-v{rev}.cif.gz'
        path.parent.mkdir(parents=True, exist_ok=True)
        digest = hashlib.sha256(body).hexdigest()
        if path.exists() and sha(path) != digest:
            raise ValueError('cached candidate original differs from download')
        if not path.exists():
            path.write_bytes(body)
        result = {'status': 'FETCHED', 'path': str(path), 'sha256': sha(path), 'bytes': path.stat().st_size,
                  'revision': rev, 'first_public_model_date': first_date, 'request': request, 'version_basis': basis}
        if not revision:
            result['metadata_request'] = info
        write_json(receipt, result)
        return result
    except DownloadBudgetExceeded as exc:
        write_json(receipt, {'status': 'FETCH_BUDGET_EXCEEDED', 'pdb': pdb, 'reason': str(exc), 'request': exc.evidence,
                             'permanent': 'size' in str(exc)})
        return None
    except Exception as exc:
        write_json(receipt, {'status': 'FETCH_FAILED', 'pdb': pdb, 'error': f'{type(exc).__name__}: {exc}'})
        return None


def required_revision(pdb, origins, config):
    pinned = config.get('pinned_candidate_revisions', {}).get(pdb)
    versions = pq.read_table(Path(config['index']) / 'entry_sources.parquet',
                             filters=[('entry_pdb_id', '=', pdb.lower())]).to_pylist() if 'PLINDER' in origins else []
    values = {f"{v['source_mmcif_major_revision']}.{v['source_mmcif_minor_revision']}" for v in versions}
    if len(values) > 1 or (pinned and values and values != {pinned}):
        raise ValueError(f'conflicting required revisions for {pdb}: {values} vs {pinned}')
    return pinned or next(iter(values), None)


def saved_original(run, pdb, wanted):
    manifest = run / 'rank_comparison' / 'supplemental_raw' / f'{pdb}.json'
    receipt = run / 'fetch' / f'{pdb}.json'
    version_receipts = sorted(receipt.parent.glob(f'{pdb}.v*.json'))
    for path in [manifest, *reversed(apo_c_pilot.attempts(receipt)), *reversed(version_receipts)]:
        if not path.is_file():
            continue
        row = read_json(path)
        if (row.get('status') in (None, 'FETCHED') and row.get('path')
                and (wanted is None or row['revision'] == wanted)):
            source = Path(row['path'])
            if not source.is_file() or sha(source) != row['sha256']:
                raise ValueError(f'fixed original unavailable or changed: {path}')
            return (str(source), row['sha256'], source.stat().st_size, row['revision']), str(path)
    return None, None


def inspect_candidate(run, case, prepared, pdb, origins, config):
    block, structure, sequence, site, accessions = prepared
    wanted = required_revision(pdb, origins, config)
    with apo_c_pilot.open_catalog(Path(read_json(run / 'run.json')['catalog'])) as db:
        source = apo_c_union_triage.existing_original(db, pdb, wanted)
        source_origin, source_receipt = 'sealed_catalog_asset_file', None
        if source is None:
            for previous in [run, *(Path(x) for x in config['reuse_runs'])]:
                source, source_receipt = saved_original(previous, pdb, wanted)
                if source is not None:
                    source_origin = 'wwpdb_fixed_version_download' if previous == run else 'verified_prior_run'
                    break
        if source is None:
            fetched = fetch_candidate(run, pdb, origins, config, wanted)
            key = pdb if wanted is None else f"{pdb}.v{wanted.replace('.', '-')}"
            history = apo_c_pilot.attempts(run / 'fetch' / f'{key}.json')
            if fetched is None:
                return {'status': 'SOURCE_UNAVAILABLE_OR_BUDGET_EXHAUSTED', 'pdb': pdb,
                        'fetch_receipt': str(history[-1]) if history else None,
                        'required_revision': wanted}
            source = (fetched['path'], fetched['sha256'], fetched['bytes'], fetched['revision'])
            source_origin, source_receipt = 'wwpdb_fixed_version_download', str(history[-1])
        path, _, size, _ = source
        if size > config['max_compressed_bytes']:
            return {'status': 'RESOURCE_LIMIT_REQUIRES_REVIEW', 'pdb': pdb,
                    'reason': 'compressed source size exceeds frozen pilot budget', 'bytes': size}
        candidate = gemmi.read_structure(path)
        atom_count = sum(len(res) for chain in candidate[0] for res in chain)
        expanded = max((sum(len(g.operators) * len(g.subchains) for g in a.generators)
                        for a in candidate.assemblies), default=0)
        if (atom_count > config['max_asymmetric_unit_atoms']
                or len(candidate.assemblies) > config['max_assemblies']
                or expanded > config['max_assembly_instances']):
            return {'status': 'RESOURCE_LIMIT_REQUIRES_REVIEW', 'pdb': pdb, 'atom_count': atom_count,
                    'assembly_count': len(candidate.assemblies), 'max_assembly_instances': expanded}
        rows = [row for v in origins.values() for row in v['source_rows']]
        # RCSB can name an entity while another source names a different author chain.
        # Include both; never silently discard the chain-only source's entity.
        authors = {r['author_chain'] for r in rows if r.get('author_chain')}
        for chain in candidate[0]:
            if chain.name in authors:
                for residue in chain:
                    if residue.entity_type == gemmi.EntityType.Polymer and residue.entity_id:
                        rows.append({'entity': pdb + '_' + str(residue.entity_id)})
        result = apo_c_union_triage.candidate_evidence(
            run, case['pdb'], pdb, rows, {'sifts_common_uniprot': accessions},
            block, structure, sequence, site, db, fixed_source=source, fixed_origin=source_origin)
        if 'candidate_source' in result:
            result['candidate_source'].update(required_revision=wanted, receipt=source_receipt)
    if any(not e.get('common_uniprot') for e in result.get('entities', [])):
        try:
            mappings, info = apo_c_pilot.uniprot_mappings(run, pdb)
            result['candidate_sifts_request'] = info
            for entity in result.get('entities', []):
                if entity.get('common_uniprot'):
                    entity['identity_basis'] = 'fixed mmCIF entity references'
                    continue
                labels = {o['label_asym'] for o in entity.get('assembly_observations', []) if o.get('label_asym')}
                matches = {acc: [m for m in detail.get('mappings', [])
                                 if m.get('struct_asym_id') in labels and
                                 str(m.get('entity_id')) == str(entity['entity_id'])]
                           for acc, detail in mappings.items() if acc in accessions}
                matches = {acc: rows for acc, rows in matches.items() if rows}
                entity['sifts_chain_mappings'] = matches
                if matches:
                    entity['common_uniprot'] = sorted(matches)
                    entity['identity_basis'] = 'Holo deposited accession plus candidate live SIFTS chain/segment mapping'
                    if entity.get('site_residues_covered_by_construct_sequence'):
                        entity['status'] = 'SITE_SEQUENCE_COVERAGE_REQUIRES_OBSERVATION_CHECK'
        except Exception as exc:
            result['candidate_sifts_error'] = f'{type(exc).__name__}: {exc}'
    for entity in result.get('entities', []):
        for observation in entity.get('assembly_observations', []):
            if observation.get('status') == 'ACTUAL_CONTACTS_MEASURED':
                observation['research_screen'] = observation_status(entity, observation)
                observation['research_flags'] = observation_flags(entity, observation)
    result['identity_evidence_note'] = 'Common accession comes from deposited references or an explicitly recorded SIFTS chain mapping; residue-level SIFTS is not established.'
    result['source_revision_note'] = 'Candidate fixed revision is explicit; match to each source release remains to be checked.'
    return result


def reuse_searches(run, previous):
    """Reuse a search only for the same frozen Holo identity and search implementation."""
    source = Path(previous)
    prior_cases = {x['id']: x for x in read_json(source / 'samples.json')}
    entries = []
    for case in read_json(run / 'samples.json'):
        old = prior_cases.get(case['id'])
        src, dest = apo_c_pilot.case_file(source, case), apo_c_pilot.case_file(run, case)
        if not old or not src.is_file():
            continue
        if old['identity'] != case['identity'] or sha(Path(old['source_uri'])) != case['identity']['source_sha256']:
            raise ValueError(f'prior Holo identity differs: {case["id"]}')
        data = read_json(src)
        if data['sample'] != case['id'] or data['code_sha256'] != sha(Path(apo_c_pilot.__file__)):
            raise ValueError(f'prior search implementation differs: {case["id"]}')
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            if sha(dest) != sha(src):
                raise ValueError(f'existing search differs: {case["id"]}')
        else:
            shutil.copy2(src, dest)
        entries.append({'system_id': case['id'], 'source': str(src), 'sha256': sha(src)})
    record = run / 'search_reuse.json'
    if record.exists():
        if read_json(record)['entries'] != entries:
            raise ValueError('prior search-reuse record differs')
    else:
        write_json(record, {'source_run': str(source), 'entries': entries})


def read_holo_role(catalog, case):
    with apo_c_pilot.open_catalog(Path(catalog)) as db:
        row = db.execute("SELECT error_message FROM run_item WHERE object_id=? AND stage='STRUCTURE_OBSERVATION' "
                         "AND output_uri IS NOT NULL ORDER BY updated_at DESC,run_id DESC LIMIT 1", (case['id'],)).fetchone()
    if row is None:
        raise ValueError('sealed Holo role row is missing')
    record = json.loads(row[0])
    if record['identity'] != case['identity']:
        raise ValueError('sealed Holo role row has a different identity')
    return {'role': record.get('role'), 'role_evidence': record.get('role_evidence')}


def has_technical_gap(evidence):
    return (evidence['scientific_decision']['status'] == 'TECHNICAL_GAP'
            or any(o.get('status') == 'TECHNICAL_FAILURE' for e in evidence.get('entities', [])
                   for o in e.get('assembly_observations', []))
            or bool(evidence.get('candidate_sifts_error')))


def source_hints(origins):
    return {(name, tuple(row.get(k) for k in ('entity', 'author_chain', 'assembly', 'chain_instance')))
            for name, route in origins.items() for row in route['source_rows']}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', required=True)
    p.add_argument('--catalog', required=True)
    p.add_argument('--index', required=True)
    p.add_argument('--samples', required=True)
    p.add_argument('--reuse-run', action='append', default=[],
                   help='prior C research Run with hash-checked fixed candidate originals')
    p.add_argument('--reuse-search-run', help='prior Run with the identical Holo search implementation and inputs')
    p.add_argument('--per-source-rank', type=int, default=3)
    p.add_argument('--max-unique-per-holo', type=int, default=6)
    p.add_argument('--retry-failed', action='store_true',
                   help='Retry missing search pages, fixed inputs, and technical observations in an existing run')
    args = p.parse_args()
    if not 1 <= args.per_source_rank <= 10 or not 1 <= args.max_unique_per_holo <= 30:
        p.error('invalid candidate inspection budget')
    run, sample_path = Path(args.run), Path(args.samples)
    deps = [Path(__file__), Path(apo_c_pilot.__file__), Path(apo_c_pair.__file__),
            Path(apo_c_rank_audit.__file__), Path(apo_c_science.__file__), Path(apo_c_triplet_regression.__file__),
            Path(apo_c_union_triage.__file__), Path(holo_source_mapping.__file__)]
    config = {'sample_spec_sha256': sha(sample_path), 'catalog': str(Path(args.catalog).resolve()),
              'index': str(Path(args.index).resolve()), 'per_source_rank': args.per_source_rank,
              'max_unique_per_holo': args.max_unique_per_holo, 'max_compressed_bytes': 20000000,
              'max_asymmetric_unit_atoms': 200000, 'max_assemblies': 6, 'max_assembly_instances': 128,
              'max_candidate_download_attempts': 10, 'max_derivative_assemblies_per_candidate': 2,
              'download_deadline_seconds': 30,
              'reuse_runs': [str(Path(x).resolve()) for x in args.reuse_run],
              'reuse_search_run': str(Path(args.reuse_search_run).resolve()) if args.reuse_search_run else None,
              'code_sha256': {x.name: sha(x) for x in deps},
              'scientific_rule': apo_c_science.RULE, 'formal_pair_acceptance': False}
    started, clock = now(), time.monotonic()
    if not run.exists():
        apo_c_pilot.initial_samples(SimpleNamespace(samples=str(sample_path), run=str(run),
                                                    catalog=args.catalog, index=args.index))
        write_json(run / 'frozen_config.json', config)
        for path in deps:
            snapshot = run / 'code' / f'{path.stem}-{sha(path)[:12]}.py'
            if not snapshot.exists():
                shutil.copy2(path, snapshot)
    elif not (run / 'frozen_config.json').exists() or read_json(run / 'frozen_config.json') != config:
        raise ValueError('existing Run does not match the frozen configuration and code')
    if args.reuse_search_run:
        reuse_searches(run, args.reuse_search_run)
    config['retry_failed'] = args.retry_failed
    apo_c_pilot.search_cases(SimpleNamespace(run=str(run), case=None, split=None, retry_failed=args.retry_failed))
    results = []
    for case in read_json(run / 'samples.json'):
        key = case['id'].replace(':', '_')
        outcome = run / 'automated' / key / 'result.json'
        search_path = apo_c_pilot.attempts(apo_c_pilot.case_file(run, case))[-1]
        old_results = apo_c_pilot.attempts(outcome)
        if old_results:
            old = read_json(old_results[-1])
            if (old['search_sha256'] == sha(search_path) and
                    (not args.retry_failed or (old['holo_evidence']['status'] != 'TECHNICAL_FAILURE'
                      and not any(has_technical_gap(read_json(Path(row['evidence_path']))) for row in old['candidates'])))):
                results.append(old)
                print('RESULT_REUSE', case['id'], flush=True)
                continue
        search = read_json(search_path)
        rank = rankings(search, case['pdb'])
        chosen = select_budget(rank, args.per_source_rank, args.max_unique_per_holo)
        result = {'system_id': case['id'], 'holo': case['pdb'], 'search_evidence': str(search_path),
                  'search_sha256': sha(search_path), 'source_status': {x: search[x]['status']
                         for x in ('plinder', 'sifts_rcsb', 'ahoj')},
                  'rankings_top10': {k: v[:10] for k, v in rank.items()}, 'selected_candidates': chosen,
                  'candidates': [], 'deliveries': [], 'formal_pair_acceptance': False}
        output_file = apo_c_pilot.next_attempt(outcome)
        try:
            case.update(read_holo_role(args.catalog, case))
            holo, prepared = prepare_holo(case)
            result['holo_evidence'] = holo
            if prepared is not None:
                derived = outcome.parent / 'structures'
                derived.mkdir(parents=True, exist_ok=True)
                holo_output = derived / 'holo_system.cif'
                source_holo = Path(case['system_uri'])
                if holo_output.exists() and sha(holo_output) != sha(source_holo):
                    raise ValueError('existing Holo system differs from sealed source')
                if not holo_output.exists():
                    shutil.copy2(source_holo, holo_output)
                if not gemmi.read_structure(str(holo_output))[0].count_atom_sites():
                    raise ValueError('Holo system readback has no atoms')
                result['holo_system'] = {'path': str(holo_output), 'sha256': sha(holo_output)}
                for item in chosen:
                    candidate_path = outcome.parent / 'candidates' / f"{item['pdb']}.json"
                    history = apo_c_pilot.attempts(candidate_path)
                    if history:
                        evidence = read_json(history[-1])
                        if (source_hints(evidence.get('search_origins', {})) == source_hints(item['sources']) and
                                (not args.retry_failed or not has_technical_gap(evidence))):
                            result['candidates'].append({'pdb': item['pdb'], 'status': evidence['status'],
                                                         'evidence_path': str(history[-1])})
                            print('CANDIDATE_REUSE', case['id'], item['pdb'], flush=True)
                            continue
                    candidate_path = apo_c_pilot.next_attempt(candidate_path)
                    print('CANDIDATE_START', case['id'], item['pdb'], flush=True)
                    try:
                        evidence = inspect_candidate(run, case, prepared, item['pdb'], item['sources'], config)
                        evidence['scientific_decision'] = apo_c_science.decide_candidate(holo, prepared, evidence)
                        selected_assemblies = sorted({o['assembly'] for e in evidence.get('entities', [])
                                for o in e.get('assembly_observations', [])
                                if o.get('status') == 'ACTUAL_CONTACTS_MEASURED'})[:config['max_derivative_assemblies_per_candidate']]
                        evidence['materialized_assemblies'] = []
                        if selected_assemblies:
                            candidate = gemmi.read_structure(evidence['candidate_source']['path'])
                            for assembly in selected_assemblies:
                                out = derived / f"{item['pdb']}_assembly_{assembly}.cif"
                                attempt = 1
                                while out.exists():
                                    attempt += 1
                                    out = derived / f"{item['pdb']}_assembly_{assembly}_attempt{attempt:03d}.cif"
                                qc = apo_c_pair.write_assembly(candidate, assembly, out)
                                evidence['materialized_assemblies'].append({'assembly': assembly, 'path': str(out), 'qc': qc})
                    except Exception as exc:
                        evidence = {'status': 'TECHNICAL_FAILURE', 'pdb': item['pdb'],
                                    'error': f'{type(exc).__name__}: {exc}',
                                    'scientific_decision': {'status': 'TECHNICAL_GAP', 'rule': apo_c_science.RULE,
                                                            'reasons': ['CANDIDATE_INSPECTION_FAILED']}}
                    evidence['search_origins'] = item['sources']
                    write_json(candidate_path, evidence)
                    result['candidates'].append({'pdb': item['pdb'], 'status': evidence['status'],
                                                 'evidence_path': str(candidate_path)})
                    print('CANDIDATE_DONE', case['id'], item['pdb'], evidence['status'], flush=True)
                for row in result['candidates']:
                    e = read_json(Path(row['evidence_path']))
                    if e['scientific_decision']['status'] != 'SUPPORTED_RESEARCH_PAIR':
                        continue
                    obs, = [o for ent in e['entities'] for o in ent['assembly_observations']
                            if o['decision']['status'] == 'SUPPORTED_RESEARCH_PAIR']
                    assembly, = [i for i in e['materialized_assemblies'] if i['assembly'] == obs['assembly']]
                    full_holo = derived / 'holo_full_assembly.cif'
                    if full_holo.exists() and not full_holo.with_suffix('.json').exists():
                        n = 2
                        while (derived / f'holo_full_assembly_attempt{n:03d}.cif').exists() and not (
                               derived / f'holo_full_assembly_attempt{n:03d}.json').exists():
                            n += 1
                        full_holo = derived / f'holo_full_assembly_attempt{n:03d}.cif'
                    if not full_holo.exists():
                        holo_qc = apo_c_pair.write_assembly(
                            prepared[1], str(case['identity']['assembly']), full_holo,
                            required_connections={edge['id'] for edge in holo['target_connections']})
                        write_json(full_holo.with_suffix('.json'), holo_qc)
                    holo_qc = read_json(full_holo.with_suffix('.json'))
                    if sha(full_holo) != holo_qc['sha256'] or sha(Path(assembly['path'])) != assembly['qc']['sha256']:
                        raise ValueError('delivered assembly differs from readback')
                    result['deliveries'].append({
                        'system_id': case['id'], 'candidate': row['pdb'], 'role_evidence': e['scientific_decision']['holo_target_role'],
                        'holo_identity': holo['identity'], 'holo_original': holo['source'],
                        'holo_assembly': {'path': str(full_holo), 'qc': holo_qc},
                        'apo_original': e['candidate_source'], 'apo_assembly': assembly,
                        'apo_observation': {k: obs[k] for k in ('model', 'instance', 'label_asym', 'author_chain', 'assembly')},
                        'holo_target_members': holo['target_members'], 'mapped_site': obs['mapped_site'],
                        'assessment': obs['decision'], 'candidate_evidence': row['evidence_path']})
            else:
                result['candidate_review'] = 'MULTICHAIN_SCOPE_REQUIRES_REVIEW'
        except Exception as exc:
            result['holo_evidence'] = {'status': 'TECHNICAL_FAILURE', 'error': f'{type(exc).__name__}: {exc}'}
            result['deliveries'] = []
        write_json(output_file, result)
        results.append(result)
        print('RESULT', case['id'], result['holo_evidence']['status'],
              'selected', len(chosen), 'examined', len(result['candidates']), flush=True)
    states = Counter(candidate['status'] for result in results for candidate in result['candidates'])
    research = {}
    for result in results:
        research[result['system_id']] = [read_json(Path(candidate['evidence_path']))['scientific_decision']['status']
                                         for candidate in result['candidates']]
    observations = Counter(o['research_screen'] for result in results for candidate in result['candidates']
                           for e in read_json(Path(candidate['evidence_path'])).get('entities', [])
                           for o in e.get('assembly_observations', [])
                           if 'research_screen' in o)
    summary = {'started_utc': started, 'finished_utc': now(), 'elapsed_s': round(time.monotonic() - clock, 3),
               'holo_dispatched': len(results), 'candidate_entry_outcomes': dict(states),
               'chain_observations': dict(observations), 'formal_pair_acceptances': 0,
               'scientific_rule': apo_c_science.RULE,
               'research_pair_outcomes': dict(Counter(s for values in research.values() for s in values)),
               'holo_with_research_pair': sum('SUPPORTED_RESEARCH_PAIR' in values for values in research.values()),
               'research_pairs_delivered': sum(len(r.get('deliveries', [])) for r in results),
               'holo_with_delivered_research_pair': sum(bool(r.get('deliveries')) for r in results),
               'limits': ['Search routes can fail independently.', 'Multi-chain targets remain review.',
                          'An observed no-contact scope is not proof of an unbound experiment.',
                          'No population accuracy, blind assessment or formal acceptance is claimed.']}
    n = len(list((run / 'receipts').glob('automation_run*.json'))) + 1 if (run / 'receipts').exists() else 1
    path = run / 'receipts' / f'automation_run{n:03d}.json'
    write_json(path, summary)
    print('SUMMARY', path, summary, flush=True)


if __name__ == '__main__':
    main()

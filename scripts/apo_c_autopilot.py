"""Bounded Holo-only candidate search and evidence pilot; no formal acceptance."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
import os
import re
from tempfile import NamedTemporaryFile
from pathlib import Path
from types import SimpleNamespace
import shutil
import time

import gemmi
import numpy as np
import pyarrow.parquet as pq
import requests
from scipy.spatial import cKDTree

import apo_c_author_documents
import apo_c_literature
import apo_c_pair
import apo_c_pilot
import apo_c_priority
import apo_c_rank_audit
import apo_c_science
import apo_c_triplet_regression
import apo_c_union_triage
import holo_source_mapping
from apo_c_pilot import now, read_json, sha, write_json


RRF_OFFSET = 60
FIRST_RUN_RANKING = 'metadata-quota-rrf-v1'


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


def select_budget(ranked, per_source, max_unique, original_rank=None):
    # Equal rank-depth first, then source order; provenance from all ranked lists is retained.
    selected = []
    for i in range(per_source):
        for name in ('PLINDER', 'AHoJ', 'RCSB'):
            items = ranked[name]
            if i < len(items) and items[i]['pdb'] not in selected and len(selected) < max_unique:
                selected.append(items[i]['pdb'])
    return [{'pdb': pdb, 'sources': {name: item for name, items in (original_rank or ranked).items()
                                    for item in items if item['pdb'] == pdb}}
            for pdb in selected]


def select_fused_budget(ranked, quota):
    candidates = {}
    for name, rows in ranked.items():
        for row in rows:
            candidates.setdefault(row['pdb'], {})[name] = row

    def priority(item):
        pdb, sources = item
        ranks = [row['entry_rank'] for row in sources.values()]
        return -sum(Fraction(1, RRF_OFFSET + rank) for rank in ranks), min(ranks), pdb

    if len(candidates) < quota:
        raise ValueError('sealed search has fewer candidates than the frozen quota')
    return [{'pdb': pdb, 'sources': sources}
            for pdb, sources in sorted(candidates.items(), key=priority)[:quota]]


def select_first_run_budget(ranked, metadata_order, per_source, max_unique):
    quota = len(select_budget(metadata_order, per_source, max_unique))
    return select_fused_budget(ranked, quota)


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


def shared_fixed_original(directory, pdb, revision, fetched=None, receipt=None):
    """Publish a receipt pointer only after validating the immutable archive and fixed identity."""
    if not directory or not revision:
        return None, None
    if not re.fullmatch(r'[0-9][A-Z0-9]{3}', pdb) or not re.fullmatch(r'\d+\.\d+', revision):
        raise ValueError('invalid archive identity')
    root = Path(directory)
    manifest = root / pdb[:2] / f'{pdb}-v{revision}.json'

    def verified(row):
        source = Path(row['path'])
        if (row.get('revision') != revision or not source.is_file() or source.stat().st_size != row['bytes']
                or sha(source) != row['sha256'] or apo_c_pair.revision(gemmi.cif.read(str(source)).sole_block(), pdb)[0] != revision):
            raise ValueError(f'invalid shared fixed original for {pdb} {revision}')
        return str(source), row['sha256'], row['bytes'], revision

    if fetched is not None:
        if fetched.get('status') != 'FETCHED' or not receipt or not Path(receipt).is_file():
            raise ValueError('shared source has no successful originating receipt')
        origin = read_json(Path(receipt))
        fields = ('path', 'revision', 'bytes', 'sha256')
        if origin.get('status') != 'FETCHED' or any(origin[k] != fetched[k] for k in fields):
            raise ValueError('original and provenance disagree')
        row = {k: fetched[k] for k in fields} | {'pdb': pdb, 'receipt': str(receipt), 'receipt_sha256': sha(Path(receipt))}
        if not manifest.exists():
            verified(row)
            manifest.parent.mkdir(parents=True, exist_ok=True)
            with NamedTemporaryFile(dir=manifest.parent, mode='w', encoding='utf8', delete=False) as tmp:
                json.dump(row, tmp, sort_keys=True)
                tmp.flush()
                os.fsync(tmp.fileno())
            try:
                try:
                    os.link(tmp.name, manifest)
                except FileExistsError:
                    pass
            finally:
                Path(tmp.name).unlink(missing_ok=True)
    if not manifest.is_file():
        return None, None
    stored = read_json(manifest)
    if (stored['pdb'] != pdb or (fetched and stored['sha256'] != fetched['sha256'])
            or sha(Path(stored['receipt'])) != stored['receipt_sha256']):
        raise ValueError('conflicting shared fixed source or changed provenance')
    return verified(stored), str(manifest)


def attempted_revision(run, pdb):
    for path in reversed(apo_c_pilot.attempts(run / 'fetch' / f'{pdb}.json')):
        row = read_json(path)
        url = row.get('request', {}).get('url', '')
        match = re.fullmatch(r'https://files-versioned\.rcsb\.org/pdb_versioned/data/entries/[a-z0-9]{2}/'
                             + f'pdb_0000{pdb.lower()}/pdb_0000{pdb.lower()}_xyz_v' + r'(\d+)-(\d+)\.cif\.gz', url)
        if match:
            return '.'.join(match.groups())
    return None


def candidate_original(run, pdb, origins, config):
    wanted = required_revision(pdb, origins, config)
    with apo_c_pilot.open_catalog(Path(read_json(run / 'run.json')['catalog'])) as db:
        source = apo_c_union_triage.existing_original(db, pdb, wanted)
    if source is not None:
        return wanted, source, 'sealed_catalog_asset_file', None
    for previous in [run, *(Path(x) for x in config['reuse_runs'])]:
        source, receipt = saved_original(previous, pdb, wanted)
        if source is not None:
            origin = 'wwpdb_fixed_version_download' if previous == run else 'verified_prior_run'
            return wanted, source, origin, receipt
    shared_revision = wanted or attempted_revision(run, pdb)
    source, receipt = shared_fixed_original(config.get('shared_originals'), pdb, shared_revision)
    if source is not None:
        return shared_revision, source, 'verified_sibling_original', receipt
    return wanted, None, None, None


def deposition_context(block, holo):
    table = block.get_mmcif_category('_pdbx_database_related.')
    rows = [(i, str(code).upper(), str(note or '')) for i, (database, code, note) in enumerate(
        zip(*(table.get(k, []) for k in ('db_name', 'db_id', 'details'))), 1) if str(database).upper() == 'PDB']
    return {'entry': holo['identity']['pdb_id'], 'sha256': holo['source_sha256'],
            'revision': holo['source_revision'], 'rows': rows}


def study_fingerprint(pdb, current, case, run, config):
    source = current.get('candidate_source') if pdb == current['pdb'] else None
    if source:
        path, digest, size, revision = (source[key] for key in ('path', 'sha256', 'bytes', 'revision'))
    elif pdb == case['pdb']:
        path, digest, size, revision = (case['source_uri'], case['identity']['source_sha256'],
                                         case['source_bytes'], case['identity']['revision'])
    else:
        search_runs = [run, *(Path(x) for x in config['reuse_runs'])]
        if previous := config.get('replay_evidence_run'):
            search_runs.append(Path(previous))
        saved = next((row for previous in search_runs if (row := saved_original(previous, pdb, None)[0])), None)
        if saved is None:
            with apo_c_pilot.open_catalog(Path(config['catalog'])) as db:
                saved = apo_c_union_triage.existing_original(db, pdb)
        if saved is None:
            return None
        path, digest, size, revision = saved
    path = Path(path)
    if path.stat().st_size != size or size > config['max_compressed_bytes'] or sha(path) != digest:
        return None
    block = gemmi.cif.read(str(path)).sole_block()
    if apo_c_pair.revision(block, pdb)[0] != revision:
        return None
    tags = {'a': '_cell.length_a', 'b': '_cell.length_b', 'c': '_cell.length_c',
            'alpha': '_cell.angle_alpha', 'beta': '_cell.angle_beta', 'gamma': '_cell.angle_gamma',
            'diffraction_high': '_reflns.d_resolution_high', 'refinement_high': '_refine.ls_d_res_high',
            'wavelength': '_diffrn_source.pdbx_wavelength', 'r_work': '_refine.ls_R_factor_R_work',
            'space_group': '_symmetry.space_group_name_H-M'}
    fields = {key: gemmi.cif.as_string(values[0]) for key, tag in tags.items()
              if len(values := block.find_values(tag)) == 1 and values[0] not in {'?', '.', ''}}
    return {'entry': pdb, 'primary_citation': apo_c_science.primary_experiment(block), 'fields': fields,
            'source_receipt': {'path': str(path), 'sha256': digest, 'bytes': size, 'revision': revision,
                               'fields': {key: tags[key] for key in fields}}}


def inspect_candidate(run, case, prepared, pdb, origins, config, original=None):
    block, structure, sequence, site, accessions = prepared
    wanted, source, source_origin, source_receipt = (
        original if original is not None else candidate_original(run, pdb, origins, config))
    if source is None:
        attempted = len(config.get('download_session', ()))
        fetched = fetch_candidate(run, pdb, origins, config, wanted)
        if len(config.get('download_session', ())) > attempted:
            print('DOWNLOAD_ATTEMPT', case['id'], pdb, wanted, flush=True)
        key = pdb if wanted is None else f"{pdb}.v{wanted.replace('.', '-')}"
        history = apo_c_pilot.attempts(run / 'fetch' / f'{key}.json')
        if fetched is None:
            reason = 'ARCHIVE_FETCH_FAILED' if history else 'BATCH_DOWNLOAD_LIMIT'
            return {'status': 'SOURCE_UNAVAILABLE_OR_BUDGET_EXHAUSTED', 'pdb': pdb,
                    'fetch_receipt': str(history[-1]) if history else None,
                    'required_revision': wanted, 'source_gap_reason': reason}
        source = (fetched['path'], fetched['sha256'], fetched['bytes'], fetched['revision'])
        source_origin, source_receipt = 'wwpdb_fixed_version_download', str(history[-1])
        shared_fixed_original(config.get('shared_originals'), pdb, fetched['revision'], fetched, source_receipt)
    path, _, size, _ = source
    if size > config['max_compressed_bytes']:
        return {'status': 'RESOURCE_LIMIT_REQUIRES_REVIEW', 'pdb': pdb,
                'reason': 'compressed source size exceeds frozen pilot budget', 'bytes': size}
    candidate = gemmi.read_structure(path)
    atom_count = sum(len(res) for chain in candidate[0] for res in chain)
    if atom_count > config['max_asymmetric_unit_atoms']:
        return {'status': 'RESOURCE_LIMIT_REQUIRES_REVIEW', 'pdb': pdb, 'atom_count': atom_count,
                'reason': 'asymmetric unit atom count exceeds frozen pilot budget'}
    rows = [row for v in origins.values() for row in v['source_rows']]
    # Include both entity and author-chain identities when sources differ.
    authors = {r['author_chain'] for r in rows if r.get('author_chain')}
    entities = {str(res.entity_id) for chain in candidate[0] if chain.name in authors
                for res in chain if res.entity_type == gemmi.EntityType.Polymer and res.entity_id}
    rows.extend({'entity': pdb + '_' + entity} for entity in sorted(entities))
    with apo_c_pilot.open_catalog(Path(read_json(run / 'run.json')['catalog'])) as db:
        result = apo_c_union_triage.candidate_evidence(
            run, case['pdb'], pdb, rows, {'sifts_common_uniprot': accessions},
            block, structure, sequence, site, db, fixed_source=source, fixed_origin=source_origin,
            max_assembly_instances=config['max_assembly_instances'])
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
    science = evidence['scientific_decision']
    literature = science.get('candidate_experiment', {}).get('literature', {})
    return (evidence['scientific_decision']['status'] == 'TECHNICAL_GAP'
            or any(o.get('status') == 'TECHNICAL_FAILURE' for e in evidence.get('entities', [])
                   for o in e.get('assembly_observations', []))
            or bool(evidence.get('candidate_sifts_error'))
            or bool((literature.get('fetch_error') or literature.get('related_deposit_fetch_error') or literature.get('publisher_source_unavailable'))
                    and not science.get('delivery_eligible')))


def source_hints(origins):
    return {(name, tuple(row.get(k) for k in ('entity', 'author_chain', 'assembly', 'chain_instance')))
            for name, route in origins.items() for row in route['source_rows']}


def materialize_candidate(evidence, pdb, directory, limit):
    observations = [o for e in evidence.get('entities', []) for o in e.get('assembly_observations', [])
                    if o.get('status') == 'ACTUAL_CONTACTS_MEASURED']
    accepted = {o['assembly'] for o in observations if o['decision']['delivery_eligible']}
    measured = {o['assembly'] for o in observations}
    selected = sorted(accepted) + sorted(measured - accepted)[:max(0, limit - len(accepted))]
    assets = []
    if selected:
        candidate = gemmi.read_structure(evidence['candidate_source']['path'])
        for assembly in selected:
            out = directory / f'{pdb}_assembly_{assembly}.cif'
            attempt = 1
            while out.exists():
                attempt += 1
                out = directory / f'{pdb}_assembly_{assembly}_attempt{attempt:03d}.cif'
            qc = apo_c_pair.write_assembly(candidate, assembly, out)
            assets.append({'assembly': assembly, 'path': str(out), 'qc': qc})
    return assets


def deliver_observations(case, holo, prepared, candidates, directory):
    deliveries = []
    holo_asset = None
    for row in candidates:
        evidence = read_json(Path(row['evidence_path']))
        scientific = evidence['scientific_decision']
        if not scientific.get('delivery_eligible'):
            continue
        assemblies = {a['assembly']: a for a in evidence['materialized_assemblies']}
        accepted = [o for e in evidence['entities'] for o in e['assembly_observations']
                    if o['decision']['delivery_eligible']]
        if holo_asset is None:
            path = directory / 'holo_full_assembly.cif'
            attempt = 1
            while path.exists() and not path.with_suffix('.json').exists():
                attempt += 1
                path = directory / f'holo_full_assembly_attempt{attempt:03d}.cif'
            if not path.exists():
                qc = apo_c_pair.write_assembly(prepared[1], str(case['identity']['assembly']), path,
                        required_connections={edge['id'] for edge in holo['target_connections']})
                write_json(path.with_suffix('.json'), qc)
            qc = read_json(path.with_suffix('.json'))
            if sha(path) != qc['sha256']:
                raise ValueError('delivered Holo assembly differs from readback')
            holo_asset = {'path': str(path), 'qc': qc}
        for assembly_id in {o['assembly'] for o in accepted}:
            assembly = assemblies[assembly_id]
            if sha(Path(assembly['path'])) != assembly['qc']['sha256']:
                raise ValueError('delivered Apo assembly differs from readback')
        for obs in accepted:
            deliveries.append({
                'system_id': case['id'], 'candidate': row['pdb'], 'role_evidence': scientific['holo_target_role'],
                'holo_identity': holo['identity'], 'holo_original': holo['source'], 'holo_assembly': holo_asset,
                'apo_original': evidence['candidate_source'], 'apo_assembly': assemblies[obs['assembly']],
                'apo_observation': {k: obs[k] for k in ('model', 'instance', 'label_asym', 'author_chain', 'assembly')},
                'holo_target_members': holo['target_members'], 'mapped_site': obs['mapped_site'],
                'assessment': obs['decision'], 'candidate_experiment': scientific['candidate_experiment'],
                'candidate_evidence': row['evidence_path']})
    return deliveries


def frozen_case(run, previous, case, limit):
    key = case['id'].replace(':', '_')
    prior_result = previous / 'automated' / key / 'result.json'
    baseline = read_json(prior_result)
    selected = baseline['selected_candidates']
    inspected = baseline['candidates']
    ids = [item['pdb'] for item in selected]
    if (baseline['system_id'] != case['id'] or baseline['holo'] != case['pdb']
            or baseline['holo_evidence']['identity'] != case['identity']
            or ids != [row['pdb'] for row in inspected] or len(set(ids)) != len(ids) or len(ids) > limit):
        raise ValueError(f'inconsistent frozen case: {prior_result}')
    holo_original = Path(case['source_uri'])
    if holo_original.stat().st_size != case['source_bytes'] or sha(holo_original) != case['identity']['source_sha256']:
        raise ValueError(f'frozen Holo original differs: {holo_original}')
    search = apo_c_pilot.case_file(run, case)
    if sha(search) != baseline['search_sha256']:
        raise ValueError(f'frozen search differs: {search}')
    holo_system = Path(case['system_uri'])
    if sha(holo_system) != baseline['holo_system']['sha256']:
        raise ValueError(f'frozen Holo system differs: {holo_system}')
    for item, row in zip(selected, inspected):
        evidence_path = Path(row['evidence_path']).resolve()
        if not evidence_path.is_relative_to(prior_result.parent.resolve()):
            raise ValueError(f'candidate evidence is outside its frozen case: {evidence_path}')
        evidence = read_json(evidence_path)
        if evidence['pdb'] != item['pdb'] or evidence['search_origins'] != item['sources']:
            raise ValueError(f'frozen candidate differs: {evidence_path}')
        if source := evidence.get('candidate_source'):
            path = Path(source['path'])
            if path.stat().st_size != source['bytes'] or sha(path) != source['sha256']:
                raise ValueError(f'frozen original differs: {path}')
    return baseline


def run_case(run, case, config, literature, baseline=None):
    """Yield before a missing original so other Holo get a turn to fetch."""
    key = case['id'].replace(':', '_')
    outcome = run / 'automated' / key / 'result.json'
    search_path = apo_c_pilot.attempts(apo_c_pilot.case_file(run, case))[-1]
    search_digest = sha(search_path)
    old_results = apo_c_pilot.attempts(outcome)
    if old_results:
        old = read_json(old_results[-1])
        if (old['search_sha256'] == search_digest and
                (not config['retry_failed'] or (old['holo_evidence']['status'] != 'TECHNICAL_FAILURE'
                  and not old.get('processing_failure') and not old.get('delivery_failures')
                  and not any(has_technical_gap(read_json(Path(row['evidence_path']))) for row in old['candidates'])))):
            print('RESULT_REUSE', case['id'], flush=True)
            return old
    search = read_json(search_path)
    rerank = 'rerank_from_run' in config
    if baseline is None or rerank:
        rank = rankings(search, case['pdb'])
        chosen = select_fused_budget(rank, len(baseline['selected_candidates'])) if rerank else []
        top10 = {k: v[:10] for k, v in rank.items()}
    else:
        chosen = baseline['selected_candidates']
        top10 = baseline['rankings_top10']
    result = {'system_id': case['id'], 'holo': case['pdb'], 'search_evidence': str(search_path),
              'search_sha256': search_digest, 'source_status': {x: search[x]['status']
                     for x in ('plinder', 'sifts_rcsb', 'ahoj')},
              'rankings_top10': top10, 'selected_candidates': chosen,
              'candidates': [], 'deliveries': [], 'formal_pair_acceptance': False}
    if baseline is None:
        result['ranking_counts'] = {name: len(rows) for name, rows in rank.items()}
        partial = 'technical_failure' in result['source_status'].values()
        if any(rank.values()):
            result['candidate_availability'] = 'CANDIDATES_WITH_SOURCE_GAP' if partial else 'CANDIDATES_FOUND'
        else:
            result['candidate_availability'] = 'SEARCH_INCOMPLETE' if partial else 'NO_CANDIDATES'
    output_file = apo_c_pilot.next_attempt(outcome)
    try:
        if baseline is None:
            case.update(read_holo_role(config['catalog'], case))
        else:
            prior_holo = baseline['holo_evidence']
            case.update({k: prior_holo.get(k) for k in ('qc', 'role', 'role_evidence')})
        holo, prepared = prepare_holo(case)
        result['holo_evidence'] = holo
        if prepared is not None:
            if baseline is None:
                plan_path = outcome.parent / f'candidate_plan-{search_digest}.json'
                if plan_path.exists():
                    plan = read_json(plan_path)
                    if plan['system_id'] != case['id'] or plan['search_sha256'] != search_digest:
                        raise ValueError('candidate plan differs from current target or search')
                else:
                    ordered, priority = apo_c_priority.prioritize(run, rank, prepared, config['max_unique_per_holo'])
                    chosen = select_first_run_budget(rank, ordered, config['per_source_rank'], config['max_unique_per_holo'])
                    plan = {'system_id': case['id'], 'search_sha256': search_digest,
                            'source_status': result['source_status'], 'ranking_counts': result['ranking_counts'],
                            'selected_candidates': chosen, 'priorities_top10': {k: v[:10] for k, v in ordered.items()},
                            'candidate_priority': {**priority, 'rule': FIRST_RUN_RANKING,
                                                   'method': 'reciprocal_rank_fusion', 'k': RRF_OFFSET,
                                                   'candidate_quota': len(chosen)}}
                    write_json(plan_path, plan)
                chosen = plan['selected_candidates']
                result.update({k: plan[k] for k in ('selected_candidates', 'priorities_top10', 'candidate_priority')},
                              candidate_plan=str(plan_path))
            elif rerank:
                result['candidate_priority'] = {'method': 'reciprocal_rank_fusion', 'k': RRF_OFFSET,
                                                'candidate_quota': len(chosen)}
            else:
                result.update({k: baseline[k] for k in ('candidate_priority', 'priorities_top10') if k in baseline})
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
            prior = {row['pdb']: row['evidence_path'] for row in baseline['candidates']} if baseline else {}

            def record_candidate(pdb, evidence, path):
                row = {'pdb': pdb, 'status': evidence['status'], 'evidence_path': str(path)}
                result['candidates'].append(row)
                if baseline is not None or not evidence['scientific_decision'].get('delivery_eligible'):
                    return
                try:
                    deliveries = deliver_observations(case, holo, prepared, [row], derived)
                    if deliveries:
                        first_path = outcome.parent / 'first_delivery.json'
                        if not first_path.exists():
                            write_json(first_path, {'system_id': case['id'], 'candidate': pdb,
                                                   'written_utc': now(), 'deliveries': deliveries})
                        result['first_delivery_receipt'] = str(first_path)
                        result['deliveries'].extend(deliveries)
                except (RuntimeError, OSError, ValueError, KeyError) as exc:
                    result.setdefault('delivery_failures', []).append({'pdb': pdb, 'error': f'{type(exc).__name__}: {exc}'})

            for item in chosen:
                candidate_path = outcome.parent / 'candidates' / f"{item['pdb']}.json"
                history = apo_c_pilot.attempts(candidate_path)
                if history:
                    evidence = read_json(history[-1])
                    if (source_hints(evidence.get('search_origins', {})) == source_hints(item['sources']) and
                            (not config['retry_failed'] or not has_technical_gap(evidence))):
                        record_candidate(item['pdb'], evidence, history[-1])
                        print('CANDIDATE_REUSE', case['id'], item['pdb'], flush=True)
                        continue
                candidate_path = apo_c_pilot.next_attempt(candidate_path)
                print('CANDIDATE_START', case['id'], item['pdb'], flush=True)
                try:
                    old_path = prior.get(item['pdb'])
                    if old_path is None:
                        original = candidate_original(run, item['pdb'], item['sources'], config)
                        if original[1] is None:
                            yield
                        evidence = inspect_candidate(run, case, prepared, item['pdb'], item['sources'], config, original)
                    else:
                        evidence = read_json(Path(old_path))
                        if evidence['search_origins'] != item['sources']:
                            raise ValueError('frozen evidence has different source hints')
                        evidence['replayed_from'] = old_path
                    def assess_experiment(experiment):
                        return literature.assess(experiment, source_deposit=deposition_context(prepared[0], holo),
                                                 fixed_fingerprint=lambda pdb: study_fingerprint(
                                                     pdb, evidence, case, run, config))

                    # Keep measurements and adjudication fixed for overlapping ranking relations.
                    if not rerank or old_path is None:
                        evidence['scientific_decision'] = apo_c_science.decide_candidate(
                            holo, prepared, evidence, assess_experiment)
                    evidence['materialized_assemblies'] = materialize_candidate(
                        evidence, item['pdb'], derived, config['max_derivative_assemblies_per_candidate'])
                except Exception as exc:
                    if baseline is not None and (not rerank or old_path is not None):
                        raise
                    evidence = {'status': 'TECHNICAL_FAILURE', 'pdb': item['pdb'],
                                'error': f'{type(exc).__name__}: {exc}',
                                'scientific_decision': {'status': 'TECHNICAL_GAP', 'rule': apo_c_science.RULE,
                                                        'reasons': ['CANDIDATE_INSPECTION_FAILED']}}
                evidence['search_origins'] = item['sources']
                write_json(candidate_path, evidence)
                record_candidate(item['pdb'], evidence, candidate_path)
                print('CANDIDATE_DONE', case['id'], item['pdb'], evidence['status'], flush=True)
            if baseline is not None:
                result['deliveries'] = deliver_observations(case, holo, prepared, result['candidates'], derived)
        else:
            result['candidate_review'] = 'MULTICHAIN_SCOPE_REQUIRES_REVIEW'
    except Exception as exc:
        failure = {'status': 'TECHNICAL_FAILURE', 'error': f'{type(exc).__name__}: {exc}'}
        if result['deliveries']:
            result['processing_failure'] = failure
        else:
            result['holo_evidence'] = failure
    write_json(output_file, result)
    print('RESULT', case['id'], result['holo_evidence']['status'],
          'selected', len(chosen), 'examined', len(result['candidates']), flush=True)
    return result


def run_in_download_order(tasks):
    pending = list(enumerate(tasks))
    results = [None] * len(pending)
    while pending:
        waiting = []
        for index, task in pending:
            try:
                next(task)
            except StopIteration as completed:
                results[index] = completed.value
            else:
                waiting.append((index, task))
        pending = waiting
    return results


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', required=True)
    p.add_argument('--catalog', required=True)
    p.add_argument('--index', required=True)
    p.add_argument('--samples', required=True)
    p.add_argument('--reuse-run', action='append', default=[],
                   help='prior C research Run with hash-checked fixed candidate originals')
    p.add_argument('--reuse-search-run', help='prior Run with the identical Holo search implementation and inputs')
    frozen_mode = p.add_mutually_exclusive_group()
    frozen_mode.add_argument('--replay-evidence-run', help='reassess experiment evidence against an unchanged Run shard')
    frozen_mode.add_argument('--rerank-from-run', help='compare RRF with the same frozen searches and each Holo quota')
    p.add_argument('--per-source-rank', type=int, default=3)
    p.add_argument('--max-unique-per-holo', type=int, default=6)
    p.add_argument('--max-candidate-download-attempts', type=int, default=10)
    p.add_argument('--author-documents', help='hash-verified manifest of raw publisher supplements')
    p.add_argument('--max-publisher-candidates', type=int, default=6,
                   help='maximum unresolved candidates per new shard to look up in publisher supplements (0 disables lookup)')
    p.add_argument('--shared-originals', help='directory of fixed-revision, hash-verified source receipts shared between shards')
    p.add_argument('--literature-cache', help='shared author-source cache, defaults to the parent of the Run directory')
    retry = p.add_mutually_exclusive_group()
    retry.add_argument('--retry-inspection', action='store_true',
                       help='retry failed inspection or writing against the recorded search and candidate order')
    retry.add_argument('--retry-failed', action='store_true',
                   help='Retry missing search pages, fixed inputs, and technical observations in an existing run')
    args = p.parse_args()
    if not 1 <= args.per_source_rank <= 10 or not 1 <= args.max_unique_per_holo <= 30:
        p.error('invalid candidate inspection budget')
    if not 0 <= args.max_publisher_candidates <= 30:
        p.error('publisher candidate budget must be between 0 and 30')
    if args.max_candidate_download_attempts < 0:
        p.error('invalid candidate download budget')
    if args.replay_evidence_run and (args.retry_failed or args.retry_inspection or args.reuse_search_run or args.reuse_run):
        p.error('evidence replay already fixes the search, measured structures and source gaps')
    if args.rerank_from_run and (args.retry_failed or args.retry_inspection or args.reuse_search_run):
        p.error('ranking replay requires a new Run and uses the frozen search')
    run, sample_path = Path(args.run), Path(args.samples)
    literature_cache = Path(args.literature_cache).resolve() if args.literature_cache else run.parent.resolve() / '.apo-literature-cache'
    deps = [Path(__file__), Path(apo_c_pilot.__file__), Path(apo_c_pair.__file__),
            Path(apo_c_rank_audit.__file__), Path(apo_c_science.__file__), Path(apo_c_triplet_regression.__file__),
            Path(apo_c_union_triage.__file__), Path(holo_source_mapping.__file__), Path(apo_c_literature.__file__),
            Path(apo_c_priority.__file__)]
    live_publisher_lookup = args.max_publisher_candidates if not (args.author_documents or args.replay_evidence_run or args.rerank_from_run) else 0
    if live_publisher_lookup:
        deps.append(Path(apo_c_author_documents.__file__))
    config = {'sample_spec_sha256': sha(sample_path), 'catalog': str(Path(args.catalog).resolve()),
              'index': str(Path(args.index).resolve()), 'per_source_rank': args.per_source_rank,
              'max_unique_per_holo': args.max_unique_per_holo, 'max_compressed_bytes': 20000000,
              'max_asymmetric_unit_atoms': 200000, 'max_assembly_instances': 128,
              'max_candidate_download_attempts': args.max_candidate_download_attempts,
              'max_derivative_assemblies_per_candidate': 0,
              'download_deadline_seconds': 30, 'download_order': 'round-robin-holo-missing-originals',
              'assembly_scope': 'mapped-receptor-site',
              'reuse_runs': [str(Path(x).resolve()) for x in args.reuse_run],
              'reuse_search_run': str(Path(args.reuse_search_run).resolve()) if args.reuse_search_run else None,
              'literature_cache': str(literature_cache),
              'code_sha256': {x.name: sha(x) for x in deps},
              'scientific_rule': apo_c_science.RULE, 'formal_pair_acceptance': False}
    if live_publisher_lookup:
        config['max_publisher_candidates'] = live_publisher_lookup
    if args.author_documents:
        config['author_documents'] = {'path': str(Path(args.author_documents).resolve()), 'sha256': sha(Path(args.author_documents))}
    if not args.replay_evidence_run:
        config['shared_originals'] = str(Path(args.shared_originals).resolve() if args.shared_originals
                                         else run.parent.resolve() / '.apo-fixed-originals')
    replay = args.rerank_from_run or args.replay_evidence_run
    previous = Path(replay).resolve() if replay else None
    replay_key = 'rerank_from_run' if args.rerank_from_run else 'replay_evidence_run'
    prior_cache = None
    if previous is None:
        config['ranking_rule'] = FIRST_RUN_RANKING
    if previous is not None:
        run = run.resolve()
        if run.exists() or run.is_relative_to(previous):
            raise ValueError('evidence replay requires a new Run outside the sealed baseline')
        frozen_config = read_json(previous / 'frozen_config.json')
        for key in ('catalog', 'index', 'per_source_rank', 'max_unique_per_holo', 'scientific_rule'):
            if config[key] != frozen_config[key]:
                raise ValueError(f'replay differs from frozen configuration: {key}')
        for name, digest in config['code_sha256'].items():
            mutable = {'apo_c_autopilot.py'} if args.rerank_from_run else {'apo_c_autopilot.py', 'apo_c_literature.py'}
            if name not in mutable and frozen_config['code_sha256'][name] != digest:
                raise ValueError(f'replay changes a frozen scientific or search dependency: {name}')
        prior_cache = Path(frozen_config['literature_cache']).resolve()
        if literature_cache.is_relative_to(prior_cache) or literature_cache.is_relative_to(previous):
            raise ValueError('replay requires a separate writable literature cache')
        config.update({replay_key: str(previous), 'prior_literature_cache': str(prior_cache),
                       'max_derivative_assemblies_per_candidate': 0})
        if args.rerank_from_run:
            source_run = frozen_config.get('replay_evidence_run')
            config['reuse_runs'] = list(dict.fromkeys([str(previous), *([source_run] if source_run else []),
                                                       *config['reuse_runs']]))
        else:
            config['max_candidate_download_attempts'] = 0
    started, clock = now(), time.monotonic()
    if not run.exists():
        if previous is None:
            apo_c_pilot.initial_samples(SimpleNamespace(samples=str(sample_path), run=str(run),
                                                        catalog=args.catalog, index=args.index))
        else:
            if (sha(sample_path) != frozen_config['sample_spec_sha256']
                    or read_json(sample_path) != read_json(previous / 'sample_spec.json')):
                raise ValueError('replay sample specification differs from the baseline')
            run.mkdir(parents=True)
            (run / 'code').mkdir()
            for name in ('sample_spec.json', 'samples.json'):
                shutil.copy2(previous / name, run / name)
            write_json(run / 'run.json', dict(read_json(previous / 'run.json'), started_utc=started,
                                             **{replay_key: str(previous)}))
        write_json(run / 'frozen_config.json', config)
        for path in deps:
            snapshot = run / 'code' / f'{path.stem}-{sha(path)[:12]}.py'
            if not snapshot.exists():
                shutil.copy2(path, snapshot)
    elif not (run / 'frozen_config.json').exists() or read_json(run / 'frozen_config.json') != config:
        raise ValueError('existing Run does not match the frozen configuration and code')
    if previous is not None:
        reuse_searches(run, previous)
    elif args.reuse_search_run:
        reuse_searches(run, args.reuse_search_run)
    config['retry_failed'] = args.retry_failed or args.retry_inspection
    scan_start = time.monotonic()
    cases = read_json(run / 'samples.json')
    if previous is None:
        apo_c_pilot.search_cases(SimpleNamespace(run=str(run), case=None, split=None, retry_failed=args.retry_failed))
        baselines = [None] * len(cases)
    else:
        baselines = [frozen_case(run, previous, case, config['max_unique_per_holo']) for case in cases]
    inspection_start = time.monotonic()
    literature = apo_c_literature.PrimaryLiterature(run, literature_cache, prior_cache=prior_cache)
    if args.author_documents:
        literature.author_documents = apo_c_literature.AuthorDocuments(args.author_documents)
    elif live_publisher_lookup:
        literature.publisher_discovery = apo_c_author_documents.CandidateDocuments(run, literature_cache, live_publisher_lookup)
    results = run_in_download_order(
        run_case(run, case, config, literature, baseline) for case, baseline in zip(cases, baselines))
    if previous is not None and any(result['holo_evidence']['status'] != 'HOLO_OBSERVED_SITE_READY' for result in results):
        raise ValueError('one or more frozen Holo cases did not complete; inspect per-target output')
    inspection_seconds = round(time.monotonic() - inspection_start, 3)
    states = Counter(candidate['status'] for result in results for candidate in result['candidates'])
    research = {}
    structure_only = Counter()
    for result in results:
        statuses = []
        for candidate in result['candidates']:
            evidence = read_json(Path(candidate['evidence_path']))
            status = evidence['scientific_decision']['status']
            statuses.append(status)
            if status != 'SUPPORTED_RESEARCH_PAIR' and any(
                o.get('decision', {}).get('site_observation', {}).get('status') == 'OBSERVED_NO_COMPETING_OCCUPANCY'
                for entity in evidence.get('entities', []) for o in entity.get('assembly_observations', [])
            ):
                structure_only[result['system_id']] += 1
        research[result['system_id']] = statuses
    observations = Counter(o['research_screen'] for result in results for candidate in result['candidates']
                           for e in read_json(Path(candidate['evidence_path'])).get('entities', [])
                           for o in e.get('assembly_observations', [])
                           if 'research_screen' in o)
    summary = {'started_utc': started, 'finished_utc': now(), 'elapsed_s': round(time.monotonic() - clock, 3),
               'search_elapsed_s': round(inspection_start - scan_start, 3), 'inspection_elapsed_s': inspection_seconds,
               'holo_dispatched': len(results), 'candidate_entry_outcomes': dict(states),
               'chain_observations': dict(observations), 'formal_pair_acceptances': 0,
               'scientific_rule': apo_c_science.RULE,
               'research_pair_outcomes': dict(Counter(s for values in research.values() for s in values)),
               'structural_candidates_with_observed_clear_site': sum(structure_only.values()),
               'holo_with_structural_candidates_with_observed_clear_site': len(structure_only),
               'holo_with_research_pair': sum('SUPPORTED_RESEARCH_PAIR' in values for values in research.values()),
               'research_pairs_delivered': sum(len(r.get('deliveries', [])) for r in results),
               'holo_with_delivered_research_pair': sum(bool(r.get('deliveries')) for r in results),
               'literature': dict(literature.stats),
               'limits': ['Search routes can fail independently.', 'Multi-chain targets remain review.',
                          'An observed no-contact scope is not proof of an unbound experiment.',
                          'No population accuracy, blind assessment or formal acceptance is claimed.']}
    n = len(list((run / 'receipts').glob('automation_run*.json'))) + 1 if (run / 'receipts').exists() else 1
    path = run / 'receipts' / f'automation_run{n:03d}.json'
    write_json(path, summary)
    print('SUMMARY', path, summary, flush=True)


if __name__ == '__main__':
    main()

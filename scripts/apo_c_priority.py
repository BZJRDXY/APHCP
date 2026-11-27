"""Use deposited receptor coverage and experiment metadata before spending the Apo inspection budget."""
from hashlib import sha256
import json

import requests

from apo_c_literature import deposited_apo
from apo_c_pair import reference_positions
from apo_c_pilot import cache_get


GRAPHQL = 'https://data.rcsb.org/graphql'
ENTRY_QUERY = '''
{ entries(entry_ids: IDS) {
    rcsb_id struct { title } exptl { method }
    polymer_entities { rcsb_polymer_entity_align {
        reference_database_name reference_database_accession
        aligned_regions { ref_beg_seq_id length }
    } }
} }
'''


def entry_priority(entry, accession, site):
    entities = entry.get('polymer_entities') or []
    coverage = None
    complete_annotations = bool(entities)
    for entity in entities:
        alignments = entity.get('rcsb_polymer_entity_align') or []
        complete_annotations &= bool(alignments)
        for alignment in alignments:
            if (alignment.get('reference_database_name') != 'UniProt'
                    or alignment.get('reference_database_accession') != accession):
                continue
            regions = alignment.get('aligned_regions') or []
            spans = [(r['ref_beg_seq_id'], r['ref_beg_seq_id'] + r['length'])
                     for r in regions if isinstance(r.get('ref_beg_seq_id'), int)
                     and isinstance(r.get('length'), int) and r['length'] > 0]
            if spans:
                count = sum(any(begin <= p < end for begin, end in spans) for p in site)
                coverage = max(coverage if coverage is not None else 0, count)
    title = (entry.get('struct') or {}).get('title') or ''
    xray = any(row.get('method') == 'X-RAY DIFFRACTION' for row in entry.get('exptl') or [])
    return {'mapped_site_positions': coverage, 'required_site_positions': len(site),
            'target_domain_absent': bool(coverage == 0 and complete_annotations and site),
            'depositor_unbound_hint': deposited_apo(title), 'xray': xray}


def order_rankings(ranked, entries, accession, site, depth):
    entries = {row['rcsb_id']: row for row in entries}
    evidence = {pdb: entry_priority(row, accession, site) for pdb, row in entries.items()}

    def priority(item):
        d = evidence.get(item['pdb'], {})
        overlap = d.get('mapped_site_positions')
        partial = overlap is not None and overlap < len(site)
        return partial, not d.get('depositor_unbound_hint'), not d.get('xray')

    adjusted = {}
    for source, items in ranked.items():
        head = [item for item in items[:depth] if not evidence.get(item['pdb'], {}).get('target_domain_absent')]
        adjusted[source] = sorted(head, key=priority) + items[depth:]
    return adjusted, evidence


def prioritize(run, ranked, prepared, depth):
    ids = sorted({item['pdb'] for items in ranked.values() for item in items[:depth]})
    if not ids:
        return ranked, {'status': 'NO_CANDIDATES'}
    block, _, _, site, accessions = prepared
    authors = {row['author_chain'] for row in site.values()}
    accession = accessions[0] if len(accessions) == 1 and len(authors) == 1 else None
    mapping = reference_positions(block, accession, next(iter(authors))) if accession else {}
    positions = {mapping[p] for p in site} if site.keys() <= mapping.keys() else set()
    query = ENTRY_QUERY.replace('IDS', json.dumps(ids))
    key = sha256(query.encode()).hexdigest()[:20]
    try:
        def validate(body, _):
            payload = json.loads(body)
            if payload.get('errors') or not isinstance((payload.get('data') or {}).get('entries'), list):
                raise ValueError('invalid candidate metadata response')

        payload, receipt = cache_get(run, 'candidate_priority/' + key, GRAPHQL,
                                     {'query': query}, timeout=10, accepted=(200,), validate=validate)
        entries = [row for row in json.loads(payload)['data']['entries'] if row and row.get('rcsb_id') in ids]
        ordered, evidence = order_rankings(ranked, entries, accession, positions, depth)
        return ordered, {'status': 'RANKED_WITH_METADATA', 'entry_count': len(entries),
                         'target_accession': accession, 'target_positions': sorted(positions),
                         'entries': evidence, 'receipt': receipt}
    except (requests.RequestException, RuntimeError, ValueError) as exc:
        return ranked, {'status': 'METADATA_UNAVAILABLE', 'error_type': type(exc).__name__}

"""Independently check preliminary holdout identities against the prior pilot inputs."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

import gemmi


def load(path):
    return json.loads(path.read_text(encoding='utf8'))


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def metadata(path, labels=None):
    b = gemmi.cif.read(str(path)).sole_block()
    asym = b.get_mmcif_category('_struct_asym.')
    label_entity = dict(zip(asym.get('id', []), asym.get('entity_id', [])))
    entities = {str(label_entity[x]) for x in labels or [] if x in label_entity}
    refs = b.get_mmcif_category('_struct_ref.')
    accessions = {str(refs['pdbx_db_accession'][i]) for i, db in enumerate(refs.get('db_name', []))
                  if str(db).upper() in {'UNP', 'UNIPROT'} and refs['pdbx_db_accession'][i]
                  and (not entities or str(refs['entity_id'][i]) in entities)}
    poly = b.get_mmcif_category('_entity_poly.')
    seqs = set()
    for i, entity in enumerate(poly.get('entity_id', [])):
        if entities and str(entity) not in entities:
            continue
        value = poly.get('pdbx_seq_one_letter_code_can', [])[i]
        if value:
            seqs.add(hashlib.sha256(''.join(str(value).split()).upper().encode('ascii')).hexdigest())
    cat = b.get_mmcif_category('_citation.')
    cites = set()
    for i, value in enumerate(cat.get('id', [])):
        if str(value).lower() != 'primary':
            continue
        for col, prefix in [('pdbx_database_id_PubMed', 'pmid:'), ('pdbx_database_id_DOI', 'doi:')]:
            vals = cat.get(col, [])
            if i < len(vals) and vals[i] not in (None, False, '', '.', '?'):
                cites.add(prefix + str(vals[i]).lower())
    return accessions, seqs, cites


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True)
    p.add_argument('--selection-name', default='provisional_nominees_v2.json')
    p.add_argument('--output-name', default='independent_group_audit.json')
    args = p.parse_args()
    run = Path(args.run)
    if Path(args.selection_name).name != args.selection_name or Path(args.output_name).name != args.output_name:
        raise ValueError('selection and output names must be simple filenames')
    selection_path = run / 'holdout_preparation' / args.selection_name
    selection = load(selection_path)
    baseline = Path(load(run / 'run.json')['baseline_run'])
    samples = load(baseline / 'samples.json')
    paths = {Path(s['source_uri']) for s in samples}
    prior_pdbs = {s['pdb'] for s in samples}
    for parent in [baseline / 'pairs', baseline / 'multichain']:
        for f in parent.glob('*/attempt*/pair.json'):
            d = load(f)
            source = d.get('candidate_source') or {}
            if source.get('uri') and Path(source['uri']).is_file():
                paths.add(Path(source['uri']))
                name = d.get('candidate_pdb') or d.get('request', {}).get('candidate')
                if name:
                    prior_pdbs.add(str(name).upper())
    for f in (baseline / 'raw').glob('*.cif.gz'):
        paths.add(f)
        prior_pdbs.add(f.name.split('-v')[0].upper())
    old_acc, old_seq, old_cit = set(), set(), set()
    for f in paths:
        a, s, c = metadata(f)
        old_acc.update(a)
        old_seq.update(s)
        old_cit.update(c)
    seen_pdb, seen_acc, seen_seq, seen_cit = set(), set(), set(), set()
    rows = []
    for sample in selection['selected']:
        path = Path(sample['source_uri'])
        labels = {x[1] for x in sample['holo_identity']['binding_instances']}
        acc, seq, cites = metadata(path, labels)
        saved = sample['protein_and_citation_metadata']
        checks = {'source_hash': digest(path) == sample['holo_identity']['source_sha256'],
                  'not_prior_or_duplicate_pdb': sample['pdb'] not in prior_pdbs and sample['pdb'] not in seen_pdb,
                  'source_accessions_match': sorted(acc) == saved['accessions'] and bool(acc),
                  'source_sequences_match': sorted(seq) == saved['sequence_sha256'] and bool(seq),
                  'source_citations_match': sorted(cites) == saved['primary_citation_ids'] and bool(cites),
                  'no_prior_or_within_accession_overlap': not (acc & (old_acc | seen_acc)),
                  'no_prior_or_within_exact_sequence_overlap': not (seq & (old_seq | seen_seq)),
                  'no_prior_or_within_primary_citation_id_overlap': not (cites & (old_cit | seen_cit))}
        rows.append({'pdb': sample['pdb'], 'checks': checks, 'pass': all(checks.values())})
        seen_pdb.add(sample['pdb'])
        seen_acc.update(acc)
        seen_seq.update(seq)
        seen_cit.update(cites)
    result = {'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
              'selection': str(selection_path), 'selection_sha256': digest(selection_path),
              'previous_v1_status': 'SUPERSEDED: a PDB reference was incorrectly included as a UniProt accession.',
              'selected': len(rows), 'multi_chain': sum(x['multi_chain'] for x in selection['selected']),
              'checks_per_nominee': 8, 'all_pass': all(x['pass'] for x in rows), 'rows': rows,
              'limitations': ['Exact identifier checks cannot prove independence of unlinked publication aliases, related experiments or protein families.',
                              'No Apo candidate search or scientific adjudication has been performed for these nominees.',
                              'Nominees and future test rules remain provisional until independently reviewed and frozen.']}
    out = run / 'holdout_preparation' / args.output_name
    with out.open('x', encoding='utf8') as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    print(json.dumps({'path': str(out), 'selected': len(rows), 'all_pass': result['all_pass'],
                      'failed': [x['pdb'] for x in rows if not x['pass']]}), flush=True)
    if not result['all_pass']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()

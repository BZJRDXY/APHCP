"""Read exactly one candidate from the frozen full M1 inventory on Pro5000."""
import json
import re
import sys
from pathlib import Path

BATCH = 'holo-m1-full-20261121a'
FIELDS = ('pdb_id', 'title', 'first_release_date', 'rcsb_latest_revision', 'rcsb_latest_revision_date',
          'selected_revision', 'methods', 'resolution_angstrom', 'model_ids', 'model_count', 'rcsb_assembly_ids',
          'protein_entity_count', 'nonpolymer_ccds', 'partner_evidence', 'discoveries', 'family_evidence_types',
          'plinder_index_presence', 'plinder_annotated_assemblies', 'plinder_candidate_system_count',
          'plinder_source_revision_if_available', 'plinder_revision_status', 'metadata_review')


def canonical_id(value):
    value = value.upper()
    return value[8:] if re.fullmatch(r'PDB_0000[A-Z0-9]{4}', value) else value


def confined(parent, relative):
    p = (parent / relative).resolve()
    try:
        p.relative_to(parent)
    except ValueError as exc:
        raise ValueError('Inventory path outside its approved directory') from exc
    return p


def find_candidate(root, pdb):
    root = root.resolve()
    pdb = canonical_id(pdb)
    if not re.fullmatch(r'(?:[A-Z0-9]{4}|PDB_[A-Z0-9]{8})', pdb):
        raise ValueError('Invalid PDB ID')
    plan = confined(root, f'batches/{BATCH}/planning/plan.json')
    plan_data = json.loads(plan.read_text())
    if plan_data.get('batch_id') != BATCH:
        raise ValueError('Unexpected inventory batch')
    windows = plan_data['windows']
    if not 0 < len(windows) <= 100:
        raise ValueError('Unexpected window count')
    unavailable = []
    for w in windows:
        run = Path(w['run_directory'])
        if run.parent != root/'runs' or not re.fullmatch(r'holo-m1-[a-z0-9-]+', run.name):
            raise ValueError('Unexpected M1 window directory')
        try:
            s = json.loads(confined(root, f'runs/{run.name}/search/summary.json').read_text())
            if not s.get('concordant') or not isinstance(s.get('entry_ids'), list):
                unavailable.append(run.name)
                continue
        except (OSError, json.JSONDecodeError):
            unavailable.append(run.name)
            continue
        if pdb not in s['entry_ids']:
            continue
        candidate = None
        path = confined(root, f'runs/{run.name}/inventory/candidate_entries.jsonl')
        with path.open(encoding='utf8') as stream:
            for line in stream:
                if pdb not in line:
                    continue
                row = json.loads(line)
                if row.get('pdb_id') == pdb:
                    candidate = {k: row[k][:40] if isinstance(row.get(k), list) else row.get(k) for k in FIELDS}
                    break
        return {'state': 'found' if candidate else 'incomplete', 'candidate': candidate, 'pdb_id': pdb,
                'window': run.name, 'bounds': [w['start_inclusive'], w['end_exclusive']], 'source': str(path),
                'other_unavailable_windows': len(unavailable)}
    return {'state': 'incomplete' if unavailable else 'absent', 'candidate': None, 'pdb_id': pdb,
            'other_unavailable_windows': len(unavailable)}


if __name__ == '__main__':
    from datetime import datetime, timezone
    try:
        result = find_candidate(Path(sys.argv[1]), sys.argv[2])
        print(json.dumps({'data': result, 'queried_at': datetime.now(timezone.utc).isoformat()}, ensure_ascii=True))
    except (OSError, ValueError, KeyError, TypeError) as exc:
        sys.exit(f'M1 read-only lookup: {exc}')

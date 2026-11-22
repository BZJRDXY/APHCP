"""Resolve immutable M3 decisions from sealed Run receipts, without using Run order."""
from collections import defaultdict


def resolve(items):
    """Items are (run_id, receipt_path, assessment) for one canonical observation."""
    by_id = {}
    for run, path, row in items:
        key = row.get('system_id') or ':'.join((str(row['batch']), row['pdb_id'], row['candidate_key']))
        decision = f'{run}/{key}'
        record = {'decision_id': decision, 'run_id': run, 'report': path, 'batch': row['batch'],
                  'candidate_key': row['candidate_key'], 'disposition': row['disposition'],
                  'files': row.get('files', {}), 'supersedes': row.get('supersedes', []),
                  'decision_basis': row.get('decision_basis')}
        if decision in by_id:
            previous = by_id[decision][0]
            if any(previous[k] != record[k] for k in ('disposition', 'supersedes', 'decision_basis')):
                raise ValueError('conflicting_decisions_in_one_run')
            # Aliases have the same adjudication. Prefer the lexicographically stable key.
            if record['candidate_key'] >= previous['candidate_key']:
                continue
        by_id[decision] = record, row
    pending = {k: set(r['supersedes']) for k, (r, _) in by_id.items()}
    if any(parent not in by_id for parents in pending.values() for parent in parents):
        raise ValueError('superseded_decision_missing_for_same_identity')
    predecessors = set().union(*pending.values())
    children = defaultdict(list)
    ready = [k for k, parents in pending.items() if not parents]
    for k, parents in pending.items():
        for parent in parents:
            children[parent].append(k)
    visited = 0
    while ready:
        parent = ready.pop()
        visited += 1
        for child in children[parent]:
            pending[child].remove(parent)
            if not pending[child]:
                ready.append(child)
    if visited != len(by_id):
        raise ValueError('decision_supersession_cycle')
    heads = sorted(set(by_id) - predecessors)
    outcomes = {by_id[k][0]['disposition'] for k in heads}
    if len(outcomes) == 1:
        shown = by_id[heads[0]][1]
        effective = heads[0] if len(heads) == 1 else None
        state = 'EFFECTIVE' if effective else 'CONCORDANT'
    else:
        prior = by_id[heads[0]][1]
        shown = {k: prior[k] for k in ('system_id', 'identity', 'pdb_id', 'batch', 'candidate_key', 'revision') if k in prior}
        common_files = {k: v for k, v in prior.get('files', {}).items()
                        if all(by_id[d][1].get('files', {}).get(k) == v for d in heads)}
        shown.update(disposition='CONFLICT', reasons=['incompatible_unsuperseded_decisions'], files=common_files)
        effective, state = None, 'CONFLICT'
    return shown | {'current_decision_id': effective, 'decision_state': state, 'current_decision_candidates': heads,
                    'provenance': [by_id[k][0] for k in sorted(by_id)]}

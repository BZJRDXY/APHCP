"""Select traceable first-round challenge cases from the frozen Holo B audit."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import random
import sqlite3

import gemmi


def read_block(path):
    with gzip.open(path, 'rt', encoding='utf-8') if str(path).endswith('.gz') else open(path, encoding='utf-8') as f:
        return gemmi.cif.read_string(f.read()).sole_block()


def value(block, tag):
    answer = block.find_value(tag)
    return None if answer in (None, '', '.', '?') else gemmi.cif.as_string(answer)


def metadata(path):
    block = read_block(path)
    citation = block.get_mmcif_category('_citation.')
    ids = citation.get('id', [])
    pos = ids.index('primary') if 'primary' in ids else None
    def get(key):
        values = citation.get(key, [])
        return values[pos] if pos is not None and len(values) > pos else None
    result = {'title': get('title'), 'doi': get('pdbx_database_id_DOI'),
              'pmid': get('pdbx_database_id_PubMed'), 'year': get('year'),
              'deposition_title': value(block, '_struct.title')}
    result['experiment_details'] = {tag: [v for v in block.find_values(tag) if v not in ('.','?')]
        for tag in ('_exptl_crystal_grow.pdbx_details', '_exptl_crystal_grow.details',
                    '_exptl.details', '_diffrn.details')}
    result['entities'] = block.get_mmcif_category('_entity.')
    return result


def publication_keys(citation):
    keys = set()
    if citation.get('doi'):
        keys.add(('doi', citation['doi'].lower().strip()))
    if citation.get('pmid') not in (None, '', '-1'):
        keys.add(('pmid', str(citation['pmid']).strip()))
    if citation.get('title'):
        keys.add(('title', ' '.join(citation['title'].lower().split())))
    return keys


def eligible_candidate(row):
    def qcs(q, fields):
        return isinstance(q, dict) and all(q.get(f) == 'PASS' for f in fields)
    return (row['disposition'] == 'REVIEW'
            and json.loads(row['reasons']) == ['target_role_evidence_insufficient_or_conflicting']
            and row['match_level'] == 'SAME_REVISION_INSTANCE_CHAIN_MATCH'
            and row['annotation_interpretation'] == 'PLINDER_PROPER_RULE_SUPPORT'
            and json.loads(row['conflicts']) == []
            and qcs(json.loads(row['old_qc']), ('coordinate_identity','experimental_fidelity','topology','partner_scope_qc'))
            and qcs(json.loads(row['new_qc']), ('coordinate_identity','experimental_fidelity','topology','partner_scope')))


def group(row):
    if eligible_candidate(row):
        return 'old_review_new_supported'
    if row['disposition'] == 'REGISTERED':
        return 'registered'
    if row['disposition'] == 'REVIEW':
        return 'review_not_supported'
    if row['disposition'] == 'PREPARATION_USE_RECORDED':
        return 'preparation_recorded'
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--database', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--per-stratum', type=int, default=8)
    p.add_argument('--seed', type=int, default=20261125)
    args = p.parse_args()
    db = sqlite3.connect(args.database.resolve().as_uri()+'?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    metadata_cache = {}
    def get_metadata(sid):
        raw = db.execute('select document from observation where system_id=?',(sid,)).fetchone()[0]
        obs = json.loads(raw)
        path = obs['files']['source_mmcif']['path']
        if path not in metadata_cache:
            metadata_cache[path] = metadata(path)
        return obs, metadata_cache[path]

    used_pdb = {'1MBN','5MYQ','8RNP','5MO0'}
    used_papers = set()
    development = {}
    for pdb in sorted(used_pdb):
        row = db.execute('select system_id from result where pdb_id=? order by system_id limit 1',(pdb,)).fetchone()
        if row:
            _, m = get_metadata(row[0])
            used_papers.update(publication_keys(m))
            development[pdb] = {k:m[k] for k in ('title','doi','pmid')}

    groups = defaultdict(lambda: defaultdict(list))
    counts = Counter()
    for row in db.execute('select * from result order by system_id'):
        g = group(row)
        if g:
            counts[g] += 1
            if row['pdb_id'] not in used_pdb:
                groups[g][row['pdb_id']].append(row['system_id'])
    rng = random.Random(args.seed)
    result = {'seed':args.seed,'per_stratum':args.per_stratum,'populations':dict(counts),
              'development':development,'strata':{}}
    for g in ('old_review_new_supported','registered','review_not_supported','preparation_recorded'):
        entries = sorted(groups[g])
        rng.shuffle(entries)
        chosen = []
        for pdb in entries:
            if pdb in used_pdb:
                continue
            systems = sorted(groups[g][pdb])
            sid = systems[rng.randrange(len(systems))]
            try:
                obs, m = get_metadata(sid)
            except (OSError, ValueError, RuntimeError) as exc:
                # An inaccessible source is a selected but unresolved case, never silently replaced.
                chosen.append({'system_id':sid,'pdb_id':pdb,'source_read_error':str(exc)})
                used_pdb.add(pdb)
                if len(chosen) == args.per_stratum:
                    break
                continue
            keys = publication_keys(m)
            if keys & used_papers:
                continue
            used_papers.update(keys)
            used_pdb.add(pdb)
            detail = dict(db.execute('select * from result where system_id=?',(sid,)).fetchone())
            chosen.append({'system_id':sid,'pdb_id':pdb,'identity':obs['identity'],
                'source_mmcif':obs['files']['source_mmcif'], 'recorded_role':obs.get('role_evidence'),
                'publication_independence': 'keys_distinct' if keys else 'unverified_missing_primary_citation',
                'citation':m, 'research_result':detail})
            if len(chosen) == args.per_stratum:
                break
        result['strata'][g] = {'eligible_pdbs_before_paper_dedup':len(entries), 'selected':chosen}
    with args.output.open('x',encoding='utf-8') as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    print(json.dumps({'populations':result['populations'],
        'selected':{g:[x['pdb_id'] for x in d['selected']] for g,d in result['strata'].items()}}, ensure_ascii=False))


if __name__ == '__main__':
    main()

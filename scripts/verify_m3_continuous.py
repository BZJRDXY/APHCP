"""Verify one completed M3 increment against its immutable per-phase evidence."""
import argparse
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from holo_m3_observations import STAGE, lookup_records
from holo_m3_register import connect, database, read, require, sha, write


def storage(paths):
    visited = set(); allocated = logical = 0
    for path in paths:
        files = [path] if path.is_file() else path.rglob('*')
        for p in files:
            if not p.is_file() or p.is_symlink():
                continue
            st = p.stat(); inode = (st.st_dev,st.st_ino)
            if inode in visited:
                continue
            visited.add(inode); logical += st.st_size; allocated += getattr(st,'st_blocks',0)*512
    return {'files':len(visited),'logical_bytes':logical,'allocated_bytes':allocated}


def verify(root, controller):
    first = controller/'first_increment_receipt.json'
    result = read(first if first.is_file() else controller/'increment_receipt.json')
    plan = read(controller/'increment_plan.json'); original = Path(plan['m2_run'])
    preflight = read(original/'inputs/preflight_receipt.json')
    require(sha(original/'frozen.json')==plan['m2_manifest_sha256'], 'm2_freeze_changed')
    frozen = {r['pdb_id']:r for r in read(original/'frozen.json')['files']}
    require(len(frozen)==result['frozen_entries']==preflight['frozen_entries'], 'frozen_denominator_differs')
    batch_cases = {}; last_entries = {}; files_received = set(); m2_sources = set(); sample = {}; stage_times = []
    for path in map(Path,result['phases']):
        phase = read(path/'phase_receipt.json'); source = Path(phase['m2_run']); m2_sources.add(source)
        evidence = read(path/'candidate_evidence.json')
        require(set(phase['selected'])==set(evidence['entries'])==set(evidence['entry_evidence']), 'phase_entry_count_differs')
        for p,row in evidence['entry_evidence'].items():
            require(p in frozen, 'entry_outside_freeze')
            last_entries[p]=row
            proof=source/'source_proofs'/f'{p}.json'
            if proof.exists():
                item=read(proof)
                if item['status']=='PASS':
                    require(item['identity']==frozen[p], 'proof_freeze_disagrees')
                    files_received.add(p)
        reports = path/'review_results.json'
        reviews = {(x['pdb_id'],x['candidate_key']):x for x in read(reports)} if reports.exists() else {}
        for row in evidence['candidates']:
            k = row['pdb_id'],row['candidate_key']; review = reviews.get(k,{})
            batch_cases[k] = row,review
            if review.get('system_id'):
                sample.setdefault((row['route'],review['disposition']), (k,review['system_id']))
        where=path/'registration_receipt.json'; reg_seconds=read(where)['analysis_seconds'] if where.is_file() else None
        stage_times.append({'phase':path.name,'selected':len(phase['selected']), 'chemistry_s':phase['chemistry']['elapsed_s'],
                            'm3_register_s':reg_seconds,'catalog_s':phase['catalog']['elapsed_s'],'downstream_total_s':phase['elapsed_s']})
    require(set(last_entries)==set(frozen), 'final_entry_coverage_differs')
    require(not (set(frozen)-files_received-set(result['pending_technical'])), 'unaccounted_for_original')
    effective=Path(result['catalog_run']); db_path=database(root/'batches'/effective.name)
    require(read(controller/'catalog_head.json')['catalog_run']==effective.name,'catalog_head_differs')
    examples=[]
    for (route,state), ((p,k),sid) in sorted(sample.items()):
        rows=lookup_records(effective,system=sid,verify_files=True)
        require(len(rows)==1 and rows[0]['status']==state and len(rows[0]['files'])>=3, 'catalog_sample_did_not_read_back')
        examples.append({'pdb_id':p,'system_id':sid,'route':route,'status':state})
    phase_labels = [Path(p).name+'-observations' for p in result['phases']]
    with connect(db_path,readonly=True) as db:
        buckets={r[0]:r[1] for r in db.execute('''SELECT status,COUNT(*) FROM run_item WHERE stage=? AND object_type='observation'
            AND run_id IN (%s) GROUP BY status''' % ','.join('?' for _ in phase_labels),(STAGE,*phase_labels))}
        catalog_count=db.execute('SELECT COUNT(DISTINCT object_id) FROM run_item WHERE stage=? AND output_uri IS NOT NULL',(STAGE,)).fetchone()[0]
        new_catalog_count=db.execute('''SELECT COUNT(DISTINCT object_id) FROM run_item WHERE stage=? AND output_uri IS NOT NULL
            AND run_id IN (%s)''' % ','.join('?' for _ in phase_labels),(STAGE,*phase_labels)).fetchone()[0]
        item_filter = 'stage=? AND output_uri IS NOT NULL AND run_id IN (%s)' % ','.join('?' for _ in phase_labels)
        params = (STAGE, *phase_labels)
        queryable_entries = db.execute('SELECT COUNT(DISTINCT pdb_id) FROM run_item WHERE '+item_filter, params).fetchone()[0]
        unknown_types = db.execute("SELECT COUNT(DISTINCT object_id) FROM run_item WHERE "+item_filter+
                                   " AND json_extract(error_message,'$.type.primary_class')='UNKNOWN'", params).fetchone()[0]
        require(db.execute('PRAGMA quick_check').fetchone()[0]=='ok' and db.execute('SELECT COUNT(*) FROM release').fetchone()[0]==0,
                'catalog_database_not_valid_unpublished_staging')
    reasons = Counter(reason for _,r in batch_cases.values() for reason in r.get('reasons',[]))
    network=sum(a['body_bytes'] for m in m2_sources for pf in (m/'source_proofs').glob('*.json') for a in read(pf).get('attempts',[]))
    old = Path(read(Path(result['phases'][0])/'phase_receipt.json')['previous_catalog_run'])
    prior_db = database(root/'batches'/old.name)
    build = read(original/'summary.json'); m2_time = build['invocation']['elapsed_s']; fetch_time=read(original/'last_fetch_invocation.json')['elapsed_s']
    statuses = Counter(x['state'] for x in last_entries.values()); structural = Counter(r['prior_status'] for r,_ in batch_cases.values())
    picks = Counter(r['disposition'] for _,r in batch_cases.values() if 'disposition' in r)
    route_counts = Counter(x['route'] or 'UNKNOWN' for x in last_entries.values())
    labels = Counter(r.get('export_reason') for r,_ in batch_cases.values() if r.get('export_reason'))
    tmax = sorted(((r.get('elapsed_seconds',0),r.get('pdb_id'),r.get('chosen_systems')) for r in
                   (read(p) for p in (original/'entry_results').glob('*.json'))),reverse=True)[:5]
    return {'status':result['status'],'verified_utc':datetime.now(timezone.utc).isoformat(), 'window':preflight['window'],
            'window_entries':preflight['window_entries'],'window_exclusions':preflight['exclusions'],
            'archive_deferred':preflight['archive_deferred'],'frozen_entries':len(frozen),'originals_verified':len(files_received),
            'entry_statuses':dict(statuses),'entry_routes':dict(route_counts),'candidates':len(batch_cases),
            'structural_states':dict(structural),'review_dispositions':dict(picks),'catalog_outcomes':buckets,
            'new_catalog_queryable':new_catalog_count,'new_queryable_entries':queryable_entries,
            'new_queryable_unknown_type':unknown_types,'catalog_queryable_total':catalog_count,'examples':examples,
            'pending_technical':result['pending_technical'],'top_non_export_reasons':dict(labels.most_common(8)),
            'top_m3_review_reasons':dict(reasons.most_common(8)), 'archive_network_bytes':network,
            'archival_batch':storage(m2_sources),'main_m2_artifacts':build['sizes'],
            'm3_storage':storage([controller, *map(Path,result['phases']),
                                  *(database(root/'batches'/name.name) for name in map(Path,result['phases']))]),
            'catalog_database_bytes':db_path.stat().st_size,
            'catalog_database_growth_bytes':db_path.stat().st_size-prior_db.stat().st_size,
            'stage_times':stage_times,'m2_seconds':m2_time,'m2_fetch_seconds':fetch_time,'preflight_seconds':preflight['elapsed_s'],
            'slow_entries':tmax,'pipeline_seconds':result['elapsed_s'], 'catalog_run':effective.name}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',required=True,type=Path)
    parser.add_argument('--run',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    a=parser.parse_args()
    res=verify(a.root.resolve(),a.run.resolve())
    write(a.output,res)
    print({k:res[k] for k in ('status','frozen_entries','originals_verified','candidates','new_catalog_queryable','catalog_queryable_total')})

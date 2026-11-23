"""Verify a bounded M3 increment and publish the rolling-production acceptance receipt."""
import argparse
from pathlib import Path
import subprocess

from holo_m3_observations import lookup_records, STAGE
from holo_m3_register import connect, database, read, require, sha, utc_now, write
from verify_m3_continuous import verify


def accept(a):
    root=a.root.resolve(strict=True); run=a.run.resolve(strict=True)
    comparison=read(a.comparison)
    require(len(comparison)>=2 and all(v['exact_result_equal_to_old_and_sealed'] for v in comparison.values()),
            'missing_readback_performance_control')
    resume_control=read(a.resume_control)
    require(resume_control['targeted_retry_verified'] and resume_control['same_input_reused'],
            'missing_targeted_retry_control')
    # The normal runner performs targeted technical retries. A new launch may continue
    # remaining technical gaps; keep the initial real run receipt before that occurs.
    first=run/'first_pipeline_attempt.json'
    if not first.is_file():
        write(first,read(run/'increment_receipt.json'))
    result=read(run/'increment_receipt.json')
    if result['pending_technical']:
        subprocess.run(['bash',str(a.launch)],check=True)
        result=read(run/'increment_receipt.json')
    require(result['status']=='STAGED_COMPLETE' and not result['pending_technical'],'small_batch_needs_technical_recovery')
    first_complete=run/'first_increment_receipt.json'
    if not first_complete.is_file():
        write(first_complete,result)
    result=verify(root,run)
    require(300 <= result['frozen_entries'] == result['originals_verified'] <= 500 and
            result['new_catalog_queryable']>0, 'small_batch_missing_verified_inventory')
    require({'REGISTERED','REVIEW'}<={x['status'] for x in result['examples']}, 'new_confirmed_or_review_not_retrievable')
    phase=root/'runs'/result['catalog_run']
    registry=database(root/'batches'/phase.name)
    seed=Path(read(run/'increment_plan.json')['seed_run'])
    before=lookup_records(seed,system=a.old_system,verify_files=True)
    after=lookup_records(phase,system=a.old_system,verify_files=True)
    require(len(before)==len(after)==1 and before==after and len(after[0]['files'])>=3,'old_catalog_record_changed')
    with connect(registry,readonly=True) as conn:
        require(conn.execute('SELECT COUNT(*) FROM release').fetchone()[0]==0,'staged_registry_contains_release')
        active=conn.execute('SELECT COUNT(DISTINCT object_id) FROM run_item WHERE stage=? AND status=? AND output_uri IS NOT NULL',
                            (STAGE,'REGISTERED')).fetchone()[0]
        unknown=conn.execute('''SELECT object_id FROM run_item WHERE run_id=? AND stage=? AND output_uri IS NOT NULL
                               AND json_extract(error_message,'$.type.primary_class')='UNKNOWN' LIMIT 1''',
                             (phase.name+'-observations',STAGE)).fetchone()
    if unknown:
        u=lookup_records(phase,system=unknown[0],verify_files=True)
        require(len(u)==1 and u[0]['type']['primary_class']=='UNKNOWN' and len(u[0]['files'])>=3,'unknown_type_not_retrievable')
    tracked=[registry]
    for folder in map(Path,read(first_complete)['phases']):
        tracked.extend(folder/name for name in ('phase_receipt.json','candidate_evidence.json','review_results.json') if (folder/name).exists())
    pre={p:sha(p) for p in tracked}
    subprocess.run(['bash',str(a.launch)],check=True)
    replay=read(run/'increment_receipt.json')
    require(not replay['executed_phases'] and replay['reused_phases']==len(replay['phases']) and
            all(sha(p)==h for p,h in pre.items()),'small_batch_replay_changed_frozen_evidence')
    result.update(accepted=True,accepted_utc=utc_now(),same_input_reused=True, replay_seconds=replay['elapsed_s'],
                  unified_lookup_verified=True,targeted_retry_verified=True,previous_targeted_retry_reused=True,previous_lookup=a.old_system,example_unknown=unknown[0] if unknown else None,
                  catalog_registered_total=active,readback_comparison=comparison)
    write(a.output,result)
    print({k:result[k] for k in ('accepted','frozen_entries','new_catalog_queryable','catalog_queryable_total','catalog_registered_total')},flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','run','comparison','resume-control','launch','output'):
        p.add_argument('--'+name,required=True,type=Path)
    p.add_argument('--old-system',required=True)
    accept(p.parse_args())
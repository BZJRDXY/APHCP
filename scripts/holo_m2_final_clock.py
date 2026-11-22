"""Record actual milestones and immutable paths for the final batch handoff."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess

from holo_m2_final_verify import read, sha


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('run',type=Path); a=ap.parse_args(); r=a.run.resolve()
    start=re.findall(r'^start (.+)$',(r/'launcher.log').read_text(),re.M)
    first=read(r/'first_pass'/'summary.json'); resume=re.findall(r'^finish start (.+)$',(r/'finish.log').read_text(),re.M)
    file=r/'handoff'/'summary.json'; hand=read(file); end=datetime.fromtimestamp(file.stat().st_mtime,timezone.utc)
    adopt=r/'parallel_tail'/'adoption.json'; continued=read(adopt) if adopt.exists() else None
    assets=subprocess.check_output(['du','-sB1','--apparent-size',str(r)],text=True).strip().split()[0]
    clock={'first_started_utc':start[0],'selection_launches_utc':start,'continuation_start_utc':resume,
           'first_structure_checkpoint_utc':first['updated_utc'],'final_results_utc':end.isoformat(),'total_apparent_bytes_unique_hardlinks':int(assets),
           'always_not_counted':'Transport attempts whose body lengths were unavailable and SSH hop bytes.',
           'main_interventions':['Selection retry after nullable M1 counts.',
                                 'Authorized exact-length alternate route filled incomplete remote transfers; one manual fallback.',
                                 'Replaced the last redundant slow transfer with local source revalidation.',
                                 'Reproduced point symmetry issue and checked only the matching data scope.'],
           'structure_tail_partition':continued,'confirmed_registry_writes':0}
    (r/'execution_clock.json').write_text(json.dumps(clock,indent=2))
    batch=r.parents[1]/'batches'/r.name; assert batch.is_dir()
    receipt={'batch':batch.name,'run':str(r),'status':'M2_COMPLETE_WITH_EXPLICIT_REVIEW','counts':hand['final'],
             'manifest_sha256':hand['frozen_input_sha256'],'first_pass_seal_sha256':hand['first_pass_seal_sha256'],
             'files':{n:{'path':str(r/n),'sha256':sha(r/n)} for n in 
                      ('execution_clock.json','handoff/summary.json','handoff/entry_results.csv','handoff/system_results.csv','m3_first_review.json','independent_audit.json')},
             'production_registration':'Not started; first-review preparation only.'}
    (batch/'m2_final_receipt.json').write_text(json.dumps(receipt,indent=2))
    print(end.isoformat(),'batch receipt written',flush=True)


if __name__=='__main__':main()

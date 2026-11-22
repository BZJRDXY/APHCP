"""Supplement a slow first-party transfer with a second client, leaving server QC unchanged."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import time


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('cache',type=Path); ap.add_argument('remote'); ap.add_argument('run'); ap.add_argument('--jobs',type=int,default=32); a=ap.parse_args(); a.cache.mkdir(exist_ok=True)
    ssh=['ssh','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes',a.remote]; t=time.monotonic()
    read=lambda n:json.loads(subprocess.check_output([*ssh,'cat',f'{a.run}/inputs/{n}']))
    manifest=read('archive_head.json'); known=lambda: set(subprocess.check_output([*ssh,f"find '{a.run}/raw' -maxdepth 1 -name '*.cif.gz' -printf '%f\\n'"]).decode().splitlines()); local_attempts=a.cache/'local_attempts.jsonl'
    def pull(pin):
        p=pin['pdb_id']; n=f"{p.lower()}-v{pin['revision']}.cif.gz"; final=a.cache/n; tmp=a.cache/(n+'.part'); row={'pdb_id':p,'revision':pin['revision'],'name':n,'bytes':0}
        if not(final.is_file() and final.stat().st_size==pin['content_length']):
            c=subprocess.run(['curl.exe','-sS','--fail','--connect-timeout','10','--max-time','55','--header','Accept-Encoding: identity','--output',str(tmp),pin['source_url']],capture_output=True)
            row['bytes']=tmp.stat().st_size if tmp.is_file() else 0
            try:
                if c.returncode!=0 or tmp.stat().st_size!=pin['content_length']:
                    raise ValueError(c.stderr.decode(errors='replace').strip() or 'incorrect content length')
                with gzip.open(tmp,'rb') as f:
                    while f.read(1<<20):pass
                tmp.replace(final)
            except Exception as exc:
                tmp.unlink(missing_ok=True); return row|{'status':'FAILED','error':str(exc)}
        return row|{'status':'CACHED','sha256':hashlib.sha256(final.read_bytes()).hexdigest()}
    # Reverse order stays mostly disjoint from the server's forward queue; recheck before each wave.
    order=list(reversed(manifest)); step=max(32,a.jobs*2); count=0
    with ThreadPoolExecutor(max_workers=a.jobs) as ex, local_attempts.open('a',encoding='utf8') as f:
        for i in range(0,len(order),step):
            remote=known(); wave=[x for x in order[i:i+step] if f"{x['pdb_id'].lower()}-v{x['revision']}.cif.gz" not in remote]
            if not wave: continue
            vals=[x.result() for x in as_completed([ex.submit(pull,x) for x in wave])]
            remote=known(); push=[x for x in vals if x['status']=='CACHED' and x['name'] not in remote]
            for v in vals: f.write(json.dumps(v)+'\n')
            f.flush()
            if push:
                tar=a.cache/'wave.tar'; receipt=a.cache/'wave.json'; receipt.write_text(json.dumps({'results':vals,'offered':len(push)})+'\n',encoding='utf8')
                with tarfile.open(tar,'w') as pack:
                    for v in push: pack.add(a.cache/v['name'],arcname=v['name'])
                    pack.add(receipt,arcname='transport.json')
                incoming=f'{a.run}/supplement_{i:04d}'; subprocess.run([*ssh,f"mkdir -p '{incoming}'"],check=True)
                subprocess.run(['scp','-q','-o','BatchMode=yes','-o','StrictHostKeyChecking=yes',str(tar),f'{a.remote}:{incoming}/wave.tar'],check=True)
                cmd=f"set -e; cd '{incoming}'; tar -xf wave.tar; for a in *.cif.gz; do if ln -- \"$a\" '{a.run}/raw/'\"$a\" 2>/dev/null; then printf '%s\\n' \"$a\"; fi; done >imported.txt"
                subprocess.run([*ssh,cmd],check=True); count+=len(push)
                print(f'local supplement: offered {count} after {time.monotonic()-t:.1f}s',flush=True)
    print('finished',count,flush=True)

if __name__=='__main__': main()

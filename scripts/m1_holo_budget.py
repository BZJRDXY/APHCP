"""APHCP M1 size planning from archive HEAD and prior measured local receipts."""
from concurrent.futures import ThreadPoolExecutor
from collections import Counter, defaultdict
from datetime import datetime, timezone
import csv, hashlib, json, math, re, sqlite3, struct, sys, time
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError

from m1_holo_config import read_config

run=Path(sys.argv[1]); parent=run.parent; inv=run/'inventory'; out=run/'budget'; out.mkdir(exist_ok=True)
config = read_config(run)
started=datetime.now(timezone.utc); t0=time.monotonic(); now=lambda:datetime.now(timezone.utc).isoformat()
read=lambda p:json.loads(p.read_text(encoding='utf8'))
rows=lambda p:[json.loads(l) for l in p.read_text(encoding='utf8').splitlines()]
def save(path,d):
    p=path.with_suffix('.tmp'); p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf8'); p.replace(path)
def write(name, data, empty_fields=()):
    with (out/name).open('w',encoding='utf8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]) if data else empty_fields); w.writeheader();w.writerows({k:json.dumps(v,ensure_ascii=False,separators=(',',':')) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in data)
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def quant(xs,q):
    ys=sorted(xs); i=(len(ys)-1)*q; n=int(i); return ys[n]+(ys[min(n+1,len(ys)-1)]-ys[n])*(i-n)
def limits(xs):return {k:round(quant(xs,v),3) for k,v in {'p10':.10,'p50':.50,'p90':.90,'max':1}.items()}
def archive(pid,ver):
    p=pid.lower(); assert re.fullmatch(r'[0-9][a-z0-9]{3}',p) and re.fullmatch(r'\d+\.\d+',ver)
    z='pdb_0000'+p; return f'https://files-versioned.rcsb.org/pdb_versioned/data/entries/{p[1:3]}/{z}/{z}_xyz_v{ver.replace(".","-")}.cif.gz'
def method(r):
    methods=r['methods']
    return next((name for key,name in [('X-RAY DIFFRACTION','XRAY'),('ELECTRON MICROSCOPY','EM'),('SOLUTION NMR','NMR')] if key in methods),'OTHER')
def scale(r):
    # deposited_atom_count for NMR is model-size-like (not sum of all coordinate models).
    return max(1,r['atom_count'] or 1)*max(1,r['model_count'] or 1)

e={r['pdb_id']:r for r in rows(inv/'candidate_entries.jsonl')}; lig=defaultdict(list)
for r in rows(inv/'plinder_ligand_rows.jsonl'):lig[r['pdb_id']].append(r)
selection_path = run / 'inputs' / 'm2_selection.json'
nominees = read(selection_path) if selection_path.is_file() else []
nomids = {row[0] for row in nominees}
prior = Path(config.get('budget_reference_run', parent / 'holo-pilot-20261120d'))
corrected = prior / 'revision_recheck'
with (prior/'inputs'/'candidates.csv').open(encoding='utf-8-sig') as f:old_ids={r['pdb_id'].upper() for r in csv.DictReader(f)}
fixed={'1M17','1B86','1HEW','3TCT'}
if selection_path.is_file() and not nominees:
    raise ValueError('The M2 selection file is empty; omit it until nominees have been selected')
if len(nomids) != len(nominees) or not nomids <= set(e) or nomids & (old_ids | fixed):
    raise ValueError('M2 nominees must be unique entries in this inventory and outside the Holo regression set')
if any(not e[p]['rcsb_latest_revision'] for p in nomids):
    raise ValueError('Resolve missing structure revisions for the selected M2 nominees')
references={}; bypid=defaultdict(list); ratios=[]; gzip_sizes=[]; pilot_files=[]
for old_base in [prior/'results',corrected]:
    for rp in sorted((old_base/'entries').glob('*.json')):
        d=read(rp); p=d.get('pdb_id','').upper(); ver=d.get('source_revision'); url=d.get('source_url'); digest=d.get('source_sha256'); src=d.get('source_relative_path')
        if not p or not ver or not url or not digest or not src or not d.get('source_revision_match') or d.get('source_parse_status')!='pass':continue
        path=(old_base/src).resolve()
        if not path.is_relative_to(old_base.resolve()) or not path.is_file() or path.stat().st_size!=d.get('source_bytes') or url!=archive(p,ver):continue
        key=(p,ver); references[key]={'path':str(path),'source_url':url,'bytes':path.stat().st_size,'sha256':digest,'receipt':str(rp)}
for (p,ver),asset in sorted(references.items()):
    gz=Path(asset['path']); n=asset['bytes'];
    with gz.open('rb') as g:g.seek(-4,2); uncompressed=struct.unpack('<I',g.read(4))[0]
    ratios.append(uncompressed/n); gzip_sizes.append(n); bypid[p].append((ver,asset)); pilot_files.append({'pdb_id':p,'revision':ver,'gzip_bytes':n,'gzip_trailer_isize':uncompressed,'receipt':asset['receipt']})
if not ratios:
    raise ValueError(f'No usable coordinate size measurements in {prior}')

db=run.parent.parent; caps=db/'registry'/'aphcp_registry.sqlite3'
with sqlite3.connect(f'file:{caps}?mode=ro',uri=True) as c:
    formal=c.execute('select count(*) from asset_file').fetchone()[0]
    formal_match=c.execute('select ev.pdb_id,ev.revision,a.uri,a.sha256,a.bytes from asset_file a join entry_version ev on ev.id=a.entry_version_id where ev.pdb_id in ('+','.join('?' for _ in e)+')',list(e)).fetchall()
# Current formal registry is empty. A future populated registry needs exact source URL/URI provenance validation.
for p,ver,uri,h,n in formal_match:
    bypid[p.upper()].append((ver,{'path':uri,'sha256':h,'bytes':n,'registered_source_not_verified':True}))

groups=defaultdict(list)
for r in e.values():groups[method(r)].append(r)
samples=set()
for m,group in groups.items():
    ordered=sorted(group,key=lambda r:(scale(r),r['pdb_id']))
    for q in (.05,.18,.35,.55,.75,.92,.99):samples.add(ordered[round(q*(len(ordered)-1))]['pdb_id'])
keys=sorted({(p,e[p]['rcsb_latest_revision']) for p in samples|nomids if e[p]['rcsb_latest_revision']} | {(p,e[p]['plinder_source_revision_if_available']) for p in nomids if e[p]['plinder_revision_status']=='DIFFERENT'})

def head(key):
    p,v=key; url=archive(p,v); t=now(); row={'pdb_id':p,'revision':v,'url':url,'method':'HEAD','started_utc':t,'status':None,'content_length':None}
    try:
        with urlopen(Request(url,method='HEAD',headers={'Accept-Encoding':'identity'}),timeout=14) as res:
            n=res.headers.get('Content-Length'); row.update(status=res.status,content_length=int(n) if n and n.isdigit() and int(n)>0 else None,final_url=res.geturl(),etag=res.headers.get('ETag'),last_modified=res.headers.get('Last-Modified'),content_encoding=res.headers.get('Content-Encoding'))
    except HTTPError as exc:row['status']=exc.code
    except Exception as exc:row['error']=repr(exc)
    row['finished_utc']=now();return row
head_file=out/'archive_head.json'; known={(r['pdb_id'],r['revision']):r for r in read(head_file)} if head_file.is_file() else {}
def cached_head(key):
    previous = known.get(key)
    if previous and previous['url'] == archive(*key) and previous.get('status') in {200, 404}:
        return previous
    return head(key)


with ThreadPoolExecutor(max_workers=5) as pool:
    known.update(zip(keys, list(pool.map(cached_head, keys))))
save(head_file, [known[key] for key in sorted(known)])
stats = defaultdict(list)
for p in sorted(samples):
    r=e[p]; h=known.get((p,r['rcsb_latest_revision'])) or {}; n=h.get('content_length')
    if n and r['atom_count']:stats[method(r)].append(n/scale(r)); stats['ALL'].append(n/scale(r))

def estimate(r,key):
    h=known.get(key,{})
    if h.get('content_length'):return h['content_length'],h['content_length'],h['content_length'],'HEAD_LENGTH'
    seq=stats[method(r)] if len(stats[method(r)])>=3 else stats['ALL']
    if not r['atom_count'] or not seq:
        return *(round(quant(gzip_sizes,q)) for q in (.1,.5,.9)),'PILOT_FILE_SIZE_FALLBACK'
    return *(round(scale(r)*quant(seq,q)) for q in (.1,.5,.9)),('METHOD_ATOM_MODEL_HEAD_RATE' if seq is not stats['ALL'] else 'POOLED_ATOM_MODEL_HEAD_RATE')

def getcase(r,ver,scope):
    p=r['pdb_id']; key=p,ver
    hit=references.get(key); external=[(v,b) for v,b in bypid[p] if b.get('registered_source_not_verified')]; cache='NO_VERIFIED_CACHE'
    if hit:
        valid=sha(Path(hit['path']))==hit['sha256']; cache='VERIFIED_REUSE' if valid else 'REVIEW_CHECKSUM_MISMATCH'
        if not valid:hit=None
    elif external or bypid[p]:cache='REVIEW_UNVERIFIED_OR_OTHER_REVISION'
    lo,mid,hi,basis=estimate(r,key)
    if hit:lo=mid=hi=hit['bytes'];basis='CACHE_BYTES'
    h=known.get(key,{})
    return {'pdb_id':p,'revision':ver,'file_type':'pdbx_mmcif_gzip','scope':scope,'version_frozen_for_m2':False,'source_url':archive(p,ver),'methods':r['methods'],'model_count':r['model_count'],'atom_count':r['atom_count'],'rcsb_assembly_dependencies':r['rcsb_assembly_ids'],'cache_status':cache,'cache_evidence':hit,'head_status':h.get('status'),'size_basis':basis,'compressed_bytes_p10':lo,'compressed_bytes_p50':mid,'compressed_bytes_p90':hi,'download_bytes_p10':0 if hit else lo,'download_bytes_p50':0 if hit else mid,'download_bytes_p90':0 if hit else hi}

plan=[getcase(r,r['rcsb_latest_revision'],'MONTH_CURRENT_REVISION_PLANNING') for r in e.values() if r['rcsb_latest_revision']]
assert len(plan)==len({(x['pdb_id'],x['revision'],x['file_type']) for x in plan})
plan_by={x['pdb_id']:x for x in plan}
optional=[]; selections=[]
for p,category,ccd,reasons in nominees:
    r=e[p]; file=plan_by[p]; optional_case=None
    if r['plinder_revision_status']=='DIFFERENT':
        old=r['plinder_source_revision_if_available']; optional_case=getcase(r,old,'M2_OPTIONAL_INDEX_REVISION_COMPARATOR'); optional.append(optional_case)
    evidence=[]
    for x in lig[p]:
        if ccd and re.search(ccd,x['ligand_ccd_code'] or '',re.I):
            evidence.append({k:x[k] for k in ('system_id','ligand_id','system_biounit_id','ligand_ccd_code','ligand_instance_chains','pli_chain_instances','ligand_covalent_linkages')})
        if len(evidence)>=4:break
    if ccd and lig[p]:assert evidence,(p,ccd)
    hist=[h for h in r['revision_history'] if h['data_content_type']=='Structure model' and f"{h['major_revision']}.{h['minor_revision']}"==r['plinder_source_revision_if_available']]
    selections.append({'pdb_id':p,'test_stratum':category,'title':r['title'],'why_and_test':reasons,'first_release_date':r['first_release_date'],'latest_revision':r['rcsb_latest_revision'],'latest_revision_date':r['rcsb_latest_revision_date'],'plinder_source_revision':r['plinder_source_revision_if_available'],'plinder_source_history_record':hist,'annotation_present':r['plinder_index_presence'],'assembly_ids':r['rcsb_assembly_ids'],'paths':r['discoveries'],'actual_index_examples':evidence,'m2_candidate_file':file,'optional_comparator_file':optional_case})
write('monthly_file_plan.csv', plan, ['pdb_id','revision','file_type'])
write('m2_nominations.csv', selections, ['pdb_id','test_stratum','why_and_test'])
save(out/'m2_nominations.json', selections)
save(out/'m2_optional_comparison_files.json', optional)

systems=sorted((prior/'results'/'systems').glob('*/*/system.cif'))
size_sys=[p.stat().st_size for p in systems]; assert size_sys
indexed_system_count = read(inv/'inventory_summary.json')['counts']['plinder_candidate_systems']
made=sum(1 for x in e.values() if x['plinder_index_presence'])
levels=(('p10',.1),('p50',.5),('p90',.9))

def scenario(files,systems_count):
    vals={}
    for tag,q in levels:
        k='compressed_bytes_'+tag; z='download_bytes_'+tag; ratio=quant(ratios,q); sd=quant(size_sys,q)
        stored=sum(x[k] for x in files); network=sum(x[z] for x in files); uncompressed=round(stored*ratio); derivative=round(systems_count*sd); scratch=round(max((x[k] for x in files),default=0)*ratio+2*sd)
        vals[tag]={'compressed_raw_retained':stored,'cache_reused_bytes':stored-network,'network':network,'raw_uncompressed_if_retained':uncompressed,'derived_cif':derivative,'serial_scratch_peak':scratch,'incremental_peak_excluding_maps_annotations':network+derivative+scratch}
    return {'file_count':len(files),'derived_systems_assumed':systems_count,'bands':vals}

n_model=Counter(method(r) for r in e.values()); sizes={m:scenario([x for x in plan if method(e[x['pdb_id']])==m],0) for m in n_model}
measurements={'pilot_raw_files_n':len(pilot_files),'pilot_raw_gzip_bytes_total':sum(gzip_sizes),'pilot_gzip_trailer_uncompressed_bytes_total':sum(x['gzip_trailer_isize'] for x in pilot_files),'pilot_raw_to_gzip_ratio':limits(ratios),'pilot_system_cif_files_n':len(size_sys),'pilot_system_cif_total_bytes':sum(size_sys),'pilot_per_system_cif_bytes':limits(size_sys),'archive_head_count':len(keys),'archive_head_status':dict(Counter(str(known[k]['status']) for k in keys)),'archive_head_numeric_length_count':sum(bool(known[k]['content_length']) for k in keys),'calibration_per_atom_model_unit':{k:{'n':len(v),'bytes_per_unit':limits(v)} for k,v in stats.items() if v}}
summary={'started_utc':started.isoformat(),'finished_utc':now(),'elapsed_seconds':round(time.monotonic()-t0,3),'formal_asset_file_count':formal,'formal_asset_file_id_matches':len(formal_match),'pilot_key_count':len(references),'pilot_entry_id_overlap':sorted(set(e)&{p for p,v in references}),'cache_status':dict(Counter(x['cache_status'] for x in plan)),'size_basis_counts':dict(Counter(x['size_basis'] for x in plan)),'files_without_current_revision':len(e)-len(plan),'file_planning_key':['pdb_id','revision','file_type'],'measurement':measurements,'month_raw_by_method':sizes,'month_all_indexed_systems_scenario':scenario(plan,indexed_system_count),'month_unannotated_native_derivatives_not_counted':len(e)-made,'m2_one_readback_per_nominee_scenario':scenario([plan_by[p] for p in sorted(nomids)],len(nomids)),'m2_full_selection_plus_optional_comparator_scenario':scenario([plan_by[p] for p in sorted(nomids)]+optional,len(nomids)+len(optional)),'m2_plus_head_verified_comparators_scenario':scenario([plan_by[p] for p in sorted(nomids)]+[o for o in optional if o['head_status']==200],len(nomids)+sum(o['head_status']==200 for o in optional)),'m2_no_annotation_count':sum(not e[p]['plinder_index_presence'] for p in nomids),'m2_gz_head_statuses':dict(Counter(str(plan_by[p]['head_status']) for p in nomids)),'official_assembly_files_planned':0,'maps_or_structure_factors_planned':0,'notes':['Planning scenarios only; revisions not frozen; no coordinates were transferred.','One process scratch assumptions; estimates are empirical scenarios, not guaranteed bounds. No indexed system inferred for unannotated entries.','Official assembly files, experimental maps and structure factors are excluded; assembly generation dependencies are retained.']}
summary['m2_nomination_count'] = len(nominees)
save(out/'pilot_size_measurements.json',pilot_files);save(out/'budget_summary.json',summary)
print(json.dumps(summary,ensure_ascii=False),flush=True)

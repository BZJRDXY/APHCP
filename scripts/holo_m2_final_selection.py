"""Select/preflight a balanced final M2 batch using the existing M1 snapshots."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from csv import DictReader
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import re
import time
from urllib.request import Request, urlopen

from holo_candidate_batch import atomic_json

TARGET = {'INDEX_MATCH':1200, 'NO_LIGAND_ANNOTATION':600, 'INDEX_VERSION_MISMATCH':200}
METHOD_TARGET = {'X-RAY DIFFRACTION':.62, 'ELECTRON MICROSCOPY':.22, 'SOLUTION NMR':.14}
COMMON = {'GOL','EDO','SO4','PO4','DMS','PEG','MPD','CL','NA','K','TRS','ACT','FMT','ACE','HEZ','EPE','BME','PG4','PGE'}
LIPIDS = {'OLA','PLM','MYR','STE','CLR','CHL','PEF','PEE','POV','BOG','LMT','LMU','LDA','DMU','PC1','P1O','POV'}
COFACTORS = {'HEM','HEC','FAD','FMN','NAD','NAP','NDP','PLP','SAM','SAH','SF4','FES','COA','B12','MG','ZN','MN'}
# Restrict sensitive structural targets outside the generic dataset engineering task.
SENSITIVE = re.compile(r'\b(virus|viral|virion|viridae|viroporin|toxin|toxins|toxic|botulinum|tetanus|anthrax|ricin|abrin|pertussis|saxitoxin|diphtheria|cholera|orthopox|vaccinia|variola|ebola|marburg|filovirus|henipavirus|nipah|hendra|coronavirus|sars|mers|influenza|hiv|simian immunodeficiency|hepatitis|enterovirus|rhinovirus|poliovirus|adenovirus|rotavirus|flavivirus|dengue|zika|chikungunya|norovirus|herpes|cytomegalovirus|rabies|bacteriophage|phage|capsid|hemagglutinin|neuraminidase|spike glycoprotein)\b',re.I)
FEATURE_GOALS = {'polymer_nucleic':.065,'polymer_peptide':.09,'macrocycle_hint':.025,'branched_entity':.04,
                 'assembly_copy':.15,'covalent_link':.10,'multi_protein_index':.12,'repeated_ccd_instance':.16,
                 'lipid':.05,'cofactor':.09,'additive_only_hint':.075,'multiple_assemblies':.12,'altloc_title_hint':.025,
                 'large':.08,'multi_model':.14,'membrane_hint':.12,'protein_single_chain':.20,'protein_multichain':.20,
                 'metal_coordination_hint':.045,'explicit_pair_hint':.045,'protein_multi_entity':.12,'native_ordinary':.19}


def get(path):
    return json.loads(path.read_text(encoding='utf-8-sig'))


def stamp():
    return datetime.now(timezone.utc).isoformat()


def head(pin):
    out = {**pin,'checked_utc':stamp()}
    for i in range(2):
        try:
            q = Request(pin['source_url'], method='HEAD', headers={'Accept-Encoding':'identity','User-Agent':'APHCP-M2-preflight/1'})
            with urlopen(q,timeout=12) as response:
                n=int(response.headers.get('Content-Length',0)); out.update(status=response.status,content_length=n,final_url=response.geturl())
                if response.status==200 and n>0 and response.geturl()==pin['source_url']:
                    return out | {'availability':'VERIFIED'}
        except Exception as e:
            out.update(status=getattr(e,'code',None),error=f'{type(e).__name__}: {e}')
    return out | {'availability':'DEFERRED_BEFORE_FREEZE'}


def compact(c):
    if c.get('plinder_index_presence'):
        group={'MATCH':'INDEX_MATCH','DIFFERENT':'INDEX_VERSION_MISMATCH'}.get(c['plinder_revision_status'])
        if not group and c['rcsb_latest_revision'] and c['plinder_source_revision_if_available'] and c['rcsb_latest_revision']!=c['plinder_source_revision_if_available']:
            group='INDEX_VERSION_MISMATCH'
    else:
        group='NO_LIGAND_ANNOTATION'
    title=c['title'] or ''; date=c['first_release_date'] or ''; rev=c['rcsb_latest_revision'] or ''
    if not group or not date or SENSITIVE.search(title) or not re.fullmatch(r'\d+\.\d+',rev):
        return
    # Preserve a workload spread while avoiding a handful of giant files dominating the blind pilot.
    atoms=c.get('atom_count') or 0; models=c.get('model_count') or 1
    if not 300<=atoms<=400000 or models>45 or not 1<=len(c['rcsb_assembly_ids'])<=18:
        return
    codes=set(c['nonpolymer_ccds']); ln=title.lower(); partner=c.get('partner_evidence',[])
    features=set(); discover=c['source_entity_ids']; m=';'.join(c['methods'])
    if discover['branched']: features.add('branched_entity')
    if len(c['rcsb_assembly_ids'])>1: features.add('multiple_assemblies')
    if models>1: features.add('multi_model')
    if atoms>=43000: features.add('large')
    if codes & LIPIDS: features.add('lipid')
    if codes & COFACTORS: features.add('cofactor')
    if codes and codes<=COMMON: features.add('additive_only_hint')
    if re.search(r'cyclic (?:peptide|inhibitor)|cyclopeptide|macrocycl',ln): features.add('macrocycle_hint')
    if re.search(r'alternate (?:conformation|state)|alternative (?:conformation|state)|dual conformation|multiple conformation',ln): features.add('altloc_title_hint')
    if re.search(r'\b(?:membrane|transmembrane|transporter|channel|gpcr)\b',ln): features.add('membrane_hint')
    if re.search(r'\b(?:dna|rna|nucleosome|duplex|oligonucleotide)\b',ln): features.add('polymer_nucleic')
    if re.search(r'\b(?:duplex|double.strand|base.pair|hairpin)\b',ln): features.add('explicit_pair_hint')
    if re.search(r'\b(?:peptide|peptidic)\b',ln): features.add('polymer_peptide')
    if re.search(r'\b(?:covalent|conjugate|cross.link|adduct)\b',ln): features.add('covalent_title_hint')
    if (c['protein_entity_count'] or 0)>1: features.add('protein_multi_entity')
    if (c['plinder_candidate_system_count'] or 0)>1: features.add('multiple_systems')
    if codes-COMMON and 'X-RAY DIFFRACTION' in m and not features&{'polymer_nucleic','polymer_peptide','branched_entity','large'}: features.add('native_ordinary')
    if group=='NO_LIGAND_ANNOTATION' and not codes and not partner and not discover['branched'] and not features&{'polymer_nucleic','polymer_peptide','multi_model','covalent_title_hint','macrocycle_hint'}: features.add('no_partner_metadata_boundary')
    if not codes and not features and not partner and group=='NO_LIGAND_ANNOTATION': return
    return {'pdb_id':c['pdb_id'],'group':group,'rcsb_revision':rev,'plinder_revision':c['plinder_source_revision_if_available'],
            'title':title,'date':date[:10],'method':m,'year_bin':int(date[:4])//4*4,'atom_count':atoms,'models':models,
            'assemblies':c['rcsb_assembly_ids'],'protein_entity_count':c['protein_entity_count'],'ccds':sorted(codes),
            'partner_evidence':partner,'pmid':c.get('pubmed_id'),'features':sorted(features),'protein':{},'index_evidence':{},
            'priority':int(sha256(c['pdb_id'].encode()).hexdigest()[:10],16)/16**10}


def shortlist(full, excluded):
    buckets=defaultdict(list)
    with (full/'candidate_entries.jsonl').open() as f:
        for line in f:
            c=json.loads(line)
            if c['pdb_id'] in excluded: continue
            p=compact(c)
            if p is not None: buckets[p['group'],p['year_bin'],p['method']].append(p)
    result={}
    for (g,_,m),bucket in buckets.items():
        # Retain sparse dates/methods and keep an enriched reserve for rare documented cases.
        bucket.sort(key=lambda x:x['priority']); cap=650 if g=='INDEX_MATCH' else 520 if g=='NO_LIGAND_ANNOTATION' else 450
        rare={'macrocycle_hint','explicit_pair_hint','altloc_title_hint','covalent_title_hint'}; kept=bucket[:cap]+[x for x in bucket[cap:] if rare.intersection(x['features'])][:170]
        for x in kept: result[x['pdb_id']]=x
    return result


def enrich_windows(run_root, pool):
    for w in sorted(run_root.glob('holo-m1-*')):
        path=w/'inventory/rcsb_entities.csv'
        if not path.is_file(): continue
        with path.open() as f:
            for r in DictReader(f):
                p=pool.get(r['pdb_id'])
                if not p or r['entity_kind']!='polymer': continue
                p['window']=w.name; tags=set(p['features']); typ=r['polymer_type'].lower(); name=r['description']; length=int(r['polymer_length'] or 0)
                if SENSITIVE.search(name): p['excluded_by_entity']=True
                if 'dna' in typ or 'rna' in typ: tags.add('polymer_nucleic')
                if 'protein' in typ or 'peptide' in typ:
                    n=len(json.loads(r['label_asym_ids'] or '[]'))
                    tags.add('protein_multichain' if n>1 else 'protein_single_chain')
                    if length<=45: tags.add('polymer_peptide')
                    if length>p['protein'].get('length',0):
                        evidence=json.loads(r['family_evidence'] or '{}'); uni=evidence.get('uniprot_ids',[]); pf=sorted({a['annotation_id'] for a in evidence.get('family_annotations',[]) if a['type'].lower()=='pfam'});
                        p['protein']={'entity':r['entity_id'],'length':length,'description':name,'chains':n,'uniprot':uni,'pfam':pf}
                p['features']=sorted(tags)


def enrich_index(full, pool):
    hits={p for p,row in pool.items() if row['group']!='NO_LIGAND_ANNOTATION'}; per=defaultdict(lambda: defaultdict(set))
    with (full/'plinder_ligand_rows.jsonl').open('rb') as f:
        for line in f:
            # The first column is entry_pdb_id in the verified M1 rollup.
            m=re.match(br'\{"entry_pdb_id":\s*"([a-z0-9]{4})"',line)
            if m and m[1].decode().upper() not in hits: continue
            row=json.loads(line); p=row['entry_pdb_id'].upper()
            if p not in hits: continue
            rec=pool[p]; tags=set(rec['features']); ev=rec['index_evidence']; chains=row.get('system_protein_chains_asym_id') or []
            if len(chains)>1: tags.add('multi_protein_index'); ev.setdefault('multi_protein_chains',chains)
            if any(str(ch).split('.',1)[0].isdigit() and int(str(ch).split('.',1)[0])>1 for ch in chains+(row.get('ligand_instance_chains') or [])):
                tags.add('assembly_copy'); ev.setdefault('copy_instance_chains',chains+(row.get('ligand_instance_chains') or []))
            for k,tag in [('ligand_is_covalent','covalent_link'),('ligand_is_cofactor','cofactor'),('ligand_is_oligopeptide','polymer_peptide'),('ligand_is_oligosaccharide','indexed_glycan'),('ligand_is_oligonucleotide','polymer_nucleic')]:
                if row.get(k): tags.add(tag); ev.setdefault(k,row.get('ligand_id'))
            if row.get('ligand_covalent_linkages'): ev.setdefault('ligand_covalent_linkages',row['ligand_covalent_linkages'])
            if row.get('ligand_is_ion') and row.get('ligand_num_interactions',0): tags.add('metal_coordination_hint')
            c=row.get('ligand_ccd_code'); instance=row.get('ligand_instance_chain')
            if c and instance and c not in COMMON:
                d=per[p][(row['system_biounit_id'],c)]; d.add(instance)
                if len(d)>1: tags.add('repeated_ccd_instance'); ev['repeated_ccd_instance']={'ccd':c,'instances':sorted(d)}
            if not rec['protein'].get('pfam') and row.get('system_pocket_Pfam'):
                ev.setdefault('pocket_pfam',row['system_pocket_Pfam'])
            rec['features']=sorted(tags)


def ranking(rows, group):
    left={x['pdb_id']:x for x in rows}; order=[]; tags=Counter(); method=Counter(); series=Counter(); year=Counter(); pm=Counter(); family=Counter(); n=TARGET[group]
    def key(x):
        p=x['protein']; u=p.get('uniprot',[]); base=re.sub(r'[^a-z]+',' ',p.get('description','').lower()).strip()
        fam=(p.get('pfam') or [x['index_evidence'].get('pocket_pfam')])[0]; return (u[0] if u else base or x['pdb_id'],fam)
    k={p:key(x) for p,x in left.items()}; t={p:{f for f in x['features'] if f in FEATURE_GOALS} for p,x in left.items()}
    # Rank reserves adaptively; accepted rows are then selected in this order after HEAD and byte checks.
    while left and len(order)<n+max(65,n//7):
        def rank(p):
            x=left[p]; ss,ff=k[p]; yr=x['year_bin']; m=x['method']; mid=x['pmid']; n0=len(order)+1
            importance=sum(1.4*max(0.04,1-tags[feat]/max(1,n*FEATURE_GOALS[feat])) for feat in t[p]);
            if {'assembly_copy','covalent_link'}<=t[p] or {'assembly_copy','polymer_nucleic'}<=t[p]: importance+=.65
            score=importance + (1.25 if method[m]<n*METHOD_TARGET.get(m,.02) else -.75) + max(-1,1.2-year[yr]/max(12,n*.055)) - 3*series[ss] - (.55*family[ff]**.6 if ff else 0) - (1.7*pm[mid] if mid and mid>0 else 0) + x['priority']*.4
            return score
        p=max(left,key=rank); x=left.pop(p); ss,ff=k[p]; tags.update(t[p]); method[x['method']]+=1; series[ss]+=1; year[x['year_bin']]+=1; family[ff]+=1; pm[x['pmid']]+=1; order.append(x)
    return order


def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument('run',type=Path); ap.add_argument('--root',type=Path,required=True); ap.add_argument('--heads',type=int,default=32); args=ap.parse_args(); run=args.run.resolve(); inp=run/'inputs'; inp.mkdir(exist_ok=True); t=time.monotonic()
    full=args.root/'batches/holo-m1-full-20261121a/planning/full_inventory'; hist=get(inp/'history_excluded.json'); excluded=set(hist['ids']); shortlist_cache=run/'selection_shortlists.json'
    if shortlist_cache.is_file(): order=get(shortlist_cache)
    else:
        pool=shortlist(full,excluded); print('M1 shortlisted',len(pool),flush=True)
        enrich_windows(args.root/'runs',pool); print('M1 entities enriched',len(pool),flush=True)
        pool={p:x for p,x in pool.items() if not x.get('excluded_by_entity') and x.get('window')}
        enrich_index(full,pool); print('M1 index evidence enriched',len(pool),flush=True)
        order={g:ranking([x for x in pool.values() if x['group']==g],g) for g in TARGET}; chosen={x['pdb_id']:x for rs in order.values() for x in rs};
        with (full/'monthly_file_plan.csv').open() as f:
            for a in DictReader(f):
                p=chosen.get(a['pdb_id'])
                if p and p['rcsb_revision']==a['revision']:
                    p['file_plan']={k:a[k] for k in ('pdb_id','revision','source_url','file_type','cache_status','cache_evidence','compressed_bytes_p10','compressed_bytes_p50','compressed_bytes_p90')}
        for g,n in TARGET.items():
            order[g]=[x for x in order[g] if x.get('file_plan')]; assert len(order[g])>=n
        atomic_json(shortlist_cache,order); print('ranked reserves', {g:len(x) for g,x in order.items()}, flush=True)
    if (inp/'blind.json').is_file():
        print('Frozen inputs already exist; selection unchanged',flush=True); return
    heads_file=run/'archive_preflight_attempts.json'; cache=get(heads_file) if heads_file.is_file() else {}; fields=('pdb_id','revision','source_url','file_type')
    accepted={g:[] for g in TARGET}; manifest=[]; evidence=[]; receipts=[]; cap=18_000_000
    with ThreadPoolExecutor(max_workers=args.heads) as pool:
        pending={g:iter(v) for g,v in order.items()}
        while any(len(accepted[g])<n for g,n in TARGET.items()):
            todo=[]
            for g,n in TARGET.items():
                needed=n-len(accepted[g]); todo.extend((g,x) for x in (next(pending[g],None) for _ in range(max(0,min(needed,96)))) if x)
            if not todo: raise RuntimeError(f'No reserve left: {[(g,len(x)) for g,x in accepted.items()]}')
            tasks=[]
            for g,x in todo:
                pin={k:x['file_plan'][k] for k in fields}; p=x['pdb_id']; old=cache.get(p)
                tasks.append((g,x,pin,None if old and old['availability']=='VERIFIED' else pool.submit(head,pin)))
            for g,x,pin,fut in tasks:
                p=x['pdb_id']; proof=fut.result() if fut else cache[p]; cache[p]=proof
                if proof['availability']=='VERIFIED' and proof['content_length']<=cap:
                    accepted[g].append(p); manifest.append(pin); receipts.append(proof); x['archive_head_bytes']=proof['content_length'];
                    x['selection_purpose']='Check '+', '.join(x['features'] or ['baseline experimental protein and fixed-source route'])+'; M1 and index evidence are candidates only.'; evidence.append(x)
            atomic_json(heads_file,cache); print('preflight', {g:len(x) for g,x in accepted.items()},'requests/results tracked',len(cache),flush=True)
    ids={x['pdb_id'] for x in manifest}; assert len(manifest)==len(ids)==2000 and not ids&excluded
    atomic_json(inp/'archive_head.json',receipts); atomic_json(run/'selection_evidence.json',evidence); atomic_json(inp/'blind.json',manifest)
    summary={'frozen_utc':stamp(),'counts':{g:len(x) for g,x in accepted.items()},'n':len(ids),'excluded_unique':len(excluded),'history_overlap':0,
        'compressed_exact_bytes':sum(a['content_length'] for a in receipts),'max_object_bytes':cap,'estimated_download_bytes':sum(a['content_length'] for a in receipts),
        'authorized_cache_reused_before_freeze':0,'source_range':[min(x['date'] for x in evidence),max(x['date'] for x in evidence)],'methods':dict(Counter(x['method'] for x in evidence)),
        'features':dict(Counter(f for x in evidence for f in x['features'])),'windows':dict(Counter(x['window'] for x in evidence)),
        'revision_mismatch':dict(Counter('major' if x['rcsb_revision'].split('.')[0]!=x['plinder_revision'].split('.')[0] else 'minor_only' for x in evidence if x['group']=='INDEX_VERSION_MISMATCH')),
        'preflight_excluded':{p:{'status':v.get('status'),'bytes':v.get('content_length'),'error':v.get('error')} for p,v in cache.items() if p not in ids},'elapsed_s':round(time.monotonic()-t,3)}
    atomic_json(run/'selection_summary.json',summary); print('FROZEN',json.dumps({k:v for k,v in summary.items() if k!='preflight_excluded'}),flush=True)

if __name__=='__main__': main()

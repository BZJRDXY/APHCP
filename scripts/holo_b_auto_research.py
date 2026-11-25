"""Small, frozen-input Holo B evidence-gathering research pilot.

All decisions are provisional research observations. Source catalogues are read-only.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import html
import json
from pathlib import Path
import random
import re
import sqlite3
import time
import urllib.parse
import urllib.request

import gemmi


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_db(path):
    db = sqlite3.connect(path.resolve().as_uri()+'?mode=ro', uri=True)
    db.row_factory = sqlite3.Row
    return db


def read_cif(path):
    path = Path(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8') as f:
        return gemmi.cif.read_string(f.read()).sole_block()


def rows(block, category):
    data = block.get_mmcif_category(category)
    if not data:
        return []
    n = len(next(iter(data.values())))
    return [{key:values[i] for key,values in data.items()} for i in range(n)]


def scalar(block, key):
    x = block.find_value(key)
    return None if x in (None, '.', '?', '') else gemmi.cif.as_string(x)


def publication(block):
    found = [r for r in rows(block,'_citation.') if r.get('id') == 'primary']
    r = found[0] if len(found) == 1 else {}
    return {'title':r.get('title'), 'doi':r.get('pdbx_database_id_DOI'),
            'pmid':r.get('pdbx_database_id_PubMed')}


def publication_keys(p):
    keys = set()
    for field in ('doi','pmid','title'):
        v = p.get(field)
        if v not in (None,'','-1','?','.'):
            keys.add((field,' '.join(str(v).lower().split())))
    return keys


def protein_keys(block):
    keys = set()
    for r in rows(block, '_entity_poly.'):
        if str(r.get('type') or '').lower().startswith('polypeptide'):
            seq = re.sub(r'\s+','',str(r.get('pdbx_seq_one_letter_code_can') or ''))
            if seq and seq not in ('?','.'):
                keys.add(('sequence',seq))
    for r in rows(block, '_struct_ref.'):
        if str(r.get('db_name') or '').upper() in ('UNP','UNIPROT'):
            value = r.get('pdbx_db_accession') or r.get('db_code')
            if value and value not in ('?','.'):
                keys.add(('uniprot',str(value).upper()))
    return keys


def aliases(block, code):
    result = []
    for r in rows(block,'_chem_comp.'):
        if r.get('id') == code:
            if r.get('name'):
                result.append(('chem_comp.name',str(r['name'])))
            for x in str(r.get('pdbx_synonyms') or '').split(';'):
                if x.strip():result.append(('chem_comp.pdbx_synonyms',x.strip()))
    for r in rows(block,'_entity.'):
        if r.get('pdbx_description'):
            # The target entity is matched by caller before adding this description.
            result.append(('entity:'+str(r.get('id')),str(r['pdbx_description'])))
    return result


def norm(text):
    return ' '.join(re.findall(r'[a-z0-9]+',str(text or '').lower()))


def title_evidence(block, target):
    title = scalar(block,'_struct.title') or ''
    title_norm = norm(title)
    code, entity = str(target[6]), str(target[3])
    names = []
    for field,name in aliases(block,code):
        if field.startswith('entity:') and field != 'entity:'+entity:
            continue
        s = norm(name)
        if len(s) >= 5 and re.search(r'(?<![a-z0-9])'+re.escape(s)+r'(?![a-z0-9])',title_norm):
            names.append({'field':field,'matched':name})
    # A CCD token can be useful, but short everyday English words are not evidence.
    if len(code) >= 3 and code.upper() not in {'AND','THE','FOR','NEW','ONE','ALL','NOT','WITH','HAS','CAN','USE','WAS'}:
        pieces = re.findall(r'[A-Za-z]+|[0-9]+',code)
        pattern = r'(?<![A-Za-z0-9])'+r'[\s-]*'.join(map(re.escape,pieces))+r'(?![A-Za-z0-9])'
        if re.search(pattern,title,re.I):
            names.append({'field':'ccd_code','matched':code})
    cue = bool(re.search(r'\b(complex(?:ed)?|bound|inhibitor|substrate|agonist|antagonist)\b',title,re.I))
    negative = bool(re.search(r'\b(without|unbound|unliganded|replaced|instead|rather than|not)\b',title,re.I))
    same_author = [r for r in rows(block,'_pdbx_nonpoly_scheme.')
                   if str(r.get('entity_id')) == entity and str(r.get('pdb_strand_id')) == str(target[4])]
    # Different copies in the same author chain cannot inherit a title assertion automatically.
    unique = len({(r.get('asym_id'),r.get('pdb_seq_num'),r.get('mon_id')) for r in same_author}) == 1
    return {'field':'_struct.title','text':title,'matches':names,'binding_cue':cue,
            'negative_cue':negative,'unique_same_author_entity':unique}


def connections(block, identity):
    targets = identity.get('target_residue_instances') or []
    proteins = identity.get('binding_instances') or []
    found = []
    for conn in rows(block, '_struct_conn.'):
        if str(conn.get('conn_type_id') or '').lower() != 'covale':
            continue
        for target in targets:
            for left,right in (('ptnr1','ptnr2'),('ptnr2','ptnr1')):
                def get(side,field):return conn.get(side+'_'+field)
                if not (str(get(left,'label_asym_id')) == str(target[1])
                    and str(get(left,'label_comp_id')) == str(target[6])
                    and str(get(left,'auth_asym_id')) == str(target[4])):
                    continue
                ins = conn.get('pdbx_'+left+'_PDB_ins_code')
                authseq = str(get(left,'auth_seq_id')) + (str(ins) if ins not in (None,False,'.','?','') else '')
                if authseq != str(target[5]):continue
                links = [p for p in proteins if str(p[1]) == str(get(right,'label_asym_id')) and str(p[2]) == str(target[2])]
                if not links:continue
                sym1 = conn.get(left+'_symmetry');sym2 = conn.get(right+'_symmetry')
                same_source = sym1 == '1_555' and sym2 == '1_555'
                found.append({'id':conn.get('id'),'type':conn.get('conn_type_id'),'pdbx_role':conn.get('pdbx_role'),
                    'target_instance':target[0], 'target_ccd':target[6], 'target_atom':get(left,'label_atom_id'),
                    'target_auth_seq':target[5], 'protein_instances':[x[0] for x in links],
                    'protein_comp':get(right,'label_comp_id'),'protein_auth_seq':get(right,'auth_seq_id'),
                    'protein_atom':get(right,'label_atom_id'),'same_source_symmetry':same_source,
                    'source_symmetry':[sym1,sym2]})
    return found


def prep_evidence(block, target):
    # Only a directly written cryoprotection statement is a preparation fact.
    if target[6] != 'GOL':return []
    out = []
    for field in ('_exptl_crystal_grow.pdbx_details','_exptl_crystal_grow.details'):
        for raw in block.find_values(field):
            value = gemmi.cif.as_string(raw)
            if re.search(r'\bglycerol\b',value,re.I) and re.search(r'\bcryo[ -]?(protect|protected|protection)',value,re.I):
                if not re.search(r'\b(without|no glycerol|not|instead)\b',value,re.I):
                    out.append({'field':field,'text':value})
    return out


def article_evidence(article, block, target, pdb):
    """Use an abstract only if one sentence names this PDB and the target."""
    text=article.get('_abstract') or ''
    code=str(target[6]);entity=str(target[3]);candidates=[]
    for field,name in aliases(block,code):
        if field.startswith('entity:') and field!='entity:'+entity:continue
        n=norm(name)
        if len(n)>=5:candidates.append((field,name,n))
    if len(code)>=3 and code.upper() not in {'AND','THE','FOR','NEW','ONE','ALL','NOT','WITH','HAS','CAN','USE','WAS'}:
        candidates.append(('ccd_code',code,norm(code)))
    for sentence in re.split(r'(?<=[.!?])\s+',text):
        if not re.search(r'(?<![A-Za-z0-9])'+re.escape(pdb)+r'(?![A-Za-z0-9])',sentence,re.I):continue
        if not re.search(r'\b(complex(?:ed)?|bound|binds?|binding|inhibitor|substrate|agonist|antagonist)\b',sentence,re.I):continue
        if re.search(r'\b(without|unbound|unliganded|replaced|instead|rather than|not)\b',sentence,re.I):continue
        normalized=norm(sentence)
        for field,name,needle in candidates:
            if re.search(r'(?<![a-z0-9])'+re.escape(needle)+r'(?![a-z0-9])',normalized):
                return {'field':'Europe PMC abstract','pmid':article.get('pmid'),'pdb_id':pdb,
                        'matched_field':field,'matched':name,'excerpt':' '.join(sentence.split()[:20])}
    return None


def qc_pass(row):
    old = json.loads(row['old_qc']) if row['old_qc'] else None
    new = json.loads(row['new_qc']) if row['new_qc'] else None
    def good(d,fields):return isinstance(d,dict) and all(d.get(k)=='PASS' for k in fields)
    return (good(old,('coordinate_identity','experimental_fidelity','topology','partner_scope_qc')) and
            good(new,('coordinate_identity','experimental_fidelity','topology','partner_scope')))


def assess(row, block, article):
    obs=json.loads(row['document']); identity=obs['identity']; targets=identity.get('target_residue_instances') or []
    if len(targets)!=1 or len(targets[0])<7:
        return {'decision':'REVIEW','reason':'target_not_one_complete_component','identity':identity}
    t=targets[0]; title=title_evidence(block,t); conn=connections(block,identity); prep=prep_evidence(block,t)
    paper=article_evidence(article,block,t,identity['pdb_id'])
    evidence={'system_id':row['system_id'],'identity':identity,
        'source_mmcif':obs['files']['source_mmcif'],'old_qc':json.loads(row['old_qc']) if row['old_qc'] else None,
        'm3_qc':json.loads(row['new_qc']) if row['new_qc'] else None,'title_evidence':title,
        'connections':conn,'preparation_records':prep,'publication':publication(block),
        'article_lookup':{k:v for k,v in article.items() if k!='_abstract'},'article_evidence':paper,
        'plinder_match':row['match_level'],'plinder_interpretation':row['annotation_interpretation'],
        'old_conflicts':json.loads(row['conflicts']) if row['conflicts'] else None,
        'prior_source_evidence':json.loads(row['source_evidence']) if row['source_evidence'] else None}
    evidence['decision']='REVIEW'; evidence['reason']='insufficient_localized_role_evidence'
    if not qc_pass(row):
        evidence['reason']='quality_or_partner_scope_not_pass_or_not_recorded';return evidence
    prior=evidence['prior_source_evidence'] or {}
    if evidence['old_conflicts'] or 'competing' in str(prior.get('reason') or '').lower():
        evidence['reason']='historical_conflict_or_competing_interpretations';return evidence
    direct=[x for x in conn if x['same_source_symmetry']]
    glyco=[x for x in direct if t[6]=='NAG' and x['target_atom']=='C1' and x['protein_comp']=='ASN' and x['protein_atom']=='ND2']
    title_support=bool(title['matches'] and title['binding_cue'] and not title['negative_cue'] and title['unique_same_author_entity'])
    if glyco:
        if title_support:
            evidence['reason']='glycosylation_link_and_explicit_title_target_both_present';return evidence
        evidence['decision']='PROTEIN_N_GLYCOSYLATION_LINK_SUPPORTED'
        evidence['reason']='fixed_source_records_target_NAG_C1_to_binding_ASN_ND2_same_assembly_operator'
    elif prep:
        # This is a preparation fact, never a claim that this particular pose has no other relevance.
        if title_support:
            evidence['reason']='preparation_record_and_explicit_title_target_both_present';return evidence
        evidence['decision']='PREPARATION_USE_DOCUMENTED'
        evidence['reason']='fixed_source_explicitly_records_glycerol_cryoprotection_not_exclusive_site_role'
    elif title_support:
        evidence['decision']='DEPOSITOR_NAMED_BINDING_PARTNER_SUPPORTED'
        evidence['reason']='fixed_entry_title_names_unique_target_entity_and_binding_with_recorded_contact'
    elif paper and title['unique_same_author_entity']:
        evidence['decision']='PAPER_NAMED_BINDING_PARTNER_SUPPORTED'
        evidence['reason']='primary_abstract_same_sentence_names_PDB_target_and_binding_with_recorded_contact'
    return evidence


def article_batch(pmids, timeout=35):
    ids=sorted({str(x) for x in pmids if x not in (None,'','-1','?','.')})
    if not ids:return {},{'state':'NO_PMIDS','ids':0}
    q='EXT_ID:('+' OR '.join(ids)+') AND SRC:MED'
    params={'query':q,'format':'json','resultType':'core','pageSize':100}
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode(params)
    started=now()
    try:
        request=urllib.request.Request(url,headers={'User-Agent':'APHCP-research-pilot/1.0'})
        with urllib.request.urlopen(request,timeout=timeout) as response:
            body=response.read();data=json.loads(body)
        found={str(v.get('pmid')):{'state':'RETRIEVED','pmid':v.get('pmid'),'doi':v.get('doi'),
               'title':v.get('title'),'abstract_available':bool(v.get('abstractText')),
               '_abstract':html.unescape(re.sub(r'<[^>]+>',' ',v.get('abstractText') or '')),
               'source':'Europe PMC MED','query_utc':started} for v in data.get('resultList',{}).get('result',[])}
        return found,{'state':'SUCCESS','query_utc':started,'ids':len(ids),'returned':len(found),'bytes':len(body),'endpoint':'Europe PMC REST'}
    except (OSError,ValueError,TimeoutError) as exc:
        return {},{'state':'SOURCE_TECHNICAL_FAILURE','query_utc':started,'ids':len(ids),'error_type':type(exc).__name__,'error':str(exc)[:400]}


def load_record(db,sid):
    row=db.execute('select result.*,observation.document from result join observation using(system_id) where system_id=?',(sid,)).fetchone()
    if row is None:raise ValueError('system not in first-round snapshot: '+sid)
    return row


def run(args):
    manifest=json.loads(args.manifest.read_text(encoding='utf-8'))
    if manifest['code_sha256'] != sha(__file__):raise ValueError('code differs from frozen manifest')
    if args.output.exists():raise ValueError('output already exists')
    started=now();clock=time.monotonic();db=read_db(args.database)
    items=[];pubs=[];fetch_errors={}
    for spec in manifest['systems']:
        sid=spec['system_id']
        try:
            row=load_record(db,sid);obs=json.loads(row['document']);path=obs['files']['source_mmcif']['path']
            if obs['files']['source_mmcif']['sha256'] != spec['source_sha256']:
                raise ValueError('source identity differs from frozen manifest')
            if sha(path) != spec['source_sha256']:
                raise ValueError('fixed-source file content differs from sealed hash')
            block=read_cif(path);pub=publication(block);pubs.append(pub.get('pmid'))
            items.append((spec,row,block,pub))
        except (OSError,ValueError,KeyError,RuntimeError) as exc:
            fetch_errors[sid]={'decision':'TECHNICAL_FAILURE','reason':type(exc).__name__, 'details':str(exc)[:400]}
    articles,article_receipt=article_batch(pubs)
    results=[]
    for spec,row,block,pub in items:
        pmid=pub.get('pmid')
        article=articles.get(str(pmid), {'state':('NO_PMID' if pmid in (None,'','-1','?','.') else
            'NOT_RETURNED' if article_receipt['state']=='SUCCESS' else 'SOURCE_TECHNICAL_FAILURE'),'pmid':pmid})
        try:
            output=assess(row,block,article)
        except (OSError,ValueError,KeyError,RuntimeError,TypeError,IndexError) as exc:
            output={'decision':'TECHNICAL_FAILURE','reason':type(exc).__name__,'details':str(exc)[:400]}
        results.append({'stratum':spec['stratum'],'system_id':spec['system_id'],**output})
    for spec in manifest['systems']:
        if spec['system_id'] in fetch_errors:
            results.append({'stratum':spec['stratum'],'system_id':spec['system_id'],**fetch_errors[spec['system_id']]})
    results.sort(key=lambda x:(x['stratum'],x['system_id']))
    counts=Counter(x['decision'] for x in results)
    receipt={'started_utc':started,'ended_utc':now(),'elapsed_s':round(time.monotonic()-clock,3),
        'inputs':len(manifest['systems']),'outputs':len(results),'decisions':dict(counts),
        'code_sha256':sha(__file__),'manifest_sha256':sha(args.manifest),
        'article_batch':article_receipt,'automated_run_manual_interventions':0,
        'manual_intervention_definition':'human intervention from launch through raw output write; independent audit excluded',
        'source':'sealed first-round research database and frozen mmCIF references','result_state':'RESEARCH_ONLY'}
    args.output.mkdir(parents=True,exist_ok=False)
    (args.output/'raw_results.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
    (args.output/'run_receipt.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(receipt,ensure_ascii=False),flush=True)


def manifest_for_ids(args, ids):
    db=read_db(args.database);systems=[]
    for stratum,sid in ids:
        row=load_record(db,sid);o=json.loads(row['document'])
        systems.append({'stratum':stratum,'system_id':sid,'pdb_id':o['identity']['pdb_id'],
                        'revision':o['identity']['revision'],'source_sha256':o['files']['source_mmcif']['sha256']})
    return {'created_utc':now(),'code_sha256':sha(__file__),'database':str(args.database.resolve()),'systems':systems}


def main():
    p=argparse.ArgumentParser(description=__doc__);s=p.add_subparsers(dest='command',required=True)
    dev=s.add_parser('dev-manifest');dev.add_argument('--database',type=Path,required=True);dev.add_argument('--output',type=Path,required=True)
    go=s.add_parser('run');go.add_argument('--database',type=Path,required=True);go.add_argument('--manifest',type=Path,required=True);go.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.command=='dev-manifest':
        ids=[('development','holo:1A6A:91edb76fe4ee34a6bdbec5307898'),
             ('development','holo:1AEC:d077ccfa00c1ee231d17344ba70f'),
             ('development','holo:4IB8:7e73180fd1bd98e0d156618f348b'),
             ('development','holo:3STH:1560341102304352bd56f4c6912e')]
        with args.output.open('x',encoding='utf-8') as f:json.dump(manifest_for_ids(args,ids),f,indent=2)
    else:run(args)


if __name__=='__main__':main()

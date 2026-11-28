"""Build a bounded, unlabelled author-supplement library from ordinary candidate PDB IDs.

The input is a JSON list of PDB accessions. Publication metadata only nominates files.
The caller still validates the fixed candidate, target site and explicit experimental group.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import re
from urllib.parse import urlsplit, urljoin
import xml.etree.ElementTree as ET

import requests
from apo_c_pilot import write_json, sha

PDBe = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/related_publications/'
PUBMED = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
STATE = re.compile(r'\b(?:unligated|apo|unliganded|unbound|(?:ligand|inhibitor)[- ]free)\b', re.I)


def discover(candidates, out, max_papers=8, max_files=12):
    out.mkdir(parents=True, exist_ok=True)
    counts, errors, papers, docs = Counter(), Counter(), {}, []

    def fetch(url, params=None, accept_pdf=False):
        key = sha256(json.dumps([url, params], sort_keys=True).encode()).hexdigest()
        p, meta = out / (key + '.body'), out / (key + '.json')
        if meta.is_file():
            rec = json.loads(meta.read_text()); data=p.read_bytes()
            if sha256(data).hexdigest()!=rec['sha256']: raise ValueError('cached author source changed')
            counts['hits']+=1
            return data,rec
        counts['requests']+=1
        with requests.get(url,params=params,timeout=(4,8),stream=True) as response:
            response.raise_for_status(); body=bytearray();limit=4_000_000 if accept_pdf else 1_500_000
            for chunk in response.iter_content(65536):
                body.extend(chunk)
                if len(body)>limit:raise ValueError('bounded author source size exceeded')
            data=bytes(body)
            rec={'url':response.url,'sha256':sha256(data).hexdigest(),'bytes':len(data),'path':str(p.resolve())}
        if accept_pdf and not data.startswith(b'%PDF-'):raise ValueError('advertised supplement is not a PDF')
        if not p.exists():p.write_bytes(data)
        if sha(p)!=rec['sha256']:raise ValueError('uncommitted response differs')
        write_json(meta,rec);return data,rec

    for code in dict.fromkeys(x.upper() for x in candidates):
        if not re.fullmatch(r'[0-9][A-Z0-9]{3}',code):raise ValueError('invalid PDB accession')
        try:
            data,_=fetch(PDBe+code.lower());entry=json.loads(data)[code.lower()]
            for group in ('appears_without_citation','cited_by'):
                for article in entry.get(group,{}).get('Articles',[]):
                    pmid=str(article.get('pubmed_id',''));title=article.get('title') or ''
                    if pmid.isdigit() and STATE.search(title) and len(papers)<max_papers:papers.setdefault(pmid,title)
        except (requests.RequestException,KeyError,ValueError) as exc:errors[type(exc).__name__]+=1
    if papers:
        xml,_=fetch(PUBMED,{'db':'pubmed','retmode':'xml','id':','.join(papers)})
        for item in ET.fromstring(xml).iter('PubmedArticle'):
            pmid=item.findtext('MedlineCitation/PMID');doi=item.findtext("PubmedData/ArticleIdList/ArticleId[@IdType='doi']",'').lower()
            if pmid not in papers or not doi.startswith('10.1107/') or len(docs)>=max_files:continue
            try:
                raw,_=fetch('https://api.crossref.org/works/'+doi);work=json.loads(raw)['message']
                if work['DOI'].lower()!=doi:continue
                uris=[x['URL'] for x in work.get('link',[]) if x.get('intended-application')=='text-mining']; links=[]
                for url in uris:
                    m=re.match(r'^https?://journals\.iucr\.org/([a-z]/issues/\d{4}/\d{2}/\d{2}/([a-z]{2}\d{4}))/[^/]+$',url,re.I)
                    if m:links.append(('https://journals.iucr.org/'+m[1]+'/',m[2]))
                if not links:continue
                directory,key=links[0];page,_=fetch('https://scripts.iucr.org/cgi-bin/sendsup?'+key);text=page.decode('utf8',errors='replace')
                if not re.search(r'<meta[^>]+name=["\']citation_doi["\'][^>]+content=["\']'+re.escape(doi)+r'["\']',text,re.I):continue
                files=sorted(set(re.findall(r'\b'+re.escape(key)+r'sup\d+\.pdf\b',text,re.I)))
                for name in files[:max_files-len(docs)]:
                    _,rec=fetch(urljoin(directory,name),accept_pdf=True)
                    if urlsplit(rec['url']).hostname=='journals.iucr.org':docs.append(rec|{'doi':doi,'pmid':pmid,'discovery':'PDBe PDB reference, PubMed DOI and publisher supplement listing'})
            except (requests.RequestException,KeyError,ValueError) as exc:errors[type(exc).__name__]+=1
    write_json(out/'documents.json',{'documents':docs});write_json(out/'collection.json',{'candidates':len(set(candidates)),'papers':len(papers),'documents':len(docs),'http':dict(counts),'unavailable':dict(errors)})
    return docs


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--candidates',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();discover(json.loads(args.candidates.read_text()),args.out)

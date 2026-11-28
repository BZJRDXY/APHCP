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
from tempfile import mkdtemp
from urllib.parse import urlsplit, urljoin
import xml.etree.ElementTree as ET

import requests
from apo_c_pilot import publish_bytes, read_json, write_json, sha

PDBe = 'https://www.ebi.ac.uk/pdbe/api/pdb/entry/related_publications/'
PUBMED = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
STATE = re.compile(r'\b(?:unligated|apo|unliganded|unbound|(?:ligand|inhibitor)[- ]free)\b', re.I)


def discover(candidates, out, max_papers=8, max_files=12, *, cache=None, max_requests=None):
    out.mkdir(parents=True, exist_ok=True)
    responses = cache if cache is not None else out
    counts, errors, papers, docs = Counter(), Counter(), {}, []

    def fetch(url, params=None, accept_pdf=False):
        key = sha256(json.dumps([url, params], sort_keys=True).encode()).hexdigest()
        meta = responses / (key + '.json')
        if meta.is_file():
            counts['hits'] += 1
        else:
            if max_requests is not None and counts['requests'] >= max_requests:
                raise RequestLimit('publisher discovery request budget reached')
            if cache is not None and url == PUBMED:
                from apo_c_literature import wait_for_ncbi
                wait_for_ncbi(cache.parent)
            counts['requests'] += 1
            with requests.get(url, params=params, timeout=(4, 8), stream=True) as response:
                response.raise_for_status()
                body = bytearray()
                limit = 4_000_000 if accept_pdf else 1_500_000
                for chunk in response.iter_content(65536):
                    body.extend(chunk)
                    if len(body) > limit:
                        raise ValueError('author source size limit exceeded')
                data = bytes(body)
                if accept_pdf and not data.startswith(b'%PDF-'):
                    raise ValueError('advertised supplement is not a PDF')
                digest = sha256(data).hexdigest()
                path = responses / 'bodies' / digest
                if not path.exists():
                    try:
                        publish_bytes(path, data)
                    except FileExistsError:
                        pass
                row = {'url': response.url, 'sha256': digest, 'bytes': len(data), 'path': str(path.resolve())}
                counts['response_bytes'] += len(data)
                try:
                    write_json(meta, row)
                except FileExistsError:
                    pass  # Concurrent workers use the first complete response for this request.
        row = read_json(meta)
        data = Path(row['path']).read_bytes()
        if sha256(data).hexdigest() != row['sha256']:
            raise RuntimeError('publisher response differs from its receipt')
        if accept_pdf and not data.startswith(b'%PDF-'):
            raise ValueError('cached response is not a PDF')
        return data, row

    for code in dict.fromkeys(x.upper() for x in candidates):
        if not re.fullmatch(r'[0-9][A-Z0-9]{3}',code):raise ValueError('invalid PDB accession')
        try:
            data,_=fetch(PDBe+code.lower());entry=json.loads(data)[code.lower()]
            for group in ('appears_without_citation','cited_by'):
                for article in entry.get(group,{}).get('Articles',[]):
                    pmid=str(article.get('pubmed_id',''));title=article.get('title') or ''
                    if pmid.isdigit() and STATE.search(title) and len(papers)<max_papers:papers.setdefault(pmid,title)
        except (requests.RequestException,KeyError,ValueError) as exc:errors[type(exc).__name__]+=1
    articles = ()
    if papers:
        try:
            xml, _ = fetch(PUBMED, {'db': 'pubmed', 'retmode': 'xml', 'id': ','.join(papers)})
            articles = ET.fromstring(xml).iter('PubmedArticle')
        except (requests.RequestException, ValueError, ET.ParseError) as exc:
            errors[type(exc).__name__] += 1
    if articles:
        for item in articles:
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


class RequestLimit(ValueError):
    """The optional lookup used its per-candidate network quota."""


class CandidateDocuments:
    """Look for a publisher supplement only when a measured site still lacks an experiment."""

    def __init__(self, run, cache, limit):
        self.root = run / 'author_discovery'
        self.cache = cache / 'publisher_responses'
        self.limit = limit
        self.answers = {}
        self.proofs = {}

    def lookup(self, entry):
        from apo_c_literature import AuthorDocuments
        if not re.fullmatch(r'[0-9][A-Z0-9]{3}', entry):
            raise ValueError('invalid candidate accession')
        if entry in self.proofs:
            return self.proofs[entry], {'status': 'INDEXED_PUBLISHER_PROOF'}
        if entry in self.answers:
            return self.answers[entry]
        path = self.root / entry
        if not path.exists() and self.root.is_dir() and sum(p.is_dir() for p in self.root.iterdir()) >= self.limit:
            return None, {'status': 'CANDIDATE_BUDGET'}
        path.mkdir(parents=True, exist_ok=True)
        done = path / 'complete.json'
        if done.is_file():
            record = read_json(done)
        else:
            attempt = Path(mkdtemp(prefix='attempt-', dir=path))
            discover([entry], attempt, max_papers=2, max_files=2, cache=self.cache, max_requests=8)
            stats = read_json(attempt / 'collection.json')
            errors = set(stats['unavailable'])
            manifest = attempt / 'documents.json'
            record = {'manifest': str(manifest.resolve()), 'sha256': sha(manifest),
                      'http': stats['http'], 'unavailable': stats['unavailable'],
                      'status': 'SOURCE_UNAVAILABLE' if errors - {'RequestLimit'} else
                                'REQUEST_BUDGET' if errors else 'COMPLETE'}
        manifest = Path(record['manifest'])
        if sha(manifest) != record['sha256']:
            raise ValueError('publisher manifest changed')
        indexed = AuthorDocuments(manifest)
        for pdb in indexed.by_entry:
            self.proofs.setdefault(pdb, indexed.support(pdb))
        if not done.is_file() and record['status'] != 'SOURCE_UNAVAILABLE':
            write_json(done, record)
        self.answers[entry] = indexed.support(entry), record
        return self.answers[entry]


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--candidates',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();discover(json.loads(args.candidates.read_text()),args.out)

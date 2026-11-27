"""Resolve experimental Apo evidence from accession-linked author sources."""
from __future__ import annotations

from base64 import b64decode, b64encode
from collections import Counter
import gzip
from hashlib import sha256
import html
import json
import os
from pathlib import Path
import re
from tempfile import NamedTemporaryFile
import time
import xml.etree.ElementTree as ET

import requests

from apo_c_pilot import cache_get


EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
PMC_OAI = 'https://pmc.ncbi.nlm.nih.gov/api/oai/v1/mh/'
APO_WORDS = r'(?:apo|unliganded|(?:ligand|inhibitor)[- ]free|unbound|uninhibited|vacant active[- ]site)'
APO_STATE = re.compile(r'\b' + APO_WORDS + r'\b', re.I)
STARTS_APO = re.compile(r'^' + APO_WORDS + r'\b', re.I)
ACTIVE_SITE_INHIBITED = re.compile(r'\b(?:active|catalytic)[- ]site[- ](?:inhibited|blocked|occupied)\b', re.I)
STRUCTURE_REPORT = re.compile(r'\b(?:x[- ]ray|crystal\w*|structures? of|structures? in|structural studies)\b', re.I)
NONEXPERIMENTAL = re.compile(r'\b(?:hypothetical|in silico|predicted|simulation|could not|not determined)\b', re.I)
PDB_ID = re.compile(r'\b(?=[0-9][A-Z0-9]{3}\b)(?=[0-9A-Z]{0,3}[A-Z])[0-9][A-Z0-9]{3}\b', re.I)


def xml_document(payload):
    root = ET.fromstring(payload)
    for node in root.iter():
        node.tag = node.tag.rsplit('}', 1)[-1]
    return root


def node_text(node):
    return ' '.join(''.join(node.itertext()).split())


def deposited_apo(title):
    return bool(APO_STATE.search(title) and not NONEXPERIMENTAL.search(title))


def reported_apo_structure(abstract):
    for sentence in re.split(r'(?<=[.!?])\s+(?=[A-Z])', abstract):
        if APO_STATE.search(sentence) and STRUCTURE_REPORT.search(sentence) and not NONEXPERIMENTAL.search(sentence):
            return sentence
    return None


def accession_state(article, entry):
    """Read exact accession associations in a table or deposition statement."""
    exact = re.escape(entry)
    for table in article.iter('table'):
        for row in table.iter('tr'):
            cells = [node_text(c) for c in row if c.tag in {'td', 'th'}]
            if ([c.upper() for c in cells if PDB_ID.fullmatch(c)] == [entry]
                    and any(STARTS_APO.match(c) for c in cells)):
                return 'ARTICLE_ACCESSION_TABLE', ' | '.join(cells)
    for container in (article.find('body'), article.find('back')):
        if container is None:
            continue
        for paragraph in container.iter('p'):
            text = node_text(paragraph)
            if not re.search(rf'\b{exact}\b', text, re.I) or not APO_STATE.search(text):
                continue
            state = re.search(rf'\b{exact}\b\s*\(([^()]{{1,140}})\)', text, re.I)
            if state and STARTS_APO.match(state.group(1)):
                return 'ARTICLE_ACCESSION_STATE', state.group()
            if not re.search(r'\brespectively\b', text, re.I):
                continue
            ordered = re.search(
                r'\bstructures? of (.+?)\s+(?:have been|were)\s+(?:submitted|deposited)\b.+?'
                r'\b(?:identifier|accession|ID)\s*(?:codes?|numbers?)?\s*(?:are|were|:)\s*([^:.]+)',
                text, re.I)
            if ordered:
                states = re.split(r'\s*,\s*|\s+and\s+', ordered.group(1))
                ids = [p.upper() for p in PDB_ID.findall(ordered.group(2))]
                if len(ids) == len(states) and ids.count(entry) == 1 and APO_STATE.search(states[ids.index(entry)]):
                    return 'ARTICLE_ORDERED_DEPOSITION', states[ids.index(entry)] + ' -> ' + entry
    for wrap in article.iter('table-wrap'):
        caption = wrap.find('caption')
        if caption is None or not re.search(r'crystallograph|diffraction', node_text(caption), re.I):
            continue
        note = next((node_text(n) for n in wrap.iter('fn') if APO_STATE.search(node_text(n))
                     and re.search(r'\b(?:collect\w*|determin\w*|diffraction|solv\w*)\b', node_text(n), re.I)), None)
        if not note:
            continue
        rows = [[node_text(c) for c in row if c.tag in {'td', 'th'}] for row in wrap.iter('tr')]
        for headers, codes in zip(rows, rows[1:]):
            if (len(headers) != len(codes) or len(codes) < 3 or not
                    re.match(r'^(?:Protein Data Bank|PDB)\s+(?:entry|code|ID|accession)\b', codes[0], re.I)):
                continue
            ids = [c.upper() for c in codes[1:]]
            if ids.count(entry) != 1 or not all(PDB_ID.fullmatch(code) for code in ids):
                continue
            bound = re.compile(r'\+|\b(?:bound|complex(?:ed)? with|liganded)\b', re.I)
            unbound = [i for i, name in enumerate(headers[1:]) if name and not bound.search(name)]
            if unbound == [ids.index(entry)] and all(headers[1:]):
                return 'ARTICLE_APO_CRYSTAL_COMPARISON', f'{entry} | {headers[unbound[0] + 1]} | {note}'
    return None


def accession_conflict(text, entry):
    """Bind an explicit inhibited-site statement to the accession in its clause."""
    pattern = rf'\b{re.escape(entry)}\b'
    for match in re.finditer(pattern, text, re.I):
        after = re.match(r'\s*\(([^()]{{1,140}})\)', text[match.end():])
        if after and ACTIVE_SITE_INHIBITED.search(after.group(1)):
            return match.group() + after.group()
        before = re.split(r'(?<=[.;])\s+', text[max(0, match.start() - 155):match.start()])[-1]
        prior_codes = list(PDB_ID.finditer(before))
        if prior_codes:
            before = before[prior_codes[-1].end():]
        if re.search(r'\b(?:PDB|Protein Data Bank)\b', before, re.I) and ACTIVE_SITE_INHIBITED.search(before):
            return before + match.group()
    return None


class PrimaryLiterature:
    """Reuse verified responses across Runs and read full text only for unresolved evidence."""

    def __init__(self, run, shared_cache=None):
        self.run = run
        self.shared_cache = Path(shared_cache) if shared_cache is not None else None
        self.documents = {}
        self.articles = {}
        self.stats = Counter()

    def get(self, name, url, params, accepted=(200,)):
        key = sha256(json.dumps([url, params], sort_keys=True).encode()).hexdigest()
        shared = self.shared_cache / (key + '.json.gz') if self.shared_cache is not None else None
        if shared is not None and shared.is_file():
            try:
                with gzip.open(shared, 'rt', encoding='utf8') as stream:
                    cached = json.load(stream)
                body, receipt = b64decode(cached['payload'], validate=True), cached['receipt']
                if receipt['status'] in accepted and sha256(body).hexdigest() == receipt['response_sha256']:
                    self.stats['shared_hits'] += 1
                    return body, receipt | {'literature_cache': 'shared'}
            except (OSError, ValueError, KeyError):
                self.stats['shared_read_errors'] += 1

        def validate(body, status):
            if status == 200:
                if name.startswith(('doi/', 'accession/')):
                    if 'esearchresult' not in json.loads(body):
                        raise ValueError('invalid literature search response')
                elif xml_document(body).tag not in {'PubmedArticleSet', 'OAI-PMH'}:
                    raise ValueError('invalid literature record response')

        local = self.run / 'http' / 'literature' / name
        reused = local.with_suffix('.gz').is_file() and local.with_suffix('.json').is_file()
        if not reused:
            time.sleep(0.35)
            self.stats['network_requests'] += 1
            if name.startswith('pmc/'):
                self.stats['network_full_text_requests'] += 1
        body, receipt = cache_get(self.run, 'literature/' + name, url, params,
                                  timeout=12, accepted=accepted, validate=validate)
        if reused:
            self.stats['run_cache_hits'] += 1
        else:
            self.stats['network_response_seconds'] += receipt['elapsed_s']
            self.stats['network_response_bytes'] += receipt['response_bytes']
        if shared is not None and receipt['status'] == 200:
            tmp = None
            try:
                shared.parent.mkdir(parents=True, exist_ok=True)
                with NamedTemporaryFile(dir=shared.parent, delete=False) as temporary:
                    tmp = Path(temporary.name)
                    contents = {'payload': b64encode(body).decode('ascii'), 'receipt': receipt}
                    temporary.write(gzip.compress(json.dumps(contents).encode()))
                os.replace(tmp, shared)
            except OSError:
                self.stats['shared_write_errors'] += 1
            finally:
                if tmp is not None:
                    tmp.unlink(missing_ok=True)
        return body, receipt

    def get_document(self, citation):
        pmid = str(citation.get('pdbx_database_id_PubMed') or '').strip()
        doi = str(citation.get('pdbx_database_id_DOI') or '').strip().lower()
        key = (pmid, doi)
        if key in self.documents:
            return self.documents[key]
        record = {'status': 'PRIMARY_CITATION_UNAVAILABLE'}
        abstract, pmc = '', None
        try:
            if not pmid.isdecimal() and doi.startswith('10.'):
                payload, _ = self.get('doi/' + sha256(doi.encode()).hexdigest(), EUTILS + '/esearch.fcgi',
                                      {'db': 'pubmed', 'term': doi + '[aid]', 'retmode': 'json', 'retmax': 2})
                ids = json.loads(payload)['esearchresult']['idlist']
                pmid = ids[0] if len(ids) == 1 else ''
            if pmid.isdecimal():
                payload, meta = self.get('pubmed/' + pmid, EUTILS + '/efetch.fcgi',
                                         {'db': 'pubmed', 'id': pmid, 'retmode': 'xml'})
                doc = xml_document(payload)
                item = next((r for r in doc.iter('PubmedArticle')
                             if r.findtext('MedlineCitation/PMID') == pmid), None)
                ids = ({i.attrib.get('IdType'): i.text for i in item.findall('PubmedData/ArticleIdList/ArticleId')}
                       if item is not None else {})
                found_doi = (ids.get('doi') or '').lower()
                if item is not None and (not doi.startswith('10.') or not found_doi or found_doi == doi):
                    abstract = ' '.join(node_text(a) for a in item.iter('AbstractText'))
                    deposits = sorted({(x.text or '').upper() for bank in item.findall(
                        'MedlineCitation/Article/DataBankList/DataBank')
                        if (bank.findtext('DataBankName') or '').upper() == 'PDB'
                        for x in bank.findall('AccessionNumberList/AccessionNumber')})
                    record = {'status': 'PRIMARY_CITATION_READ', 'pmid': pmid, 'doi': found_doi or doi,
                              'abstract_source': meta, 'pdb_deposits': deposits}
                    pmc = ids.get('pmc') if re.fullmatch(r'PMC\d+', ids.get('pmc') or '') else None
                else:
                    record = {'status': 'PRIMARY_CITATION_IDENTITY_UNRESOLVED', 'pmid': pmid, 'doi': doi}
        except (requests.RequestException, RuntimeError, ValueError, ET.ParseError, KeyError) as exc:
            record |= {'fetch_error': type(exc).__name__}
        if 'fetch_error' not in record or abstract:
            self.documents[key] = record, abstract, pmc
        return record, abstract, pmc

    def full_article(self, pmc):
        if pmc in self.articles:
            return self.articles[pmc]
        payload, receipt = self.get('pmc/' + pmc, PMC_OAI,
            {'verb': 'GetRecord', 'identifier': 'oai:pubmedcentral.nih.gov:' + pmc[3:],
             'metadataPrefix': 'pmc'}, accepted=(200, 400, 404))
        article = xml_document(payload).find('.//record/metadata/article') if receipt['status'] == 200 else None
        if article is not None:
            ids = {i.text.removeprefix('PMC') for i in article.findall('front/article-meta/article-id')
                   if i.text and i.get('pub-id-type') in {'pmc', 'pmcid'}}
            if pmc[3:] not in ids:
                article = None
        self.articles[pmc] = article, receipt
        return article, receipt

    def candidate_sources(self, experiment, entry):
        citation = experiment['primary_citation']
        original = self.get_document(citation)
        yield original, True
        if original[0].get('status') == 'PRIMARY_CITATION_UNAVAILABLE' and PDB_ID.fullmatch(entry):
            payload, receipt = self.get('accession/' + entry, EUTILS + '/esearch.fcgi',
                {'db': 'pubmed', 'term': f'"PDB/{entry}"[si]', 'retmode': 'json', 'retmax': 3})
            for pmid in json.loads(payload)['esearchresult']['idlist']:
                record, abstract, pmc = self.get_document({'pdbx_database_id_PubMed': pmid})
                if entry in record.get('pdb_deposits', ()):
                    yield (record | {'accession_source': receipt}, abstract, pmc), True
        holo = experiment.get('holo_primary_citation')
        if holo and holo != citation:
            yield self.get_document(holo), False

    def assess(self, experiment):
        entry = experiment['source_entry'].upper()
        outcome = {'status': 'PRIMARY_CITATION_UNAVAILABLE'}
        checked = set()
        try:
            for (record, abstract, pmc), candidate_linked in self.candidate_sources(experiment, entry):
                if not record.get('pmid'):
                    if record.get('fetch_error'):
                        outcome['fetch_error'] = record['fetch_error']
                    continue
                if record['pmid'] in checked:
                    continue
                checked.add(record['pmid'])
                outcome = dict(record)
                text = html.unescape(re.sub(r'<[^>]+>', '', abstract))
                conflict = accession_conflict(text, entry)
                sentence = reported_apo_structure(text)
                source = {'entry': entry, 'pmid': record['pmid'], 'doi': record.get('doi')}
                if conflict:
                    return outcome | {'conflict': source | {'basis': 'ACCESSION_ACTIVE_SITE_INHIBITED', 'excerpt': conflict}}
                if candidate_linked and deposited_apo(experiment['title']) and sentence:
                    return outcome | {'support': source | {'basis': 'MATCHED_APO_DEPOSITION_AND_PRIMARY_ABSTRACT',
                                                           'excerpt': sentence[:600]}}
                if not pmc:
                    continue
                article, receipt = self.full_article(pmc)
                outcome['full_text_source'] = receipt
                if article is None:
                    continue
                source['pmcid'] = pmc
                conflict = next((issue for p in article.iter('p')
                                 if (issue := accession_conflict(node_text(p), entry))), None)
                if conflict:
                    return outcome | {'conflict': source | {'basis': 'ACCESSION_ACTIVE_SITE_INHIBITED', 'excerpt': conflict}}
                precise = accession_state(article, entry)
                if precise:
                    basis, excerpt = precise
                    return outcome | {'support': source | {'basis': basis, 'excerpt': excerpt[:600]}}
        except (requests.RequestException, RuntimeError, ValueError, ET.ParseError, KeyError) as exc:
            outcome['fetch_error'] = type(exc).__name__
        return outcome

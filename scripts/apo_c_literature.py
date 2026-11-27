"""Resolve candidate-deposition Apo evidence from the primary experiment."""
from __future__ import annotations

from hashlib import sha256
import html
import json
import re
import time
import xml.etree.ElementTree as ET

import requests

from apo_c_pilot import cache_get


EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
PMC_OAI = 'https://pmc.ncbi.nlm.nih.gov/api/oai/v1/mh/'
APO_STATE = re.compile(r'\b(?:apo|unliganded|ligand[- ]free|unbound)\b', re.I)
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
                    and any(re.match(r'^(?:apo|unliganded|ligand[- ]free|unbound)\b', c, re.I) for c in cells)):
                return 'ARTICLE_ACCESSION_TABLE', ' | '.join(cells)
    for container in (article.find('body'), article.find('back')):
        if container is None:
            continue
        for paragraph in container.iter('p'):
            text = node_text(paragraph)
            if not re.search(rf'\b{exact}\b', text, re.I) or not APO_STATE.search(text):
                continue
            state = re.search(rf'\b{exact}\b\s*\(([^()]{{1,140}})\)', text, re.I)
            if state and re.match(r'^(?:apo|unliganded|ligand[- ]free|unbound)\b', state.group(1), re.I):
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
    return None


class PrimaryLiterature:
    """Share downloaded primary records within a Run and retry transient failures on resume."""

    def __init__(self, run):
        self.run = run
        self.documents = {}

    def get(self, name, url, params, accepted=(200,)):
        # NCBI public E-utilities permit at most three requests per second without a key.
        time.sleep(0.35)
        return cache_get(self.run, 'literature/' + name, url, params,
                         timeout=12, accepted=accepted)

    def get_document(self, citation):
        pmid = str(citation.get('pdbx_database_id_PubMed') or '').strip()
        doi = str(citation.get('pdbx_database_id_DOI') or '').strip().lower()
        key = pmid if pmid.isdecimal() else doi
        if key in self.documents:
            return self.documents[key]
        record = {'status': 'PRIMARY_CITATION_UNAVAILABLE'}
        abstract, article = '', None
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
                    record = {'status': 'PRIMARY_CITATION_READ', 'pmid': pmid, 'doi': found_doi or doi,
                              'abstract_source': meta}
                    pmc = ids.get('pmc') or ''
                    if re.fullmatch(r'PMC\d+', pmc):
                        payload, meta = self.get('pmc/' + pmc, PMC_OAI,
                            {'verb': 'GetRecord', 'identifier': 'oai:pubmedcentral.nih.gov:' + pmc[3:],
                             'metadataPrefix': 'pmc'}, accepted=(200, 400, 404))
                        record['full_text_source'] = meta
                        if meta['status'] == 200:
                            full = xml_document(payload)
                            article = full.find('.//record/metadata/article')
                            if article is not None:
                                actual_ids = {i.text.removeprefix('PMC') for i in article.findall(
                                    'front/article-meta/article-id') if i.text and i.get('pub-id-type') in {'pmc', 'pmcid'}}
                                if pmc[3:] not in actual_ids:
                                    article = None
                else:
                    record = {'status': 'PRIMARY_CITATION_IDENTITY_UNRESOLVED', 'pmid': pmid, 'doi': doi}
        except (requests.RequestException, RuntimeError, ValueError, ET.ParseError, KeyError) as exc:
            record |= {'fetch_error': type(exc).__name__}
        if 'fetch_error' not in record or abstract:
            self.documents[key] = record, abstract, article
        return record, abstract, article

    def assess(self, experiment):
        record, abstract, article = self.get_document(experiment['primary_citation'])
        entry = experiment['source_entry'].upper()
        precise = accession_state(article, entry) if article is not None else None
        sentence = reported_apo_structure(html.unescape(re.sub(r'<[^>]+>', '', abstract)))
        if precise:
            basis, excerpt = precise
        elif deposited_apo(experiment['title']) and sentence:
            basis, excerpt = 'MATCHED_APO_DEPOSITION_AND_PRIMARY_ABSTRACT', sentence
        else:
            return dict(record)
        return record | {'support': {'basis': basis, 'entry': entry, 'pmid': record['pmid'],
                                     'doi': record['doi'], 'excerpt': excerpt[:600]}}

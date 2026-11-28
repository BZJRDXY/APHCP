"""Resolve experimental Apo evidence from accession-linked author sources."""
from __future__ import annotations

from base64 import b64decode, b64encode
from collections import Counter
from decimal import Decimal, InvalidOperation
import gzip
from hashlib import sha256
import html
from html.parser import HTMLParser
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from tempfile import NamedTemporaryFile
import time
from urllib.parse import parse_qs, quote, urljoin, urlsplit, urlunsplit
import xml.etree.ElementTree as ET

import requests

from apo_c_pilot import cache_get


EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
PMC_OAI = 'https://pmc.ncbi.nlm.nih.gov/api/oai/v1/mh/'
CROSSREF_WORK = 'https://api.crossref.org/works/'
BIOSTUDIES = 'https://www.ebi.ac.uk/biostudies/'
ARCHIVE_AVAILABILITY = 'https://archive.org/wayback/available'
PUBLISHER_TEXT_HOSTS = {'journals.iucr.org', 'pnas.org', 'www.pnas.org'}
PDB_ENTRY = 'https://data.rcsb.org/rest/v1/core/entry/'
APO_WORDS = r'(?:apo|unliganded|(?:ligand|inhibitor)[- ]free|unbound|uninhibited|vacant active[- ]site)'
APO_STATE = re.compile(r'\b' + APO_WORDS + r'\b', re.I)
STARTS_APO = re.compile(r'^' + APO_WORDS + r'\b', re.I)
ACTIVE_SITE_INHIBITED = re.compile(r'\b(?:active|catalytic)[- ]site[- ](?:inhibited|blocked|occupied)\b', re.I)
BOUND_STATE = re.compile(r'\b(?:bound to|ligand[- ]bound|inhibitor[- ]bound|complex(?:ed)? with|'
                         r'molecular replacement|template|starting model)\b', re.I)
STRUCTURE_REPORT = re.compile(r'\b(?:x[- ]ray|crystal\w*|structures? of|structures? in|structural studies)\b', re.I)
NONEXPERIMENTAL = re.compile(r'\b(?:hypothetical|in silico|predicted|simulation|could not|not determined)\b', re.I)
PDB_ID = re.compile(r'\b(?=[0-9][A-Z0-9]{3}\b)(?=[0-9A-Z]{0,3}[A-Z])[0-9][A-Z0-9]{3}\b', re.I)
PDB_LIST = re.compile(r'\b(?:PDB|Protein Data Bank)\b.{0,28}\b(?:codes?|references?|accessions?)\b|'
                      r'\bData deposition\b|\bdeposited\b.*\b(?:PDB|Protein Data Bank)\b', re.I)
ACCESSION_LABEL = re.compile(r'\(\s*(?:(?:PDB\s+)?(?:accession|id|code)(?:\s+number)?\s*#?\s*|#\s*)'
                             r'(' + PDB_ID.pattern + r')\s*\)', re.I)
CRYSTAL_FIELDS = ('a', 'b', 'c', 'alpha', 'beta', 'gamma', 'diffraction_high', 'refinement_high')
PDF_DECIMAL = r'\d+(?:\.\d+)?'


def pdf_layout(payload):
    converter = shutil.which('pdftotext')
    if converter is None:
        return None
    result = subprocess.run([converter, '-layout', '-enc', 'UTF-8', '-', '-'], input=payload,
                            capture_output=True, check=True, timeout=20)
    return result.stdout.decode('utf8')


def crystal_apo_groups(text):
    """Read crystallographic measurements from the first, explicitly labelled group."""
    for page_number, page in enumerate(text.split('\f'), 1):
        lines = page.splitlines()
        for index, line in enumerate(lines):
            if not re.match(r'^\s*Table\s+\w+[. :]\s*.*\b(?:data collection|refinement)\b', line, re.I):
                continue
            rows = lines[index + 1:]
            stop = next((i for i, row in enumerate(rows) if re.match(r'^\s*Table\s+\w+[. :]', row, re.I)), len(rows))
            rows = rows[:stop]
            for head, row in enumerate(rows[:10]):
                if not re.match(r'^\s*Protein\s{2,}', row, re.I) or head == 0:
                    continue
                labels = list(re.finditer(r'\S(?:.*?\S)?(?=\s{2,}|$)', row.rstrip()))
                groups = list(re.finditer(r'\S(?:.*?\S)?(?=\s{2,}|$)', rows[head - 1].rstrip()))
                if len(labels) < 3 or len(groups) != 2:
                    continue
                state, comparison = groups
                first = labels[1]
                if (not explicitly_unbound(state.group()) or not re.search(r'\b(?:bound|complex|crystalli[sz]ed with)\b', comparison.group(), re.I)
                        or not labels[0].end() < first.start() <= state.end() < comparison.start()
                        or first.end() >= comparison.start()):
                    continue
                # Require each number to align with the first data column.
                center = (first.start() + first.end()) / 2
                column_start = center - 9
                column_end = (center + (labels[2].start() + labels[2].end()) / 2) / 2
                rules = {
                    'a': r'a\s*,.*?', 'b': r'b\s*,.*?', 'c': r'c\s*,.*?',
                    'alpha': r'(?:alpha|\u03b1)\s*,.*?', 'beta': r'(?:beta|\u03b2)\s*,.*?',
                    'gamma': r'(?:gamma|\u03b3)\s*,.*?', 'wavelength': r'wavelength\s*,.*?',
                    'diffraction_high': r'resolution\s*,.*?', 'refinement_high': r'resolution range\s*,.*?',
                }
                measurements = {}
                positions = {}
                for field, label in rules.items():
                    value = rf'{PDF_DECIMAL}\s*[-\u2013\u2212]\s*({PDF_DECIMAL})' if field.endswith('_high') else rf'({PDF_DECIMAL})'
                    pattern = re.compile(rf'^\s*{label}\s{{2,}}(?=\S)', re.I)
                    matching_rows = [(i, m.end()) for i, record in enumerate(rows) if (m := pattern.match(record))]
                    if len(matching_rows) != 1:
                        continue
                    position, start = matching_rows[0]
                    if field in {'a', 'b', 'c', 'alpha', 'beta', 'gamma'} and not any(
                            re.match(r'^\s*Unit cell\s*$', previous, re.I) for previous in rows[max(0, position - 8):position]):
                        continue
                    found = re.match(value + r'(?=\s|$)', rows[position][start:])
                    if found and column_start <= start and start + found.end() <= column_end:
                        measurements[field] = found.group(1)
                        positions[field] = index + position + 2
                if all(key in measurements for key in CRYSTAL_FIELDS):
                    yield measurements, {'page': page_number, 'table_line': index + 1,
                                         'state_line': index + head + 1, 'measurements': positions,
                                         'group': first.group(), 'state': state.group()}


def same_crystal(author, deposited):
    """Compare the eight required values at the precision reported by the article."""
    missing = False
    for field, value in author.items():
        try:
            written, actual = Decimal(value), Decimal(str(deposited.get(field)))
            if not actual.is_finite():
                raise InvalidOperation
        except InvalidOperation:
            missing |= field in CRYSTAL_FIELDS
            continue
        precision = Decimal(5).scaleb(written.as_tuple().exponent - 1)
        if abs(written - actual) > precision:
            return False
    return None if missing else True


def biostudies_links(block):
    links = block.get('links') or []
    for group in links:
        for link in (group if isinstance(group, list) else [group]):
            for attribute in link.get('attributes', []):
                if attribute.get('name') == 'Type' and link.get('url'):
                    yield attribute.get('value', '').upper(), link['url']


def biostudies_article(study, pmc, doi, pmid):
    if study.get('accno') != 'S-E' + pmc:
        return None
    section = study.get('section', {})
    for index in section.get('subsections', []):
        for subsection in (index if isinstance(index, list) else [index]):
            if subsection.get('type') != 'Publication':
                continue
            links = {(t, url.lower()) for t, url in biostudies_links(subsection)}
            if {('PMC', pmc.lower()), ('DOI', doi.lower()), ('PMID', pmid)} <= links:
                return section
    return None


def xml_document(payload):
    root = ET.fromstring(payload)
    for node in root.iter():
        node.tag = node.tag.rsplit('}', 1)[-1]
    return root


def node_text(node):
    return ' '.join(''.join(node.itertext()).split())


def deposited_apo(title):
    return bool(APO_STATE.search(title) and not NONEXPERIMENTAL.search(title))


def explicitly_unbound(description):
    negated = re.search(r'\b(?:not|non[- ]?)\s*' + APO_WORDS + r'\b', description, re.I)
    return (deposited_apo(description) and not negated and not BOUND_STATE.search(description)
            and not ACTIVE_SITE_INHIBITED.search(description)
            and not re.search(r'\+|\bor\b|\b(?:liganded|holo)\b', description, re.I))


BACKGROUND_STATE = re.compile(
    r'\b(?:no (?:\w+ ){0,4}(?:structures?|data)|neither|never|previous(?:ly)?|yet to be|remains? unknown|'
    r'not (?:yet |previously )?(?:available|reported|known|solved|determined)|lacking|unavailable|'
    r'would|might|predicted|template|used (?:as|for) (?:a )?(?:model|molecular replacement))\b', re.I)
PRESENT_RESULT = re.compile(
    r'\b(?:we (?:also |have |here |now )*(?:report|present|determined?|describe|obtained?|solved?|show|reveal)\w*|'
    r'here[, ]|this (?:work|study)|we crystall[iy][sz]ed|(?:crystalli[sz]es? only at|determined at .{0,30}resolution|'
    r'(?:apo states?|structural studies|crystal structures?|co-crystallization).{0,110}?reveal(?:s|ed)?\b)|'
    r'(?:was|were|have been|has been) (?:also )?(?:determined|solved|obtained|reported)|'
    r'(?:the|a) (?:first )?(?:crystal |x[- ]ray )?structures? (?:of|in) .{0,80} (?:is|are) (?:presented|reported))\b', re.I)


def reported_apo_structure(abstract):
    for sentence in re.split(r'(?<=[.!?])\s+(?=[A-Z])', abstract):
        if (APO_STATE.search(sentence) and STRUCTURE_REPORT.search(sentence) and PRESENT_RESULT.search(sentence)
                and not NONEXPERIMENTAL.search(sentence) and not BACKGROUND_STATE.search(sentence)
                and not re.search(r'\b(?:not|never|no)\b.{0,34}\b' + APO_WORDS + r'\b', sentence, re.I)):
            return sentence
    return None


def fixed_relation_support(experiment, source):
    """A relation labels the accession on its own row, never the structure containing it."""
    if not source or experiment.get('preparation_conflicts_or_unresolved_additions'):
        return None
    entry = experiment['source_entry'].upper()
    if entry == source['entry'] or BOUND_STATE.search(experiment['title']) or ACTIVE_SITE_INHIBITED.search(experiment['title']):
        return None
    rows = [(i, note) for i, pdb, note in source['rows'] if pdb == entry]
    if not rows or not all(note and explicitly_unbound(note) and not BACKGROUND_STATE.search(note)
                           and not re.search(r'\b(?:hypothes\w*|predicted|cofactor[- ]free|metal[- ]free|'
                                             r'apolipoprotein|apo[- ](?:heme|zinc|calcium|metal)|by analogy)\b', note, re.I)
                           for _, note in rows):
        return None
    return {'basis': 'FIXED_HOLO_DEPOSITION_NAMES_CANDIDATE_STATE', 'entry': entry,
            'source_entry': source['entry'], 'source_sha256': source['sha256'], 'source_revision': source['revision'],
            'field': '_pdbx_database_related.details', 'row_numbers': [i for i, _ in rows], 'excerpt': rows[0][1]}


def author_deposition_state(text, entry, paragraphs):
    if not PDB_LIST.search(text):
        return None
    records = []
    previous = 0
    for m in ACCESSION_LABEL.finditer(text):
        description = text[previous:m.start()].split(':')[-1].strip(' ,;.')
        description = re.sub(r'^and\s+', '', description, flags=re.I)
        records.append((m.group(1).upper(), description))
        previous = m.end()
    selected = [description for code, description in records if code == entry]
    if len(selected) != 1 or len(selected[0]) > 150:
        return None
    description = selected[0]
    if explicitly_unbound(description):
        return 'ARTICLE_STATE_BEFORE_ACCESSION', f'{description} ({entry})', None
    # Some deposition lists label the apo arm with the plain construct name.
    # Require a matched methods statement and a bound arm of that same construct.
    if not re.fullmatch(r'[\w/. -]{3,75}', description) or re.search(r'\b(?:bound|complex|mutant|model)\b', description, re.I):
        return None
    if sum(d.casefold() == description.casefold() for _, d in records) != 1:
        return None
    name = re.escape(description)
    if not any(re.match(name + r'\s+bound to\b', d, re.I) for _, d in records):
        return None
    method = re.compile(r'\b' + name + r'\s+(?:was|were)\s+crystalli[sz]ed\s+in\s+(?:its|the)\s+'
                        r'apo(?:[ -]form)?\s+and\s+in\s+complex', re.I)
    for paragraph, content in paragraphs:
        if match := method.search(content):
            return 'ARTICLE_METHODS_WITH_DEPOSITION_GROUP', f'{description} ({entry}) | {match.group()}', paragraph
    return None


def accession_state(article, entry):
    """Read exact accession associations in a table or deposition statement."""
    exact = re.escape(entry)
    locations = {node: i for i, node in enumerate(article.iter(), 1)}
    parents = {child: parent for parent in locations for child in parent}

    def location(node):
        address = {'element': node.tag, 'document_element_number': locations[node]}
        while node is not None:
            if node.get('id'):
                address['anchor_id'] = node.get('id')
                break
            node = parents.get(node)
        return address

    for table in article.iter('table'):
        for row in table.iter('tr'):
            cells = [node_text(c) for c in row if c.tag in {'td', 'th'}]
            if ([c.upper() for c in cells if PDB_ID.fullmatch(c)] == [entry]
                    and any(STARTS_APO.match(c) and explicitly_unbound(c) for c in cells)
                    and not any(BOUND_STATE.search(c) or ACTIVE_SITE_INHIBITED.search(c) for c in cells)
                    and len(PDB_ID.findall(' '.join(cells))) == 1):
                return 'ARTICLE_ACCESSION_TABLE', ' | '.join(cells), location(row)
        if proof := accession_table_column(table, entry):
            excerpt, label, row = proof
            return 'ARTICLE_ACCESSION_COLUMN', excerpt, {'state': location(label), 'accession': location(row)}
    paragraphs = [(p, node_text(p)) for container in (article.find('body'), article.find('back'))
                  if container is not None for p in container.iter('p')]
    for paragraph, text in paragraphs:
        if proof := author_deposition_state(text, entry, paragraphs):
            basis, excerpt, method = proof
            return basis, excerpt, ({'deposition': location(paragraph), 'methods': location(method)}
                                    if method is not None else location(paragraph))
        if not re.search(rf'\b{exact}\b', text, re.I) or not APO_STATE.search(text):
            continue
        state = re.search(rf'\b{exact}\b\s*\(([^()]{{1,140}})\)', text, re.I)
        if state and STARTS_APO.match(state.group(1)) and explicitly_unbound(state.group(1)):
            return 'ARTICLE_ACCESSION_STATE', state.group(), location(paragraph)
        if not re.search(r'\brespectively\b', text, re.I):
            continue
        ordered = re.search(
            r'\bstructures? of (.+?)\s+(?:have been|were)\s+(?:submitted|deposited)\b.+?'
            r'\b(?:identifier|accession|ID)\s*(?:codes?|numbers?)?\s*(?:are|were|:)\s*([^:.]+)',
            text, re.I)
        if ordered:
            states = re.split(r'\s*,\s*|\s+and\s+', ordered.group(1))
            ids = [p.upper() for p in PDB_ID.findall(ordered.group(2))]
            if len(ids) == len(states) and ids.count(entry) == 1 and explicitly_unbound(states[ids.index(entry)]):
                return 'ARTICLE_ORDERED_DEPOSITION', states[ids.index(entry)] + ' -> ' + entry, location(paragraph)
    for wrap in article.iter('table-wrap'):
        caption = wrap.find('caption')
        if caption is None or not re.search(r'crystallograph|diffraction', node_text(caption), re.I):
            continue
        note = next((n for n in wrap.iter('fn') if explicitly_unbound(node_text(n))
                     and re.search(r'\b(?:collect\w*|determin\w*|diffraction|solv\w*)\b', node_text(n), re.I)), None)
        if note is None:
            continue
        rows = list(wrap.iter('tr'))
        for head, row in zip(rows, rows[1:]):
            cells = [c for c in head if c.tag in {'td', 'th'}], [c for c in row if c.tag in {'td', 'th'}]
            if any(cell.get(a, '1') != '1' for group in cells for cell in group for a in ('colspan', 'rowspan')):
                continue
            headers, codes = [[node_text(c) for c in group] for group in cells]
            if (len(headers) != len(codes) or len(codes) < 3 or not
                    re.match(r'^(?:Protein Data Bank|PDB)\s+(?:entry|code|ID|accession)\b', codes[0], re.I)):
                continue
            ids = [c.upper() for c in codes[1:]]
            if ids.count(entry) != 1 or not all(PDB_ID.fullmatch(code) for code in ids):
                continue
            bound = re.compile(r'\+|\b(?:bound|complex(?:ed)? with|liganded)\b', re.I)
            unbound = [i for i, name in enumerate(headers[1:]) if name and not bound.search(name)]
            if unbound == [ids.index(entry)] and all(headers[1:]):
                excerpt = f'{entry} | {headers[unbound[0] + 1]} | {node_text(note)}'
                return 'ARTICLE_APO_CRYSTAL_COMPARISON', excerpt, {'accession': location(row), 'note': location(note)}
    return None


def accession_table_column(table, entry):
    """Resolve an explicit experimental label in the same table column as the accession."""
    rows = [[cell for cell in row if cell.tag in {'td', 'th'}] for row in table.iter('tr')]
    for headers, codes in zip(rows, rows[1:]):
        if (len(headers) != len(codes) or len(codes) < 3 or
                any(cell.get(a, '1') != '1' for cell in headers + codes for a in ('colspan', 'rowspan'))):
            continue
        values = [node_text(cell) for cell in codes]
        if not re.fullmatch(r'(?:Protein Data Bank|PDB)(?:\s+(?:entry|code|ID|accession)s?)?\s*[:.]?',
                            values[0], re.I):
            continue
        ids = [code.upper() for code in values[1:]]
        if ids.count(entry) != 1 or not all(PDB_ID.fullmatch(code) for code in ids):
            continue
        label = node_text(headers[1 + ids.index(entry)])
        if explicitly_unbound(label):
            return f'{values[0]}: {entry} | {label}', headers[1 + ids.index(entry)], codes[1 + ids.index(entry)]
    return None


def publisher_text_link(work):
    for link in work.get('link', []):
        if link.get('intended-application') != 'text-mining':
            continue
        parts = urlsplit(link.get('URL') or '')
        if parts.scheme in {'http', 'https'} and parts.netloc.lower() in PUBLISHER_TEXT_HOSTS:
            return urlunsplit(parts._replace(scheme='https', fragment=''))
    return None


def iucr_issue_url(link):
    if link is None:
        return None
    parts = urlsplit(link)
    if parts.netloc != 'journals.iucr.org':
        return None
    issue = re.fullmatch(r'(/[a-z]+/issues/\d{4}/\d{2}/\d{2}/)[a-z0-9]+/[a-z0-9]+\.pdf', parts.path, re.I)
    return urljoin(link, issue.group(1) + 'index.html') if issue else None


def iucr_article_metadata_url(registered):
    if iucr_issue_url(registered) is None:
        return None
    article_code, filename = urlsplit(registered).path.rsplit('/', 2)[-2:]
    if filename.casefold() != (article_code + '.pdf').casefold():
        return None
    return urljoin(registered, '/paper?buy=yes&cnor=' + quote(article_code))


def archive_capture(availability, original):
    snapshot = availability.get('archived_snapshots', {}).get('closest') or {}
    timestamp = snapshot.get('timestamp') or ''
    if snapshot.get('available') is not True or str(snapshot.get('status')) != '200' or not re.fullmatch(r'\d{14}', timestamp):
        return None
    archive = urlsplit(snapshot.get('url', ''))
    prefix = '/web/' + timestamp + '/'
    if archive.hostname != 'web.archive.org' or not archive.path.startswith(prefix):
        return None
    original_capture = archive.path[len(prefix):]
    reference = urlsplit(original)
    captured = urlsplit(original_capture)
    if (captured.scheme not in {'http', 'https'} or captured.hostname != reference.hostname or captured.path != reference.path
            or captured.port not in {None, 80 if captured.scheme == 'http' else 443}
            or captured.username or archive.query or archive.fragment):
        return None
    archive_url = 'https://web.archive.org/web/' + timestamp + 'id_/' + original_capture
    return archive_url, timestamp


class IssueArticleLinks(HTMLParser):
    """Retain complete links and the adjacent text needed for PDB labels."""

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.links = []
        self.text = []
        self.href = None
        self.label = []

    def handle_starttag(self, tag, attrs):
        if tag == 'a':
            self.href = dict(attrs).get('href')
            self.label = []
        elif tag in {'article', 'h1', 'h2', 'h3', 'h4', 'h5', 'h6'} and self.href is None:
            self.links.append(('', '', '', self.getpos()))
            self.text = []

    def handle_data(self, data):
        (self.label if self.href is not None else self.text).append(data)

    def handle_endtag(self, tag):
        if tag == 'a' and self.href is not None:
            self.links.append((''.join(self.text), self.href, ''.join(self.label), self.getpos()))
            self.text, self.href = [], None


def issue_accession_state(payload, issue_url, work, registered_text, entry, metadata=False):
    parser = IssueArticleLinks()
    parser.feed(payload.decode('utf8'))
    primary = urlsplit(work.get('resource', {}).get('primary', {}).get('URL', ''))
    registered_article = urlsplit(registered_text).path.rsplit('/', 1)[0] + '/'
    doi = work['DOI'].lower()
    is_article, pdb_references = False, False
    matches = []
    if metadata:
        # Require the primary article identifier on the publisher's own standalone metadata page.
        linked_dois = {urlsplit(href).path[1:].lower() for _, href, _, _ in parser.links
                       if urlsplit(href).netloc in {'doi.org', 'dx.doi.org'} and urlsplit(href).path.startswith('/10.')}
        if linked_dois != {doi}:
            return None
        is_article = True
    for preceding, href, label, (line, _) in parser.links:
        if not href:
            pdb_references = False
            if not metadata:
                is_article = False
            continue
        link = urlsplit(urljoin(issue_url, href))
        if link.netloc in {'doi.org', 'dx.doi.org'} and link.path.startswith('/10.'):
            is_article = link.path[1:].lower() == doi
        elif not metadata and link.netloc == 'journals.iucr.org' and link.path == '/paper' and link.query:
            is_article = link.path == primary.path and link.query.lower() == primary.query.lower()
        elif not metadata and link.netloc == 'journals.iucr.org' and '/issues/' in link.path and distinct_article_dir(link.path, issue_url):
            is_article = link.path.startswith(registered_article)
        else:
            marker = re.search(r'\bPDB references\s*:\s*(.*)$', ' '.join(preceding.split()), re.I)
            if marker:
                pdb_references = is_article
                group = marker.group(1)
            else:
                group = preceding
            accession = label.strip()
            state, _, linked_code = accession.rpartition(',')
            if state and PDB_ID.fullmatch(linked_code.strip()):
                accession = linked_code.strip()
                group = state
            text = ' '.join(group.split()).strip(' ,;')
            if not pdb_references or not pdb_archive_link(link, accession):
                pdb_references = False
                continue
            if accession.upper() == entry:
                matches.append((text, {'line': line, 'article_doi': doi, 'link': link.geturl(),
                                       'heading': 'PDB references'}))
            continue
        pdb_references = False
    if len(matches) == 1:
        text, location = matches[0]
        if len(text) <= 140 and explicitly_unbound(text) and not PDB_ID.search(text):
            basis = 'PUBLISHER_ARTICLE_PDB_REFERENCE' if metadata else 'PUBLISHER_ISSUE_PDB_REFERENCE'
            return basis, f'{text} | {entry}', location
    return None


def distinct_article_dir(path, issue_url):
    directory = urlsplit(issue_url).path.rsplit('/', 1)[0]
    return bool(re.fullmatch(re.escape(directory) + r'/[a-z]{2}\d+[a-z0-9]*/[^/]*', path, re.I))


def pdb_archive_link(link, label):
    if not PDB_ID.fullmatch(label):
        return False
    host, path = link.netloc.lower(), link.path.rstrip('/').upper()
    if host == 'scripts.iucr.org' and path == '/CGI-BIN/CR.CGI':
        query = parse_qs(link.query)
        return [v.upper() for v in query.get('pdbId', [])] == [label.upper()] and query.get('rm') == ['pdb']
    if host in {'www.rcsb.org', 'rcsb.org'}:
        return path == '/STRUCTURE/' + label.upper() or any(
            value.upper() == label.upper() for key, values in parse_qs(link.query).items()
            if key.lower() in {'structureid', 'pdbid'} for value in values)
    if host in {'www.ebi.ac.uk', 'ebi.ac.uk', 'pdbe.org'}:
        return path in {'/PDBE/ENTRY/PDB/' + label.upper(), '/' + label.upper()}
    if host in {'pdbj.org', 'www.pdbj.org'}:
        return path == '/MINE/SUMMARY/' + label.upper()
    return False


def pdf_accession_state(text, entry):
    """Read exact IDs in aligned table rows or in explicit deposition statements."""
    for page, page_text in enumerate(text.split('\f'), 1):
        rows = page_text.splitlines()
        for line_number, line in enumerate(rows, 1):
            if not re.search(rf'\b{re.escape(entry)}\b', line, re.I):
                continue
            cells = re.split(r'\s{2,}', line.strip())
            if ([c.upper() for c in cells if PDB_ID.fullmatch(c)] == [entry]
                    and any(explicitly_unbound(c) for c in cells)
                    and not any(BOUND_STATE.search(c) or ACTIVE_SITE_INHIBITED.search(c) for c in cells)
                    and any(re.search(r'\b(?:PDB|Protein Data Bank)\b', row, re.I)
                            for row in rows[max(0, line_number - 16):line_number])):
                return ('PUBLISHER_ACCESSION_TABLE_ROW', ' | '.join(cells),
                        f'page {page}, row {line_number}')
            if not re.match(r'^\s*(?:Protein Data Bank|PDB)\s+(?:entry|code|ID|accession)s?\b', line, re.I):
                continue
            ids = list(PDB_ID.finditer(line))
            target = [i for i, match in enumerate(ids) if match.group().upper() == entry]
            if len(ids) < 2 or len(target) != 1 or line_number == 1:
                continue
            previous = rows[line_number - 2]
            headers = list(re.finditer(r'\S(?:.*?\S)?(?=\s{2,}|$)', previous.rstrip()))
            if len(headers) != len(ids):
                continue
            if any(abs((m.start() + m.end()) - (n.start() + n.end())) > 5 for m, n in zip(ids, headers)):
                continue
            label = headers[target[0]].group()
            if explicitly_unbound(label):
                return ('PUBLISHER_ACCESSION_TABLE_COLUMN', f'{entry} | {label}',
                        f'page {page}, row {line_number}')
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


def same_primary_study(left, right):
    fields = ('pdbx_database_id_PubMed', 'pdbx_database_id_DOI')
    matches = [str(left[k]).strip().casefold() == str(right[k]).strip().casefold()
               for k in fields if left.get(k) and right.get(k)]
    return bool(matches) and all(matches)


def wait_for_ncbi(shared_cache):
    if shared_cache is None or os.name != 'posix':
        time.sleep(0.35)
        return
    import fcntl

    shared_cache.mkdir(parents=True, exist_ok=True)
    with (shared_cache / 'ncbi-request-rate').open('a+', encoding='ascii') as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.seek(0)
        previous = float(handle.read() or 0)
        elapsed = time.monotonic() - previous
        if 0 <= elapsed < 0.35:
            time.sleep(0.35 - elapsed)
        handle.seek(0)
        handle.truncate()
        handle.write(str(time.monotonic()))
        handle.flush()


class AuthorDocuments:
    """Index explicit experimental groups in archived publisher files once per process."""
    def __init__(self, manifest):
        self.by_entry = {}
        manifest = Path(manifest)
        for doc in json.loads(manifest.read_text(encoding='utf8'))['documents']:
            path = Path(doc['path'])
            payload = path.read_bytes()
            url = urlsplit(doc['url'])
            if (url.scheme != 'https' or url.hostname != 'journals.iucr.org' or not doc['doi'].startswith('10.1107/')
                    or sha256(payload).hexdigest() != doc['sha256'] or not payload.startswith(b'%PDF-')):
                raise ValueError('unverified author supplement')
            layout = pdf_layout(payload)
            if layout is None:
                raise ValueError('author supplement could not be read')
            for page, text in enumerate(layout.split('\f'), 1):
                for paragraph in re.split(r'(?=^\s*(?:Figure|Table)\s+S?\d+\b)', text, flags=re.M):
                    caption = ' '.join(paragraph.split())
                    figure = re.match(r'^(Figure\s+S?\d+)\b', caption, re.I)
                    if not figure or not re.search(r'(?<![a-z])in (?:the |an? )?(?:unligated|unliganded|apo|unbound) (?:state|form) for\b', caption, re.I):
                        continue
                    group = re.search(r'\bstructures (?:represented|shown) are PDB (?:codes|entries)\s+([^.;]{4,330})(?:\.|;)', caption, re.I)
                    if not group or BACKGROUND_STATE.search(caption[:group.end()]) or BOUND_STATE.search(caption[:group.end()]):
                        continue
                    codes = [m.group().upper() for m in PDB_ID.finditer(group[1])]
                    if not 2 <= len(codes) <= 20 or len(codes) != len(set(codes)):
                        continue
                    # Only one unbound group, with accession/color labels and no competing arm.
                    stripped = re.sub(r'\([A-Za-z -]{2,24}\)', '', PDB_ID.sub('', group[1]))
                    if not re.fullmatch(r'[\s,]*(?:and[\s,]*)?', stripped, re.I):
                        continue
                    for code in codes:
                        self.by_entry.setdefault(code, []).append({'basis':'AUTHOR_SUPPLEMENT_ENUMERATES_UNBOUND_STRUCTURES',
                            'entry':code, 'doi':doc['doi'], 'source_receipt':dict(doc), 'excerpt':caption[:group.end()],
                            'location':{'page':page,'figure':figure[1], 'accessions':codes}})

    def support(self, entry):
        return next(iter(self.by_entry.get(entry, [])), None)


class PrimaryLiterature:
    """Reuse verified responses across Runs and read full text only for unresolved evidence."""

    def __init__(self, run, shared_cache=None, prior_cache=None):
        self.run = run
        self.shared_cache = Path(shared_cache) if shared_cache is not None else None
        self.prior_cache = Path(prior_cache) if prior_cache is not None else None
        self.documents = {}
        self.articles = {}
        self.stats = Counter()

    def get(self, name, url, params, accepted=(200,)):
        key = sha256(json.dumps([url, params], sort_keys=True).encode()).hexdigest()
        shared = self.shared_cache / (key + '.json.gz') if self.shared_cache is not None else None
        for directory, source in ((self.shared_cache, 'shared'), (self.prior_cache, 'prior_cache')):
            if directory is None:
                continue
            cache = directory / (key + '.json.gz')
            if not cache.is_file():
                continue
            try:
                with gzip.open(cache, 'rt', encoding='utf8') as stream:
                    cached = json.load(stream)
                body, receipt = b64decode(cached['payload'], validate=True), cached['receipt']
                if receipt['status'] in accepted and sha256(body).hexdigest() == receipt['response_sha256']:
                    self.stats[source + '_hits'] += 1
                    return body, receipt | {'literature_cache': source}
            except (OSError, ValueError, KeyError):
                self.stats[source + '_read_errors'] += 1

        def validate(body, status):
            if status == 200:
                if name.startswith(('doi/', 'accession/')):
                    if 'esearchresult' not in json.loads(body):
                        raise ValueError('invalid literature search response')
                elif name.startswith('related/'):
                    if json.loads(body).get('rcsb_id') != name.split('/', 1)[1]:
                        raise ValueError('related deposit identity mismatch')
                elif name.startswith('crossref/'):
                    if 'message' not in json.loads(body):
                        raise ValueError('invalid Crossref work response')
                elif name.startswith('publisher/'):
                    if not body.startswith(b'%PDF-') and xml_document(body).tag not in {'article', 'article-set'}:
                        raise ValueError('publisher did not provide PDF or JATS')
                elif name.startswith(('publisher_issue/', 'publisher_metadata/', 'archive_capture/')):
                    return  # The parser requires an identity-bound article section and archive links.
                elif name.startswith('archive_available/'):
                    if 'archived_snapshots' not in json.loads(body):
                        raise ValueError('invalid archive availability response')
                elif name.startswith('biostudies/'):
                    if json.loads(body).get('accno') != name.split('/', 1)[1]:
                        raise ValueError('BioStudies accession mismatch')
                elif name.startswith('supplement/'):
                    if not body.startswith(b'%PDF-'):
                        raise ValueError('BioStudies supplement is not PDF')
                elif xml_document(body).tag not in {'PubmedArticleSet', 'OAI-PMH'}:
                    raise ValueError('invalid literature record response')

        local = self.run / 'http' / 'literature' / name
        reused = local.with_suffix('.gz').is_file() and local.with_suffix('.json').is_file()
        if not reused:
            if url.startswith((EUTILS, PMC_OAI)):
                wait_for_ncbi(self.shared_cache)
            self.stats['network_requests'] += 1
            if name.startswith(('pmc/', 'publisher/', 'supplement/')):
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
        record = {'status': 'PRIMARY_CITATION_UNAVAILABLE', 'pmid': pmid, 'doi': doi}
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
        root = xml_document(payload)
        error = root.find('./error')
        if error is not None and error.get('code'):
            receipt = receipt | {'oai_error': error.get('code')}
        article = root.find('.//record/metadata/article') if receipt['status'] == 200 else None
        if article is not None:
            ids = {i.text.removeprefix('PMC') for i in article.findall('front/article-meta/article-id')
                   if i.text and i.get('pub-id-type') in {'pmc', 'pmcid'}}
            if pmc[3:] not in ids:
                article = None
        self.articles[pmc] = article, receipt
        return article, receipt

    def archived_issue(self, page_url, work):
        receipts, tried = [], set()
        date = (work.get('issued', {}).get('date-parts') or [[]])[0]
        timestamp_hints = [None]
        if len(date) >= 2 and all(isinstance(part, int) for part in date[:3]):
            timestamp_hints.insert(0, f'{date[0]:04d}{date[1]:02d}{(date[2] if len(date) >= 3 else 1):02d}')
        for hint in timestamp_hints:
            params = {'url': page_url} | ({'timestamp': hint} if hint else {})
            key = sha256(json.dumps(params, sort_keys=True).encode()).hexdigest()
            try:
                payload, receipt = self.get('archive_available/' + key, ARCHIVE_AVAILABILITY, params,
                                            accepted=(200, 403, 404, 429))
                receipts.append({'source': 'WEB_ARCHIVE_AVAILABILITY', **receipt})
                capture = archive_capture(json.loads(payload), page_url) if receipt['status'] == 200 else None
                if capture is None or capture[0] in tried:
                    continue
                url, timestamp = capture
                tried.add(url)
                key = sha256(url.encode()).hexdigest()
                body, receipt = self.get('archive_capture/' + key, url, {}, accepted=(200, 403, 404, 429))
                receipt = {'source': 'ARCHIVED_PUBLISHER_ISSUE_PDB_REFERENCE', 'publisher_url': page_url,
                           'capture_timestamp': timestamp, 'capture_url': url, **receipt}
                receipts.append(receipt)
                if receipt['status'] == 200 and receipt.get('url') == url:
                    return body, receipt, receipts
            except (requests.RequestException, RuntimeError, ValueError) as exc:
                receipts.append({'source': 'ARCHIVED_PUBLISHER_ISSUE_PDB_REFERENCE', 'publisher_url': page_url,
                                 'fetch_error': type(exc).__name__})
        return None, None, receipts

    def publisher_article(self, doi):
        if doi in self.articles:
            return self.articles[doi]
        key = sha256(doi.encode()).hexdigest()
        receipts = []
        article, text, issue_pages = None, None, []
        try:
            payload, receipt = self.get('crossref/' + key, CROSSREF_WORK + quote(doi, safe='/'), {},
                                        accepted=(200, 404, 429))
            receipts.append({'source': 'CROSSREF_WORK', **receipt})
            if receipt['status'] == 200:
                work = json.loads(payload)['message']
                link = publisher_text_link(work) if str(work.get('DOI', '')).lower() == doi else None
                if link:
                    for page_url, source, metadata in (
                        (iucr_issue_url(link), 'PUBLISHER_ISSUE_PDB_REFERENCE', False),
                        (iucr_article_metadata_url(link), 'PUBLISHER_ARTICLE_PDB_REFERENCE', True),
                    ):
                        if page_url is None:
                            continue
                        try:
                            cache_kind = 'publisher_metadata/' if metadata else 'publisher_issue/'
                            issue_body, issue_receipt = self.get(cache_kind + sha256(page_url.encode()).hexdigest(),
                                page_url, {}, accepted=(200, 401, 403, 404, 429, 451))
                            issue_receipt = {'source': source, **issue_receipt}
                            if issue_receipt['status'] == 200:
                                issue_pages.append((issue_body, page_url, work, link, issue_receipt, metadata))
                            else:
                                issue_receipt['document_error'] = 'PUBLISHER_METADATA_HTTP_UNAVAILABLE'
                        except (requests.RequestException, RuntimeError, ValueError) as exc:
                            issue_receipt = {'source': source, 'url': page_url,
                                             'fetch_error': type(exc).__name__}
                        receipts.append(issue_receipt)
                        if not metadata and issue_receipt.get('status') in {403, 404}:
                            archive_body, archive_source, archive_attempts = self.archived_issue(page_url, work)
                            receipts.extend(archive_attempts)
                            if archive_body is not None:
                                issue_pages.append((archive_body, page_url, work, link, archive_source, False))

                    payload, receipt = self.get('publisher/' + key, link, {},
                                                accepted=(200, 401, 403, 404, 429, 451))
                    receipts.append({'source': 'PUBLISHER_REGISTERED_TEXT_MINING', **receipt})
                    if receipt['status'] == 200:
                        if payload.startswith(b'%PDF-'):
                            if converted := pdf_layout(payload):
                                front = '\n'.join(converted.split('\f', 1)[0].splitlines()[:40])
                                if re.search(r'\b(?:doi\s*:\s*|doi\.org/)' + re.escape(doi) + r'(?=\s|$|[.,])', front, re.I):
                                    text = converted
                                else:
                                    receipts[-1]['document_error'] = 'PUBLISHER_PDF_IDENTITY_UNRESOLVED'
                            else:
                                receipts[-1]['document_error'] = 'PDFTOTEXT_UNAVAILABLE'
                        else:
                            root = xml_document(payload)
                            matches = [a for a in root.iter('article') if any(
                                (i.text or '').strip().lower() == doi and i.get('pub-id-type') == 'doi'
                                for i in a.findall('front/article-meta/article-id'))]
                            if len(matches) == 1:
                                article = matches[0]
                            else:
                                receipts[-1]['document_error'] = 'PUBLISHER_JATS_IDENTITY_UNRESOLVED'
                    else:
                        receipts[-1]['document_error'] = 'PUBLISHER_HTTP_UNAVAILABLE'
                else:
                    receipts[-1]['document_error'] = 'REGISTERED_TEXT_MINING_LINK_UNAVAILABLE'
        except (requests.RequestException, RuntimeError, ValueError, ET.ParseError, KeyError, OSError,
                subprocess.SubprocessError) as exc:
            receipts.append({'source': 'PUBLISHER_REGISTERED_TEXT_MINING', 'fetch_error': type(exc).__name__})
        self.articles[doi] = (article, text, issue_pages, receipts)
        return self.articles[doi]

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

    def supplementary_groups(self, pmc, doi, pmid):
        key = 'supplement/' + pmc
        if key in self.articles:
            return self.articles[key]
        reports, groups = [], []
        accession = 'S-E' + pmc
        payload, receipt = self.get('biostudies/' + accession, BIOSTUDIES + 'api/v1/studies/' + accession, {},
                                    accepted=(200, 403, 404, 429, 502, 503))
        reports.append({'kind': 'BIOSTUDIES_STUDY', 'source_receipt': receipt})
        section = biostudies_article(json.loads(payload), pmc, doi, pmid) if receipt['status'] == 200 else None
        deposits = sorted({url.upper() for t, url in biostudies_links(section or {})
                           if t == 'PDB' and PDB_ID.fullmatch(url)})
        if section is not None and 2 <= len(deposits) <= 30:
            files = [f for f in section.get('files', []) if re.fullmatch(r'[a-z0-9_][a-z0-9_.-]*\.pdf',
                      f.get('path', ''), re.I) and 0 < int(f.get('size') or 0) <= 4_000_000]
            for file in files if len(files) <= 8 else []:
                name = file['path']
                url = BIOSTUDIES + 'files/' + accession + '/' + quote(name)
                key_name = 'supplement/' + sha256(url.encode()).hexdigest()
                try:
                    pdf, source = self.get(key_name, url, {'download': 'true'}, accepted=(200, 403, 404, 429))
                    reports.append({'kind': 'BIOSTUDIES_AUTHOR_SUPPLEMENT', 'source_receipt': source})
                    if source['status'] == 200 and len(pdf) == int(file['size']):
                        text = pdf_layout(pdf)
                        if text is None:
                            reports[-1]['document_error'] = 'PDFTOTEXT_UNAVAILABLE'
                            break
                        groups.extend((fields, location, source) for fields, location in crystal_apo_groups(text))
                except (requests.RequestException, RuntimeError, ValueError, OSError, subprocess.SubprocessError) as exc:
                    reports.append({'kind': 'BIOSTUDIES_AUTHOR_SUPPLEMENT', 'url': url,
                                    'fetch_error': type(exc).__name__})
        self.articles[key] = groups, deposits, reports
        return self.articles[key]

    def supplement_support(self, source, experiment, pmc, fixed_fingerprint, attempts):
        if not all((pmc, source.get('pmid'), source.get('doi'), fixed_fingerprint)):
            return None
        groups, deposits, reports = self.supplementary_groups(pmc, source['doi'], str(source['pmid']))
        attempts.extend({**source, **report} for report in reports)
        entry = source['entry']
        named_related = {link['entry'] for link in experiment.get('related_deposits') or []}
        if not groups or entry not in deposits or not named_related.issubset(deposits):
            return None
        candidates = {accession: fixed_fingerprint(accession) for accession in deposits}
        if not all(data and data['entry'] == accession and
                   same_primary_study(experiment['primary_citation'], data['primary_citation'])
                   for accession, data in candidates.items()):
            return None
        proofs = []
        for fields, location, receipt in groups:
            verdicts = {accession: same_crystal(fields, data['fields']) for accession, data in candidates.items()}
            if verdicts[entry] is True and all(verdicts[peer] is False for peer in deposits if peer != entry):
                proofs.append({'basis': 'AUTHOR_SUPPLEMENT_UNIQUE_CRYSTAL_FINGERPRINT',
                               'location': location, 'source_receipt': receipt, 'article_measurements': fields,
                               'fixed_depositions': candidates, 'study_source': reports[0]['source_receipt'],
                               'excerpt': location['state'] + ' ' + location['group'] + ' crystallographic statistics'})
        if len(proofs) == 1:
            return source | proofs[0]
        return None

    def related_deposit_support(self, experiment, entry):
        citation = experiment['primary_citation']
        related = experiment.get('related_deposits') or []
        if not related or not any(citation.get(k) for k in ('pdbx_database_id_PubMed', 'pdbx_database_id_DOI')):
            return None
        for link in sorted(related, key=lambda r: not deposited_apo(r['details']))[:5]:
            source = link['entry']
            if source == entry or not PDB_ID.fullmatch(source):
                continue
            payload, receipt = self.get('related/' + source, PDB_ENTRY + source, {}, accepted=(200, 404))
            if receipt['status'] != 200:
                continue
            observed = json.loads(payload)
            if not same_primary_study(citation, observed.get('rcsb_primary_citation') or {}):
                continue
            notes = [str(row.get('details') or '') for row in observed.get('pdbx_database_related') or []
                     if str(row.get('db_name') or '').upper() == 'PDB'
                     and str(row.get('db_id') or '').upper() == entry and row.get('details')]
            if notes and all(explicitly_unbound(note) for note in notes):
                return {'basis': 'RECIPROCAL_PDB_DEPOSITION_SAME_PRIMARY_STUDY', 'entry': entry,
                        'related_pdb': source, 'pmid': citation.get('pdbx_database_id_PubMed'),
                        'doi': citation.get('pdbx_database_id_DOI'), 'excerpt': notes[0][:600],
                        'source_receipt': receipt}
        return None

    def assess(self, experiment, *, fixed_fingerprint=None, source_deposit=None):
        entry = experiment['source_entry'].upper()
        attempts = []
        outcome = {'status': 'PRIMARY_CITATION_UNAVAILABLE', 'source_attempts': attempts}
        checked = set()
        issue_support, supplement_source = None, None
        try:
            for (record, abstract, pmc), candidate_linked in self.candidate_sources(experiment, entry):
                identity = record.get('pmid'), record.get('doi')
                if identity in checked or not any(identity):
                    if identity not in checked:
                        attempts.append({'kind': 'CITATION', 'candidate_linked': candidate_linked, **record})
                    continue
                checked.add(identity)
                attempts.append({'kind': 'CITATION', 'candidate_linked': candidate_linked, **record})
                if candidate_linked:
                    outcome.update(record)
                text = html.unescape(re.sub(r'<[^>]+>', '', abstract))
                conflict = accession_conflict(text, entry)
                sentence = reported_apo_structure(text)
                source = {'entry': entry, 'pmid': record.get('pmid'), 'doi': record.get('doi'),
                          'candidate_linked': candidate_linked}
                if conflict:
                    return outcome | {'conflict': source | {'basis': 'ACCESSION_ACTIVE_SITE_INHIBITED', 'excerpt': conflict}}
                if candidate_linked and explicitly_unbound(experiment['title']) and sentence:
                    return outcome | {'support': source | {'basis': 'MATCHED_APO_DEPOSITION_AND_PRIMARY_ABSTRACT',
                                                           'location': 'PubMed abstract + deposition title',
                                                           'excerpt': sentence[:600]}}
                article, pdf = None, None
                full_source = {}
                if pmc:
                    full_source = {'pmcid': pmc}
                    try:
                        article, receipt = self.full_article(pmc)
                        full_source['source_receipt'] = receipt
                        doi = source['doi']
                        if article is not None and doi:
                            article_dois = {(i.text or '').strip().lower() for i in article.findall(
                                'front/article-meta/article-id') if i.get('pub-id-type') == 'doi'}
                            if article_dois and article_dois != {doi}:
                                article = None
                                receipt = receipt | {'document_error': 'PMC_ARTICLE_DOI_MISMATCH'}
                                full_source['source_receipt'] = receipt
                        if candidate_linked:
                            outcome['full_text_source'] = receipt
                    except (requests.RequestException, RuntimeError, ValueError, ET.ParseError, KeyError) as exc:
                        full_source['fetch_error'] = type(exc).__name__
                    attempts.append({'kind': 'PMC_OAI', **source, **full_source})
                doi = source['doi']
                if candidate_linked and pmc and doi:
                    supplement_source = source, pmc
                if article is None and candidate_linked and doi and doi.startswith('10.'):
                    article, pdf, issue_pages, receipts = self.publisher_article(doi)
                    attempts.extend({'kind': attempt['source'], **source, 'source_receipt': attempt}
                                    for attempt in receipts)
                    full_source = {'source_receipt': receipts[-1]}
                    for body, page_url, work, link, receipt, metadata in issue_pages or []:
                        if precise := issue_accession_state(body, page_url, work, link, entry, metadata=metadata):
                            basis, excerpt, location = precise
                            issue_support = source | {'source_receipt': receipt, 'basis': basis,
                                                    'excerpt': excerpt, 'location': location}
                if pdf:
                    conflict = next((issue for line in pdf.splitlines()
                                     if (issue := accession_conflict(line, entry))), None)
                    if conflict:
                        return outcome | {'conflict': source | full_source | {
                            'basis': 'ACCESSION_ACTIVE_SITE_INHIBITED', 'excerpt': conflict}}
                    if precise := pdf_accession_state(pdf, entry):
                        basis, excerpt, location = precise
                        return outcome | {'support': source | full_source | {
                            'basis': basis, 'excerpt': excerpt[:600], 'location': location}}
                if article is None:
                    continue
                conflict = next((issue for p in article.iter('p')
                                 if (issue := accession_conflict(node_text(p), entry))), None)
                if conflict:
                    return outcome | {'conflict': source | full_source | {
                        'basis': 'ACCESSION_ACTIVE_SITE_INHIBITED', 'excerpt': conflict}}
                precise = accession_state(article, entry)
                if precise:
                    basis, excerpt, location = precise
                    return outcome | {'support': source | full_source | {
                        'basis': basis, 'excerpt': excerpt[:600], 'location': location}}
        except (requests.RequestException, RuntimeError, ValueError, ET.ParseError, KeyError) as exc:
            outcome['fetch_error'] = type(exc).__name__
        if support := fixed_relation_support(experiment, source_deposit):
            return outcome | {'support': support}
        documents = getattr(self, 'author_documents', None)
        if (documents and not BOUND_STATE.search(experiment['title']) and not ACTIVE_SITE_INHIBITED.search(experiment['title'])
                and not experiment.get('preparation_conflicts_or_unresolved_additions') and (support := documents.support(entry))):
            return outcome | {'support': support}
        if issue_support is not None:
            return outcome | {'support': issue_support}
        if supplement_source is not None and fixed_fingerprint is not None:
            source, pmc = supplement_source
            try:
                if support := self.supplement_support(source, experiment, pmc, fixed_fingerprint, attempts):
                    return outcome | {'support': support}
            except (requests.RequestException, RuntimeError, ValueError, KeyError, TypeError, OSError) as exc:
                attempts.append({'kind': 'BIOSTUDIES_STUDY', **source, 'fetch_error': type(exc).__name__})
        try:
            if support := self.related_deposit_support(experiment, entry):
                return outcome | {'support': support}
        except (requests.RequestException, RuntimeError, ValueError, KeyError) as exc:
            outcome['related_deposit_fetch_error'] = type(exc).__name__
        return outcome

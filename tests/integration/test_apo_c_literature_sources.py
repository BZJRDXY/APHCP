from pathlib import Path
import sys
import unittest
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from apo_c_literature import (PrimaryLiterature, accession_state, archive_capture, crystal_apo_groups, issue_accession_state,
                              iucr_article_metadata_url, iucr_issue_url, same_crystal)


class LiteratureSourceTests(unittest.TestCase):
    def test_xml_receipt_identifies_supporting_row(self):
        root = ET.fromstring('''<article><body>
          <p id='incidental'>Previously submitted 9XYZ</p>
          <table-wrap id='wrong'><table><tr><td>Apo or liganded receptor</td><td>9XYZ</td></tr></table></table-wrap>
          <table-wrap id='experiment'><table><tr><td>Unliganded receptor</td><td>9XYZ</td></tr></table></table-wrap>
          </body></article>''')
        basis, _, location = accession_state(root, '9XYZ')
        self.assertEqual(basis, 'ARTICLE_ACCESSION_TABLE')
        self.assertEqual(location['anchor_id'], 'experiment')

    def test_issue_requires_the_candidate_article_group(self):
        journal = 'https://journals.iucr.org'
        registered = journal + '/d/issues/2020/05/00/ab1234/ab1234.pdf'
        issue_url = iucr_issue_url(registered)
        work = {'DOI': '10.0000/article'}
        other = '<h3><a href=cd5678/index.html>Other paper</a></h3>'
        paper = '<h3><a href=ab1234/index.html>Primary paper</a></h3>'
        archive = 'https://www.rcsb.org' + '/structure/9XYZ'
        reference = f'<p>PDB references: apo receptor, <a href={archive}>9xyz</a>.</p>'

        def parse(value):
            return issue_accession_state(value.encode(), issue_url, work, registered, '9XYZ')

        self.assertIsNone(parse(other + reference + paper))
        self.assertIsNone(parse(paper + reference.replace('/9XYZ', '/9XYB')))
        self.assertIsNone(parse(paper + reference.replace('apo receptor', 'apo or liganded receptor')))
        self.assertIsNone(parse(paper + '<h3>Next paper</h3>' + reference))
        evidence = parse(other + reference + paper + reference)
        self.assertEqual(evidence[0], 'PUBLISHER_ISSUE_PDB_REFERENCE')
        self.assertEqual(evidence[2]['article_doi'], work['DOI'])

    def test_publisher_inline_label_requires_doi_and_exact_link(self):
        journal = 'https://journals.iucr.org'
        registered = journal + '/d/issues/2020/05/00/ab1234/ab1234.pdf'
        metadata = iucr_article_metadata_url(registered)
        work = {'DOI': '10.0000/article'}
        doi = '<a href=https://doi.org/10.0000/article>DOI</a>'
        reference = ('<p>PDB references: <a href=https://scripts.iucr.org/cgi-bin/cr.cgi?pdbId=9XYZ&amp;rm=pdb>'
                     'apo receptor, 9xyz</a></p>')

        def parse(value):
            return issue_accession_state(value.encode(), metadata, work, registered, '9XYZ', metadata=True)

        self.assertEqual(parse(doi + reference)[0], 'PUBLISHER_ARTICLE_PDB_REFERENCE')
        self.assertIsNone(parse(doi.replace('/article', '/another') + reference))
        self.assertIsNone(parse(doi + reference.replace('pdbId=9XYZ', 'pdbId=9XYA')))

    def test_public_issue_archive_requires_matching_original(self):
        original = 'https://journals.iucr.org/d/issues/2020/05/00/index.html'
        captured = 'http://journals.iucr.org:80/d/issues/2020/05/00/index.html'
        timestamp = '20200515091213'
        prefix = 'https://web.archive.org/web/' + timestamp
        snapshot = {'available': True, 'status': '200', 'timestamp': timestamp, 'url': prefix + '/' + captured}

        def locate(url):
            return archive_capture({'archived_snapshots': {'closest': snapshot | {'url': url}}}, original)

        self.assertEqual(locate(snapshot['url']), (prefix + 'id_/' + captured, timestamp))
        self.assertIsNone(locate(snapshot['url'].replace('/05/', '/06/')))
        self.assertIsNone(locate(snapshot['url'].replace('journals.iucr.org', 'example.invalid')))
        self.assertIsNone(locate(snapshot['url'].replace(':80/', ':91/')))
        self.assertIsNone(locate(snapshot['url'].replace(timestamp, '20200514105000')))

    def test_crystal_group_comparison_requires_one_same_study_match(self):
        title = 'Table 1. Data collection and refinement statistics'
        rows = [('a, A', '41.0'), ('b, A', '75.4'), ('c, A', '42.4'), ('alpha, d', '90.0'),
                ('beta, d', '108.5'), ('gamma, d', '90.0'), ('Wavelength, A', '1.00000'),
                ('Resolution, A', '50 - 2.15'), ('Resolution range, A', '30 - 2.2')]
        lines = [title, ' ' * 48 + '(apo-form)' + ' ' * 16 + 'co-crystallized with ligand',
                 f'{"Protein":30}{"Native":^39}{"SeMet":^13}{"Native":^20}', 'Unit cell']
        lines += [f'{key:47}{value:>9}{"99.9":>20}' for key, value in rows]
        groups = list(crystal_apo_groups('\n'.join(lines)))
        fields, location = groups[0]
        exact = fields | {'a': '40.969'}
        self.assertTrue(same_crystal(fields, exact))
        self.assertFalse(same_crystal(fields, exact | {'wavelength': '0.97919'}))
        self.assertIsNone(same_crystal(fields, {key: value for key, value in exact.items() if key != 'b'}))
        self.assertEqual(location['group'], 'Native')

        reader = PrimaryLiterature(None)
        reader.supplementary_groups = lambda *_: ([(fields, location, {'status': 200})], ['9XYZ', '9XYA'],
                                                   [{'source_receipt': {'status': 200}}])
        citation = {'pdbx_database_id_DOI': '10.0000/article', 'pdbx_database_id_PubMed': '123'}
        source = {'entry': '9XYZ', 'doi': citation['pdbx_database_id_DOI'], 'pmid': '123'}
        experiment = {'primary_citation': citation, 'related_deposits': [{'entry': '9XYA'}]}
        fixed = {entry: {'entry': entry, 'fields': dict(exact), 'primary_citation': citation,
                         'source_receipt': {'sha256': 'fixture'}} for entry in ('9XYZ', '9XYA')}
        self.assertIsNone(reader.supplement_support(source, experiment, 'PMC123', fixed.get, []))
        fixed['9XYA']['fields']['refinement_high'] = '1.80'
        self.assertEqual(reader.supplement_support(source, experiment, 'PMC123', fixed.get, [])['basis'],
                         'AUTHOR_SUPPLEMENT_UNIQUE_CRYSTAL_FINGERPRINT')
        fixed.pop('9XYA')
        self.assertIsNone(reader.supplement_support(source, experiment, 'PMC123', fixed.get, []))

    def test_candidate_format_failure_survives_distinct_holo_paper(self):
        reader = PrimaryLiterature(None)
        candidate = {'pdbx_database_id_PubMed': '123'}
        holo = {'pdbx_database_id_PubMed': '456'}
        reader.get_document = lambda cite: ({'pmid': cite['pdbx_database_id_PubMed']}, '',
                                            'PMC' + cite['pdbx_database_id_PubMed'])
        candidate_receipt = {'status': 200, 'oai_error': 'cannotDisseminateFormat'}
        holo_article = ET.fromstring('<article><body><p>We studied an apo receptor.</p></body></article>')
        reader.full_article = lambda pmc: ((None, candidate_receipt) if pmc == 'PMC123'
                                          else (holo_article, {'status': 200}))
        experiment = {'source_entry': '9XYZ', 'title': 'Receptor', 'primary_citation': candidate,
                      'holo_primary_citation': holo}
        outcome = reader.assess(experiment)
        self.assertEqual(outcome['pmid'], '123')
        self.assertEqual(outcome['full_text_source'], candidate_receipt)
        self.assertNotIn('support', outcome)
        attempts = [r for r in outcome['source_attempts'] if r['kind'] == 'PMC_OAI']
        self.assertEqual([(r['pmid'], r['candidate_linked']) for r in attempts], [('123', True), ('456', False)])

    def test_conflicting_pmc_doi_rejects_article(self):
        reader = PrimaryLiterature(None)
        doi = '10.0000/primary'
        reader.get_document = lambda _: ({'pmid': '123', 'doi': doi}, '', 'PMC123')
        article = ET.fromstring('''<article><front><article-meta>
          <article-id pub-id-type='doi'>10.0000/wrong</article-id></article-meta></front>
          <body><p>Coordinates: 9XYZ (apo receptor).</p></body></article>''')
        reader.full_article = lambda _: (article, {'status': 200})
        reader.publisher_article = lambda _: (None, None, None, [{'source': 'CROSSREF_WORK', 'status': 404}])
        outcome = reader.assess({'source_entry': '9XYZ', 'title': 'Receptor', 'primary_citation': {}})
        self.assertEqual(outcome['full_text_source']['document_error'], 'PMC_ARTICLE_DOI_MISMATCH')
        self.assertNotIn('support', outcome)

"""Check when a deposition supports a target-specific unbound state."""
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from apo_c_literature import PrimaryLiterature, accession_conflict, accession_state
from apo_c_pair import anchored_end_alignment, tagged_end_alignment
from apo_c_priority import order_rankings
from apo_c_science import complete_altlocs, experiment_record


def deposition(title, citation=None, preparation=None, doi='10.0000/control'):
    block = gemmi.cif.Document().add_new_block('test')
    block.set_pair('_entry.id', '9XYZ')
    block.set_pair('_struct.title', gemmi.cif.quote(title))
    if citation:
        row = block.init_loop('_citation.', ['id', 'title', 'pdbx_database_id_DOI'])
        row.add_row(['primary', gemmi.cif.quote(citation), doi])
    if preparation:
        block.set_pair('_exptl_crystal_grow.pdbx_details', gemmi.cif.quote(preparation))
    return block


class ExperimentScopeTests(unittest.TestCase):
    def setUp(self):
        self.target = [{'component': 'LIG'}]
        self.term = {'matched_term': 'LIG'}
        self.citation = 'Experimental structures of the apo and LIG-bound states of the receptor'
        self.holo = deposition('CRYSTAL STRUCTURE OF THE LIG COMPLEX', self.citation)

    def test_native_and_a_title_without_an_experiment_do_not_prove_absence(self):
        for title, state in [('Native receptor', 'DEPOSITION_CONTEXT'), ('Receptor structure', 'NOT_ESTABLISHED'),
                             ('Apo receptor', 'DEPOSITION_CONTEXT')]:
            with self.subTest(title=title):
                row = experiment_record(deposition(title), self.target, self.holo, self.term)
                self.assertEqual(row['status'], state)
                self.assertEqual(row['evidence'], [])

    def test_scoped_primary_experiment_and_conflicting_preparation(self):
        apo = deposition('Apo receptor', self.citation)
        outcome = experiment_record(apo, self.target, self.holo, self.term)
        self.assertEqual(outcome['status'], 'EXPLICIT_UNBOUND')
        self.assertEqual(outcome['evidence'][0]['bound_target_term'], 'LIG')
        different = deposition('Apo receptor', self.citation, doi='10.0000/unrelated')
        self.assertEqual(experiment_record(different, self.target, self.holo, self.term)['status'],
                         'DEPOSITION_CONTEXT')
        conflict = deposition('Apo receptor', self.citation, 'Co-crystallized with LIG')
        self.assertEqual(experiment_record(conflict, self.target, self.holo, self.term)['status'], 'CONFLICT')

    def test_historical_role_can_use_the_deposited_ligand_identity(self):
        apo = deposition('Apo receptor', self.citation)
        outcome = experiment_record(apo, self.target, self.holo, {'status': 'HISTORICAL_EVIDENCE'})
        self.assertEqual(outcome['status'], 'EXPLICIT_UNBOUND')

    def test_article_assigns_apo_to_the_specific_accession(self):
        doc = ET.fromstring("""
            <article><body><p>Structures of the receptor Apo, receptor-inhibitor and mutant-inhibitor
            have been deposited in the Protein Data Bank and their PDB ID codes are 9XYA, 9XYB and
            9XYC: 9XYD, 9XYE, respectively.</p></body></article>
        """)
        self.assertEqual(accession_state(doc, '9XYA')[0], 'ARTICLE_ORDERED_DEPOSITION')
        self.assertIsNone(accession_state(doc, '9XYB'))
        self.assertIsNone(accession_state(doc, '9XYE'))

    def test_article_table_with_two_accessions_does_not_infer_state(self):
        doc = ET.fromstring("""
            <article><body><table><tr><td>Apo receptor</td><td>9XYA</td><td>9XYB</td></tr>
            <tr><td>Apo receptor</td><td>9XYC</td></tr></table></body></article>
        """)
        self.assertIsNone(accession_state(doc, '9XYA'))
        self.assertEqual(accession_state(doc, '9XYC')[0], 'ARTICLE_ACCESSION_TABLE')

    def test_abstract_support_requires_the_matching_deposition_description(self):
        reader = PrimaryLiterature(None)
        reader.get_document = lambda _: (
            {'pmid': '123', 'doi': '10.0000/control'},
            'We report crystal structures of the receptor in the apo form and in complex with a ligand.', None)
        data = {'primary_citation': {}, 'source_entry': '9XYZ', 'title': 'Crystal structure of apo receptor'}
        self.assertEqual(reader.assess(data)['support']['basis'], 'MATCHED_APO_DEPOSITION_AND_PRIMARY_ABSTRACT')
        self.assertNotIn('support', reader.assess(data | {'title': 'Receptor structure'}))

    def test_pmc_oai_identifier_matches_the_primary_paper(self):
        reader = PrimaryLiterature(None)
        pubmed = b"""<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>123</PMID></MedlineCitation>
            <PubmedData><ArticleIdList><ArticleId IdType="pmc">PMC456</ArticleId></ArticleIdList>
            <ReferenceList><Reference><ArticleIdList><ArticleId IdType="pmc">PMC999</ArticleId>
            </ArticleIdList></Reference></ReferenceList></PubmedData>
            </PubmedArticle></PubmedArticleSet>"""
        full = b"""<OAI-PMH xmlns="http://www.openarchives.org/OAI/2.0/"><GetRecord><record><metadata>
            <article xmlns=""><front><article-meta><article-id pub-id-type="pmcid">PMC456</article-id>
            </article-meta></front><body><p>Coordinates: 9XYZ (apo receptor).</p></body></article>
            </metadata></record></GetRecord></OAI-PMH>"""
        def get(name, *args, **kwargs):
            self.assertIn(name, {'pubmed/123', 'pmc/PMC456'})
            return (pubmed if name.startswith('pubmed/') else full), {'status': 200}

        reader.get = get
        data = {'primary_citation': {'pdbx_database_id_PubMed': '123'},
                'source_entry': '9XYZ', 'title': 'Receptor structure'}
        self.assertEqual(reader.assess(data)['support']['basis'], 'ARTICLE_ACCESSION_STATE')

    def test_altlocs_do_not_fill_missing_atoms_from_another_conformation(self):
        atoms = [{'name': 'CA', 'altloc': ''}, {'name': 'CB', 'altloc': 'A'},
                 {'name': 'CB', 'altloc': 'B'}, {'name': 'CG', 'altloc': 'A'}]
        self.assertFalse(complete_altlocs(atoms, {'CA', 'CB', 'CG'}))
        atoms.append({'name': 'CG', 'altloc': 'B'})
        self.assertTrue(complete_altlocs(atoms, {'CA', 'CB', 'CG'}))

    def test_uninhibited_deposition_and_experimental_abstract_skip_full_text(self):
        reader = PrimaryLiterature(None)
        reader.get_document = lambda _: (
            {'status': 'PRIMARY_CITATION_READ', 'pmid': '123'},
            'We determined the crystal structure of the inhibitor-free receptor.', 'PMC456')
        reader.full_article = lambda _: self.fail('sufficient abstract should skip full text')
        evidence = reader.assess({'primary_citation': {}, 'source_entry': '9XYZ',
                                   'title': 'Crystal structure of uninhibited receptor'})
        self.assertEqual(evidence['support']['basis'], 'MATCHED_APO_DEPOSITION_AND_PRIMARY_ABSTRACT')

    def test_explicit_inhibition_belongs_to_the_accession_in_its_clause(self):
        clause = ('We solved the crystal structures of the active site-inhibited receptor (PDB 9XYZ) '
                  'and the other complex (PDB 9XYA).')
        self.assertIsNotNone(accession_conflict(clause, '9XYZ'))
        self.assertIsNone(accession_conflict(clause, '9XYA'))
        reader = PrimaryLiterature(None)
        reader.get_document = lambda _: ({'status': 'PRIMARY_CITATION_READ', 'pmid': '123'}, '', 'PMC456')
        reader.full_article = lambda _: (ET.fromstring(f'<article><body><p>{clause}</p></body></article>'),
                                          {'status': 200})
        evidence = reader.assess({'primary_citation': {}, 'source_entry': '9XYZ', 'title': 'Receptor complex'})
        self.assertIn('conflict', evidence)
        self.assertNotIn('support', evidence)

    def test_accession_link_and_crystal_table_identify_the_apo_column(self):
        reader = PrimaryLiterature(None)
        article = ET.fromstring("""
            <article><body><table-wrap><caption>Crystallographic data collection</caption>
            <table><tr><th></th><th>Receptor</th><th>Receptor + peptide</th></tr>
            <tr><td>Protein Data Bank code</td><td>9XYZ</td><td>9XYA</td></tr></table>
            <table-wrap-foot><fn>Apo-structure diffraction was collected experimentally.</fn></table-wrap-foot>
            </table-wrap></body></article>
        """)
        self.assertIsNone(accession_state(article, '9XYA'))
        reader.get = lambda *args, **kwargs: (b'{"esearchresult":{"idlist":["123"]}}', {'status': 200})
        reader.get_document = lambda citation: (
            ({'status': 'PRIMARY_CITATION_READ', 'pmid': '123', 'pdb_deposits': ['9XYZ']}, '', 'PMC456')
            if citation else ({'status': 'PRIMARY_CITATION_UNAVAILABLE'}, '', None))
        reader.full_article = lambda _: (article, {'status': 200})
        evidence = reader.assess({'primary_citation': {}, 'source_entry': '9XYZ', 'title': 'Receptor'})
        self.assertEqual(evidence['support']['basis'], 'ARTICLE_APO_CRYSTAL_COMPARISON')

    def test_document_response_is_reused_across_different_runs(self):
        from hashlib import sha256

        body = b'<PubmedArticleSet />'
        receipt = {'status': 200, 'response_sha256': sha256(body).hexdigest(),
                   'response_bytes': len(body), 'elapsed_s': 0.1}
        with TemporaryDirectory() as root, patch('apo_c_literature.cache_get', return_value=(body, receipt)) as fetch:
            first = PrimaryLiterature(Path(root) / 'first', Path(root) / 'cache')
            second = PrimaryLiterature(Path(root) / 'second', Path(root) / 'cache')
            args = ('pubmed/123', 'https://example.invalid/source', {'id': '123'})
            first.get(*args)
            data, meta = second.get(*args)
            self.assertEqual(data, body)
            self.assertEqual(meta['literature_cache'], 'shared')
            self.assertEqual(second.stats['shared_hits'], 1)
            self.assertEqual(fetch.call_count, 1)

    def test_reference_anchors_resolve_equal_terminal_repeat_alignments(self):
        query, candidate = list('ACDEPH'), list('ACDEPHPH')
        old = {1: 1, 2: 2, 3: 3, 4: 4, 5: 7, 6: 8}
        q_ref = {key: 100 + key for key in range(1, 7)}
        c_ref = {key: 100 + key for key in range(1, 9)}
        resolved, alignment = anchored_end_alignment(query, candidate, old, {'mismatches': []},
                                                       q_ref, c_ref, {3, 4})
        self.assertEqual(resolved[6], 6)
        self.assertEqual(alignment['unaligned_candidate_positions'], [7, 8])

    def test_tagged_ends_allow_only_an_identical_deposited_core(self):
        left, right = list('HHACDE'), list('MSGACDE')
        left_ref = dict(zip(range(3, 7), range(10, 14)))
        right_ref = dict(zip(range(4, 8), range(10, 14)))
        mapping, identity = tagged_end_alignment(
            left, right, {'mismatches': [1]}, left_ref, right_ref, {4, 5}, {1, 2}, {1, 2, 3})
        self.assertEqual(mapping, {3: 4, 4: 5, 5: 6, 6: 7})
        self.assertEqual(identity['unaligned_holo_positions'], [1, 2])
        self.assertEqual(identity['mismatches'], [])
        self.assertIsNone(tagged_end_alignment(
            left, list('MSGACNE'), {}, left_ref, right_ref, {4, 5}, {1, 2}, {1, 2, 3}))
        self.assertIsNone(tagged_end_alignment(
            left, right, {}, left_ref, right_ref, {4, 5}, {1, 2}, {1, 3}))
        self.assertIsNone(tagged_end_alignment(
            left, right, {}, left_ref, right_ref, {1, 5}, {1, 2}, {1, 2, 3}))

    def test_author_state_by_accession_and_matching_crystallization_methods(self):
        doc = ET.fromstring("""
            <article><body><sec><p>Human enzyme-A was crystallized in its apo form and in complexes with phosphate.</p></sec>
            <p>PDB codes: Human enzyme-A (#9XYZ), Human enzyme-A bound to phosphate (#9XYA).</p>
            <p>Data deposition: Inhibitor-bound enzyme (accession number 9XYB) and apo-form enzyme-C
            (accession number 9XYC) were deposited in the Protein Data Bank.</p></body></article>
        """)
        self.assertEqual(accession_state(doc, '9XYZ')[0], 'ARTICLE_METHODS_WITH_DEPOSITION_GROUP')
        self.assertEqual(accession_state(doc, '9XYC')[0], 'ARTICLE_STATE_BEFORE_ACCESSION')
        self.assertIsNone(accession_state(doc, '9XYA'))
        self.assertIsNone(accession_state(doc, '9XYB'))
        self.assertIsNone(accession_state(ET.fromstring(
            '<article><body><p>PDB codes: Human enzyme-A (#9XYZ), Human enzyme-A bound to phosphate (#9XYA).'
            '</p></body></article>'), '9XYZ'))

    def test_related_deposition_requires_same_primary_study_and_incoming_apo_label(self):
        import json

        reader = PrimaryLiterature(None)
        citation = {'pdbx_database_id_PubMed': '123', 'pdbx_database_id_DOI': '10.0000/control'}
        reader.get_document = lambda _: ({'status': 'PRIMARY_CITATION_READ', 'pmid': '123'}, '', None)
        deposited = {'rcsb_id': '9XYA', 'rcsb_primary_citation': dict(citation),
                     'pdbx_database_related': [{'db_name': 'PDB', 'db_id': '9XYZ', 'details': 'unliganded protein'}]}
        reader.get = lambda *args, **kwargs: (json.dumps(deposited).encode(), {'status': 200})
        info = {'primary_citation': citation, 'source_entry': '9XYZ', 'title': 'Crystal structure',
                'related_deposits': [{'entry': '9XYA', 'details': 'related crystal form'}]}
        self.assertEqual(reader.assess(info)['support']['basis'], 'RECIPROCAL_PDB_DEPOSITION_SAME_PRIMARY_STUDY')
        deposited['rcsb_primary_citation']['pdbx_database_id_DOI'] = '10.0000/unrelated'
        self.assertNotIn('support', reader.assess(info))
        deposited['rcsb_primary_citation'] = citation
        deposited['pdbx_database_related'][0]['details'] = 'not unliganded'
        self.assertNotIn('support', reader.assess(info))

    def test_rank_order_removes_wrong_domain_and_promotes_unbound_hint(self):
        def entry(pdb, title, begin):
            return {'rcsb_id': pdb, 'struct': {'title': title}, 'exptl': [{'method': 'X-RAY DIFFRACTION'}],
                    'polymer_entities': [{'rcsb_polymer_entity_align': [{'reference_database_name': 'UniProt',
                         'reference_database_accession': 'P00001',
                         'aligned_regions': [{'ref_beg_seq_id': begin, 'length': 20}]}]}]}

        ranked = {'one': [{'pdb': pdb} for pdb in ('9XYA', '9XYB', '9XYZ')]}
        rows = [entry('9XYA', 'Extracellular domain', 1), entry('9XYB', 'Bound receptor', 100),
                entry('9XYZ', 'Uninhibited receptor', 100)]
        ordered, evidence = order_rankings(ranked, rows, 'P00001', {105, 106}, 3)
        self.assertEqual([r['pdb'] for r in ordered['one']], ['9XYZ', '9XYB'])
        self.assertTrue(evidence['9XYA']['target_domain_absent'])


if __name__ == '__main__':
    unittest.main()

"""Check when a deposition supports a target-specific unbound state."""
from pathlib import Path
import sys
import unittest
import xml.etree.ElementTree as ET

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from apo_c_literature import PrimaryLiterature, accession_state
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


if __name__ == '__main__':
    unittest.main()

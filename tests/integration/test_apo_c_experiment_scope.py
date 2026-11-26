"""Check when a deposition supports a target-specific unbound state."""
from pathlib import Path
import sys
import unittest

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from apo_c_science import experiment_record


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


if __name__ == '__main__':
    unittest.main()

"""Checks for mapping and source parsing used by the exploratory Apo pilot."""
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from apo_c_pair import align_coordinates, sequence_map
from apo_c_pilot import parse_ahoj_table, relevant_uniprot


class ApoPilotCoreTests(unittest.TestCase):
    def test_alignment_retains_insertion_and_deletion_positions(self):
        mapping, metrics = sequence_map(list('ACD'), list('AXCD'))
        self.assertEqual(mapping, {1: 1, 2: 3, 3: 4})
        self.assertEqual(metrics['paired_positions'], 3)
        inverse, _ = sequence_map(list('AXCD'), list('ACD'))
        self.assertEqual(inverse, {1: 1, 3: 2, 4: 3})

    def test_rigid_alignment_recovers_rotation_and_translation(self):
        source = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
        expected_rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        target = source @ expected_rotation + np.array([7., -4., 2.])
        rotation, translation, residual = align_coordinates(source, target)
        np.testing.assert_allclose(source @ rotation + translation, target, atol=1e-10)
        self.assertLess(residual, 1e-10)

    def test_ahoj_rows_preserve_entry_chain_and_order(self):
        page = '''<table id="apo-results"><tbody>
        <tr data-row-id="apo_1"><td><a href="https://www.rcsb.org/structure/2xtt">2xtt</a>
        <span class="chain-id">B</span></td><td>1.5</td></tr>
        <tr data-row-id="apo_2"><td><a href="https://www.rcsb.org/structure/1s0q">1s0q</a>
        <span class="chain-id">A</span></td><td>1.1</td></tr>
        </tbody></table>'''
        rows = parse_ahoj_table(page, 'apo')
        self.assertEqual([(r['pdb'], r['auth_chain'], r['source_order']) for r in rows],
                         [('2XTT', 'B', 1), ('1S0Q', 'A', 2)])

    def test_sifts_filter_uses_labels_of_all_required_chains(self):
        case = {'identity': {'binding_instances': [['1.A', 'A', '1'], ['2.A', 'A', '2'], ['1.C', 'C', '1']]}}
        mappings = {'U1': {'mappings': [{'struct_asym_id': 'A'}, {'struct_asym_id': 'C'}]},
                    'U2': {'mappings': [{'struct_asym_id': 'B'}]}}
        result = relevant_uniprot(case, mappings)
        self.assertEqual(set(result), {'U1'})
        self.assertEqual(len(result['U1']['mappings']), 2)


if __name__ == '__main__':
    unittest.main()

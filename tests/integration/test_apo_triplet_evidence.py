"""Regression checks for identity and uncertainty in fixed Apo/Holo evidence."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from apo_c_triplet_regression import (atom_alt_compatibility, ccd_heavy_atoms,
                                      connections_for_target, unobserved_records)
from apo_c_chain_equivalence import validate_assignment
from apo_c_rank_audit import unique_entry_ranking
from apo_c_holdout_prepare import source_metadata


class TripletEvidenceTests(unittest.TestCase):
    def test_connection_uses_label_author_number_and_component(self):
        target = {'target_residue_instances': [['1.B', 'B', '1', '2', 'B', '107', 'PHE']]}
        right = {'kind': 'covale', 'ends': [
            {'label_asym': 'E', 'auth_asym': 'B', 'auth_seq': '101', 'ccd': '112', 'insertion': None},
            {'label_asym': 'B', 'auth_asym': 'B', 'auth_seq': '107', 'ccd': 'PHE', 'insertion': None}]}
        wrong = {'kind': 'covale', 'ends': [
            {'label_asym': 'Q', 'auth_asym': 'B', 'auth_seq': '107', 'ccd': 'PHE', 'insertion': None},
            {'label_asym': 'E', 'auth_asym': 'B', 'auth_seq': '101', 'ccd': '112', 'insertion': None}]}
        with patch('apo_c_triplet_regression.connection_edges', return_value=[right, wrong]):
            found = connections_for_target(None, target)
        self.assertEqual(len(found), 1)
        self.assertEqual(found[0]['target_endpoint_membership'], [False, True])

    def test_ccd_heavy_atom_count_excludes_hydrogen_and_deuterium(self):
        block = gemmi.cif.read_string('''data_x
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
RCB C1 C
RCB O1 O
RCB H1 H
RCB D1 D
OTHER C2 C
''').sole_block()
        self.assertEqual(ccd_heavy_atoms(block, 'RCB'), ['C1', 'O1'])

    def test_unobserved_records_match_selected_model_and_label(self):
        block = gemmi.cif.read_string('''data_x
loop_
_pdbx_unobs_or_zero_occ_residues.PDB_model_num
_pdbx_unobs_or_zero_occ_residues.label_asym_id
_pdbx_unobs_or_zero_occ_residues.label_seq_id
_pdbx_unobs_or_zero_occ_residues.occupancy_flag
1 A 190 1
2 A 191 1
1 B 192 0
1 A 193 0
''').sole_block()
        result = unobserved_records(block, '_pdbx_unobs_or_zero_occ_residues.', 'A', 1)
        self.assertEqual(set(result), {190, 193})
        self.assertEqual(result[190][0]['deposited_observation_status'], 'UNOBSERVED')
        self.assertEqual(result[193][0]['deposited_observation_status'], 'ZERO_OCCUPANCY')

    def test_equal_nonblank_altlocs_are_not_assumed_correlated(self):
        self.assertEqual(atom_alt_compatibility({'altloc': 'A'}, {'altloc': 'A'}),
                         'cross_residue_altloc_compatibility_unverified')
        self.assertEqual(atom_alt_compatibility({'altloc': ''}, {'altloc': 'A'}), 'one_or_both_blank_altloc')

    def test_multichain_assignment_rejects_missing_or_reused_copy(self):
        required = ['1.A', '1.C']
        validate_assignment([['1.A', '1.B'], ['1.C', '2.B']], required)
        with self.assertRaises(ValueError):
            validate_assignment([['1.A', '1.B']], required)
        with self.assertRaises(ValueError):
            validate_assignment([['1.A', '1.B'], ['1.C', '1.B']], required)

    def test_entry_ranking_retains_distinct_chain_rows(self):
        rows = [{'pdb': '1AKE', 'author_chain': 'A'}, {'pdb': '4AKE', 'author_chain': 'A'},
                {'pdb': '4AKE', 'author_chain': 'B'}, {'pdb': '4X8H', 'author_chain': 'A'}]
        ranked = unique_entry_ranking(rows, '1AKE')
        self.assertEqual([(x['pdb'], x['entry_rank']) for x in ranked], [('4AKE', 1), ('4X8H', 2)])
        self.assertEqual([x['author_chain'] for x in ranked[0]['source_rows']], ['A', 'B'])

    def test_holdout_metadata_distinguishes_uniprot_from_pdb_references(self):
        content = '''data_x
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_struct_ref.db_name
_struct_ref.pdbx_db_accession
_struct_ref.entity_id
PDB 1E2E 1
UNP P23919 1
loop_
_entity_poly.entity_id
_entity_poly.pdbx_seq_one_letter_code_can
1 ACDE
loop_
_citation.id
_citation.pdbx_database_id_PubMed
_citation.pdbx_database_id_DOI
primary 12345 10.1234/example
'''
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'example.cif'
            path.write_text(content, encoding='ascii')
            result = source_metadata(path, {'A'})
        self.assertEqual(result['accessions'], ['P23919'])
        self.assertEqual(len(result['sequence_sha256']), 1)


if __name__ == '__main__':
    unittest.main()

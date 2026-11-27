"""Checks for bounded candidate selection and conservative evidence states."""
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from apo_c_autopilot import DownloadBudgetExceeded, download_bounded, observation_flags, observation_status, select_budget
from apo_c_pilot import sha
from apo_c_union_triage import candidate_evidence


class ApoAutopilotTests(unittest.TestCase):
    def test_round_robin_budget_preserves_cross_source_provenance(self):
        ranked = {
            'PLINDER': [{'pdb': 'AAAA', 'source_rows': [{'assembly': '1'}]},
                        {'pdb': 'BBBB', 'source_rows': []}],
            'AHoJ': [{'pdb': 'BBBB', 'source_rows': [{'author_chain': 'A'}]}],
            'RCSB': [{'pdb': 'CCCC', 'source_rows': [{'entity': 'CCCC_1'}]}],
        }
        selected = select_budget(ranked, per_source=2, max_unique=2)
        self.assertEqual([x['pdb'] for x in selected], ['AAAA', 'BBBB'])
        self.assertEqual(set(selected[1]['sources']), {'PLINDER', 'AHoJ'})

    def test_empty_observed_site_never_claims_unbound_or_acceptance(self):
        entity = {'common_uniprot': ['P12345'], 'holo_site_residues': 2}
        observation = {'site_observed_residues': 2, 'site_complete_relative_holo_residues': 2,
                       'actual_component_contacts': {'status': 'MEASURED', 'nonwater': []}}
        self.assertEqual(observation_status(entity, observation),
                         'OBSERVED_SCOPE_NO_NONWATER_CONTACT_UNBOUND_UNPROVEN')
        observation['actual_component_contacts']['nonwater'] = [{
            'component': 'MG', 'min_distance_to_mapped_site_atoms_A': 9.65, 'atoms_within_4A_of_mapped_site': 0}]
        self.assertEqual(observation_status(entity, observation),
                         'OBSERVED_SCOPE_NO_NONWATER_CONTACT_UNBOUND_UNPROVEN')
        observation['site_observed_residues'] = 1
        self.assertEqual(observation_status(entity, observation), 'SITE_NOT_FULLY_OBSERVED')
        observation['actual_component_contacts']['nonwater'].append({
            'component': 'ACP', 'min_distance_to_mapped_site_atoms_A': 2.66, 'atoms_within_4A_of_mapped_site': 20})
        self.assertEqual([x['component'] for x in observation_flags(entity, observation)['nonwater_within_4A']], ['ACP'])
        observation['site_observed_residues'] = 2
        self.assertEqual(observation_status(entity, observation), 'NONWATER_CONTACT_REQUIRES_REVIEW')
        entity['common_uniprot'] = []
        self.assertEqual(observation_status(entity, observation), 'IDENTITY_REQUIRES_REVIEW')

    def test_oversized_source_is_rejected_before_get(self):
        head = SimpleNamespace(status_code=200, headers={'Content-Length': '69723946'}, close=lambda: None)
        with patch('apo_c_autopilot.requests.head', return_value=head), patch('apo_c_autopilot.requests.get') as get:
            with self.assertRaises(DownloadBudgetExceeded) as caught:
                download_bounded('https://example.invalid/example.cif.gz', 20000000, 30)
        self.assertEqual(caught.exception.evidence['head_declared_bytes'], 69723946)
        self.assertEqual(caught.exception.evidence['bytes_received'], 0)
        get.assert_not_called()

    def test_expand_only_assemblies_with_mapped_receptor_and_keep_companions(self):
        def chain(name, entity, position, occupancy):
            receptor = gemmi.Chain(name)
            residue = gemmi.Residue()
            residue.name = 'ALA'
            residue.entity_type = gemmi.EntityType.Polymer
            residue.entity_id, residue.subchain = entity, name
            residue.label_seq, residue.seqid = 1, gemmi.SeqId(1, ' ')
            atom = gemmi.Atom()
            atom.name, atom.element = 'CA', gemmi.Element('C')
            atom.pos, atom.occ = gemmi.Position(position, 0, 0), occupancy
            residue.add_atom(atom)
            receptor.add_residue(residue)
            return receptor

        class Model(list):
            num = 1

        class Structure(list):
            pass

        structure = Structure([Model([chain('A', '1', 0, 1), chain('B', '2', 2, 1), chain('C', '1', 5, 0)])])
        structure.entities = [SimpleNamespace(name='1', full_sequence=['ALA'])]

        def assembly(name, labels, count=1):
            return SimpleNamespace(name=name, generators=[SimpleNamespace(subchains=labels, operators=[None] * count)])

        structure.assemblies = [assembly(str(n), ['B'], 200) for n in range(1, 8)]
        structure.assemblies += [assembly('8', ['A', 'B']), assembly('9', ['C']), assembly('10', ['A'], 200)]
        operator = SimpleNamespace(transform=gemmi.Transform())
        instances = {f'1.{c.name}': (operator, c.name, [(c.name, c[0])]) for c in structure[0][:2]}
        with TemporaryDirectory() as tmp:
            original = Path(tmp) / 'source.cif'
            original.write_text("""data_9XYZ
_entry.id 9XYZ
_struct.title 'Synthetic receptor and neighbor'
loop_
_pdbx_audit_revision_history.data_content_type
_pdbx_audit_revision_history.major_revision
_pdbx_audit_revision_history.minor_revision
_pdbx_audit_revision_history.revision_date
'Structure model' 1 1 2020-01-01
""", encoding='utf8')
            site = {1: {'author_seq': '1', 'author_chain': 'A', 'component': 'ALA', 'atom_names': {'CA'}}}
            native = str(original), sha(original), original.stat().st_size, '1.1'
            with patch('apo_c_union_triage.gemmi.make_structure_from_block', return_value=structure), \
                    patch('apo_c_union_triage.assembly_instances', return_value=instances) as expand:
                result = candidate_evidence(Path(tmp), '9XYA', '9XYZ', [{'entity': '9XYZ_1'}],
                                            {'sifts_common_uniprot': []}, None, None, ['A'], site, None,
                                            fixed_source=native, max_assembly_instances=128)
        self.assertEqual([call.args[1] for call in expand.call_args_list], ['8'])
        entity, = result['entities']
        self.assertEqual(entity['assembly_scope']['without_observed_site'], ['9'])
        self.assertEqual(entity['assembly_scope']['size_limited'], [{'assembly': '10', 'instances': 200}])
        measured = next(o for o in entity['assembly_observations'] if o['status'] == 'ACTUAL_CONTACTS_MEASURED')
        self.assertEqual(measured['actual_component_contacts']['nonwater'][0]['instance'], '1.B')


if __name__ == '__main__':
    unittest.main()

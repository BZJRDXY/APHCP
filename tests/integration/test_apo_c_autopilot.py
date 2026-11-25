"""Checks for bounded candidate selection and conservative evidence states."""
from pathlib import Path
import sys
import unittest
from types import SimpleNamespace
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from apo_c_autopilot import DownloadBudgetExceeded, download_bounded, observation_flags, observation_status, select_budget


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


if __name__ == '__main__':
    unittest.main()

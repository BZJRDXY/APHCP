"""Exercise scientific decisions separately from coordinate and role assessment."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from holo_m3_decisions import resolve


def result(run, outcome, parents=()):
    return (run, f'{run}.json', {'system_id':'example-system', 'batch':'pilot', 'pdb_id':'TEST',
                               'candidate_key':'example', 'disposition':outcome, 'supersedes':list(parents),
                               'files':{'system_cif':{'path':'example.cif','sha256':'fixture'}}})


class DecisionTests(unittest.TestCase):
    def test_supersession_preserves_order_independent_current_and_history(self):
        hold=result('z-prior','REVIEW'); accept=result('a-next','REGISTERED',['z-prior/example-system'])
        chosen=resolve([hold,accept])
        self.assertEqual(chosen,resolve([accept,hold]))
        self.assertEqual(chosen['current_decision_id'],'a-next/example-system')
        self.assertEqual(chosen['disposition'],'REGISTERED')
        self.assertEqual({x['disposition'] for x in chosen['provenance']},{'REVIEW','REGISTERED'})
        self.assertEqual(chosen['files'],hold[2]['files'])

    def test_independent_disagreement_and_agreement_both_retain_provenance(self):
        a=result('a','REGISTERED'); z=result('z','REVIEW')
        for seq in ([a,z],[z,a]):
            value=resolve(seq)
            self.assertEqual(value['disposition'],'CONFLICT')
            self.assertIsNone(value['current_decision_id'])
            self.assertEqual(len(value['provenance']),2)
            self.assertEqual(value['files'],a[2]['files'])
        resolved=resolve([a,result('b','REGISTERED')])
        self.assertEqual(resolved['decision_state'],'CONCORDANT')
        self.assertEqual(resolved['disposition'],'REGISTERED')
        self.assertIsNone(resolved['current_decision_id'])

    def test_orphan_and_circular_replacements_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'missing_for_same_identity'):
            resolve([result('a','REVIEW',['unknown/another-system'])])
        with self.assertRaisesRegex(ValueError,'supersession_cycle'):
            resolve([result('a','REVIEW',['b/example-system']),result('b','REGISTERED',['a/example-system'])])


if __name__=='__main__':
    unittest.main()

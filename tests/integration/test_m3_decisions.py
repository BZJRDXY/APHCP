"""Exercise scientific decisions separately from coordinate and role assessment."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from holo_m3_decisions import decision_view, resolve


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

    def test_observation_query_uses_current_evidence_and_clears_conflicts(self):
        old=result('old','REVIEW')
        new=result('new','REGISTERED',['old/example-system'])
        new[2].update(role={'status':'PASS','source_evidence':{'reason':'new evidence'}},
                      qc={'topology':'PASS'})
        frozen={'status':'REVIEW','role':'REVIEW','role_evidence':{'reason':'old evidence'},
                'm3_readback_qc':{'topology':'REVIEW'},'qc':{'coordinate_identity':'PASS'}}
        shown=decision_view(frozen,[old,new])
        self.assertEqual((shown['status'],shown['role']),('REGISTERED','PASS'))
        self.assertEqual(shown['role_evidence'],new[2]['role']['source_evidence'])
        self.assertEqual(shown['m3_readback_qc'],new[2]['qc'])
        self.assertEqual(shown['qc'],frozen['qc'])
        unresolved=decision_view(frozen,[old,new,result('other','REVIEW')])
        self.assertEqual((unresolved['status'],unresolved['role']),('CONFLICT','REVIEW'))
        self.assertIsNone(unresolved['role_evidence'])
        self.assertIsNone(unresolved['m3_readback_qc'])

    def test_orphan_and_circular_replacements_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'missing_for_same_identity'):
            resolve([result('a','REVIEW',['unknown/another-system'])])
        with self.assertRaisesRegex(ValueError,'supersession_cycle'):
            resolve([result('a','REVIEW',['b/example-system']),result('b','REGISTERED',['a/example-system'])])


if __name__=='__main__':
    unittest.main()

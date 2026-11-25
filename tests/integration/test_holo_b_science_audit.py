"""Focused checks for instance linkage and immutable decision resolution."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from holo_b_science_audit import array_stream, assess_match, resolve


class ResearchAuditTests(unittest.TestCase):
    def setUp(self):
        self.identity = {'pdb_id':'1ABC','revision':'1.2','assembly':'1','model_scope':[1],
            'binding_instances':[['1.A','A','1']],
            'target_residue_instances':[['1.B','B','1','2','A','203','LIG']]}
        self.row = {'system_biounit_id':'1','system_id':'s','ligand_id':'l',
            'ligand_ccd_code':'LIG','ligand_asym_id':'B','ligand_auth_id':'A',
            'ligand_instance_chain':'1.B','ligand_instance_chains':['1.B'],
            'ligand_member_asym_ids':['B'],'ligand_residue_numbers':[1],
            'ligand__is_multi_residue':False,'ligand_protein_chains_asym_id':['1.A'],
            'ligand_is_proper':True,'ligand_is_artifact':False,'ligand_is_ion':False}
        self.chains = {('1abc','1','1.A'),('1abc','1','1.B')}

    def match(self, rows, rev='1.2'):
        return assess_match(self.identity,rows,{'1abc':rev},self.chains)[:2]

    def test_copies_and_versions_do_not_inherit_support(self):
        self.assertEqual(self.match([self.row]),('SAME_REVISION_INSTANCE_CHAIN_MATCH','PLINDER_PROPER_RULE_SUPPORT'))
        self.assertEqual(self.match([self.row],rev='1.1'),('DIFFERENT_REVISION_UNVERIFIED','UNKNOWN'))
        other_copy = dict(self.row,ligand_instance_chain='2.B')
        self.assertEqual(self.match([other_copy]),('CCD_ONLY','UNKNOWN'))
        other_pocket = dict(self.row,ligand_protein_chains_asym_id=['2.A'])
        self.assertEqual(self.match([other_pocket]),('SAME_REVISION_SCOPE_OR_SITE_AMBIGUOUS','UNKNOWN'))

    def test_inconsistent_labels_are_not_silent_support(self):
        contradiction = dict(self.row,ligand_is_proper=False,ligand_is_artifact=True)
        self.assertEqual(self.match([self.row,contradiction]),('SAME_REVISION_INSTANCE_CHAIN_MATCH','CONFLICT'))

    def test_immutable_decisions_resolve_supersession(self):
        older = {'system_id':'s','disposition':'REVIEW','supersedes':None}
        newer = {'system_id':'s','disposition':'REGISTERED','supersedes':['old/s']}
        row, state = resolve([('old',older),('new',newer)])
        self.assertEqual((row['disposition'],state),('REGISTERED','EFFECTIVE'))
        with self.assertRaises(ValueError):
            resolve([('new',newer)])

    def test_sealed_json_is_streamed_and_verified(self):
        data = json.dumps([{'key':'配体'}, {'key':'value'}],ensure_ascii=False).encode('utf-8')
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'records.json'; path.write_bytes(data)
            self.assertEqual(len(list(array_stream(path,hashlib.sha256(data).hexdigest()))),2)
            with self.assertRaises(ValueError):
                list(array_stream(path,'0'*64))


if __name__ == '__main__':
    unittest.main()

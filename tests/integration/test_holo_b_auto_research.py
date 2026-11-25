"""Safeguards for the research-only automatic evidence pilot."""
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import gemmi

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import holo_b_auto_research as pilot


def block(title='protein',code='NAG',second_copy=False,symmetry='1_555',target_atom='C1'):
    scheme='B 2 '+code+' A 200\n'
    if second_copy:scheme+='C 2 '+code+' A 201\n'
    return gemmi.cif.read_string('''data_test
_struct.title 'TITLE'
loop_
_chem_comp.id
_chem_comp.name
_chem_comp.pdbx_synonyms
CODE 'target compound' ?
loop_
_entity.id
_entity.pdbx_description
1 'protein'
2 'target compound'
loop_
_pdbx_nonpoly_scheme.asym_id
_pdbx_nonpoly_scheme.entity_id
_pdbx_nonpoly_scheme.mon_id
_pdbx_nonpoly_scheme.pdb_strand_id
_pdbx_nonpoly_scheme.pdb_seq_num
SCHEME
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.ptnr2_symmetry
c1 covale A ASN ND2 A 118 1_555 B CODE ATOM A 200 SYM
'''.replace('TITLE',title).replace('CODE',code).replace('SCHEME',scheme)
   .replace('ATOM',target_atom).replace('SYM',symmetry)).sole_block()


def row(code='NAG',operator='1',new_qc=True,conflicts=None):
    identity={'pdb_id':'1ABC','revision':'1.0','assembly':'1','model_scope':[1],
        'binding_instances':[['1.A','A','1']],
        'target_residue_instances':[['1.B','B',operator,'2','A','200',code]]}
    obs={'identity':identity,'files':{'source_mmcif':{'path':'example.cif.gz','sha256':'abc'}}}
    old={'coordinate_identity':'PASS','experimental_fidelity':'PASS','topology':'PASS','partner_scope_qc':'PASS'}
    new={'coordinate_identity':'PASS','experimental_fidelity':'PASS','topology':'PASS','partner_scope':'PASS'}
    return {'system_id':'s','document':json.dumps(obs),'old_qc':json.dumps(old),'new_qc':json.dumps(new) if new_qc else 'null',
        'conflicts':json.dumps(conflicts or []),'source_evidence':json.dumps({'reason':'no explicit statement'}),
        'match_level':'SAME_REVISION_INSTANCE_CHAIN_MATCH','annotation_interpretation':'PLINDER_PROPER_RULE_SUPPORT'}


class AutomaticEvidenceTests(unittest.TestCase):
    def test_exact_glycosylation_site_only(self):
        self.assertEqual(pilot.assess(row(),block(),{})['decision'],'PROTEIN_N_GLYCOSYLATION_LINK_SUPPORTED')
        self.assertEqual(pilot.assess(row(operator='2'),block(),{})['decision'],'REVIEW')
        self.assertEqual(pilot.assess(row(),block(symmetry='2_555'),{})['decision'],'REVIEW')
        self.assertEqual(pilot.assess(row(),block(target_atom='C4'),{})['decision'],'REVIEW')

    def test_covalency_and_upstream_proper_are_not_enough(self):
        self.assertEqual(pilot.assess(row(code='LIG'),block(code='LIG'),{})['decision'],'REVIEW')
        self.assertEqual(pilot.assess(row(code='LIG'),block(code='LIG',title='complex with LIG',second_copy=True),{})['decision'],'REVIEW')

    def test_qc_and_explicit_conflicts_override_title(self):
        cif=block(code='LIG',title='complex with LIG')
        self.assertEqual(pilot.assess(row(code='LIG'),cif,{})['decision'],'DEPOSITOR_NAMED_BINDING_PARTNER_SUPPORTED')
        self.assertEqual(pilot.assess(row(code='LIG',new_qc=False),cif,{})['decision'],'REVIEW')
        self.assertEqual(pilot.assess(row(code='LIG',conflicts=['conflict']),cif,{})['decision'],'REVIEW')

    def test_article_failure_is_not_a_scientific_negative(self):
        with patch.object(pilot.urllib.request,'urlopen',side_effect=TimeoutError('timeout')):
            articles,receipt=pilot.article_batch(['123'],timeout=1)
        self.assertEqual(articles,{})
        self.assertEqual(receipt['state'],'SOURCE_TECHNICAL_FAILURE')
        self.assertEqual(pilot.assess(row(),block(),{'state':'SOURCE_TECHNICAL_FAILURE'})['decision'],
                         'PROTEIN_N_GLYCOSYLATION_LINK_SUPPORTED')


if __name__=='__main__':unittest.main()

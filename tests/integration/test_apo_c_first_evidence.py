"""Specific regressions for accession binding, misleading extracts and fixed-source reuse."""
from pathlib import Path
from tempfile import TemporaryDirectory
import gzip
import sys
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[2] / 'scripts'))
from apo_c_autopilot import shared_fixed_original
from apo_c_literature import reported_apo_structure, fixed_relation_support
from apo_c_pilot import write_json, sha
from verify_apo_c_evidence_replay import check_experimental_citation


class FirstApoEvidence(unittest.TestCase):
    def test_result_and_absence_in_same_abstract(self):
        missing = 'Crystal structures have been described, but currently no structures of the apo form are available.'
        result = 'Here, the apo structure of this receptor was solved at a resolution of 2.1 Å.'
        self.assertEqual(reported_apo_structure(missing + ' ' + result), result)
        self.assertIsNone(reported_apo_structure(missing))
        with self.assertRaises((ValueError, AssertionError)):
            check_experimental_citation({'basis':'MATCHED_APO_DEPOSITION_AND_PRIMARY_ABSTRACT', 'excerpt':missing}, {'candidate':'1ABC'})

    def test_depositor_must_name_this_candidate_without_conflicting_rows(self):
        source={'entry':'1BCD','rows':[(2,'1ABC','Apo form of the receptor')], 'revision':'1.1','sha256':'0'*64}
        ex={'source_entry':'1ABC','title':'X-ray study of receptor'}
        self.assertEqual(fixed_relation_support(ex,source)['row_numbers'], [2])
        self.assertIsNone(fixed_relation_support(ex|{'source_entry':'2ABC'},source))
        for bad in ['not apo', 'apo or ligand-bound', 'apo template']:
            self.assertIsNone(fixed_relation_support(ex,source|{'rows':source['rows']+[(3,'1ABC',bad)]}))

    def test_optional_discovery_bounds_and_reuses_unresolved_accessions(self):
        from apo_c_author_documents import CandidateDocuments

        def unavailable(candidates, out, **options):
            self.assertEqual(candidates, ['1ABC'])
            self.assertLessEqual(options['max_requests'], 8)
            write_json(out / 'documents.json', {'documents': []})
            write_json(out / 'collection.json', {'http': {'requests': 1}, 'unavailable': {'Timeout': 1}})

        def no_documents(candidates, out, **options):
            write_json(out / 'documents.json', {'documents': []})
            write_json(out / 'collection.json', {'http': {'requests': 1}, 'unavailable': {}})

        with TemporaryDirectory() as td:
            root = Path(td)
            resolver = CandidateDocuments(root, root / 'cache', 1)
            with patch('apo_c_author_documents.discover', side_effect=unavailable) as discover:
                self.assertEqual(resolver.lookup('1ABC')[1]['status'], 'SOURCE_UNAVAILABLE')
                self.assertEqual(resolver.lookup('1ABC')[1]['status'], 'SOURCE_UNAVAILABLE')
                self.assertEqual(resolver.lookup('2ABC')[1]['status'], 'CANDIDATE_BUDGET')
                self.assertEqual(discover.call_count, 1)
            with patch('apo_c_author_documents.discover', side_effect=no_documents) as discover:
                resumed = CandidateDocuments(root, root / 'cache', 1)
                self.assertEqual(resumed.lookup('1ABC')[1]['status'], 'COMPLETE')
                second = CandidateDocuments(root, root / 'cache', 1)
                self.assertEqual(second.lookup('1ABC')[1]['status'], 'COMPLETE')
                self.assertEqual(discover.call_count, 1)

    def test_fixed_revision_and_provenance_survive_shard_reuse(self):
        from test_apo_c_recovery import fixed_cif
        with TemporaryDirectory() as td:
            root = Path(td); shared = root/'shared'; source = root/'original.cif.gz'; source.write_bytes(fixed_cif(1))
            receipt = root/'origin.json'; row={'status':'FETCHED','path':str(source),'revision':'1.1','sha256':sha(source),'bytes':source.stat().st_size};write_json(receipt,row)
            first = shared_fixed_original(shared,'1ABC','1.1',row,receipt)
            self.assertEqual(shared_fixed_original(shared,'1ABC','1.1'),first)
            self.assertEqual(shared_fixed_original(shared,'1ABC','1.0'),(None,None))
            source.write_bytes(gzip.compress(b'not the checked file'))
            with self.assertRaises(ValueError): shared_fixed_original(shared,'1ABC','1.1')

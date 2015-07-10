import os
import unittest
import shutil

from mod_sbml.constraint_based_analysis.efm.reaction_classification_by_efm import get_important_reactions
from mod_sbml.cobra_tests.SBMLTestCase import DATA_DIR, TEST_SBML, create_test_sbml
from mod_sbml.constraint_based_analysis.efm.efm_analyser import get_efms
from mod_sbml.utils.path_manager import create_dirs

__author__ = 'anna'

EFM_DIR = os.path.join(DATA_DIR, 'efms')

class EFMTestCase(unittest.TestCase):

    def setUp(self):
        """
                     r2-> m2 <-r3->
                  <-      ^
        <-r1-> m1         r6
                  <-      |
                     r4-> m3 --r5->


        Expect to find 2 EFMs that include r3:
            10 r1	10 r2	10 r3
            10 r1	10 r3	10 r4	10 r6
        """
        create_test_sbml()
        create_dirs(EFM_DIR)

    def tearDown(self):
        if os.path.exists(TEST_SBML):
            os.remove(TEST_SBML)
        if os.path.exists(EFM_DIR):
            shutil.rmtree(EFM_DIR, True)

    def test_efm_num(self):
        id2efm = get_efms(target_r_id='r3', target_r_reversed=False, sbml=TEST_SBML, directory=EFM_DIR, r_id2rev={})
        self.assertEqual(2, len(id2efm), 'EFM number was supposed to be 2, got %g instead.' % len(id2efm))

    def test_efm_content(self):
        id2efm = get_efms(target_r_id='r3', target_r_reversed=False, sbml=TEST_SBML, directory=EFM_DIR, r_id2rev={})
        self.assertIn({'r1': 10, 'r2': 10, 'r3': 10}, [efm.to_r_id2coeff() for efm in id2efm.itervalues()],
                      'Failed to detect EFM: 10 r1, 10 r2, 10 r3')

    def test_important_reactions(self):
        id2efm = get_efms(target_r_id='r3', target_r_reversed=False, sbml=TEST_SBML, directory=EFM_DIR, r_id2rev={})
        _, important_r_ids = get_important_reactions(id2efm, imp_rn_threshold=0)
        imp_rs = {'r1', 'r2', 'r3', 'r4', 'r6'}
        self.assertEqual(imp_rs, important_r_ids, "Important reactions were supposed to be %s, got %s"
                         % (imp_rs, important_r_ids))

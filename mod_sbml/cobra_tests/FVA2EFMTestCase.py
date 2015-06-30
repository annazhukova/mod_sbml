import os
import unittest
import shutil

from cobra.io.sbml import create_cobra_model_from_sbml_file

from mod_sbml.cobra_tests.SBMLTestCase import DATA_DIR, TEST_SBML, create_test_sbml
from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.fva_analyser import analyse_by_fva
from mod_sbml.constraint_based_analysis.efm.efm_analyser import perform_efma
from mod_sbml.utils.path_manager import create_dirs

__author__ = 'anna'

EFM_DIR = os.path.join(DATA_DIR, 'efms')
FVA_DIR = os.path.join(DATA_DIR, 'fva')

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
        create_dirs(FVA_DIR)

    def tearDown(self):
        if os.path.exists(TEST_SBML):
            os.remove(TEST_SBML)
        if os.path.exists(EFM_DIR):
            shutil.rmtree(EFM_DIR, True)
        if os.path.exists(FVA_DIR):
            shutil.rmtree(FVA_DIR, True)

    def test_important_reaction_num_to_fva_reaction_num_correspondence(self):
        """
        FVA finds all the reactions that can participate in the optimal solution,
        i.e. all the reactions for which v_fva_min != 0 or V_fva_max != 0.

        Important reactions that participate in more than 0 EFMs are those
        that can participate in any solution that gives a non-zero flux through the objective reaction.

        Thus fva_reactions \subseteq important_reactions.
        """
        cobra_model = create_cobra_model_from_sbml_file(TEST_SBML)
        r_id2bounds, _, _ = analyse_by_fva(cobra_model, bm_r_id='r3', directory=FVA_DIR, objective_sense='maximize')

        _, important_r_ids = \
            perform_efma(target_r_id='r3', target_r_reversed=False, sbml=TEST_SBML, directory=EFM_DIR, r_id2rev={},
                         acom_path=None, calculate_patterns=False, calculate_important_reactions=True,
                         imp_rn_threshold=0, rewrite=True)
        for r_id in r_id2bounds.iterkeys():
            self.assertIn(r_id, important_r_ids, "%s was supposed to be important." % r_id)

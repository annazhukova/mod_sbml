import os
import shutil
import unittest

from cobra.io.sbml import create_cobra_model_from_sbml_file

from constraint_based_analysis.cobra_constraint_based_analysis.fva_analyser import analyse_by_fva
from cobra_tests.SBMLTestCase import DATA_DIR, TEST_SBML, create_test_sbml
from utils.path_manager import create_dirs

__author__ = 'anna'

FVA_DIR = os.path.join(DATA_DIR, 'fva')

class FVATestCase(unittest.TestCase):

    def setUp(self):
        """
                     r2-> m2 <-r3->
                  <-      ^
        <-r1-> m1         r6
                  <-      |
                     r4-> m3 --r5->
        """
        create_test_sbml()
        create_dirs(FVA_DIR)

    def tearDown(self):
        if os.path.exists(TEST_SBML):
            os.remove(TEST_SBML)
        if os.path.exists(FVA_DIR):
            shutil.rmtree(FVA_DIR, True)

    def test_fva_objective(self):
        cobra_model = create_cobra_model_from_sbml_file(TEST_SBML)
        r_id2bounds, _ = analyse_by_fva(cobra_model, bm_r_id='r3', directory=FVA_DIR, objective_sense='maximize')
        r3_flux = r_id2bounds['r3'][0]
        self.assertEqual(10, r3_flux, 'Optimal flux through r3 was supposed to be 10, got %g instead.' % r3_flux)

    def test_fva_ess_rxns_num(self):
        cobra_model = create_cobra_model_from_sbml_file(TEST_SBML)
        r_id2bounds, _ = analyse_by_fva(cobra_model, bm_r_id='r3', directory=FVA_DIR, objective_sense='maximize')
        essential_rs_num = len([1 for (l_b, u_b) in r_id2bounds.itervalues() if l_b * u_b > 0])
        self.assertEqual(2, essential_rs_num, 'Number of essential reactions was supposed to be 2, got %g instead.'
                         % essential_rs_num)

    def test_fva_ess_rxn_in(self):
        cobra_model = create_cobra_model_from_sbml_file(TEST_SBML)
        r_id2bounds, _ = analyse_by_fva(cobra_model, bm_r_id='r3', directory=FVA_DIR, objective_sense='maximize')
        essential_rs = {r_id for (r_id, (l_b, u_b)) in r_id2bounds.iteritems() if l_b * u_b > 0}
        self.assertIn('r1', essential_rs, 'Reaction r1 should have been essential')

    def test_fva_ess_rxn_out(self):
        cobra_model = create_cobra_model_from_sbml_file(TEST_SBML)
        r_id2bounds, _ = analyse_by_fva(cobra_model, bm_r_id='r3', directory=FVA_DIR, objective_sense='maximize')
        essential_rs = {r_id for (r_id, (l_b, u_b)) in r_id2bounds.iteritems() if l_b * u_b > 0}
        self.assertNotIn('r5', essential_rs, 'Reaction r5 should not have been essential')

    def test_fva_var_rxn_num(self):
        cobra_model = create_cobra_model_from_sbml_file(TEST_SBML)
        r_id2bounds, _ = analyse_by_fva(cobra_model, bm_r_id='r3', directory=FVA_DIR, objective_sense='maximize')
        variable_rs = {r_id for (r_id, (l_b, u_b)) in r_id2bounds.iteritems() if l_b * u_b <= 0}
        self.assertEqual(3, len(variable_rs), 'Number of variable reactions was supposed to be 3, got %g instead.'
                         % len(variable_rs))

    def test_fva_var_rxn_in(self):
        cobra_model = create_cobra_model_from_sbml_file(TEST_SBML)
        r_id2bounds, _ = analyse_by_fva(cobra_model, bm_r_id='r3', directory=FVA_DIR, objective_sense='maximize')
        variable_rs = {r_id for (r_id, (l_b, u_b)) in r_id2bounds.iteritems() if l_b * u_b <= 0}
        self.assertIn('r4', variable_rs, 'Reaction r4 should have been variable')

    def test_fva_var_rxn_out(self):
        cobra_model = create_cobra_model_from_sbml_file(TEST_SBML)
        r_id2bounds, _ = analyse_by_fva(cobra_model, bm_r_id='r3', directory=FVA_DIR, objective_sense='maximize')
        variable_rs = {r_id for (r_id, (l_b, u_b)) in r_id2bounds.iteritems() if l_b * u_b <= 0}
        self.assertNotIn('r5', variable_rs, 'Reaction r5 should not have been variable')
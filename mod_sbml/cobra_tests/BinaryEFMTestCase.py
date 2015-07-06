import unittest
from mod_sbml.constraint_based_analysis.efm.EFM import EFM

__author__ = 'anna'


class BinaryEFMTestCase(unittest.TestCase):

    def setUp(self):
        """
        Model:

        <-r0                 r2-> m2 <-r3->
            ->             <-      ^
              m0 <-r1-> m1       r6
            ->             <-      |
        -r00                 r4-> m3 --r5->
        """
        self.r_id2coeff = {'r0': 1, 'r1': 1, 'r2': 1, 'r3': 1}
        self.r_ids = ['r00', 'r0', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6']
        self.rev_r_ids = {'r0', 'r1', 'r2', 'r3', 'r4'}
        self.efm = EFM(r_ids=self.r_ids,
                       rev_r_ids=self.rev_r_ids,
                       r_id2coeff=self.r_id2coeff)

    def tearDown(self):
        self.efm = None
        self.r_id2coeff = None
        self.r_ids = None
        self.rev_r_ids = None

    def test_efm_len(self):
        efm_len = len(self.efm)
        self.assertEqual(4, efm_len, 'EFM length was supposed to be 4, got %d instead' % efm_len)

    def test_efm_binary(self):
        binary_efm = self.efm.binary_efm
        shift = 0
        expected_result = 0
        for r_i in self.r_ids:
            expected_result += ((1 if r_i in self.r_id2coeff and self.r_id2coeff[r_i] > 0 else 0) << shift)
            shift += 1
            if r_i in self.rev_r_ids:
                expected_result += ((1 if r_i in self.r_id2coeff and self.r_id2coeff[r_i] < 0 else 0) << shift)
                shift += 1
        self.assertEqual((expected_result,), binary_efm,
                         "EFM binary representation was supposed to be %s, got %s instead"
                         % (bin(expected_result), bin(binary_efm[0])))

    def test_efm_r_id2coeff(self):
        r_id2coeff = self.efm.to_r_id2coeff()
        self.assertEqual(self.r_id2coeff, r_id2coeff,
                         "EFM was supposed to be %s, got %s instead" % (self.r_id2coeff, r_id2coeff))

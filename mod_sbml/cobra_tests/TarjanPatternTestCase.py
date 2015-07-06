import unittest
from mod_sbml.constraint_based_analysis.efm.pattern_detection import detect_patterns
from mod_sbml.constraint_based_analysis.efm.EFM import EFM

__author__ = 'anna'


class TarjanPatternEFMTestCase(unittest.TestCase):

    def setUp(self):
        """
        Model:

        <-r0                 r2-> m2 <-r3->
            ->             <-      ^
              m0 <-r1-> m1       r6
            ->             <-      |
        -r00                 r4-> m3 --r5->
        """
        self.r_ids = ['r00', 'r0', 'r1', 'r2', 'r3', 'r4', 'r5', 'r6']
        self.rev_r_ids = {'r0', 'r1', 'r2', 'r3', 'r4'}

    def tearDown(self):
        self.r_ids = None
        self.rev_r_ids = None

    def test_pattern_len(self):
        r_id2coeffs = [{'r0': 1, 'r1': 1, 'r2': 1, 'r3': 1},
                       {'r0': 1, 'r1': 1, 'r4': 1, 'r5': 1},
                       {'r00': 1, 'r0': -1},
                       {'r00': 1, 'r1': 1, 'r4': 1, 'r6': 1, 'r3': 1},
                       {'r00': 1, 'r1': 1, 'r4': 1, 'r5': 1}]
        id2efm = dict(zip(xrange(1, len(r_id2coeffs) + 1),
                          (EFM(r_ids=self.r_ids, rev_r_ids=self.rev_r_ids, r_id2coeff=it) for it in r_id2coeffs)))
        p_id2efm_ids, id2pattern, _ = detect_patterns(id2efm, 3, 2)
        # Expect to detect the following patterns:
        # {'r1': 1, 'r4': 1, 'r5': 1}
        # {'r1': 1, 'r4': 1, 'r00': 1}
        p_len = len(id2pattern)
        self.assertEqual(p_len, 2, "Expected to find 3 patterns, found %s" % [it.to_r_id2coeff() for it in id2pattern.itervalues()])

    def test_patterns(self):
        r_id2coeffs = [{'r0': 1, 'r1': 1, 'r2': 1, 'r3': 1},
                       {'r0': 1, 'r1': 1, 'r4': 1, 'r5': 1},
                       {'r00': 1, 'r0': -1},
                       {'r00': 1, 'r1': 1, 'r4': 1, 'r6': 1, 'r3': 1},
                       {'r00': 1, 'r1': 1, 'r4': 1, 'r5': 1}]
        id2efm = dict(zip(xrange(1, len(r_id2coeffs) + 1),
                          (EFM(r_ids=self.r_ids, rev_r_ids=self.rev_r_ids, r_id2coeff=it) for it in r_id2coeffs)))
        p_id2efm_ids, id2pattern, _ = detect_patterns(id2efm, 3, 2)
        # Expect to detect the following patterns:
        # {'r1': 1, 'r4': 1, 'r5': 1}
        # {'r1': 1, 'r4': 1, 'r00': 1}
        self.assertIn({'r1': 1, 'r4': 1, 'r5': 1}, [it.to_r_id2coeff() for it in id2pattern.itervalues()],
                      "Did not find pattern {'r1' : 1, 'r4' : 1, 'r5': 1}")
import unittest
from mod_sbml.constraint_based_analysis.efm.EFM import EFM
from efm.efm_classification import classify_efms

__author__ = 'anna'


class PatternTestCase(unittest.TestCase):

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

    def test_efm_intersection1(self):
        r_id2coeff1 = {'r0': 1, 'r1': 1, 'r2': 1, 'r3': 1}
        r_id2coeff2 = {'r0': 1, 'r1': 1, 'r4': 1, 'r5': 1}
        efm1 = EFM(r_ids=self.r_ids, rev_r_ids=self.rev_r_ids, r_id2coeff=r_id2coeff1)
        efm2 = EFM(r_ids=self.r_ids, rev_r_ids=self.rev_r_ids, r_id2coeff=r_id2coeff2)
        pattern = efm2.intersection(efm1)
        intersection_r_id2coeff = {r_id: (1 if coeff > 0 else -1) for (r_id, coeff) in r_id2coeff2.iteritems() if
                                   r_id in r_id2coeff1 and r_id2coeff1[r_id] * coeff > 0}
        self.assertEqual(pattern.to_r_id2coeff(), intersection_r_id2coeff,
                         "Pattern was supposed to be %s, got %s instead"
                         % (intersection_r_id2coeff, pattern.to_r_id2coeff()))

    def test_efm_intersection2(self):
        r_id2coeff1 = {'r00': 1, 'r0': -1}
        r_id2coeff2 = {'r0': 1, 'r1': 1, 'r4': 1, 'r5': 1}
        efm1 = EFM(r_ids=self.r_ids, rev_r_ids=self.rev_r_ids, r_id2coeff=r_id2coeff1)
        efm2 = EFM(r_ids=self.r_ids, rev_r_ids=self.rev_r_ids, r_id2coeff=r_id2coeff2)
        pattern = efm1.intersection(efm2)
        intersection_r_id2coeff = {r_id: (1 if coeff > 0 else -1) for (r_id, coeff) in r_id2coeff2.iteritems() if
                                   r_id in r_id2coeff1 and r_id2coeff1[r_id] * coeff > 0}
        self.assertEqual(pattern.to_r_id2coeff(), intersection_r_id2coeff,
                         "Pattern was supposed to be %s, got %s instead"
                         % (intersection_r_id2coeff, pattern.to_r_id2coeff()))

    def test_efm_classification_pattern(self):
        r_id2coeffs = [{'r0': 1, 'r1': 1, 'r2': 1, 'r3': 1},
                       {'r0': 1, 'r1': 1, 'r4': 1, 'r5': 1},
                       {'r00': 1, 'r0': -1},
                       {'r00': 1, 'r1': 1, 'r4': 1, 'r6': 1, 'r3': 1},
                       {'r00': 1, 'r1': 1, 'r4': 1, 'r5': 1}]
        id2efm = dict(zip(xrange(1, len(r_id2coeffs) + 1),
                          (EFM(r_ids=self.r_ids, rev_r_ids=self.rev_r_ids, r_id2coeff=it) for it in r_id2coeffs)))
        p_id2efm_ids, id2pattern = classify_efms(id2efm, min_pattern_len=2, max_pattern_num=None)
        patterns = id2pattern.values()
        p = {'r0': 1, 'r1': 1}
        self.assertIn(p, [it.to_r_id2coeff() for it in patterns], "Pattern %s was not found" % p)

    def test_efm_classification_len(self):
        r_id2coeffs = [{'r0': 1, 'r1': 1, 'r2': 1, 'r3': 1},
                       {'r0': 1, 'r1': 1, 'r4': 1, 'r5': 1},
                       {'r00': 1, 'r0': -1},
                       {'r00': 1, 'r1': 1, 'r4': 1, 'r6': 1, 'r3': 1},
                       {'r00': 1, 'r1': 1, 'r4': 1, 'r5': 1}]
        id2efm = dict(zip(xrange(1, len(r_id2coeffs) + 1),
                          (EFM(r_ids=self.r_ids, rev_r_ids=self.rev_r_ids, r_id2coeff=it) for it in r_id2coeffs)))
        p_id2efm_ids, id2pattern = classify_efms(id2efm, min_pattern_len=2, max_pattern_num=None)
        # Expect to detect the following patterns:
        # {'r0': 1, 'r1': 1}
        # {'r1': 1, 'r3': 1}
        # {'r1': 1, 'r4': 1}
        p_num = len(id2pattern)
        self.assertEqual(p_num, 3, "Expected to find 3 patterns, found %s" % p_num)

    def test_efm_classification_with_min_efm_num(self):
        r_id2coeffs = [{'r0': 1, 'r1': 1, 'r2': 1, 'r3': 1},
                       {'r0': 1, 'r1': 1, 'r4': 1, 'r5': 1},
                       {'r00': 1, 'r0': -1},
                       {'r00': 1, 'r1': 1, 'r4': 1, 'r6': 1, 'r3': 1},
                       {'r00': 1, 'r1': 1, 'r4': 1, 'r5': 1}]
        id2efm = dict(zip(xrange(1, len(r_id2coeffs) + 1),
                          (EFM(r_ids=self.r_ids, rev_r_ids=self.rev_r_ids, r_id2coeff=it) for it in r_id2coeffs)))
        p_id2efm_ids, id2pattern = classify_efms(id2efm, min_pattern_len=2, max_pattern_num=None, min_efm_num=3)
        # Expect to detect the following pattern:
        # {'r1': 1, 'r4': 1}
        patterns = [p.to_r_id2coeff() for p in id2pattern.itervalues()]
        expected_patterns = [{'r1': 1, 'r4': 1}]
        self.assertEqual(expected_patterns, patterns, "Expected to find %s, found %s" % (expected_patterns, patterns))

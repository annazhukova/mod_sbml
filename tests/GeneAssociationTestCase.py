import unittest

from mod_sbml.sbml.sbml_manager import _filter, parse_gene_association


class GATestCase(unittest.TestCase):

    def test_simple_in(self):
        res = _filter(parse_gene_association('A', flatten=False), ['A', 'B'], flatten=True)
        self.assertEqual(res, 'A', 'Was expecting A, got %s' % res)

    def test_simple_out(self):
        res = _filter(parse_gene_association('C', flatten=False), ['A', 'B'], flatten=True)
        self.assertIsNone(res, 'Was expecting None, got %s' % res)

    def test_or_all_in(self):
        res = _filter(parse_gene_association('A or B', flatten=False), ['A', 'B'], flatten=True)
        self.assertIn(res, ['(A or B)', '(B or A)'], 'Was expecting "A or B", got %s' % res)

    def test_or_some_in(self):
        res = _filter(parse_gene_association('A or C', flatten=False), ['A', 'B'], flatten=True)
        self.assertEqual('A', res, 'Was expecting "A", got %s' % res)

    def test_or_none_in(self):
        res = _filter(parse_gene_association('D or C', flatten=False), ['A', 'B'], flatten=True)
        self.assertIsNone(res, 'Was expecting None, got %s' % res)

    def test_and_all_in(self):
        res = _filter(parse_gene_association('A and B', flatten=False), ['A', 'B'], flatten=True)
        self.assertIn(res, ['(A and B)', '(B and A)'], 'Was expecting "A and B", got %s' % res)

    def test_and_some_in(self):
        res = _filter(parse_gene_association('A and C', flatten=False), ['A', 'B'], flatten=True)
        self.assertIsNone(res, 'Was expecting None, got %s' % res)

    def test_and_none_in(self):
        res = _filter(parse_gene_association('D and C', flatten=False), ['A', 'B'], flatten=True)
        self.assertIsNone(res, 'Was expecting None, got %s' % res)

    def test_complicated_expression_in(self):
        res = _filter(parse_gene_association('A or B and C', flatten=False), ['A', 'B'], flatten=True)
        self.assertEqual('A', res, 'Was expecting "A", got %s' % res)

    def test_complicated_expression_out(self):
        res = _filter(parse_gene_association('A and D or B and C', flatten=False), ['A', 'B'], flatten=True)
        self.assertIsNone(res, 'Was expecting None, got %s' % res)

    def test_complicated_expression(self):
        res = _filter(parse_gene_association('(A or B or C) and D', flatten=False), ['A', 'B', 'D'], flatten=True)
        self.assertIn(res, ['((A or B) and D)', '(D and (A or B))', '((B or A) and D)', '(D and (B or A))'],
                      'Was expecting "((A or B) and D)", got %s' % res)
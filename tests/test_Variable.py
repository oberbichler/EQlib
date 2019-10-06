import unittest
import EQlib as eq
import numpy as np
from numpy.testing import assert_almost_equal

class TestVariable(unittest.TestCase):

    def test_constructor(self):
        variable = eq.Variable(value=5)

        self.assertEqual(variable.ref_value, 5)
        self.assertEqual(variable.act_value, 5)
        self.assertEqual(variable.lower_bound, float('-inf'))
        self.assertEqual(variable.upper_bound, float('inf'))
        self.assertTrue(variable.is_active)
        self.assertEqual(variable.multiplier, 1)
        self.assertEqual(variable.name, '')

    def test_constructor_with_inactive(self):
        variable = eq.Variable(value=5, is_active=False)

        self.assertEqual(variable.ref_value, 5)
        self.assertEqual(variable.act_value, 5)
        self.assertEqual(variable.lower_bound, float('-inf'))
        self.assertEqual(variable.upper_bound, float('inf'))
        self.assertFalse(variable.is_active)
        self.assertEqual(variable.multiplier, 1)
        self.assertEqual(variable.name, '')

    def test_constructor_with_boundaries(self):
        variable = eq.Variable(value=5, lower_bound=-1, upper_bound=8)

        self.assertEqual(variable.ref_value, 5)
        self.assertEqual(variable.act_value, 5)
        self.assertEqual(variable.lower_bound, -1)
        self.assertEqual(variable.upper_bound, 8)
        self.assertTrue(variable.is_active)
        self.assertEqual(variable.multiplier, 1)
        self.assertEqual(variable.name, '')

    def test_constructor_with_name(self):
        variable = eq.Variable(value=5, name='test')

        self.assertEqual(variable.ref_value, 5)
        self.assertEqual(variable.act_value, 5)
        self.assertEqual(variable.lower_bound, float('-inf'))
        self.assertEqual(variable.upper_bound, float('inf'))
        self.assertTrue(variable.is_active)
        self.assertEqual(variable.multiplier, 1)
        self.assertEqual(variable.name, 'test')

if __name__ == '__main__':
    unittest.main()
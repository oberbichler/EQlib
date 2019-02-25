import unittest
import EQlib as eq
from numpy.testing import assert_almost_equal

class TestElement(unittest.TestCase):

    def test_element__dynamic_attributes(self):
        element = eq.Element()

        element.test_attribute = 'test value'

        assert(hasattr(element, 'test_attribute'))
        self.assertEqual(element.test_attribute, 'test value')

    def test_element_override(self):
        class Element(eq.Element):
            def __init__(self, a, b):
                self.a = a
                self.b = b

            def dofs(self):
                return [a.dof, b.dof]

            def compute(self):
                return [[1, 2], [3, 4]], [5, 6]

        a = eq.Parameter(0)
        b = eq.Parameter(0)

        element = Element(a, b)

        self.assertSequenceEqual(element.dofs(), [a.dof, b.dof])

        lhs, rhs = element.compute()

        assert_almost_equal(lhs, [[1, 2], [3, 4]])
        assert_almost_equal(rhs, [5, 6])

    def test_element_pass_options(self):
        class Element(eq.Element):
            def compute(self, options):
                return options['lhs'], options['rhs']

        element = Element()

        options = {
            'lhs': [[1, 2], [3, 4]],
            'rhs': [5, 6],
        }

        lhs, rhs = element.compute(options)

        assert_almost_equal(lhs, [[1, 2], [3, 4]])
        assert_almost_equal(rhs, [5, 6])

if __name__ == '__main__':
    unittest.main()
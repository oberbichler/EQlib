import unittest
import EQlib as eq
from numpy.testing import assert_almost_equal, assert_equal

class TestSystem(unittest.TestCase):

    def test_dofs(self):
        class Element(eq.Element):
            def __init__(self, a, b):
                eq.Element.__init__(self)
                self.a = a
                self.b = b

            def dofs(self):
                return [self.a, self.b]

        dof_1 = eq.Parameter()
        dof_2 = eq.Parameter()
        dof_3 = eq.Parameter()

        element_1 = Element(dof_1, dof_2)
        element_2 = Element(dof_2, dof_3)

        system = eq.System([element_1, element_2])

        self.assertSequenceEqual(system.dofs, [dof_1, dof_2, dof_3])

if __name__ == '__main__':
    unittest.main()
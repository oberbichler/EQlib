import unittest
import EQlib as eq
from numpy.testing import assert_almost_equal, assert_equal

class TestNode(unittest.TestCase):

    def test_dynamic_attributes(self):
        node = eq.Node()

        node.test_attribute = 'test value'

        assert(hasattr(node, 'test_attribute'))
        self.assertEqual(node.test_attribute, 'test value')

    def test_constructor(self):
        node = eq.Node(1, 2, 3)

        self.assertEqual(node.ref_x.value, 1)
        self.assertEqual(node.ref_y.value, 2)
        self.assertEqual(node.ref_z.value, 3)

        self.assertEqual(node.x.value, 1)
        self.assertEqual(node.y.value, 2)
        self.assertEqual(node.z.value, 3)

    def test_getters_and_setters(self):
        node = eq.Node(1, 2, 3)

        node.displacements = [2, 3, 4]

        assert_equal(node.ref_location, [1, 2, 3])
        assert_equal(node.act_location, [3, 5, 7])
        assert_equal(node.displacements, [2, 3, 4])

        node.forces = [5, 7, 9]

        assert_equal(node.forces, [5, 7, 9])

if __name__ == '__main__':
    unittest.main()
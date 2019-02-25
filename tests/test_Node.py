import unittest
import EQlib as eq
from numpy.testing import assert_almost_equal, assert_equal

class TestNode(unittest.TestCase):

    def test_node__dynamic_attributes(self):
        node = eq.Node()

        node.test_attribute = 'test value'

        assert(hasattr(node, 'test_attribute'))
        self.assertEqual(node.test_attribute, 'test value')

    def test_node_constructor(self):
        node = eq.Node(1, 2, 3)

        self.assertEqual(node.x.ref_value, 1)
        self.assertEqual(node.y.ref_value, 2)
        self.assertEqual(node.z.ref_value, 3)

        self.assertEqual(node.x.act_value, 1)
        self.assertEqual(node.y.act_value, 2)
        self.assertEqual(node.z.act_value, 3)

        self.assertEqual(node.x.delta, 0)
        self.assertEqual(node.y.delta, 0)
        self.assertEqual(node.z.delta, 0)

        self.assertEqual(node.x.target, 0)
        self.assertEqual(node.y.target, 0)
        self.assertEqual(node.z.target, 0)

        self.assertEqual(node.x.result, 0)
        self.assertEqual(node.y.result, 0)
        self.assertEqual(node.z.result, 0)

        self.assertEqual(node.x.residual, 0)
        self.assertEqual(node.y.residual, 0)
        self.assertEqual(node.z.residual, 0)

    def test_node_getters_and_setters(self):
        node = eq.Node(1, 2, 3)

        node.displacements = [2, 3, 4]

        assert_equal(node.ref_location, [1, 2, 3])
        assert_equal(node.act_location, [3, 5, 7])
        assert_equal(node.displacements, [2, 3, 4])

        node.forces = [5, 7, 9]

        assert_equal(node.forces, [5, 7, 9])

if __name__ == '__main__':
    unittest.main()
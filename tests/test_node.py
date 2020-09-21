import eqlib as eq

import pytest

from numpy.testing import assert_equal

if __name__ == '__main__':
    import sys
    import os
    print(f'pid: {os.getpid()}')
    pytest.main(sys.argv)


@pytest.fixture
def node():
    return eq.Node(1, 2, 3)


def test_init(node):
    assert_equal(node.ref_x.value, 1)
    assert_equal(node.ref_y.value, 2)
    assert_equal(node.ref_z.value, 3)

    assert_equal(node.x.value, 1)
    assert_equal(node.y.value, 2)
    assert_equal(node.z.value, 3)


def test_dynamic_attributes(node):
    node.test_attribute = 'test value'

    assert_equal(hasattr(node, 'test_attribute'), True)
    assert_equal(node.test_attribute, 'test value')


def test_getters_and_setters(node):
    node.displacements = [2, 3, 4]

    assert_equal(node.ref_location, [1, 2, 3])
    assert_equal(node.act_location, [3, 5, 7])
    assert_equal(node.displacements, [2, 3, 4])

    node.forces = [5, 7, 9]

    assert_equal(node.forces, [5, 7, 9])

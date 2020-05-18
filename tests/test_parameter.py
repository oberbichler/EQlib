import eqlib as eq

import pytest

from numpy.testing import assert_equal

if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)


def test_init():
    parameter = eq.Variable(value=5)

    assert_equal(parameter.value, 5)
    assert_equal(parameter.name, '')


def test_init_with_name():
    parameter = eq.Variable(value=5, name='test')

    assert_equal(parameter.value, 5)
    assert_equal(parameter.name, 'test')

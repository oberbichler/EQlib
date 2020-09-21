import eqlib as eq

import pytest

from numpy.testing import assert_equal

if __name__ == '__main__':
    import sys
    import os
    print(f'pid: {os.getpid()}')
    pytest.main(sys.argv)


def test_init():
    variable = eq.Variable(value=5)

    assert_equal(variable.value, 5)
    assert_equal(variable.lower_bound, float('-inf'))
    assert_equal(variable.upper_bound, float('inf'))
    assert_equal(variable.is_active, True)
    assert_equal(variable.multiplier, 1)
    assert_equal(variable.name, '')


def test_init_with_active():
    variable = eq.Variable(value=5, is_active=False)

    assert_equal(variable.value, 5)
    assert_equal(variable.lower_bound, float('-inf'))
    assert_equal(variable.upper_bound, float('inf'))
    assert_equal(variable.is_active, False)
    assert_equal(variable.multiplier, 1)
    assert_equal(variable.name, '')


def test_init_with_boundaries():
    variable = eq.Variable(value=5, lower_bound=-1, upper_bound=8)

    assert_equal(variable.value, 5)
    assert_equal(variable.lower_bound, -1)
    assert_equal(variable.upper_bound, 8)
    assert_equal(variable.is_active, True)
    assert_equal(variable.multiplier, 1)
    assert_equal(variable.name, '')


def test_init_with_name():
    variable = eq.Variable(value=5, name='test')

    assert_equal(variable.value, 5)
    assert_equal(variable.lower_bound, float('-inf'))
    assert_equal(variable.upper_bound, float('inf'))
    assert_equal(variable.is_active, True)
    assert_equal(variable.multiplier, 1)
    assert_equal(variable.name, 'test')

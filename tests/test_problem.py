import pytest

import eqlib as eq
import numpy as np

from numpy.testing import assert_equal

if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)


class ZeroObjective(eq.Objective):
    def __init__(self, variables):
        eq.Objective.__init__(self)
        self.variables = variables

    def compute(self, g, h):
        g.fill(0)
        h.fill(0)
        return 0


class ConstantObjective(eq.Objective):
    def __init__(self, variables, f, g, h):
        eq.Objective.__init__(self)
        self.variables = variables
        self.f = f
        self.g = g
        self.h = h

    def compute(self, g, h):
        g[:] = self.g
        h[:] = self.h
        return self.f


@pytest.fixture
def zero_problem():
    x1 = eq.Variable(name='x1', value=2.0)
    x2 = eq.Variable(name='x2', value=7.0)
    x3 = eq.Variable(name='x3', value=9.0)

    elements = [
        ZeroObjective([x1, x2]),
        ZeroObjective([x2, x3]),
    ]

    problem = eq.Problem(elements)

    return (elements, problem)


@pytest.fixture
def problem():
    x1 = eq.Variable(name='x1', value=2.0)
    x2 = eq.Variable(name='x2', value=7.0)
    x3 = eq.Variable(name='x3', value=9.0)

    elements = [
        ConstantObjective([x1, x2], 1, [2, 3], [[4, -5.6], [-5.6, 6.2]]),
        ConstantObjective([x2, x3], 3, [5.1, 7], [[-2, 6], [6, 11.3]]),
    ]

    problem = eq.Problem(elements)

    return (elements, problem)


def test_hm_add_diagonal(zero_problem):
    elements, problem = zero_problem

    problem.compute()

    assert_equal(problem.f, 0)
    assert_equal(problem.df, np.zeros(3))
    assert_equal(problem.hm.toarray(), np.zeros((3, 3)))

    problem.hm_add_diagonal(1)

    assert_equal(problem.f, 0)
    assert_equal(problem.df, np.zeros(3))
    assert_equal(problem.hm.toarray(), np.eye(3))


def test_hm_norm_inf(problem):
    elements, problem = problem

    problem.compute()

    assert_equal(problem.f, 4)
    assert_equal(problem.df, [2, 8.1, 7])
    assert_equal(problem.hm.toarray(), [[4, -5.6, 0], [0, 4.2, 6], [0, 0, 11.3]])

    norm = problem.hm_norm_inf

    assert_equal(norm, 17.3)

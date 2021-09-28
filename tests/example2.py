import eqlib as eq

import numpy as np
import pytest

from numpy.testing import assert_equal

if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)


class ConstantObjective(eq.Objective):
    def __init__(self, variables, f, df, hm):
        eq.Objective.__init__(self)
        self.variables = variables
        self.f = f
        self.df = df
        self.hm = hm

    def compute(self, df, hm, info):
        if len(df) != 0:
            df[:] = self.df
        if len(hm) != 0:
            hm[:] = self.hm
        return self.f


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

    return problem


def test_init(problem):
    assert_equal(problem.nb_variables, 3)
    assert_equal(problem.nb_equations, 0)

    assert_equal(problem.structure_hm.rows, 3)
    assert_equal(problem.structure_hm.cols, 3)
    assert_equal(problem.structure_hm.ia, [0, 2, 4, 5])
    assert_equal(problem.structure_hm.ja, [0, 1, 1, 2, 2])
    assert_equal(problem.structure_hm.nb_nonzeros, 5)

    assert_equal(len(problem.hm_values), 5)


def test_variable_index(problem):
    for i, variable in enumerate(problem.variables):
        assert_equal(problem.variable_index(variable), i)


def test_nb_objectives(problem):
    assert_equal(problem.nb_objectives, 2)


def test_nb_constraints(problem):
    assert_equal(problem.nb_constraints, 0)


def test_eval(problem):
    problem.eval(problem.x)

    assert_equal(problem.f, 4)
    assert_equal(problem.df, [2, 8.1, 7])
    assert_equal(problem.hm.toarray(), [[4, -5.6, 0], [0, 4.2, 6], [0, 0, 11.3]])


def test_hm_diagonal(problem):
    problem.eval(problem.x)

    assert_equal(problem.hm_diagonal, [4, 4.2, 11.3])

    problem.hm_diagonal = [1, 2, 3]

    assert_equal(problem.hm.toarray(), [[1, -5.6, 0], [0, 2, 6], [0, 0, 3]])


def test_hm_add_diagonal(problem):
    problem.eval(problem.x)

    problem.hm_add_diagonal(1)

    assert_equal(problem.f, 4)
    assert_equal(problem.df, [2, 8.1, 7])
    assert_equal(problem.hm.toarray(),
                 [[5, -5.6, 0], [0, 5.2, 6], [0, 0, 12.3]])


def test_hm_norm_inf(problem):
    problem.eval(problem.x)

    norm = problem.hm_norm_inf

    assert_equal(norm, 17.3)


def test_scale_hm(problem):
    problem.eval(problem.x)

    problem.scale_hm(1.5)

    assert_equal(problem.hm.toarray(), np.multiply([[4, -5.6, 0], [0, 4.2, 6], [0, 0, 11.3]], 1.5))

import unittest
import eqlib as eq
import hyperjet as hj
import numpy as np
from numpy.testing import assert_equal


def explode(value, g, h):
    g[:] = value.g
    h[:] = value.h
    return value.f


class TestLambdaElements(unittest.TestCase):

    def test_lambda_constraint(self):
        g1 = eq.Equation()
        g2 = eq.Equation()

        x1 = eq.Variable(1)
        x2 = eq.Variable(2)

        def compute(equations, variables, fs, gs, hs):
            x1, x2 = hj.HyperJet.variables(variables)
            f1 = x1**2 + x1 * x2 + x2**3
            f2 = x1**3 + x1 * x2 + x2**2
            fs[0] = explode(f1, gs[0], hs[0])
            fs[1] = explode(f2, gs[1], hs[1])

        element = eq.LambdaConstraint([g1, g2], [x1, x2], compute)

        fs = np.empty(2, float)
        gs = np.empty((2, 2), float)
        hs = np.empty((2, 2, 2), float)

        element.compute(fs, gs, hs)

        assert_equal(fs, [11, 7])
        assert_equal(gs, [[4, 13], [5, 5]])
        assert_equal(hs, [[[2, 1], [1, 12]], [[6, 1], [1, 2]]])

    def test_lambda_objective(self):
        x1 = eq.Variable(1)
        x2 = eq.Variable(2)

        def compute(variables, g, h):
            x1, x2 = hj.HyperJet.variables(variables)
            f = x1**2 + x1 * x2 + x2**3
            return explode(f, g, h)

        element = eq.LambdaObjective([x1, x2], compute)

        g = np.empty(2)
        h = np.empty((2, 2))

        f = element.compute(g, h)

        assert_equal(f, 11)
        assert_equal(g, [4, 13])
        assert_equal(h, [[2, 1], [1, 12]])


if __name__ == '__main__':
    unittest.main()

import unittest
import EQlib as eq
import numpy as np
import hyperjet as hj
from numpy.testing import assert_almost_equal, assert_equal

class F1(eq.Objective):
    def __init__(self, x1, x2):
        eq.Objective.__init__(self)
        self.x1 = x1
        self.x2 = x2

    def variables(self):
        return [self.x1, self.x2]

    def compute(self, g, h):
        x1, x2 = hj.HyperJet.variables([self.x1.act_value, self.x2.act_value])
        r = x1**2 + x2**2
        g[:] = r.g
        h[:] = r.h
        return r.f

class F2(eq.Objective):
    def __init__(self, x2, x3):
        eq.Objective.__init__(self)
        self.x2 = x2
        self.x3 = x3

    def variables(self):
        return [self.x2, self.x3]

    def compute(self, g, h):
        x2, x3 = hj.HyperJet.variables([self.x2.act_value, self.x3.act_value])
        r = x2**2 - x3
        g[:] = r.g
        h[:] = r.h
        return r.f

class C1(eq.Constraint):
    def __init__(self, g1, x1, x3):
        eq.Constraint.__init__(self)
        self.g1 = g1
        self.x1 = x1
        self.x3 = x3

    def equations(self):
        return [self.g1]

    def variables(self):
        return [self.x1, self.x3]

    def compute(self, fs, gs, hs):
        x1, x3 = hj.HyperJet.variables([self.x1.act_value, self.x3.act_value])
        rs = [x1**2 + x1 * x3]
        for k in range(len(rs)):
            fs[k] = rs[k].f
            gs[k][:] = rs[k].g
            hs[k][:] = rs[k].h

class C2(eq.Constraint):
    def __init__(self, g1, g2, x3):
        eq.Constraint.__init__(self)
        self.g1 = g1
        self.g2 = g2
        self.x3 = x3

    def equations(self):
        return [self.g1, self.g2]

    def variables(self):
        return [self.x3]

    def compute(self, fs, gs, hs):
        x3, = hj.HyperJet.variables([self.x3.act_value])
        rs = [x3**2, x3]
        for k in range(len(rs)):
            fs[k] = rs[k].f
            gs[k][:] = rs[k].g
            hs[k][:] = rs[k].h

class C3(eq.Constraint):
    def __init__(self, g2, g3, x4):
        eq.Constraint.__init__(self)
        self.g2 = g2
        self.g3 = g3
        self.x4 = x4

    def equations(self):
        return [self.g2, self.g3]

    def variables(self):
        return [self.x4]

    def compute(self, fs, gs, hs):
        x4, = hj.HyperJet.variables([self.x4.act_value])
        rs = [-x4, x4]
        for k in range(len(rs)):
            fs[k] = rs[k].f
            gs[k][:] = rs[k].g
            hs[k][:] = rs[k].h

class C4(eq.Constraint):
    def __init__(self, g3, x2):
        eq.Constraint.__init__(self)
        self.g3 = g3
        self.x2 = x2

    def equations(self):
        return [self.g3]

    def variables(self):
        return [self.x2]

    def compute(self, fs, gs, hs):
        x2, = hj.HyperJet.variables([self.x2.act_value])
        rs = [x2]
        for k in range(len(rs)):
            fs[k] = rs[k].f
            gs[k][:] = rs[k].g
            hs[k][:] = rs[k].h

class TestProblem(unittest.TestCase):

    def test_example(self):
        g1 = eq.Equation(name='g1', lower_bound=1, upper_bound=1, multiplier=3.2)
        g2 = eq.Equation(name='g2', lower_bound=float('-inf'), upper_bound=-1, multiplier=9.3)
        g3 = eq.Equation(name='g3', lower_bound=2.5, upper_bound=5, multiplier=11.6)

        x1 = eq.Variable(name='x1', value=2.0, lower_bound=-0.5, upper_bound=float('inf'))
        x2 = eq.Variable(name='x2', value=7.0, lower_bound=-2.0, upper_bound=float('inf'))
        x3 = eq.Variable(name='x3', value=1.0, lower_bound= 0.0, upper_bound=2)
        x4 = eq.Variable(name='x4', value=5.0, lower_bound=-2.0, upper_bound=2)

        objectives = [F1(x1, x2), F2(x2, x3)]
        constraints = [C1(g1, x1, x3), C2(g1, g2, x3), C3(g2, g3, x4), C4(g3, x2)]

        problem = eq.Problem(objectives, constraints)

        problem.sigma = 1.5

        assert_equal(problem.nb_equations, 3)
        assert_equal(problem.nb_variables, 4)

        assert_equal(problem.dg.nonzero(), [[0, 0, 1, 1, 2, 2], [0, 2, 2, 3, 1, 3]])
        assert_equal(problem.hl.nonzero(), [[0, 1, 1, 2, 2, 2, 3], [0, 0, 1, 0, 1, 2, 3]])

        assert_almost_equal(problem.equation_multipliers, [3.2, 9.3, 11.6])
        assert_almost_equal(problem.variable_multipliers, [1, 1, 1, 1])

        problem.compute()

        assert_almost_equal(problem.f, 151.5)
        assert_almost_equal(problem.g, [7, -4, 12])
        assert_almost_equal(problem.df, [6, 42, -1.5, 0])
        assert_almost_equal(problem.dg.todense(), [[5, 0, 4, 0], [0, 0, 1, -1], [0, 1, 0, 1]])
        assert_almost_equal(problem.hl.todense(), [[9.4, 0, 0, 0], [0, 6, 0, 0], [3.2, 0, 6.4, 0], [0, 0, 0, 0]])

if __name__ == '__main__':
    unittest.main()
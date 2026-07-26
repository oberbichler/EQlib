"""A LinearSolver subclassed in Python must have its ``solve`` dispatched.

``Problem.hm_inv_v`` calls ``factorize(ia, ja, a)`` and then
``solve(ia, ja, a, b, x)`` on its linear solver. For a solver implemented in
Python this goes through the pybind11 trampoline ``PyLinearSolver``. The
trampoline's ``solve`` override used to dispatch to the Python method named
``factorize`` (with 5 arguments) instead of ``solve``, so a Python ``solve``
implementation was never invoked.
"""

import eqlib as eq
import numpy as np
import pytest
from numpy.testing import assert_almost_equal

if __name__ == "__main__":
    import os
    import sys

    print(f"pid: {os.getpid()}")
    pytest.main(sys.argv)


class ConstantObjective(eq.Objective):
    def __init__(self, variables, f, g, h):
        eq.Objective.__init__(self)
        self.variables = variables
        self._f = f
        self._g = g
        self._h = h

    def compute(self, g, h):
        if len(g) != 0:
            g[:] = self._g
        if len(h) != 0:
            h[:] = self._h
        return self._f


class RecordingSolver(eq.LinearSolver):
    """Records how it is called and applies a trivial, matrix-independent
    transform (``x = b + 1``) so the result is easy to check."""

    def __init__(self):
        eq.LinearSolver.__init__(self)
        self.factorize_args = []
        self.solve_args = []

    def analyze(self, ia, ja, a):
        return False

    # Accept extra positional args so a misrouted 5-argument call does not raise
    # a TypeError but is recorded instead, giving a clean assertion failure.
    def factorize(self, ia, ja, a, *extra):
        self.factorize_args.append(3 + len(extra))
        return False

    def solve(self, ia, ja, a, b, x):
        self.solve_args.append(np.asarray(b).copy())
        x[:] = np.asarray(b) + 1.0
        return False


def test_python_linear_solver_solve_is_dispatched():
    x1 = eq.Variable(name="x1", value=1.0)
    x2 = eq.Variable(name="x2", value=2.0)

    problem = eq.Problem(
        [ConstantObjective([x1, x2], 0.0, [0.0, 0.0], [[2.0, 0.0], [0.0, 3.0]])]
    )
    problem.compute()

    solver = RecordingSolver()
    problem.linear_solver = solver

    b = np.array([1.0, 2.0])
    result = problem.hm_inv_v(b)

    # solve must be dispatched exactly once, with the right rhs.
    assert len(solver.solve_args) == 1
    assert_almost_equal(solver.solve_args[0], b)

    # factorize is called once by hm_inv_v itself, with its 3-argument signature
    # (i.e. never as a misrouted 5-argument solve call).
    assert solver.factorize_args == [3]

    # the x written by solve flows back to the caller.
    assert_almost_equal(result, b + 1.0)

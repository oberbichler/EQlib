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


class ConstantConstraint(eq.Constraint):
    def __init__(self, equations, variables, fs, gs, hs):
        eq.Constraint.__init__(self)
        self.equations = equations
        self.variables = variables
        self._fs = fs
        self._gs = gs
        self._hs = hs

    def compute(self, fs, gs, hs):
        for k in range(len(fs)):
            fs[k] = self._fs[k]
            if len(gs[k]) != 0:
                gs[k][:] = self._gs[k]
            if len(hs[k]) != 0:
                hs[k][:] = self._hs[k]


@pytest.fixture
def problem():
    g1 = eq.Equation(name="g1", lower_bound=1, upper_bound=1, multiplier=2.0)

    x1 = eq.Variable(name="x1", value=2.0)
    x2 = eq.Variable(name="x2", value=3.0)

    objectives = [
        ConstantObjective([x1, x2], 1.0, [2.0, 3.0], [[4.0, 0.0], [0.0, 5.0]])
    ]
    constraints = [
        ConstantConstraint(
            [g1], [x1, x2], [10.0], [[2.0, 3.0]], [[[6.0, 0.0], [0.0, 0.0]]]
        )
    ]

    return eq.Problem(objectives, constraints)


def test_remove_inactive_elements_preserves_active_constraints(problem):
    """When every element is active, ``remove_inactive_elements`` must be a
    no-op: the assembled ``g``/``dg``/``hm`` must be identical before and after.

    A copy-paste defect in ``remove_inactive_constraints`` used to reset the
    per-element equation/variable index maps to empty, silently dropping every
    constraint's contribution from ``g`` and ``dg`` after the call.
    """
    problem.compute()

    g_before = problem.g.copy()
    dg_before = problem.dg.toarray()
    hm_before = problem.hm.toarray()
    df_before = problem.df.copy()

    # Guard against a vacuous test: the constraint must actually contribute.
    assert np.any(g_before != 0.0)
    assert np.any(dg_before != 0.0)

    problem.remove_inactive_elements()

    problem.compute()

    assert_almost_equal(problem.g, g_before)
    assert_almost_equal(problem.dg.toarray(), dg_before)
    assert_almost_equal(problem.hm.toarray(), hm_before)
    assert_almost_equal(problem.df, df_before)

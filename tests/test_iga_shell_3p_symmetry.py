"""Independent correctness oracle for IgaShell3PAD's curvature ordering.

The existing IgaShell3PAD test asserts against baked reference values produced
by the element itself, and its reference configuration is flat (the reference
curvature b vanishes), so it cannot reveal a wrong curvature. A
finite-difference check of g against f is no help either: the AD differentiates
whatever energy expression is written, so g stays consistent with f even if the
expression is wrong.

This test uses a value-level physical invariant on a *curved* reference: for an
isotropic material the Kirchhoff-Love shell energy is invariant under renaming
the two parametric directions (u1 <-> u2). Relabelling only permutes the
shape-function rows (1<->2 for the first derivatives, 3<->5 for the second
derivatives); the energy must be unchanged.

The membrane metric and the transformation matrix are built in Voigt order
[11, 22, 12], and act_b uses that order, so the reference curvature ref_b must
use it too. If ref_b instead uses [11, 12, 22], the relabel symmetry breaks on
a curved reference.
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


# 6x6 identity shape functions: row k selects node k, so rows 1..5 map directly
# to a1, a2, a1_1, a1_2, a2_2 of nodes 1..5 (row 0 = values, unused by the shell).
SHAPE_FUNCTIONS = np.eye(6)

# Curved reference: the reference second-derivative vectors have a normal
# component with distinct b11=0.5, b22=0.9, b12=0.2, so the ref_b ordering
# actually matters (a flat reference would hide the bug).
REF_LOCATIONS = [
    [0.0, 0.0, 0.0],  # node 0 (unused by geometry)
    [2.0, 0.0, 0.0],  # a1
    [0.0, 3.0, 0.0],  # a2
    [0.0, 0.0, 0.5],  # a1_1
    [0.0, 0.0, 0.2],  # a1_2
    [0.0, 0.0, 0.9],  # a2_2
]
ACT_LOCATIONS = [
    [0.0, 0.0, 0.0],
    [2.1, 0.0, 0.0],
    [0.0, 3.05, 0.0],
    [0.0, 0.0, 0.6],
    [0.0, 0.0, 0.25],
    [0.0, 0.0, 1.0],
]

# relabel u1 <-> u2: rows [val, d/du1, d/du2, d2/du1^2, d2/du1du2, d2/du2^2]
PERM = [0, 2, 1, 5, 4, 3]


def _shell_energy(shape_functions):
    nodes = []
    for ref, act in zip(REF_LOCATIONS, ACT_LOCATIONS):
        node = eq.Node()
        node.ref_location = ref
        node.act_location = act
        nodes.append(node)

    element = eq.IgaShell3PAD(nodes, 1, 100, 0)  # thickness=1, E=100, nu=0 (isotropic)
    element.add(np.asarray(shape_functions), 1.5)
    return element.compute_all()[0]


def test_shell_energy_invariant_under_parametric_relabel():
    f_orig = _shell_energy(SHAPE_FUNCTIONS)
    f_swapped = _shell_energy(SHAPE_FUNCTIONS[PERM])

    # non-vacuous: the curved reference must carry bending energy
    assert abs(f_orig) > 1e-6

    assert_almost_equal(f_swapped, f_orig)

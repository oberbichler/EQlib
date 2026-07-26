"""Multiple integration points per AD objective must accumulate.

Each ``Iga*AD`` objective can hold several integration points (repeated
``add()`` calls). The assembled gradient ``g`` and Hessian ``h`` of an element
with N points must equal the sum of the contributions of N single-point
elements over the same nodes -- exactly as the scalar value ``f`` does.

The AD objectives used to write ``g``/``h`` with ``=`` inside the loop over the
integration points, so only the *last* point survived in ``g``/``h`` while
``f`` was summed over all points.
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


REF_A = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
ACT_A = [[0.0, 0.0, 0.05], [1.1, 0.0, 0.1], [0.0, 1.05, 0.2]]
REF_B = [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [0.0, 1.0, 1.0]]
ACT_B = [[0.02, 0.0, 1.0], [1.0, 0.03, 1.1], [0.0, 1.0, 0.95]]

# 3-row shape functions (values, d/du, d/dv) for the surface-derivative elements
SF1 = [[1 / 3, 1 / 3, 1 / 3], [-1.0, 1.0, 0.0], [-1.0, 0.0, 1.0]]
SF2 = [[0.3, 0.3, 0.4], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]]


def _make_nodes(ref_locations, act_locations):
    nodes = []
    for ref, act in zip(ref_locations, act_locations):
        node = eq.Node()
        node.ref_location = ref
        node.act_location = act
        nodes.append(node)
    return nodes


def _assert_multipoint_is_sum(element_both, element_0, element_1):
    f_both, g_both, h_both = element_both.compute_all()
    f0, g0, h0 = element_0.compute_all()
    f1, g1, h1 = element_1.compute_all()

    # non-vacuous: each single point must actually contribute a gradient
    assert np.any(g0 != 0.0)
    assert np.any(g1 != 0.0)

    assert_almost_equal(f_both, f0 + f1)
    assert_almost_equal(g_both, g0 + g1)
    assert_almost_equal(np.triu(h_both), np.triu(h0) + np.triu(h1))


def test_point_distance_ad_accumulates_multiple_points():
    nodes_a = _make_nodes(REF_A, ACT_A)
    nodes_b = _make_nodes(REF_B, ACT_B)

    sf_a1, sf_b1, w1 = [[0.2, 0.3, 0.5]], [[0.1, 0.6, 0.3]], 1.3
    sf_a2, sf_b2, w2 = [[0.5, 0.1, 0.4]], [[0.25, 0.25, 0.5]], 0.7

    both = eq.IgaPointDistanceAD(nodes_a, nodes_b)
    both.add(sf_a1, sf_b1, w1)
    both.add(sf_a2, sf_b2, w2)

    e0 = eq.IgaPointDistanceAD(nodes_a, nodes_b)
    e0.add(sf_a1, sf_b1, w1)
    e1 = eq.IgaPointDistanceAD(nodes_a, nodes_b)
    e1.add(sf_a2, sf_b2, w2)

    _assert_multipoint_is_sum(both, e0, e1)


def test_normal_distance_ad_accumulates_multiple_points():
    nodes_a = _make_nodes(REF_A, ACT_A)
    nodes_b = _make_nodes(REF_B, ACT_B)

    both = eq.IgaNormalDistanceAD(nodes_a, nodes_b)
    both.add(SF1, SF1, 1.3)
    both.add(SF2, SF2, 0.7)

    e0 = eq.IgaNormalDistanceAD(nodes_a, nodes_b)
    e0.add(SF1, SF1, 1.3)
    e1 = eq.IgaNormalDistanceAD(nodes_a, nodes_b)
    e1.add(SF2, SF2, 0.7)

    _assert_multipoint_is_sum(both, e0, e1)


def test_rotation_coupling_ad_accumulates_multiple_points():
    nodes_a = _make_nodes(REF_A, ACT_A)
    nodes_b = _make_nodes(REF_B, ACT_B)

    # axis must not be aligned with the surface normal a3, otherwise
    # omega . axis vanishes and the coupling has no gradient.
    axis = [1.0, 0.0, 0.0]

    both = eq.IgaRotationCouplingAD(nodes_a, nodes_b)
    both.add(SF1, SF1, axis, 1.3)
    both.add(SF2, SF2, axis, 0.7)

    e0 = eq.IgaRotationCouplingAD(nodes_a, nodes_b)
    e0.add(SF1, SF1, axis, 1.3)
    e1 = eq.IgaRotationCouplingAD(nodes_a, nodes_b)
    e1.add(SF2, SF2, axis, 0.7)

    _assert_multipoint_is_sum(both, e0, e1)

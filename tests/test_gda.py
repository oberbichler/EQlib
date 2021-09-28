from numpy.testing import assert_allclose, assert_equal
from gda import gradient_descent_akin, gradient_descent_akin_nesterov
from ieee import ieee_g4
import pytest

if __name__ == '__main__':
    import sys
    pytest.main(sys.argv)


def test_gda_ieee_g4():
    problem = ieee_g4()

    x0 = problem.x
    maxiter = 1000
    step_size = 0.2
    zeta = 0.999

    solution = gradient_descent_akin(problem, x0, maxiter, step_size, zeta)

    assert_allclose(
        solution.x,
        [78.07778778, 33.16535679, 30.1024419, 44.61602068, 36.65566796])
    assert_allclose(solution.f, -30633.57307)
    assert_allclose(solution.g, 0.0003849975963845509)
    assert_equal(solution.iterations, 412)

    print(solution)


def test_gdan_ieee_g4():
    problem = ieee_g4()

    x0 = problem.x
    maxiter = 1000
    step_size = 0.2
    zeta = 0.999
    momentum = 0.95

    solution = gradient_descent_akin_nesterov(problem, x0, maxiter, step_size, zeta, momentum)

    assert_allclose(solution.x, [79.58736123, 33.76108145, 30.99552582, 41.86309924, 35.81270578])
    assert_allclose(solution.f, -30398.686754857117)
    assert_allclose(solution.g, 0.39778516857504087)
    assert_equal(solution.iterations, 19)

    print(solution)

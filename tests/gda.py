from __future__ import annotations
import numpy as np


class SolutionInfo:
    def __init__(self, x: np.ndarray[float], f: float, g: np.ndarray[float],
                 iterations: int):
        self.x = x
        self.f = f
        self.g = g
        self.iterations = iterations

    def __repr__(self) -> str:
        return '\n'.join([
            f'Solution:',
            f'> x          = {self.x}',
            f'> f(x)       = {self.f}',
            f'> g(x)       = {self.g}',
            f'> iterations = {self.iterations}',
        ])


def compute_searchdirection(problem, x, zeta):
    problem.x = x

    problem.eval(x)

    g_lb = problem.x_lower_bounds - problem.x
    g_ub = problem.x - problem.x_upper_bounds

    df = problem.df

    g = problem.g

    weighted_dg = -problem.dg.toarray() / g[:, None]

    dg = np.sum(weighted_dg, axis=0) + 1 / g_lb - 1 / g_ub

    s = -df / np.linalg.norm(df) - zeta * dg / np.linalg.norm(dg)

    return s


def gradient_descent_akin(problem, x0, maxiter, step_size, zeta):
    # set initial point
    x = np.copy(x0)

    # track all solutions
    solutions = list()

    # run gradient descent akin method
    for iteration in range(maxiter):
        # calculate the chen's direction
        s = compute_searchdirection(problem, x, zeta)

        # normalize the search direction
        s = s / np.linalg.norm(s)

        # build a new solution
        x += step_size * s
        solutions.append((problem.f, np.copy(x)))

        # evaluate candidate point
        problem.eval(x)

        g_lb = problem.x_lower_bounds - x
        g_ub = x - problem.x_upper_bounds

        max_g_lb = np.amax(g_lb)
        max_g_ub = np.amax(g_ub)
        g = np.amax([*problem.g, max_g_lb, max_g_ub])

        if g < 0:
            continue

        f, _ = solutions[-1]
        _, x = solutions[-2]

        return SolutionInfo(x, f, g, iteration)

    raise RuntimeError('Iteration limit exceeded')


def gradient_descent_akin_nesterov(problem, x0, maxiter, step_size, zeta, momentum):
    # set initial point
    x = np.copy(x0)

    # track all solutions
    solutions = list()

    change = np.zeros_like(x)

    # run the gradient descent akin method with Nesterov momemtum
    for iteration in range(maxiter):
        # calculate the projected solution
        projected = x + momentum * change

        s = compute_searchdirection(problem, projected, zeta)

        # normalize the search direction
        s = s / np.linalg.norm(s)

        # build a new solution
        change[:] = momentum * change + step_size * s
        x += change
        solutions.append((problem.f, np.copy(x)))

        # evaluate candidate point
        problem.eval(x)

        g_lb = problem.x_lower_bounds - x
        g_ub = x - problem.x_upper_bounds

        max_g_lb = np.amax(g_lb)
        max_g_ub = np.amax(g_ub)
        g = np.amax([*problem.g, max_g_lb, max_g_ub])

        if g < 0:
            continue

        f, _ = solutions[-1]
        _, x = solutions[-2]

        return SolutionInfo(x, f, g, iteration)

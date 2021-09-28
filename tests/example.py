from __future__ import annotations
import eqlib as eq
import numpy as np
from numpy.testing import assert_allclose, assert_equal
from rich.console import Console

console = Console()

# ---

class SolutionInfo:
    def __init__(self, x: np.ndarray[float], f: float, g: np.ndarray[float],
                 iterations: int):
        self.x = x
        self.f = f
        self.g = g
        self.iterations = iterations

    def __repr__(self) -> str:
        return '\n'.join([
            f'> x          = {self.x}',
            f'> f(x)       = {self.f}',
            f'> g(x)       = {self.g}',
            f'> iterations = {self.iterations}',
        ])


class SearchDirection:
    def __init__(self, problem, zeta):
        self.problem = problem
        self.zeta = zeta
    
    def __call__(self, x):
        problem, zeta = self.problem, self.zeta

        problem.eval(x)

        g_lb = problem.x_lower_bounds - problem.x
        g_ub = problem.x - problem.x_upper_bounds

        df = problem.df

        g = problem.g

        weighted_dg = -problem.dg.toarray() / g[:, None]

        dg = np.sum(weighted_dg, axis=0) + 1 / g_lb - 1 / g_ub

        s = df / np.linalg.norm(df) + zeta * dg / np.linalg.norm(dg)
        s /= np.linalg.norm(s)

        return s


def gradient_descent_akin(problem, x0, zeta, method, maxiter):
    x = x0

    solutions = list()

    fn = SearchDirection(problem, zeta)

    # run gradient descent akin method
    for iteration in range(maxiter + 1):
        problem.eval(x)

        solutions.append((problem.f, np.copy(x)))

        g_lb = problem.x_lower_bounds - x
        g_ub = x - problem.x_upper_bounds

        max_g_lb = np.amax(g_lb)
        max_g_ub = np.amax(g_ub)
        g = np.amax([*problem.g, max_g_lb, max_g_ub])

        if g >= 0:
            break
        if iteration == maxiter:
            raise RuntimeError('Iteration limit exceeded')

        x = method.step(fn, x)

    f, x = solutions[-1]

    return SolutionInfo(x, f, g, iteration)


def solve(method, maxiter=4000):
    try:
        solution = gradient_descent_akin(problem, x0, zeta, method, maxiter=maxiter)

        console.print(f'{method.__class__.__name__}:', style='bold green')
        print(solution)
        print()

        return solution
    except Exception:
        console.print(f'{method.__class__.__name__}:', style='bold red')
        console.print_exception(show_locals=True)
        print()

        return None

# ---

from ieee import ieee_g4 as create_problem

problem = create_problem()

x0 = problem.x

zeta = 0.999

# ---

class GradientDescent:
    def __init__(self, problem, alpha):
        self.alpha = alpha

    def step(self, fn, x):
        alpha, g = self.alpha, fn(x)
        return x - alpha * g

method = GradientDescent(problem, alpha=0.2)

solution = solve(method)

assert_allclose(solution.x, [78.007875, 32.999615, 30.088128, 44.695283, 36.689646]) # [78.07778778, 33.16535679, 30.1024419, 44.61602068, 36.65566796])
assert_allclose(solution.f, -30640.722925) # -30633.57307)
assert_allclose(solution.g, 0.0003849975963845509)
assert_equal(solution.iterations, 413) # 412)

# ---

class Momentum:
    def __init__(self, problem, alpha, beta):
        self.alpha = alpha
        self.beta = beta
        self.v = np.zeros(problem.nb_variables)

    def step(self, fn, x):
        alpha, beta, v, g = self.alpha, self.beta, self.v, fn(x)
        v[:] = beta * v - alpha * g
        return x + v

method = Momentum(problem, alpha=0.2, beta=0.95)

solution = solve(method)

# ---

class NesterovMomentum:
    def __init__(self, problem, alpha, beta):
        self.alpha = alpha
        self.beta = beta
        self.v = np.zeros(problem.nb_variables)

    def step(self, fn, x):
        alpha, beta, v, g = self.alpha, self.beta, self.v, fn(x)
        v[:] = beta * v - alpha * fn(x + beta * v)
        return x + v

method = NesterovMomentum(problem, alpha=0.2, beta=0.95)

solution = solve(method)

assert_allclose(solution.x, [77.602215, 33.754973, 30.950484, 42.593188, 35.975437]) # [79.58736123, 33.76108145, 30.99552582, 41.86309924, 35.81270578])
assert_allclose(solution.f, -30432.585447) # -30398.686754857117)
assert_allclose(solution.g, 0.39778516857504087)
assert_equal(solution.iterations, 20) # 19)

# ---

class Adagrad:
    def __init__(self, problem, alpha, epsilon):
        self.alpha = alpha
        self.epsilon = epsilon
        self.s = np.zeros(problem.nb_variables)

    def step(self, fn, x):
        alpha, epsilon, s, g = self.alpha, self.epsilon, self.s, fn(x)
        s[:] += g * g
        return x - alpha * g / (np.sqrt(s) + epsilon)

method = Adagrad(problem, alpha=0.2, epsilon=1e-8)

solution = solve(method)

# ---

class RMSProp:
    def __init__(self, problem, alpha, gamma, epsilon):
        self.alpha = alpha
        self.gamma = gamma
        self.epsilon = epsilon
        self.s = np.zeros(problem.nb_variables)

    def step(self, fn, x):
        alpha, gamma, epsilon, s, g = self.alpha, self.gamma, self.epsilon, self.s, fn(x)
        s[:] = gamma * s + (1 - gamma) * (g * g)
        return x - alpha * g / (np.sqrt(s) + epsilon)

method = RMSProp(problem, alpha=0.2, gamma=0.9, epsilon=1e-8)

solution = solve(method)

# ---

class Adadelta:
    def __init__(self, problem, gamma_s, gamma_x, epsilon):
        self.gamma_s = gamma_s
        self.gamma_x = gamma_x
        self.epsilon = epsilon
        self.s = np.zeros(problem.nb_variables)
        self.u = np.zeros(problem.nb_variables)

    def step(self, fn, x):
        gamma_s, gamma_x, epsilon, s, u, g = self.gamma_s, self.gamma_x, self.epsilon, self.s, self.u, fn(x)
        s[:] = gamma_s * s + (1 - gamma_s) * (g * g)
        delta_x = -(np.sqrt(u) + epsilon) / (np.sqrt(s) + epsilon) * g
        u[:] = gamma_x * u + (1 - gamma_x) * delta_x * delta_x
        return x + delta_x

method = Adadelta(problem, gamma_s=0.999, gamma_x=0.9, epsilon=1e-8)

solution = solve(method)

# ---

class Adam:
    def __init__(self, problem, alpha, gamma_v, gamma_s, epsilon):
        self.alpha = alpha
        self.gamma_v = gamma_v
        self.gamma_s = gamma_s
        self.epsilon = epsilon
        self.k = 0
        self.v = np.zeros(problem.nb_variables)
        self.s = np.zeros(problem.nb_variables)

    def step(self, fn, x):
        alpha, gamma_v, gamma_s, epsilon, k = self.alpha, self.gamma_v, self.gamma_s, self.epsilon, self.k
        s, v, g = self.s, self.v, fn(x)
        v[:] = gamma_v * v + (1 - gamma_v) * g
        s[:] = gamma_s * s + (1 - gamma_s) * g * g
        self.k = k = k + 1
        v_hat = v / (1 - gamma_v**k)
        s_hat = s / (1 - gamma_s**k)
        return x - alpha * v_hat / (np.sqrt(s_hat) + epsilon)

method = Adam(problem, alpha=0.2, gamma_v=0.9, gamma_s=0.999, epsilon=1e-8)

solution = solve(method)

# ---

class HyperGradientDescent:
    def __init__(self, problem, alpha_0, mu):
        self.mu = mu
        self.alpha = alpha_0
        self.g_prev = np.zeros(problem.nb_variables)

    def step(self, fn, x):
        alpha, mu, g, g_prev = self.alpha, self.mu, fn(x), self.g_prev
        alpha = alpha + mu * np.dot(g, g_prev)
        self.g_prev, self.alpha = g, alpha
        return x - alpha * g        

method = HyperGradientDescent(problem, alpha_0=0.2, mu=0.2)

solution = solve(method)

# ---

class HyperNesterovMomentum:
    def __init__(self, problem, alpha_0, mu, beta):
        self.mu = mu
        self.beta = beta
        self.v = np.zeros(problem.nb_variables)
        self.alpha = alpha_0
        self.g_prev = np.zeros(problem.nb_variables)

    def step(self, fn, x):
        alpha, beta, mu = self.alpha, self.beta, self.mu
        v, g, g_prev = self.v, fn(x), self.g_prev
        alpha = alpha - mu * np.dot(g, -g_prev - beta * v)
        v[:] = beta * v + g
        self.g_prev, self.alpha = g, alpha
        return x - alpha * (g + beta * v)        

method = HyperNesterovMomentum(problem, alpha_0=0.2, mu=0.2, beta=0.95)

solution = solve(method)

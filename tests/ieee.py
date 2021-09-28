import eqlib as eq
import numpy as np
import hyperjet as hj

def ieee_g4():
    variables = [
        eq.Variable(value=100, bounds=(78, 102), name='x1'),
        eq.Variable(value= 38, bounds=(33,  45), name='x2'),
        eq.Variable(value= 34, bounds=(27,  45), name='x3'),
        eq.Variable(value= 35, bounds=(27,  45), name='x4'),
        eq.Variable(value= 28, bounds=(27,  45), name='x5'),
    ]

    equations = [eq.Equation(name=f'g{i+1}') for i in range(6)]

    @eq.hj_objective_1(variables)
    def eval_f(x, info):
        return 5.3578547 * x[2]**2 + 0.8356891 * x[0] * x[4] + 37.293239 * x[0] - 40792.141

    @eq.hj_constraint_1(equations, variables)
    def eval_g(x, info):
        return [
            +85.334407 + 0.0056858 * x[1] * x[4] + 0.0006262 * x[0] * x[3] - 0.0022053 * x[2] * x[4] - 92,
            -85.334407 - 0.0056858 * x[1] * x[4] - 0.0006262 * x[0] * x[3] + 0.0022053 * x[2] * x[4],
            +80.512490 + 0.0071317 * x[1] * x[4] + 0.0029955 * x[0] * x[1] + 0.0021813 * x[2]**2 - 110,
            -80.512490 - 0.0071317 * x[1] * x[4] - 0.0029955 * x[0] * x[1] - 0.0021813 * x[2]**2 + 90,
            +9.3009610 + 0.0047026 * x[2] * x[4] + 0.0012547 * x[0] * x[2] + 0.0019085 * x[2] * x[3] - 25,
            -9.3009610 - 0.0047026 * x[2] * x[4] - 0.0012547 * x[0] * x[2] - 0.0019085 * x[2] * x[3] + 20,
        ]

    problem = eq.Problem([eval_f], [eval_g])

    return problem

def ieee_g9():
    variables = [
        eq.Variable(value=0.5, bounds=(-10, 10), name='x1'),
        eq.Variable(value=0.5, bounds=(-10, 10), name='x2'),
        eq.Variable(value=0.5, bounds=(-10, 10), name='x3'),
        eq.Variable(value=0.5, bounds=(-10, 10), name='x4'),
        eq.Variable(value=0.5, bounds=(-10, 10), name='x5'),
        eq.Variable(value=0.5, bounds=(-10, 10), name='x6'),
        eq.Variable(value=0.5, bounds=(-10, 10), name='x7'),
    ]

    equations = [eq.Equation(name=f'g{i+1}') for i in range(4)]

    @eq.hj_objective_1(variables)
    def f(x, info):
        return (x[0] - 10)**2 + 5 * (x[1] - 12)**2 + x[2]**4 + 3 * (x[3] - 11)**2 + 10 * x[4]**6 + 7 * x[5]**2 + x[6]**4 - 4 * x[5] * x[6] - 10 * x[5] - 8 * x[6]

    @eq.hj_constraint_1(equations, variables)
    def g(x, info):
        return [
            -127 + 2 * x[0]**2 + 3 * x[1]**4 + x[2] + 4 * x[3]**2 + 5 * x[4],
            -282 + 7 * x[0] + 3 * x[1] + 10 * x[2]**2 + x[3] - x[4],
            -196 + 23 * x[0] + x[1]**2 + 6 * x[5]**2 - 8 * x[6],
            4 * x[0]**2 + x[1]**2 - 3 * x[0] * x[1] + 2 * x[2]**2 + 5 * x[5] - 11 * x[6],
        ]

    problem = eq.Problem([f], [g])

    return problem

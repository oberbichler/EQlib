from eqlib.eqlib_ext import *
import hyperjet as hj
import numpy as np


def hj_objective_1(variables):
    def factory(fn):
        def fn_(x, df, hm, info):
            x = np.array(hj.variables(x, order=1))
            
            result = fn(x, info)

            df[:] = hj.d(result)

            return hj.f(result)
        return objective(variables)(fn_)
    return factory

def hj_constraint_1(equations, variables):
    def factory(fn):
        def fn_(x, g, dg, hm, info):
            x = np.array(hj.variables(x, order=1))
            
            result = fn(x, info)

            g[:] = hj.f(result)
            dg[:,:] = hj.d(result)

            return hj.f(result)
        return constraint(equations, variables)(fn_)
    return factory

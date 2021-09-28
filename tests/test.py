import eqlib as eq
import numpy as np
import hyperjet as hj
from numpy.testing import assert_allclose

def ieee_g1():
    x = [eq.Variable(name=f'x{i+1}', value=value) for i, value in enumerate([0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])]
    g = [eq.Equation(name=f'g{i+1}') for i in range(9)]

    lower_bounds = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], float)
    upper_bounds = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 100, 100, 100, 1], float)

    for i, xi in enumerate(x):
      xi.lower_bound = -x[i].value + lower_bounds[i]
      xi.upper_bound = x[i].value - upper_bounds[i]

    @eq.objective(x)
    def eval_f(x, df, hm, info):
        x = np.array(hj.variables(x, order=1))
        
        p = 5 * x[:4].sum() - 5 * x[:4] @ x[:4] - x[4:].sum()

        df[:] = hj.d(p)
        # hm[:,:] = hj.dd(p)

        return hj.f(p)

    @eq.constraint(g, x)
    def eval_g(x, g, dg, hm, info):
        x = hj.variables(x, order=1)

        gs = [
            2 * x[0] + 2 * x[1] + x[9] + x[10] - 10,
            2 * x[0] + 2 * x[2] + x[9] + x[11] - 10,
            2 * x[1] + 2 * x[2] + x[10] + x[11] - 10,
            -8 * x[0] + x[9],
            -8 * x[1] + x[10],
            -8 * x[2] + x[11],
            -2 * x[3] - x[4] + x[9],
            -2 * x[5] - x[6] + x[10],
            -2 * x[7] - x[8] + x[11],
        ]

        g[:] = hj.f(gs)
        dg[:,:] = hj.d(gs)
        # hm[:,:] = 0.0   # not used

    problem = eq.Problem([eval_f], [eval_g])

    # problem.x_lower_bounds = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    # problem.x_upper_bounds = [1, 1, 1, 1, 1, 1, 1, 1, 1, 100, 100, 100, 1]

    return problem

problem = ieee_g1()

problem.eval(problem.x)

problem.nb_equations, problem.nb_variables

assert(problem.f == 0.5)

assert_allclose(problem.df, [0, 0, 0, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1])
print(problem.f)





def compute_searchdirection(problem, solution, zeta):
    problem.x = solution
    
    problem.eval(solution)

    dg = problem.dg.toarray()

    weighted_grad_gs = [-1 / problem.g[i] * dg[i] for i in range(problem.nb_equations)]

    glb = problem.x_lower_bounds
    gub = problem.x_upper_bounds
    dglb = np.eye(13)
    dgub = np.eye(13)
    weighted_grad_glbs = [-1 / glb[i] * dglb[i] for i in range(glb.shape[0])]
    weighted_grad_gubs = [-1 / gub[i] * dgub[i] for i in range(gub.shape[0])]

    gradient_f = problem.df
    gradient_g = np.sum(weighted_grad_gs, axis=0) + np.sum(weighted_grad_glbs, axis=0) + np.sum(weighted_grad_gubs, axis=0)

    s = -gradient_f / np.linalg.norm(gradient_f) - zeta * gradient_g / np.linalg.norm(gradient_g)

    return s

def gradient_descent_akin(problem, n_iter, step_size, zeta):
    # set initial point
    initial_point = problem.x
    
    # track all solutions
    solutions = list()
    
    # generate and initial point
    solution = initial_point

    nb_variables = problem.nb_variables
    
    # list of changes mode to each variable
    change = np.zeros(nb_variables)
    
    # run gradient descent akin method
    for it in range(n_iter):
        # calculate the chen's direction
        s = compute_searchdirection(problem, solution, zeta)

        # normalize the search direction
        s = s / np.linalg.norm(s)

        # build a solution
        new_solution = solution + step_size * s

        # save the last solution
        last_solution = solution

        # evaluate candidate point 
        solution = new_solution
        solutions.append(solution)

        problem.eval(last_solution)
        solution_eval = problem.f

        problem.eval(solution)
        max_glb = np.amax(problem.x_lower_bounds)
        max_gub = np.amax(problem.x_upper_bounds)
        solution_feas = np.amax([*problem.g, max_glb, max_gub])

        if solution_feas >= 0:
            break
            
        # report progress
        print(solution_eval)
        # print('>%d f(%s) = %.5f' % (it, solution, solution_eval))

    print('> n = %d ' %(it))
    print('> solution = %s' %(last_solution))
    print('> f(x_n-1) = %.5f,'% (solution_eval))
    print('> g(x_n) = %.5f' % (solution_feas))
    return last_solution

def gradient_descent_akin_nesterov(problem, n_iter, step_size, zeta, momentum):
    # set initial point
    initial_point = get_initialization(problem)    
    # track all solutions
    solutions = list()
    # initialization
    solution = initial_point
    # list of changes made to each variable
    change = [0.0 for _ in range(len(initial_point))]
    # run the gradient descent akin method with Nesterov momemtum
    for it in range(n_iter):
        # calculate the projected solution
        projected = [solution[i] + momentum * change[i] for i in range(solution.shape[0])]
        # calculate the search direction s for the projection
        s = compute_searchdirection(problem, asarray(projected), zeta)
        # normalize the search direction
        s = s/np.linalg.norm(s)
        # build a solution one variable at a time
        new_solution = list()
        for i in range(solution.shape[0]):
            # calculate the change
            change[i] = (momentum * change[i]) + step_size * s[i]
            # calculate the new position in this variable
            value = solution[i] + change[i]
            # store this variable
            new_solution.append(value)
        # save the last solution
        last_solution = solution            
        # store the new solution
        solution = asarray(new_solution)
        solutions.append(solution)
        # evaluate candidate point
        solution_eval = problem.Objective(last_solution)
        solution_feas = problem.Constraint(solution)
        if solution_feas >= 0:
            #print('> n = %d ' %(it))
            break            
        # report progress
        #print('>%d f(%s) = %.5f' % (it, solution, solution_eval))
    print('> n = %d ' %(it))    
    print('> solution = %s' %(last_solution))
    print('> f(x_n-1) = %.5f,'% (solution_eval))
    print('> g(x_n) = %.5f' % (solution_feas))
    return solutions


# define the total iterations
n_iter = 4000
# define the step size
step_size = 0.002
# define zeta
zeta = 0.99
# define momentum
momentum = 0.95 

gradient_descent_akin(problem, n_iter, step_size, zeta)

print('done!')
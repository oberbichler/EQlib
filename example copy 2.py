import bin.EQlib  as eq
import numpy as np

import os
print(os.getpid())

node_1 = eq.Node(x= 0, y= 0, z=0)
node_2 = eq.Node(x= 5, y= 5, z=0)
node_3 = eq.Node(x= 8, y= 5, z=0)
node_4 = eq.Node(x=15, y=-2, z=0)
node_4 = eq.Node(x=14, y= 0, z=0) #

node_5 = eq.Node(x= 2, y= 7, z=0)
node_6 = eq.Node(x= 5, y= 7, z=0)
node_7 = eq.Node(x= 8, y= 7, z=0)
node_8 = eq.Node(x=11, y= 7, z=0)

node_1.x.is_fixed = True
node_1.y.is_fixed = True
node_1.z.is_fixed = True

node_2.x.is_fixed = True
node_2.y.is_fixed = True
node_2.z.is_fixed = True

node_3.x.is_fixed = True
# node_3.y.is_fixed = True #
node_3.z.is_fixed = True

node_4.x.is_fixed = True
node_4.y.is_fixed = True
node_4.z.is_fixed = True

node_5.x.is_fixed = True
node_5.y.is_fixed = True
node_5.z.is_fixed = True

node_6.x.is_fixed = True
node_6.y.is_fixed = True
node_6.z.is_fixed = True

node_7.x.is_fixed = True
node_7.y.is_fixed = True
node_7.z.is_fixed = True

node_8.x.is_fixed = True
node_8.y.is_fixed = True
node_8.z.is_fixed = True

truss_12 = eq.Parameter(value=0) # name='s_12', 
truss_23 = eq.Parameter(value=0) # name='s_23', 
truss_34 = eq.Parameter(value=0) # name='s_34', 
truss_56 = eq.Parameter(value=0) # name='s_56', 
truss_67 = eq.Parameter(value=0) # name='s_67', 
truss_78 = eq.Parameter(value=0) # name='s_78', 
truss_26 = eq.Parameter(value=0) # name='s_26', 
truss_37 = eq.Parameter(value=0) # name='s_37', 

equilibrium_2 = eq.NodalEquilibrium(
    node=node_2,
    connections=[
        (truss_12, node_1),
        (truss_23, node_3),
        (truss_26, node_6),
    ],
    loads=[
    ],
)

equilibrium_3 = eq.NodalEquilibrium(
    node=node_3,
    connections=[
        (truss_23, node_2),
        (truss_34, node_4),
        (truss_37, node_7),
    ],
    loads=[
    ],
)

equilibrium_6 = eq.NodalEquilibrium(
    node=node_6,
    connections=[
        (truss_26, node_2),
        (truss_56, node_5),
        (truss_67, node_7),
    ],
    loads=[
        (10, [0, 1, 0]),
    ],
)

equilibrium_7 = eq.NodalEquilibrium(
    node=node_7,
    connections=[
        (truss_37, node_3),
        (truss_67, node_6),
        (truss_78, node_8),
    ],
    loads=[
        (10, [0, 1, 0]),
    ],
)

elements = [
    equilibrium_2,
    equilibrium_3,
    equilibrium_6,
    equilibrium_7,
]

# 


eq.Log.info_level = 10
system = eq.SymmetricSystem(elements)#, linear_solver={'type': 'sparse_lu'})



# from scipy.optimize import differential_evolution
# from scipy.optimize import minimize

# iterations = 0
# iterations_g = 0

# def f(x):
#     global iterations
#     # x = [10*np.sqrt(2), 10, 10*np.sqrt(2), 10, 10, 0, 0, 0, 5.2]
#     system.x = x
#     system.compute()
#     iterations += 1
#     return system.f

# def g(x):
#     global iterations_g
#     system.x = x
#     system.compute()
#     iterations_g += 1
#     return system.g

# res = differential_evolution(
#     func=f,
#     bounds=[
#         (-50, 50),
#         (-50, 50),
#         (-50, 50),
#         (-50, 50),
#         (-50, 50),
#         (-50, 50),
#         (-50, 50),
#         (-50, 50),
#         (-50, 50), #
#     ],
#     atol=1e-6,
#     # callback=lambda xk, convergence: print(convergence),
#     # strategy='currenttobest1bin',
#     disp=True,
# )

# # res = minimize(fun=f, x0=system.x, jac=g, tol=1e-7, options={'disp': True})

# print(iterations)
# print(iterations_g)

# print(res.x)

# # # system.compute()

# # # print('f =', system.f)
# # # print('g =', system.g)

for element in elements:
    print(element.compute())

# 


eq.Log.info_level = 10

solver = eq.Worhp(system)
solver.minimize()

print(node_3.act_location)
print(system.x)

print('Done!')
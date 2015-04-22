import numpy as np

indep_vec = np.linspace(0,1,101)
# indep_vec is a discrete representation 
# of the closed interval [0,1]

# Build a LINEAR Function
f = lambda x: x
#  We have created a map
#  f: REALS to REALS by
#  f(x) = x

#  Build a Quadratic Function
g = lambda x: x**2

# Discrete Functions
# f: [0,1] to REALS
linear_func = f(indep_vec)
# g: [0,1] to REALS
quad_func = g(indep_vec)

print linear_func

print quad_func
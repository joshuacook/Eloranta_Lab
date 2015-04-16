# Finite Difference Research Paper

## Prelimiaries
Consdider the derivative:

f'(x)=lim(h to 0) (f(x+h)-f(x))/h

### Issues with computational representations

1.  Can't have a function f: REALS to REALS.  Can't have an continuous function. Computers require discrete representation.
1.  What does it mean h goes to 0?

##  Discrete Representation of Functions
We will use `numpy' and it's build-in array functionality to represent functions as vectors.

~~~

import numpy as np
indep_vec = np.linespace(0,1,100)
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
guad_func = g(indep_vec)

~~~

##  Next
1.  Build arrays to represent functions
1.  Calculate first and second derivatives of functions using finite differences
1.  Plot first and second derivatives
1.  Solve Time-Dependent Schrodinger using finite differencing techniques
1.  Crank-Nicolson Real-Time
1.  Take a guassian and momentum to get a traveling wave
1.  Add in harmonic potential
1.  Try to get tunneling

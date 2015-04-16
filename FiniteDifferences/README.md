# Finite Difference Research Paper

## Prelimiaries
Consdider the derivative:

    f'(x)=lim(h to 0) (f(x+h)-f(x))/h

### Issues with computational representations

1.  Can't have a function f: REALS to REALS.  Can't have an continuous function. Computers require discrete representation.
1.  What does it mean h goes to 0?

##  Discrete Representation of Functions
We will use `numpy` and it's build-in array functionality to represent functions as vectors.

## Build arrays to represent functions

~~~

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

~~~

## Calculate first and second derivatives of functions using finite differences

###  Approximate First Derivative using Finite Differences

~~~

indep_vec = [0,0.1,0.2,...,1]
indep_vec = [0,0.1,0.2,...,1]
quad_vec = [0,0.01,0.02,...,2]
cub_vec = [0,0.001,0.008,...,3]
const_vec = [1,1,1,...,1]

~~~

WE SEEK:

~~~

lin_prm_vec = [1,1,1,...,1]
qud_prm_vec = [0,0.2,0.4,...,2]
cub_prm_vec = [0,0.03,0.12,...,3]
con_prm_vec = [0,0,0,...,0]

~~~

###  Finite Differences

Consider: `f'(x) APPROX [f(x+h)-f(x)]/h`

Also: `f(x) APPROX [f(x)-f(x-h)]/h`

BEST: `f(x) APPROX [f(x+h)-f(x-h)]/2h`

###  Linear First Difference

look at ith element of `lin_prm_vec` with `h+0.1`.
`lin_prm_vec[i] = [linear_vec[i+1]]-linear_vec[i-1]/2*h`

e.g. `lin_prm_vec[1] = [linear_vec[2]-linear_vec[0]]/2*h`
	 `lin_prm_vec[2] = [linear_vec[3]-linear_vec[1]]/2*h`

###  Toward First Difference Operator

Let `D = d/dx`, Then `Df = f'`.
Here `D` is an operator, namely the first derivative.  Seek `FD` a matrix operator that approximates the first derviative of a vector using first finite difference,
i.e. `FD*linear_vec = lin_prm_vec`

~~~
         | -2 2 0 0 0 0 0 0 0 0 0 | |0.0| = |1|
         | -1 0 1 0 0 0 0 0 0 0 0 | |0.1| = |1| 
         | 0 -1 0 1 0 0 0 0 0 0 0 | |0.2| = |1|
         | 0 0 -1 0 1 0 0 0 0 0 0 | |0.3| = |1|
         | 0 0 0 -1 0 1 0 0 0 0 0 | |0.4| = |1|
1/(2*h)  | 0 0 0 0 -1 0 1 0 0 0 0 | |0.5| = |1|
         | 0 0 0 0 0 -1 0 1 0 0 0 | |0.6| = |1|
         | 0 0 0 0 0 0 -1 0 1 0 0 | |0.7| = |1|
         | 0 0 0 0 0 0 0 -1 0 1 0 | |0.8| = |1|
         | 0 0 0 0 0 0 0 0 -1 0 1 | |0.9| = |1|
         | 0 0 0 0 0 0 0 0 0 -2 2 | |1.0| = |1|
~~~

###  Toward Second Difference Operator

Let `D2 = d^2/dx^2`

Note to Josh, Add Second Difference Operator


1.  
1.  
1.  Plot first and second derivatives
1.  Solve Time-Dependent Schrodinger using finite differencing techniques
1.  Crank-Nicolson Real-Time
1.  Take a guassian and momentum to get a traveling wave
1.  Add in harmonic potential
1.  Try to get tunneling 

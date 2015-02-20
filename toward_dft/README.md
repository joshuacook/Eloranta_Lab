README.md

# Toward DFT

## Basic Power Method Iteration
An eigenvector algorithm. 

1. Start with a random vector
2. At each iteration, multiply by the matrix and normalize
3. The iteration converges on the largest eigenvector

### in IPython
    In [1]:
    import numpy.linalg, numpy.random, numpy as np, math
    from random import random as rand
    from numpy.linalg import eig
     
    In [2]:
    B = np.array([[2,-12],[1,-5]])

    In [3]:
    eig(B)[1]
    array([[ 0.9701425 ,  0.9486833 ],
           [ 0.24253563,  0.31622777]])

    In [4]:
    eig(B)[0]
    array([-1., -2.])

Generalize to larger matrices.


Algorithm in C converges to 0.948683, 0.316228 in 2 steps. 

## Finding the eigenvalue
To find the eigenvalue we need an algorithm for solving the matrix equation Ax=b. We have the matrix

## Finding the next largest eigenvector. 
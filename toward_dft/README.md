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
    print numpy.transpose(eig(B)[1])
    [[ 0.9701425   0.24253563]
     [ 0.9486833   0.31622777]]

## Next
Finding the eigenvalue.

Finding the next largest eigenvector. 
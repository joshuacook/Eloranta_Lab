import numpy as np
from time import clock
from scipy.linalg import eigh as largest_eigh
from scipy.sparse.linalg.eigen.arpack import eigsh as largest_eigsh

np.set_printoptions(suppress=True)
np.random.seed(0)
N=500
k=10
X = np.random.random((N,N)) - 0.5
X = np.dot(X, X.T) #create a symmetric matrix

# Benchmark the dense routine
start = clock()
evals_large, evecs_large = largest_eigh(X, eigvals=(N-k,N-1))
elapsed = (clock() - start)
v

# Benchmark the sparse routine
start = clock()
evals_large_sparse, evecs_large_sparse = largest_eigsh(X, k, which='LM')
elapsed = (clock() - start)
print "eigsh elapsed time: ", elapsed
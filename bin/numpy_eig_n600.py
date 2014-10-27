import numpy as np
from time import clock
import math
import numpy.linalg as la
from scipy.linalg import eigh as largest_eigh
from scipy.sparse.linalg.eigen.arpack import eigsh as largest_eigsh
import matplotlib.pyplot as plt

'''Scipy Dense v Sparse Comparison'''


np.set_printoptions(suppress=True)
np.random.seed(0)

N=155
k=1
eps = 10E-6
times = np.zeros((N,4))
H = np.random.random((N,N)) - 0.5
H = np.dot(H, H.T) #create a symmetric matrix

for i in range(N-50):
	print i
	i = i+15
	n = i
	
	conv = 1

	times[i-2][0] = i

	# Benchmark the dense routine
	start = clock()
	evals_large, evecs_large = largest_eigh(H[0:n,0:n], eigvals=(n-k,n-1))
	elapsed = (clock() - start)
	times[i-2][1] = elapsed
# 	print "eigh elapsed time: ", elapsed

	# Benchmark the sparse routine
	start = clock()
	evals_large_sparse, evecs_large_sparse = largest_eigsh(H[0:n,0:n], k, which='LM')
	elapsed = (clock() - start)
	times[i-2][2] = elapsed
# 	print "eigsh elapsed time: ", elapsed

	start = clock()

	phi0 = np.random.rand(n)
	# print la.eig(H)[1].T
	CayleyN = (np.identity(n)-0.5*H[0:n,0:n])
	CayleyP = (np.identity(n)+0.5*H[0:n,0:n])

	while(conv > eps):
		phi1 = la.solve(CayleyP,CayleyN.dot(phi0))
		mu = math.sqrt(phi1.dot(phi1))
		phi1 = phi1/mu  
		conv = math.sqrt((np.abs(phi1)-np.abs(phi0)).dot(np.abs(phi1)-np.abs(phi0)))
		phi0 = phi1

	end = clock()

	times[i-2][3] = end-start
	
p1, = plt.plot(times[:k-50,0],times[:k-50,1])
p2, = plt.plot(times[:k-50,0],times[:k-50,2])
p3, = plt.plot(times[:k-50,0],times[:k-50,3])
plt.legend([p1, p2,p3], ["scipy.linalg.eigh", "scipy.sparse.linalg.eigsh", "Cayley Method"])
plt.show()
np.savetxt("numpy_n600.csv",times,fmt='%.4e')
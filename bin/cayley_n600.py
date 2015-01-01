import numpy as np
import numpy.linalg as la
import math, time
import matplotlib.pyplot as plt
from sys import argv

'''This script runs a Cayley Expansion and Power Iteration method to find the first eigenvector of a random matrix.'''

k = int(argv[1])
print k 

eps = 10E-6
times = np.zeros((k,2))
H = np.random.rand(k+1,k+1)
H = H.T.dot(H)

for i in range(k-50):
	print i
	i = i+2
	n = i
	err = 1


	start = time.clock()

	phi0 = np.random.rand(n)
	# print la.eig(H)[1].T
	CayleyN = (np.identity(n)-0.5*H[0:n,0:n])
	CayleyP = (np.identity(n)+0.5*H[0:n,0:n])

	while(err > eps):
		phi1 = la.solve(CayleyP,CayleyN.dot(phi0))
		mu = math.sqrt(phi1.dot(phi1))
		phi1 = phi1/mu  
		err = math.sqrt(2)*math.sqrt(abs(phi1.dot(int_H.dot(int_H)).dot(phi1)- (phi1.dot(int_H).dot(phi1))**2))
		phi0 = phi1

	end = time.clock()

	times[i-2][0] = i
	times[i-2][1] = end-start

plt.plot(times[:k-50,0],times[:k-50,1])

plt.show()

np.savetxt("cayle_n500.csv",times,fmt='%.4e')

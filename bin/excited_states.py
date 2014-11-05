import numpy as np
import numpy.linalg as la
import math, time
import matplotlib.pyplot as plt
from sys import argv
import datetime

'''This script runs a Cayley Expansion and Power Iteration method to find the first eigenvector of a random matrix.'''

k = int(argv[1])
print k 

eps = 10E-6
times = np.array([[0.,0.]]) 
temp_times = times
H = np.random.rand(k+1,k+1)
H = H.T.dot(H)
file = datetime.datetime.now().strftime("%Y%m%d%H%M%S")

for i in range(k-50):
	print i
	i = i+2
	n = i
	err = 1
	conv = 1
	int_H = H[0:n,0:n]

	start = time.clock()

	phi0 = np.random.rand(n)
	# print la.eig(H)[1].T
	CayleyN = (np.identity(n)-0.5*int_H)
	CayleyP = (np.identity(n)+0.5*int_H)

	while(conv > eps):
		phi1 = la.solve(CayleyP,CayleyN.dot(phi0))
		mu = math.sqrt(phi1.dot(phi1))
		phi1 = phi1/mu  
		conv = math.sqrt((np.abs(phi1)-np.abs(phi0)).dot(np.abs(phi1)-np.abs(phi0)))
		# err = math.sqrt(2)*math.sqrt(abs(phi1.dot(int_H.dot(int_H)).dot(phi1)- (phi1.dot(int_H.dot(int_H)).dot(phi1))**2))
		# print err
		phi0 = phi1

	end = time.clock()
	delta_t = end-start
	temp_times[0][0] = i
	temp_times[0][1] = delta_t
        times = np.concatenate((times,temp_times),axis=0)
	np.savetxt(file,times,fmt='%.4e')

plt.plot(times[:k-50,0],times[:k-50,1])

plt.show()


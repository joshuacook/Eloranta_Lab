#! /usr/bin/python

import numpy as np
import numpy.linalg as la
import math, time
import matplotlib.pyplot as plt
from scipy.linalg import eigh as largest_eigh
import datetime

k = 10000
eps = 10E-6
times = np.array([[0.,0.]]) 
temp_times = times
H = np.random.rand(k+200,k+200)
H = H.T.dot(H)
file = datetime.datetime.now().strftime("%Y%m%d%H%M%S")+".csv"

for i in range(k):

	i = i+2
	n = i
	conv = 1

	start = time.clock()
	evals_large, evecs_large = largest_eigh(H[0:n,0:n])
	elapsed = (time.clock() - start)
	times[i-2][1] = elapsed
	

	temp_times[0][0] = i
	temp_times[0][1] = elapsed
	times = np.concatenate((times,temp_times),axis=0)
	np.savetxt(file,times,fmt='%.4e')
	print i

plt.plot(times[:,0],times[:,1])

plt.show()

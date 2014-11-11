import numpy as np
#! /usr/bin/python

import numpy.linalg as la
import math, time
import matplotlib.pyplot as plt
from sys import argv
import scipy.sparse.linalg as spsla
import scipy.linalg as spla
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg.dsolve import linsolve as linsolve
import datetime

k = 1000

eps = 10E-6
times = np.array([[0.,0.]]) 
temp_times = times
H = np.random.rand(k+200,k+200)
H = H.T.dot(H)
file = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-numpy")+".csv"

for i in range(k):
	print i
	i = i+3
	n = i

	int_H = H[0:n,0:n]
	
	# Numpy Solver
	err = 1
	start = time.clock()
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*int_H)
	CayleyP = (np.identity(n)+0.5*int_H)
	while(err > eps):
			phi1 = la.solve(CayleyP,CayleyN.dot(phi0))
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(int_H.dot(int_H)).dot(phi1)- (phi1.dot(int_H).dot(phi1))**2))
			phi0 = phi1
	elapsed_numpysolver = (time.clock() - start)

	

	temp_times[0][0] = i
	temp_times[0][1] = elapsed_numpysolver
	times = np.concatenate((times,temp_times),axis=0)
	np.savetxt(file,times,fmt='%.4e')



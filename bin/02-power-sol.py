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


k = 600

eps = 10E-6
times = np.array([[0.,0.,0.,0.,0.,0.,0.,0.]]) 
temp_times = times
H = np.random.rand(k+200,k+200)
H = H.T.dot(H)
file = datetime.datetime.now().strftime("%Y%m%d%H%M%S")+".csv"

for i in range(k):
	print "%d " % i
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
	
	

	# PreInvert
	err = 1
	start = time.clock()
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*int_H)
	CayleyP = (np.identity(n)+0.5*int_H)
	CayleyP_inv = la.inv(CayleyP)
	while(err > eps):
			phi1 = CayleyP_inv.dot(CayleyN.dot(phi0))
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(int_H.dot(int_H)).dot(phi1)- (phi1.dot(int_H).dot(phi1))**2))
			phi0 = phi1
	
	elapsed_preinvert = (time.clock() - start)
	
	# Scipy Solver
	err = 1
	start = time.clock()
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*int_H)
	CayleyP = (np.identity(n)+0.5*int_H)
	while(err > eps):
			phi1 = spla.solve(CayleyP,CayleyN.dot(phi0))
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(int_H.dot(int_H)).dot(phi1)- (phi1.dot(int_H).dot(phi1))**2))
			phi0 = phi1
	elapsed_scipysolver = (time.clock() - start)
	
	
	# Conjugate Gradient Squared
	err = 1
	start = time.clock()
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*int_H)
	CayleyP = (np.identity(n)+0.5*int_H)
	while(err > eps):
			phi1 = spsla.cgs(CayleyP,CayleyN.dot(phi0))
			phi1=phi1[0]
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(int_H.dot(int_H)).dot(phi1)- (phi1.dot(int_H).dot(phi1))**2))
			phi0 = phi1
	elapsed_cgs = (time.clock() - start)
	
	# Cholesky Factorization
	err = 1
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*int_H)
	CayleyP = (np.identity(n)+0.5*int_H)
	cho = spla.cho_factor(CayleyP)
	while(err > eps):
			phi1 = spla.cho_solve(cho, CayleyN.dot(phi0))
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(int_H.dot(int_H)).dot(phi1)- (phi1.dot(int_H).dot(phi1))**2))
			phi0 = phi1	
	elapsed_cho = (time.clock() - start)
	
	start = time.clock()
	eigsh(int_H,1)
	elapsed_scipyeig = (time.clock() - start)
	

	temp_times[0][0] = i
	temp_times[0][1] = elapsed_numpysolver
	temp_times[0][2] = elapsed_cgs
	temp_times[0][3] = elapsed_cho
	temp_times[0][4] = elapsed_preinvert
	temp_times[0][5] = elapsed_scipyeig
	times = np.concatenate((times,temp_times),axis=0)
	np.savetxt(file,times,fmt='%.4e')



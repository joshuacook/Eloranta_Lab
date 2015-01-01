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

k = 25

eps = 10E-6
times = np.array([[0.,0.,0.,0.,0.,0.,0.,0.,0.]]) 
temp_times = times
file = datetime.datetime.now().strftime("%Y%m%d%H%M%S")+".csv"

for i in range(k):
	print i
	i = i+3
	n = 2**i
	H = np.random.rand(2**i,2**i)
	H = H.T.dot(H)
	
	# Numpy Solver
	err = 1
	start = time.clock()
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*H)
	CayleyP = (np.identity(n)+0.5*H)
	while(err > eps):
			phi1 = la.solve(CayleyP,CayleyN.dot(phi0))
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
			phi0 = phi1
	elapsed_numpysolver = (time.clock() - start)
	
	

	# PreInvert
	err = 1
	start = time.clock()
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*H)
	CayleyP = (np.identity(n)+0.5*H)
	CayleyP_inv = la.inv(CayleyP)
	while(err > eps):
			phi1 = CayleyP_inv.dot(CayleyN.dot(phi0))
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
			phi0 = phi1
	
	elapsed_preinvert = (time.clock() - start)
	
	# Scipy Solver
	err = 1
	start = time.clock()
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*H)
	CayleyP = (np.identity(n)+0.5*H)
	while(err > eps):
			phi1 = spla.solve(CayleyP,CayleyN.dot(phi0))
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
			phi0 = phi1
	elapsed_scipysolver = (time.clock() - start)
	
	# Scipy Solver
	err = 1
	start = time.clock()
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*H)
	CayleyP = (np.identity(n)+0.5*H)
	while(err > eps):
			phi1 = spsla.spsolve(CayleyP,CayleyN.dot(phi0))
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
			phi0 = phi1
	elapsed_scipysparsesolver = (time.clock() - start)
	
	# Scipy DSolve 
	err = 1
	start = time.clock()
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*H)
	CayleyP = (np.identity(n)+0.5*H)
	while(err > eps):
			phi1 = linsolve.spsolve(CayleyP,CayleyN.dot(phi0))
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
			phi0 = phi1
	elapsed_linsolve = (time.clock() - start)
	
	# Conjugate Gradient
	err = 1
	start = time.clock()
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*H)
	CayleyP = (np.identity(n)+0.5*H)
	while(err > eps):
			phi1 = spsla.cgs(CayleyP,CayleyN.dot(phi0))
			phi1=phi1[0]
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
			phi0 = phi1
	elapsed_cgs = (time.clock() - start)
	
	# Cholesky Factorization
	err = 1
	phi0 = np.random.rand(n)
	CayleyN = (np.identity(n)-0.5*H)
	CayleyP = (np.identity(n)+0.5*H)
	cho = spla.cho_factor(CayleyP)
	while(err > eps):
			phi1 = spla.cho_solve(cho, CayleyN.dot(phi0))
			mu = math.sqrt(phi1.dot(phi1))
			phi1 = phi1/mu  
			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
			phi0 = phi1	
	elapsed_cho = (time.clock() - start)
	
	start = time.clock()
	eigsh(H,1)
	elapsed_scipyeig = (time.clock() - start)
	

	temp_times[0][0] = i
	temp_times[0][1] = elapsed_numpysolver
	temp_times[0][2] = elapsed_scipysolver
	temp_times[0][3] = elapsed_cgs
	temp_times[0][4] = elapsed_cho
	temp_times[0][5] = elapsed_linsolve
	temp_times[0][6] = elapsed_scipysparsesolver
	temp_times[0][7] = elapsed_preinvert
	temp_times[0][8] = elapsed_scipyeig
	times = np.concatenate((times,temp_times),axis=0)
	np.savetxt(file,times,fmt='%.4e')



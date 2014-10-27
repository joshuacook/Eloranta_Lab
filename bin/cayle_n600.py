import numpy as np
import numpy.linalg as la
import math, time
import matplotlib.pyplot as plt

k = 600
eps = 10E-6
times = np.zeros((k,2))
H = np.random.rand(k+1,k+1)
H = H.T.dot(H)

for i in range(k-50):
	print i
	i = i+2
	n = i
	conv = 1


	start = time.clock()

	phi0 = np.random.rand(n)
	# print la.eig(H)[1].T
	CayleN = (np.identity(n)-0.5*H[0:n,0:n])
	CayleP = (np.identity(n)+0.5*H[0:n,0:n])

	while(conv > eps):
		phi1 = la.solve(CayleP,CayleN.dot(phi0))
		mu = math.sqrt(phi1.dot(phi1))
		phi1 = phi1/mu  
		conv = math.sqrt((np.abs(phi1)-np.abs(phi0)).dot(np.abs(phi1)-np.abs(phi0)))
		phi0 = phi1

	end = time.clock()

	times[i-2][0] = i
	times[i-2][1] = end-start

plt.plot(times[:,0],times[:,1])

plt.show()

np.savetxt("cayle_n500.csv",times,fmt='%.4e')
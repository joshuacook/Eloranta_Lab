

    import numpy as np
    import numpy.linalg as la
    import math, time
    import matplotlib.pyplot as plt
    %matplotlib inline

Algorithm

1. generate random matrix, $H$
1. take $H^TH$ to make positive definite
1. run a for loop for the desired number of tests
    1. start the clock
    1. generate a random vector
    1. do a Cayley transform on H based on $e^{-H\Delta \tau}=
\bigg(1+\frac{1}{2}H\Delta\tau\bigg)^{-1}\Big(1-\frac{1}{2}H\Delta\tau\Big)$,
    to obtain $\Big(1+\frac{1}{2}H\Delta\tau\Big)$ and
$\Big(1-\frac{1}{2}H\Delta\tau\Big)$
    1. use an iterative power method to solve
    $\Big(1+\frac{1}{2}H\Delta\tau\Big)\mathbf{\phi}=\Big(1-\frac{1}{2}H\Delta\t
au\Big)\mathbf{\phi}$
    1. stop clock
    1. use stopping criteria checking difference between iterated vectors
1. plot and export data


    k = 150
    eps = 10E-6
    times = np.zeros((k,2))
    H = np.random.rand(k+1,k+1)
    H = H.T.dot(H)
    
    for i in range(k-50):
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




    [<matplotlib.lines.Line2D at 0x1131823d0>]




![png](Cayley%20Expansion%20and%20Iterative%20Solution_files/Cayley%20Expansion%20and%20Iterative%20Solution_2_1.png)



    k = 350
    eps = 10E-6
    times = np.zeros((k,2))
    H = np.random.rand(k+1,k+1)
    H = H.T.dot(H)
    
    for i in range(k-50):
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
    
        
    np.savetxt("cayle_n500.csv",times,fmt='%.4e')    
        
    plt.plot(times[:,0],times[:,1])




    [<matplotlib.lines.Line2D at 0x1101a7390>]




![png](Cayley%20Expansion%20and%20Iterative%20Solution_files/Cayley%20Expansion%20and%20Iterative%20Solution_3_1.png)



    k = 550
    eps = 10E-6
    times = np.zeros((k,2))
    H = np.random.rand(k+1,k+1)
    H = H.T.dot(H)
    
    for i in range(k-50):
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




    [<matplotlib.lines.Line2D at 0x1048a7e90>]




![png](Cayley%20Expansion%20and%20Iterative%20Solution_files/Cayley%20Expansion%20and%20Iterative%20Solution_4_1.png)



    k = 750
    eps = 10E-6
    times = np.zeros((k,2))
    H = np.random.rand(k+1,k+1)
    H = H.T.dot(H)
    
    for i in range(k-50):
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




    [<matplotlib.lines.Line2D at 0x108ef1e90>]




![png](Cayley%20Expansion%20and%20Iterative%20Solution_files/Cayley%20Expansion%20and%20Iterative%20Solution_5_1.png)



    

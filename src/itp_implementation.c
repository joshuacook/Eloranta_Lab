int itp_method_test(int n){
	// declare necessary variables
	int i, err;
	double H[n*n];
	double CayleyN[n*n];
	double CayleyP[n*n];
	double CayleyP_inv[n*n];
	double alpha = 1.0;
  double beta = 1.0;					 
	double mu;
  double mu_inv;
  int one = 1;
  int two = 2;
  char no_trans='N';	
  double phi0[n];
	double phi1[n];
	
	// run test on a matrix of size n
	printf("%d", n);

  // generate a random matrix of size n x n
	// H = np.random.rand(2**i,2**i)
	// multiply H by its transpose to make it symmetric and thus Hermitian
	// H = H.T.dot(H)

	// PreInvert
	err = 1

	// randomize phi
	// phi0 = np.random.rand(n)
	
	// take the Cayley form of H
	// CayleyN = (np.identity(n)-0.5*H)
	// CayleyP = (np.identity(n)+0.5*H)
	
	// invert Cayley P
	// CayleyP_inv = la.inv(CayleyP)

  // iterate 
	/* while(err > eps):
	phi1 = CayleyP_inv.dot(CayleyN.dot(phi0))
	mu = math.sqrt(phi1.dot(phi1))
	phi1 = phi1/mu  
	err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
	phi0 = phi1 */

}

/* Alternative methods for solving the matrix problem
	
	// Numpy Solver
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
	
	


	
	// Scipy Solver
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
	
	// Scipy Solver
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
	
	// Scipy DSolve 
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
	
	// Conjugate Gradient
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
	
	// Cholesky Factorization
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
	
*/


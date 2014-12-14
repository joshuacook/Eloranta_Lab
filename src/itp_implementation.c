#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "blas_fortran_double.h"
#include "lapack_fortran_double.h"
#include "blas_utilities.h"
#include "itp_method.h"

int itp_method_test(double * H, int n, int print_mat, int DEBUG);


int main(int argc, char* argv[]){

  // Performance Testing
	double time_elapsed_in_seconds;

	// Matrix parameters
	int n = 100;
	double A[n*n];

  // create matrix to be used for all
  random_matrix(A,n,n);


  //--------------------- ITP Method Test  
  
  // start timer
	clock_t start = clock();
  // call method
  // get first eigenvector
  itp_method_test(A,2,1,0);  

  // end timer
	clock_t end = clock();  
  // store time
	printf("Found eigenvector in %.2g seconds\n", time_elapsed_in_seconds = (end - start)/(double)CLOCKS_PER_SEC);
  
  //--------------------- ARPACK Method
  // start timer
  // call method
  // get first eigenvector
  // end timer
  // store time
  
  //--------------------- LAPACK Method
  // start timer
  // call method
  // get first eigenvector
  // end timer
  // store time
  
	return 0;


}




// 	
// 
//  // generate a random matrix of size n x n
// 	// H = np.random.rand(2**i,2**i)
// 	// multiply H by its transpose to make it symmetric and thus Hermitian
// 	// H = H.T.dot(H)
// 
// 	// PreInvert
// 	err = 1
// 
// 	// randomize phi
// 	// phi0 = np.random.rand(n)
// 	
// 	// take the Cayley form of H
// 	// CayleyN = (np.identity(n)-0.5*H)
// 	// CayleyP = (np.identity(n)+0.5*H)
// 	
// 	// invert Cayley P
// 	// CayleyP_inv = la.inv(CayleyP)
// 
//  // iterate 
// 	/* while(err > eps):
// 	phi1 = CayleyP_inv.dot(CayleyN.dot(phi0))
// 	mu = math.sqrt(phi1.dot(phi1))
// 	phi1 = phi1/mu  
// 	err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
// 	phi0 = phi1 */
// 
// }

// Alternative methods for solving the matrix problem
// 	
// 	// Numpy Solver
// 	err = 1
// 	start = time.clock()
// 	phi0 = np.random.rand(n)
// 	CayleyN = (np.identity(n)-0.5*H)
// 	CayleyP = (np.identity(n)+0.5*H)
// 	while(err > eps):
// 			phi1 = la.solve(CayleyP,CayleyN.dot(phi0))
// 			mu = math.sqrt(phi1.dot(phi1))
// 			phi1 = phi1/mu  
// 			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
// 			phi0 = phi1
// 	elapsed_numpysolver = (time.clock() - start)
// 	
// 	
// 
// 
// 	
// 	// Scipy Solver
// 	err = 1
// 	start = time.clock()
// 	phi0 = np.random.rand(n)
// 	CayleyN = (np.identity(n)-0.5*H)
// 	CayleyP = (np.identity(n)+0.5*H)
// 	while(err > eps):
// 			phi1 = spla.solve(CayleyP,CayleyN.dot(phi0))
// 			mu = math.sqrt(phi1.dot(phi1))
// 			phi1 = phi1/mu  
// 			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
// 			phi0 = phi1
// 	elapsed_scipysolver = (time.clock() - start)
// 	
// 	// Scipy Solver
// 	err = 1
// 	start = time.clock()
// 	phi0 = np.random.rand(n)
// 	CayleyN = (np.identity(n)-0.5*H)
// 	CayleyP = (np.identity(n)+0.5*H)
// 	while(err > eps):
// 			phi1 = spsla.spsolve(CayleyP,CayleyN.dot(phi0))
// 			mu = math.sqrt(phi1.dot(phi1))
// 			phi1 = phi1/mu  
// 			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
// 			phi0 = phi1
// 	elapsed_scipysparsesolver = (time.clock() - start)
// 	
// 	// Scipy DSolve 
// 	err = 1
// 	start = time.clock()
// 	phi0 = np.random.rand(n)
// 	CayleyN = (np.identity(n)-0.5*H)
// 	CayleyP = (np.identity(n)+0.5*H)
// 	while(err > eps):
// 			phi1 = linsolve.spsolve(CayleyP,CayleyN.dot(phi0))
// 			mu = math.sqrt(phi1.dot(phi1))
// 			phi1 = phi1/mu  
// 			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
// 			phi0 = phi1
// 	elapsed_linsolve = (time.clock() - start)
// 	
// 	// Conjugate Gradient
// 	err = 1
// 	start = time.clock()
// 	phi0 = np.random.rand(n)
// 	CayleyN = (np.identity(n)-0.5*H)
// 	CayleyP = (np.identity(n)+0.5*H)
// 	while(err > eps):
// 			phi1 = spsla.cgs(CayleyP,CayleyN.dot(phi0))
// 			phi1=phi1[0]
// 			mu = math.sqrt(phi1.dot(phi1))
// 			phi1 = phi1/mu  
// 			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
// 			phi0 = phi1
// 	elapsed_cgs = (time.clock() - start)
// 	
// 	// Cholesky Factorization
// 	err = 1
// 	phi0 = np.random.rand(n)
// 	CayleyN = (np.identity(n)-0.5*H)
// 	CayleyP = (np.identity(n)+0.5*H)
// 	cho = spla.cho_factor(CayleyP)
// 	while(err > eps):
// 			phi1 = spla.cho_solve(cho, CayleyN.dot(phi0))
// 			mu = math.sqrt(phi1.dot(phi1))
// 			phi1 = phi1/mu  
// 			err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
// 			phi0 = phi1	
// 	elapsed_cho = (time.clock() - start)

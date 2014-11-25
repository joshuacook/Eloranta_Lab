#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void print_matrix(double * A, int n, int m);
void create_matrix(double * A, int n, int m);

int main(int argc, char* argv[]){
	int n = 4;
	double A[n*n];

	create_matrix(A,n,n);

	print_matrix(A,n,n);

	return 0;

}

void create_matrix(double * A, int m, int n){
	time_t timer;
	int i,j;

	srand48(time(&timer));
	
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			A[i*m+j] = (i+1)*(j+1);
			printf("element %d, %d at A[%d] is %f\n", j,i,i*m+j,A[i*m+j]);
		}
	}
}

void print_matrix(double * A, int n, int m){
	int i,j;
	for (i=0; i < n*n; i++){
		printf("%f\n",A[i]);
}
	for (i = 0; i < n; i++){
		for (j = 0; j < m; j++){
			printf("element %d,%d: %f\n", i,j,(float)A[j*(n+i)]);
		}
	}
}








// int itp_method_test(int n){
// 	
// 	// declare necessary variables
// 	int i;
//  int err;
// 	double H[n*n];
// 	double CayleyN[n*n];
// 	double CayleyP[n*n];
// 	double CayleyP_inv[n*n];
// 	double alpha = 1.0;
//  double beta = 1.0;					 
// 	double mu;
//  double mu_inv;
//  int one = 1;
//  int two = 2;
//  char no_trans='N';	
//  double phi0[n];
// 	double phi1[n];*/
// 	
// 	// run test on a matrix of size n
// 	printf("Find an eigenvector for an %d by %d matrix", n,n);
// 	for(int i = 0; i < n*n ; i++ ){
// 		H[i]=rand();
// 	}
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

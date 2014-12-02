
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "blas_fortran_double.h"

#define NEWLINE printf("\n");

void print_matrix(double * A, int n, int m);
void print_vector(double * A, int n);
void random_matrix(double * A, int n, int m);
void random_vector(double *A, int n);
void identity_matrix(double * I, int n);
int itp_method_test(int n, int print_mat);


int main(int argc, char* argv[]){
	int n = 6;
	double A[n*n];

  itp_method_test(4,1);
	return 0;

}

void random_matrix(double * A, int m, int n){
	time_t timer;
	int i;

	srand48(time(&timer));
	
	for (i = 0; i < m*n; i++){
		A[i] = drand48();
	}
}

void identity_matrix(double * I, int n){
  int i;

  for(i = 0; i < n + 1; i++){
    I[i*(n+1)] = 1;
  } 
}

void random_vector(double * A, int n){
  time_t timer;
  int i;

  srand48(time(&timer));

  for (i = 0; i < n; i++){
    A[i] = drand48();
  } 
}

void print_vector(double * A, int n){
  int i;
  printf("<");
  for (i = 0; i < n-1; i++){
    printf("%f, ", (float)A[i]);
  }
  printf("%f", (float)A[i]);
  printf(">");
  NEWLINE;
  NEWLINE;
}

void print_matrix(double * A, int n, int m){
	int i,j;
	for (i = 0; i < m; i++){
		for (j = 0; j < n; j++){
			printf("%f ", (float)A[j*m+i]);
		}
		NEWLINE;	
	}
  NEWLINE;
}



int itp_method_test(int n,int print_mat){

	// declare necessary variables
 	int i;
  int err;
 	double H[n*n];
  double I[n*n];
 	double CayleyN[n*n];
 	double CayleyP[n*n];
 	double CayleyP_inv[n*n];
 	double d_zero = 0.0;
  double d_neghalf = -0.5;
  double d_poshalf = 0.5;
  double d_one = 1.0;					 
 	double mu;
  double mu_inv;
  int one = 1;
  int two = 2;
  char no_trans='N';
  char trans='T';	
  double phi0[n];
 	double phi1[n];



 	// run test on a matrix of size n
 	printf("Find an eigenvector for an %d by %d matrix\n", n,n);
  random_matrix(H,n,n);
  if(print_mat==1){
		print_matrix(H,n,n);
	}

  // multiply H by its transpose to make it symmetric and thus Hermitian
  // one*H*H'+zero*CayleyN
  dgemm_(&no_trans,&trans,&n,&n,&n,&d_one,H,&n,H,&n,&d_zero,CayleyN,&n);
  // one*H*H'+zero*CayleyP
  dgemm_(&no_trans,&trans,&n,&n,&n,&d_one,H,&n,H,&n,&d_zero,CayleyP,&n);
  print_matrix(CayleyN,n,n);
  print_matrix(CayleyP,n,n);

  // preInvert

  err = 1;

 	// randomize phi

  random_vector(phi0,n);
  print_vector(phi0,n);

 	// take the Cayley form of H
  identity_matrix(I,n);
  print_matrix(I,n,n);

  // CayleyN = (one*I*I-0.5*CayleyN)
  dgemm_(&no_trans,&no_trans,&n,&n,&n,&d_one,I,&n,I,&n,&d_neghalf,CayleyN,&n);
  print_matrix(CayleyN,n,n);
 	// CayleyP = (one*I*I+0.5*CayleyP)
  dgemm_(&no_trans,&no_trans,&n,&n,&n,&d_one,I,&n,I,&n,&d_poshalf,CayleyP,&n);
  print_matrix(CayleyP,n,n);

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

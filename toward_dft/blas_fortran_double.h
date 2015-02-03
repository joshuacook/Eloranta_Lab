// Variables
// alpha (double)	a scalar multiple
// beta (double)	a scalar multiple
// incX (int) 		the increment for the vector elements in X
// incY (int)		the increment for the vector elements in Y
// lda (int)		the leading dimension of a
// N (int) 			the number of elements in X and  Y
// Trans (char)		transpose specification 
//					'N' = none
//					'T' = transpose
//					'C' = conjugate transpose
// X (double) 		a vector of double-type elements
// Y (double)		a vector of double-type elements

//------------------------------------------------------ Level 1 Routines -- //
// BLAS Level 1 routines and functions perform vector-vector operations. 

// ?copy - Copies vector to another vector.
extern void dcopy_(int* N,double* X, int* incX, double* Y, int* incY);

// ?dot - Computes a vector-vector dot product
extern double ddot_(int* N, double* X, int* incX, double* Y, int* incY);

// ?scal - Computes the product of a vector by a scalar.
extern void dscal_(int* N,double* alpha,double* X, int* incX); 

// ?nrm2 - Computes the Euclidean norm of a vector.
extern double dnrm2_(int* N, double* X, int* incX);

//------------------------------------------------------ Level 2 Routines -- //
// BLAS Level 2 routines perform matrix-vector operations.

// ?gemv
// Computes a matrix-vector product using a general matrix
extern void dgemv_(char* Trans, int* M, int* N, double* alpha, 
					double* A, int* lda, double* X, int* incX, 
					double* beta, double* Y, int* incY );

//------------------------------------------------------ Level 3 Routines -- //
// BLAS Level 3 routines perform matrix-matrix operations.

// ?gemm
// Computes a scalar-matrix-matrix product and adds the result to a scalar-matrix product.
extern void dgemm_(char * transa, char* transb, int* m, int* n, int* k, double* alpha, double* A, int* lda, double* B, int* ldb, double* beta, double* C, int* ldc);

// LAPACK Routines
// Variables 
// m - the number of rows in a
// n - the number of columns in a
// A - a matrix
// lda - the leading dimension of A
// ipiv - an output array, the pivot variables
// work - a workspace array 
// lwork - the size of work, lwork geq n
// info - result (0, successful,-i,the ith variablehas an illegal parameter, i, matrix is singular and u_ii is zero

// LU decomoposition of a general matrix                                                                                                    
extern  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);                                                            
// generate inverse of a matrix given its LU decomposition                                                                                  
extern  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);                                          

void inverse(double* A, int N)                                                                                                              
{                                                                                                                                           
   int *IPIV;                                                                                                                              
   IPIV = (int *)malloc(sizeof(int)*(N+1));                                                                                                
   int LWORK = N*N;                                                                                                                        
   double *WORK;                                                                                                                           
   WORK = (double *)malloc(sizeof(double)*LWORK);                                                                                          
   int INFO;                                                                                                                               
                                                                                                                                            
   dgetrf_(&N,&N,A,&N,IPIV,&INFO);                                                                                                         
   dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);                                                                                                
                                                                                                                                           
   free(IPIV);                                                                                                                             
   free(WORK);                                                                                                                             
}                                  

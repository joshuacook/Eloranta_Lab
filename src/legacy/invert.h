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

int main(){

    double A [2*2] = {
        1,2,
        3,4
    };

    inverse(A, 2);

    printf("%f %f\n", A[0], A[1]);
    printf("%f %f\n", A[2], A[3]);

    return 0;
}     

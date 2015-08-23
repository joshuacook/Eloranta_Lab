#include "grid.h"
#include "private.h"

void zhegv_(int *, char *, char *, int *, double complex *, int *, double complex *, int *, double *, double complex *, int *, double *, int *);
void zheev_(char *, char *, int *, double complex *, int *, double *, double complex *, int *, double *, int *);

/* NOTE: These are not parallel versions !!! */

EXPORT int grid_generalized_hermitian_eigenvalue_problem(double *eigenvalue, double complex *hamiltonian, double complex *overlap, int states) {

  int itype, n, lda, ldb, lwork, info;
  char jobz, uplo;
  double complex *a, *b, *work;
  double *w, *rwork;
  
  itype = 1;
  jobz = 'V';
  uplo = 'U';
  n = states;
  lda = n;
  ldb = n;
  a = hamiltonian;
  b = overlap;
  w = eigenvalue;
  lwork = 2*n;
  work = (double complex *) malloc(lwork * sizeof(double complex));
  rwork = (double *) malloc((3*n - 2) * sizeof(double));
  
  if (!work || !rwork) return -999;
  
  zhegv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &info);
  
  return info;
}

EXPORT int grid_hermitian_eigenvalue_problem(double *eigenvalue, double complex *overlap, int states) {

  int n, lda, lwork, info;
  char jobz, uplo;
  double complex *a, *work;
  double *w, *rwork;
  
  jobz = 'V';
  uplo = 'U';
  n = states;
  lda = n;
  a = overlap;
  w = eigenvalue;
  lwork = 2*n;
  work = (double complex *) malloc(lwork * sizeof(double complex));
  rwork = (double *) malloc((3*n - 2) * sizeof(double));
  
  if (!work || !rwork) return -999;
  
  zheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info);
  
  return info;
}

/*
 * Solve tridiagonal matrix equation A x = b (Thomas algorithm).
 *
 * n = dimensionality (number of equations).
 * a = subdiagonal (indexed 1 ... n-1). 0 unused.
 * b = diagonal (overwritten).
 * c = supdiagonal (indexed 0 ... n-2). n-1 unused.
 * v = right hand side vector (overwritten).
 * x = solution on exit (double complex *).
 *
 * Note: a and x may be the same array.
 *
 */

EXPORT void grid_solve_tridiagonal_system(long n, double complex *a, double complex *b, double complex *c, double complex *v, double complex *x) {

  double complex m;
  long i;

  for (i = 1; i < n; i++) {
    m = a[i] / b[i-1];
    b[i] = b[i] - m * c[i-1];
    v[i] = v[i] - m * v[i-1];
  }
 
  /* a not used anymore (it is possible to have a = x) */
  x[n-1] = v[n-1] / b[n-1];
 
  for (i = n - 2; i >= 0; --i)
    x[i] = (v[i] - c[i] * x[i+1]) / b[i];
}

/*
 * Solve tridiagonal matrix equation A x = b (Sherman-Morrison).
 * Alpha and beta specify the extreme non-zero elements in A.
 * (this arises from finite diff. and periodic boundary conds)
 *
 * n = dimensionality (number of equations).
 * a = subdiagonal (indexed 1 ... n-1). 0 unused.
 * b = diagonal (overwritten).
 * c = supdiagonal (indexed 0 ... n-2). n-1 unused.
 * v = right hand side vector (overwritten).
 * alpha = left corner matrix element for the last row of A (double complex).
 * beta  = right corner matrix element for the first rotw of A (double complex).
 * x = solution on exit (double complex *).
 * bb = worksapce1 (0, ..., n) (double complex).
 *
 * Notes:  - with finite difference, use alpha = beta = 1.0 for periodic boundary
 *           and -1.0 for both if anti periodc.
 *
 */

EXPORT void grid_solve_tridiagonal_system_cyclic(long n, double complex *a, double complex *b, double complex *c, double complex *v, double complex alpha, double complex beta, double complex *x, double complex *bb) {

  double complex gamma, fact;
  long i;

  gamma = -b[0];
  bb[0] = b[0] - gamma;
  bb[n-1] = b[n-1] - alpha * beta / gamma;
  for (i = 1; i < n-1; i++) bb[i] = b[i];
  grid_solve_tridiagonal_system(n, a, bb, c, v, x); /* both b and v overwritten */
  v[0] = gamma;
  v[n-1] = alpha;
  bb[0] = b[0] - gamma;       /* bb was overwritten - reconstruct */
  bb[n-1] = b[n-1] - alpha * beta / gamma;
  for (i = 1; i < n-1; i++) {
    v[i] = 0.0;
    bb[i] = b[i];
  }
  grid_solve_tridiagonal_system(n, a, bb, c, v, a); /* note: a and x may be the same arrays */
  fact = (x[0] + beta * x[n-1] / gamma) / (1.0 + a[0] + beta * a[n-1] / gamma);
  for(i = 0; i < n; i++)
    x[i] -= fact * a[i];
}

/*
 * Creates Cholesky decomposition for tridiagonal matrix with super- and sub-diagonals equal to one.
 *
 * a = Diagonal elements (double complex *). On output replaced by the decomposition.
 * n = Dimension of a (long).
 * w = 0: Dirichlet problem, 1: Neumann (int).
 *
 * For periodic boundaries, use the Sherman-Morrison formula - see grid2d_wf.c, for example.
 *
 * TODO: It seems that this routine does not work?
 *
 */

EXPORT void grid_cholesky_decomposition(double complex *a, long n, int w) {

  long i;
  double complex lsm, *ld = a;
  
  /* A   x = L.L^T x = b
   * L   y = b
   * L^T x = y
   * (Note: ^T is transpose not hermitian conjugate)
   * A is tridiagonal matrix with super and sub diagonals equal to 1
   * Diagonal of L, ld_i = sqrt( a_i - ls_(i-1) * ls_(i-1) ) = sqrt( a_i -  1 / ( ld_(i-1) * ld_(i-1) ) )
   * Subdiagonal of L, ls_i = 1 / ld_i 
   */

  lsm = 0.0;
  for(i = 0; i < n; i++) {
    /* ld_i = sqrt(a_i - ls_(i-1) * ls_(i-1)) */
    // added conj()
    if(w && i == n-1) 
      ld[i] = csqrt(a[i] - 4.0 * conj(lsm) * lsm);  /* last row: ... 0 2 X */
    else
      ld[i] = csqrt(a[i] - conj(lsm) * lsm);  /* last row: ... 0 1 X */
    /* ls_(i-1) = 1 / ld_(i-1) */
    lsm = 1.0 / ld[i];
  }
}

/*
 * Solves A x = b, where A is tridiagonal matrix with super and sub diagonals equal to one using
 * given Cholesky decomposition.
 *
 * ld = Cholsky decomposition computed by grid_cholesky_decomposition() (double complex).
 * b  = RHS of Ax = b (double complex *). On exit, contains the solution x.
 * n  = Dimension of ld and b (long).
 * w  = 0: Dirichlet, 1: Neumann (int).
 *
 * TODO: It seems that Neumann does not work!? Use grid_solve_tridiagonal_system for now!
 *
 */

EXPORT void grid_cholesky_substitute(double complex *ld, double complex *b, long n, int w) {

  long i;
  double complex *y = b, *x = b;
  
  /*
   * L y = b
   * y_0 = b_0 / ld_0 
   *
   *        b_i - ls_(i-1) y_(i-1)  
   * y_i = ------------------------ 
   *                ld_i           
   */
  y[0] = b[0] / ld[0];
  for(i = 1; i < n; i++)
    /*
     *        b_i - ls_(i-1) y_(i-1)     b_i ld_(i-1) - y_(i-1)
     * y_i = ------------------------ = ------------------------
     *                ld_i                    ld_i  ld_(i-1)
     */
    // debug  !!!! what's this?
    if(w && i == n-1) 
      y[i] = (b[i] * ld[i-1] - 2.0 * y[i-1]) / (ld[i] * ld[i-1]);
    else
      y[i] = (b[i] * ld[i-1] - y[i-1]) / (ld[i] * ld[i-1]);
 
  /*
   * L^T x = y
   * x_(n-1) = y_(n-1) / ld_(n-1)
   *
   *        y_i - ls_(i+1) x_(i+1)  
   * x_i = ------------------------ 
   *               ld_i             
   */
  // added conj()
  x[n-1] = y[n-1] / conj(ld[n-1]);

  for(i = n - 2; i >= 0; i--)
    /*
     *        y_i - ls_i x_(i+1)     y_i - x_(i+1) / ld_i     y_i ld_i - x_(i+1)
     * x_i = -------------------- = ---------------------- = --------------------
     *               ld_i                     ld_i                ld_i  ld_i
     */
    // added conj()
    if(w && i == 0)
      x[i] = (y[i] * conj(ld[i]) - 2.0 * x[i+1]) / (conj(ld[i]) * conj(ld[i]));
    else
      x[i] = (y[i] * conj(ld[i]) - x[i+1]) / (conj(ld[i]) * conj(ld[i]));
}

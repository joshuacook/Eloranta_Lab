int itp_method_test(double * H, int n,int print_mat, int DEBUG){

  // declare necessary variables
  int i;
  double err, err1, err2, mu, mu_inv;
  double eps = 10E-6;
  double I[n*n];
  double Htemp[n*n];
  double HdotH[n*n];
  double CayleyN[n*n];
  double CayleyP_inv[n*n];
  double d_zero = 0.0;
  double d_neghalf = -0.5;
  double d_poshalf = 0.5;
  double d_one = 1.0;					 
  int one = 1;
  int two = 2;
  char no_trans='N';
  char trans='T';	
  double phi0[n];
  double phi1[n];
  double phitemp[n];

  identity_matrix(I,n);

 	// run test on a matrix of size n
 	if (DEBUG) printf("Find an eigenvector for an %d by %d matrix\n", n,n);

  // multiply H by its transpose to make it symmetric and thus Hermitian
  dgemm_(&no_trans,&trans,&n,&n,&n,&d_one,H,&n,H,&n,&d_zero,Htemp,&n);
  dgemm_(&no_trans,&trans,&n,&n,&n,&d_one,Htemp,&n,I,&n,&d_zero,H,&n);
  if (DEBUG) printf("Generated matrix\n");
 	if (DEBUG) print_matrix(H,n,n);

  err = 1;

 	// randomize phi

  random_vector(phi0,n);
  if (DEBUG) printf("random vector\n");
  if (DEBUG) print_vector(phi0,n);

 	// take the Cayley form of H

  // CayleyN = (one*I*I-0.5*CayleyN)
  if (DEBUG) printf("CayleyN of generated matrix\n");
  dgemm_(&no_trans,&no_trans,&n,&n,&n,&d_one,H,&n,I,&n,&d_zero,CayleyN,&n);
  dgemm_(&no_trans,&no_trans,&n,&n,&n,&d_one,I,&n,I,&n,&d_neghalf,CayleyN,&n);
  if (DEBUG) print_matrix(CayleyN,n,n);
 	// CayleyP = (one*I*I+0.5*CayleyP)
 	if (DEBUG) printf("CayleyP of generated matrix\n");
  dgemm_(&no_trans,&no_trans,&n,&n,&n,&d_one,H,&n,I,&n,&d_zero,CayleyP_inv,&n);
  dgemm_(&no_trans,&no_trans,&n,&n,&n,&d_one,I,&n,I,&n,&d_poshalf,CayleyP_inv,&n);
  if (DEBUG) print_matrix(CayleyP_inv,n,n);

  // preInvert
  if (DEBUG) printf("Inverse of CayleyP\n");
  inverse(CayleyP_inv,n);
  if (DEBUG) print_matrix(CayleyP_inv,n,n); 

  // iterate 
 	// repeatedly solve phi1 = CayleyP_inv*CayleyN*phi0 
  while(err > eps){
  dgemv_(&no_trans,&n,&n, &d_one,CayleyN,&n,phi0,&one,&d_zero,phitemp,&one);
  dgemv_(&no_trans,&n,&n, &d_one,CayleyP_inv,&n,phitemp,&one,&d_zero,phi1,&one);
  if (DEBUG) printf("phi1: ");
  if (DEBUG) print_vector(phi1,n);
  mu = dnrm2_(&n,phi1,&one);
  mu_inv = 1/mu;
  
  
  // should probably write separate method for finding the error this logic is pretty complicated and hard to follow
  if (DEBUG) printf("%f, %f\n\n", mu, mu_inv);
  dscal_(&n,&mu_inv,phi1,&one);
  if (DEBUG) print_vector(phi1,n);
	// 	err = math.sqrt(2)*math.sqrt(abs(phi1.dot(H.dot(H)).dot(phi1)- (phi1.dot(H).dot(phi1))**2))
  dgemv_(&no_trans,&n,&n, &d_one,HdotH,&n,phi1,&one,&d_zero,phitemp,&one);
  if (DEBUG) printf("phitemp: ");
  if (DEBUG) print_vector(phitemp,n);
  err1 = ddot_(&n,phi1,&one,phitemp,&one);
  if (DEBUG) printf("err1: %f\n", err1);
  dgemv_(&no_trans,&n,&n, &d_one,H,&n,phi1,&one,&d_zero,phitemp,&one);
  err2 = ddot_(&n,phi1,&one,phitemp,&one);
  err2 = err2 * err2;
  if (DEBUG) printf("err2: %f\n", err2);
  err = sqrt(abs(err1)-err2);
  if (DEBUG) printf("err: %f\n", err);
  dcopy_(&n,phi1,&one,phi0,&one);
  }

  print_vector(phi1, n); 

	return 0; 	
}


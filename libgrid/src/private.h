/* linear_weight */
static inline double linear_weight(double x) {

  if (x < 0.0) x = -x;
  return (1.0 - x);
}

static inline double sqnorm(double complex x) {

  return creal(x) * creal(x) + cimag(x) * cimag(x);
}

static inline double rms(double *data, int n) {

  int i; 
  double sum = 0;

  for(i = 0; i < n; i++)
    sum += data[i] * data[i];
  return sqrt(sum / n);
}

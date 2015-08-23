#ifndef __GRID1DPRIVATE__
#define __GRID1DPRIVATE__

typedef struct sShiftParametersc1d_struct {
  double x;
  const cgrid1d *grid;
} sShiftParametersc1d;

typedef struct sShiftParametersr1d_struct {
  double x;
  const rgrid1d *grid;
} sShiftParametersr1d;


/* linearly weighted integral
 *
 * in 1d:
 *      x-step        x         x+step
 *                    1
 * ns = 1 |-----------|-----------| 
 *             1/2         1/2
 * ns = 2 |-----+-----|-----+-----|
 *           1/3     1/3     1/3
 * ns = 3 |---+-------|-------+---| 
 */

static inline double complex linearly_weighted_integralc1d(double complex (*func)(void *arg, double x), void *farg, double x, double step, int ns) {

  int i;
  double xs;
  double wx;
  double substep = 2.0 * step / ns;
  double complex sum;
  double wsum;
  
  sum = 0.0; wsum = 0.0;
  
  for( i = 0; i < ns; i++ ) {
    xs = -step + (i + .5) * substep;
    wx = linear_weight(xs / step);
    /* sum += w * func(x,y,z) */
    sum += wx * func(farg, x + xs);
    wsum += wx;
  }
  
  return sum / wsum;
}

/* linearly weighted integral
 *
 * in 1d:
 *      x-step        x         x+step
 *                    1
 * ns = 1 |-----------|-----------| 
 *             1/2         1/2
 * ns = 2 |-----+-----|-----+-----|
 *           1/3     1/3     1/3
 * ns = 3 |---+-------|-------+---| 
 */

static inline double linearly_weighted_integralr1d(double (*func)(void *arg, double x), void *farg, double x, double step, int ns) {

  int i;
  double xs;
  double wx;
  double substep = 2.0 * step / ns;
  double sum;
  double wsum;
  
  sum = 0.0; wsum = 0.0;
  
  for( i = 0; i < ns; i++ ) {
    xs = -step + (i + .5) * substep;
    wx = linear_weight(xs / step);
    /* sum += w * func(x,y,z) */
    sum += wx * func(farg, x + xs);
    wsum += wx;
  }
  
  return sum / wsum;
}

static inline double complex shift_cgrid1d(void *arg, double x) {

  sShiftParametersc1d *params = (sShiftParametersc1d *) arg;

  return cgrid1d_value(params->grid, x - params->x);
}

static inline double shift_rgrid1d(void *arg, double x) {

  sShiftParametersr1d *params = (sShiftParametersr1d *) arg;

  return rgrid1d_value(params->grid, x - params->x);
}

#endif /* __GRID1DPRIVATE__ */

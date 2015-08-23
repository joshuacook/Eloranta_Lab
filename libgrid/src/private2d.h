#ifndef __GRID2DPRIVATE__
#define  __GRID2DPRIVATE__

typedef struct sShiftParametersc2d_struct {
  double x,y;
  const cgrid2d *grid;
} sShiftParametersc2d;

typedef struct sShiftParametersr2d_struct {
  double x,y;
  const rgrid2d *grid;
} sShiftParametersr2d;


/* linearly weighted integral
 *
 * in 1D:
 *      x-step        x         x+step
 *                    1
 * ns = 1 |-----------|-----------| 
 *             1/2         1/2
 * ns = 2 |-----+-----|-----+-----|
 *           1/3     1/3     1/3
 * ns = 3 |---+-------|-------+---| 
 */
static inline double complex linearly_weighted_integralc2d(double complex (*func)(void *arg, double x, double y), void *farg, double x, double y, double step, int ns) {

  int i,j;
  double xs, ys;
  double w, wx, wy;
  double substep = 2.0 * step / ns;
  double complex sum;
  double wsum;
  
  sum = 0.0;
  wsum = 0.0;
  
  for(i = 0; i < ns; i++) {
    xs = -step + (i + 0.5) * substep;
    wx = linear_weight(xs / step);
    
    for(j = 0; j < ns; j++) {
      ys = -step + (j + 0.5) * substep;
      wy = linear_weight(ys / step);

      w = wx * wy;
      sum += w * func(farg, x + xs, y + ys);
      wsum += w;
    }
  }
  
  return sum / wsum;
}

/* linearly weighted integral
 *
 * in 1D:
 *      x-step        x         x+step
 *                    1
 * ns = 1 |-----------|-----------| 
 *             1/2         1/2
 * ns = 2 |-----+-----|-----+-----|
 *           1/3     1/3     1/3
 * ns = 3 |---+-------|-------+---| 
 */
static inline double linearly_weighted_integralr2d(double (*func)(void *arg, double x, double y), void *farg, double x, double y, double step, int ns) {

  int i,j;
  double xs, ys;
  double w, wx, wy;
  double substep = 2.0 * step / ns;
  double sum;
  double wsum;
  
  sum = 0.0;
  wsum = 0.0;
  
  for(i = 0; i < ns; i++) {
    xs = -step + (i + 0.5) * substep;
    wx = linear_weight(xs / step);
    
    for(j = 0; j < ns; j++) {
      ys = -step + (j + 0.5) * substep;
      wy = linear_weight(ys / step);

      w = wx * wy;
      sum += w * func(farg, x + xs, y + ys);
      wsum += w;
    }
  }
  
  return sum / wsum;
}

static inline double complex shift_cgrid2d(void *arg, double x, double y) {

  sShiftParametersc2d *params = (sShiftParametersc2d *) arg;

  return cgrid2d_value(params->grid, x - params->x, y - params->y);
}

static inline double shift_rgrid2d(void *arg, double x, double y) {

  sShiftParametersr2d *params = (sShiftParametersr2d *) arg;

  return rgrid2d_value(params->grid, x - params->x, y - params->y);
}

#endif /*  __GRID2DPRIVATE__ */

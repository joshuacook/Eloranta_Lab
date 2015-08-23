#ifndef __GRID3DPRIVATE__
#define  __GRID3DPRIVATE__

typedef struct sShiftParametersc3d_struct {
  double x, y, z;
  const cgrid3d *grid;
} sShiftParametersc3d;

typedef struct sShiftParametersr3d_struct {
  double x, y, z;
  const rgrid3d *grid;
} sShiftParametersr3d;

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
static inline double complex linearly_weighted_integralc3d(double complex (*func)(void *arg, double x, double y, double z), void *farg, double x, double y, double z, double step, int ns ) {

  int i,j,k;
  double xs, ys, zs;
  double w, wx, wy, wz;
  double substep = 2 * step / ns;
  double complex sum;
  double wsum;
  
  sum = 0.0;
  wsum = 0.0;
  
  for(i = 0; i < ns; i++) {
    xs = -step + (i + .5) * substep;
    wx = linear_weight(xs / step);
    
    for(j = 0; j < ns; j++) {
      ys = -step + (j + .5) * substep;
      wy = linear_weight(ys / step);
      
      for(k = 0; k < ns; k++) {
        zs = -step + (k + .5) * substep;
        wz = linear_weight(zs / step);
        
        w = wx * wy * wz;
	
	sum += w * func(farg, x + xs, y + ys, z + zs);
        wsum += w;
      }
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
static inline double linearly_weighted_integralr3d(double (*func)(void *arg, double x, double y, double z), void *farg, double x, double y, double z, double step, int ns) {

  int i,j,k;
  double xs, ys, zs;
  double w, wx, wy, wz;
  double substep = 2 * step / ns;
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
      
      for(k = 0; k < ns; k++) {
        zs = -step + (k + 0.5) * substep;
        wz = linear_weight(zs / step);
        
        w = wx * wy * wz;
	
	sum += w * func(farg, x + xs, y + ys, z + zs);
        wsum += w;
      }
    }
  }
  
  return sum / wsum;
}

static inline double complex shift_cgrid3d(void *arg, double x, double y, double z) {

  sShiftParametersc3d *params = (sShiftParametersc3d *) arg;

  return cgrid3d_value(params->grid, x - params->x, y - params->y, z - params->z);
}

static inline double shift_rgrid3d(void *arg, double x, double y, double z) {

  sShiftParametersr3d *params = (sShiftParametersr3d *) arg;

  return rgrid3d_value(params->grid, x - params->x, y - params->y, z - params->z);
}

#endif /*  __GRID3DPRIVATE__ */

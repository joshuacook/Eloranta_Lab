#define GRID_EPS 1.0E-12

#define CGRID1D_DIRICHLET_BOUNDARY cgrid1d_value_outside_constantdirichlet
#define CGRID1D_NEUMANN_BOUNDARY cgrid1d_value_outside_neumann
#define CGRID1D_PERIODIC_BOUNDARY cgrid1d_value_outside_periodic

#define RGRID1D_DIRICHLET_BOUNDARY rgrid1d_value_outside_constantdirichlet
#define RGRID1D_NEUMANN_BOUNDARY rgrid1d_value_outside_neumann
#define RGRID1D_PERIODIC_BOUNDARY rgrid1d_value_outside_periodic

#define CGRID2D_DIRICHLET_BOUNDARY cgrid2d_value_outside_constantdirichlet
#define CGRID2D_NEUMANN_BOUNDARY cgrid2d_value_outside_neumann
#define CGRID2D_PERIODIC_BOUNDARY cgrid2d_value_outside_periodic

#define RGRID2D_DIRICHLET_BOUNDARY rgrid2d_value_outside_constantdirichlet
#define RGRID2D_NEUMANN_BOUNDARY rgrid2d_value_outside_neumann
#define RGRID2D_PERIODIC_BOUNDARY rgrid2d_value_outside_periodic

#define CGRID3D_DIRICHLET_BOUNDARY cgrid3d_value_outside_constantdirichlet
#define CGRID3D_NEUMANN_BOUNDARY cgrid3d_value_outside_neumann
#define CGRID3D_PERIODIC_BOUNDARY cgrid3d_value_outside_periodic

#define RGRID3D_DIRICHLET_BOUNDARY rgrid3d_value_outside_constantdirichlet
#define RGRID3D_NEUMANN_BOUNDARY rgrid3d_value_outside_neumann
#define RGRID3D_PERIODIC_BOUNDARY rgrid3d_value_outside_periodic

/* Special boundaries for vortex solutions in superfluid helium */
#define CGRID3D_VORTEX_X_BOUNDARY cgrid3d_value_outside_vortex_x
#define CGRID3D_VORTEX_Y_BOUNDARY cgrid3d_value_outside_vortex_y
#define CGRID3D_VORTEX_Z_BOUNDARY cgrid3d_value_outside_vortex_z
/* Warning: The above DO NOT give the usual periodic boundaries */

#define WF1D_DIRICHLET_BOUNDARY 1
#define WF1D_NEUMANN_BOUNDARY   2
#define WF1D_PERIODIC_BOUNDARY  3

#define WF1D_2ND_ORDER_PROPAGATOR 2
#define WF1D_4TH_ORDER_PROPAGATOR 4

#define WF2D_DIRICHLET_BOUNDARY 1
#define WF2D_NEUMANN_BOUNDARY   2
#define WF2D_PERIODIC_BOUNDARY  3

#define WF2D_2ND_ORDER_PROPAGATOR 2
#define WF2D_4TH_ORDER_PROPAGATOR 4

#define WF3D_DIRICHLET_BOUNDARY     1
#define WF3D_NEUMANN_BOUNDARY       2
#define WF3D_PERIODIC_BOUNDARY      3
#define WF3D_VORTEX_X_BOUNDARY      4
#define WF3D_VORTEX_Y_BOUNDARY      5
#define WF3D_VORTEX_Z_BOUNDARY      6

#define WF3D_2ND_ORDER_PROPAGATOR 2
#define WF3D_4TH_ORDER_PROPAGATOR 4

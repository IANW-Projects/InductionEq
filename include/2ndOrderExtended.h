//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

//------------------------------------------------------------------------------
// Derivative and Dissipation Operators
//------------------------------------------------------------------------------

#define ORDER 2
#define NUM_BOUNDS 3
#define STENCIL_WIDTH 5
#define STENCIL_WIDTH_HOD 2 //Stencil width for high order dissipation (HOD)

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff[2*NUM_BOUNDS+1][STENCIL_WIDTH] = {
  // left boundary coefficients
  {0.0, 0.0, -6.0/5.0, 7.0/5.0, -1.0/5.0},
  {0.0, -1.0/2.0, 0.0, 1.0/2.0, 0.0},
  {1.0/11.0, -7.0/11.0, 0.0, 6.0/11.0, 0.0},
  // central coefficients
  {0.0, -1.0/2.0, 0.0, 1.0/2.0, 0.0},
  // right boundary coefficients
  {0.0, -6.0/11.0, 0.0, 7.0/11.0, -1.0/11.0},
  {0.0, -1.0/2.0, 0.0, 1.0/2.0, 0.0},
  {1.0/5.0, -7.0/5.0, 6.0/5.0, 0.0, 0.0}
};


/*
Coefficients of the inverse mass/norm matrix.


The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV[2*NUM_BOUNDS+1] = {
  // left boundary coefficients
  12.0/5.0,
  6.0/7.0,
  12.0/11.0,
  // central coefficients
  1.0,
  // right boundary coefficients
  12.0/11.0,
  6.0/7.0,
  12.0/5.0
};


/*
Coeffcients of the high-order dissipation operator.

TODO: Explanation...
*/
REAL constant D_HO[2*NUM_BOUNDS+1][STENCIL_WIDTH_HOD] = {
  // left boundary coefficients
  { 0.0, 1.0},
  {-1.0, 1.0},
  {-1.0, 1.0},
  // central coefficients
  {-1.0, 1.0},
  // right boundary coefficients
  {-1.0, 1.0},
  {-1.0, 1.0},
  {-1.0, 0.0}
};



//------------------------------------------------------------------------------
// Divergence Cleaning Operators
//------------------------------------------------------------------------------

// coefficients for the wide stencil operator using homogeneous Dirichlet boundary conditions
#define NUM_BOUNDS_LAPLACE_WS_D0 4
#define STENCIL_WIDTH_LAPLACE_WS_D0 5

REAL constant SBP_minus_laplace_WS_D0[2*NUM_BOUNDS_LAPLACE_WS_D0+1][STENCIL_WIDTH_LAPLACE_WS_D0] = {
  // left boundary coefficients
  {0.0, 0.0, 0.0, 0.0, 0.0},
  {0.0, 0.0, 56.0/55.0, -1.0/10.0, -3.0/11.0},
  {0.0, -7.0/55.0, 67.0/110.0, 0.0, -3.0/11.0},
  {-7.0/22.0, 0.0, 23.0/44.0, 0.0, -1.0/4.0},
  // central coefficients
  {-1.0/4.0, 0.0, 1.0/2.0, 0.0, -1.0/4.0},
  // right boundary coefficients
  {-1.0/4.0, 0.0, 23.0/44.0, 0.0, -7.0/22.0},
  {-3.0/11.0, 0.0, 67.0/110.0, -7.0/55.0, 0.0},
  {-3.0/11.0, -1.0/10.0, 56.0/55.0, 0.0, 0.0},
  {0.0, 0.0, 0.0, 0.0, 0.0}
};

// coefficients for the narrow stencil operator using homogeneous Dirichlet boundary conditions
// have not been derived
#ifdef USE_LAPLACE_NARROW_STENCIL_DIRICHLET
#error "Error in include/2ndOrderExtended.h: Coefficients of the narrow stencil operator using homogeneous Dirichlet boundary conditions have not been derived."
#endif

// coefficients for the wide stencil operator using the least norm solution
#define NUM_BOUNDS_LAPLACE_WS_LN 4
#define STENCIL_WIDTH_LAPLACE_WS_LN 7

REAL constant SBP_minus_laplace_WS_LN[2*NUM_BOUNDS_LAPLACE_WS_LN+1][STENCIL_WIDTH_LAPLACE_WS_LN] = {
  // left boundary coefficients
  {0.0, 0.0, 0.0, 1187.0/550.0, 427.0/275.0, -47.0/50.0, 6.0/55.0},
  {0.0, 0.0, 61.0/110.0, 56.0/55.0, -1.0/10.0, -3.0/11.0, 0.0},
  {0.0, -47.0/110.0, -7.0/55.0, 67.0/110.0, 0.0, -3.0/11.0, 0.0},
  {1.0/22.0, -7.0/22.0, 0.0, 23.0/44.0, 0.0, -1.0/4.0, 0.0},
  // central coefficients
  {0.0, -1.0/4.0, 0.0, 1.0/2.0, 0.0, -1.0/4.0, 0.0},
  // right boundary coefficients
  {0.0, -1.0/4.0, 0.0, 23.0/44.0, 0.0, -7.0/22.0, 1.0/22.0},
  {0.0, -3.0/11.0, 0.0, 67.0/110.0, -7.0/55.0, -47.0/110.0, 0.0},
  {0.0, -3.0/11.0, -1.0/10.0, 56.0/55.0, 61.0/110.0, 0.0, 0.0},
  {6.0/55.0, -47.0/50.0, 427.0/275.0, 1187.0/550.0, 0.0, 0.0, 0.0}
};

REAL constant SBP_diff_adjoint[2*NUM_BOUNDS+1][STENCIL_WIDTH] = {
  // left boundary coefficients
  {0.0, 0.0, -6.0/5.0, -7.0/5.0, 1.0/5.0},
  {0.0, 1.0/2.0, 0.0, -1.0/2.0, 0.0},
  {-1.0/11.0, 7.0/11.0, 0.0, -6.0/11.0, 0.0},
  // central coefficients
  {0.0, 1.0/2.0, 0.0, -1.0/2.0, 0.0},
  // right boundary coefficients
  {0.0, 6.0/11.0, 0.0, -7.0/11.0, 1.0/11.0},
  {0.0, 1.0/2.0, 0.0, -1.0/2.0, 0.0},
  {-1.0/5.0, 7.0/5.0, 6.0/5.0, 0.0, 0.0}
};

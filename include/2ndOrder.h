//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

//Contains the coefficients of second order derivative, dissipation and divergence cleaning operators

//------------------------------------------------------------------------------
// Derivative and Dissipation Operators
//------------------------------------------------------------------------------

#define ORDER 2
#define NUM_BOUNDS 1
#define STENCIL_WIDTH 3
#define STENCIL_WIDTH_HOD 2 //Stencil width for high order dissipation (HOD)

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff[2*NUM_BOUNDS+1][STENCIL_WIDTH] = {
  // left boundary coefficients
  {0.0, -1.0, 1.0},
  // central coefficients
  {-1.0/2.0, 0.0, 1.0/2.0},
  // right boundary coefficients
  {-1.0, 1.0, 0.0}
};


/*
Coefficients of the inverse mass/norm matrix.

The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV[2*NUM_BOUNDS+1] = {
  // left boundary coefficients
  2.0,
  // central coefficients
  1.0,
  // right boundary coefficients
  2.0
};


/*
Coeffcients of the high-order dissipation operator.

TODO: Explanation...
*/
REAL constant D_HO[2*NUM_BOUNDS+1][STENCIL_WIDTH_HOD] = {
  // left boundary coefficients
  { 0.0, 1.0},
  // central coefficients
  {-1.0, 1.0},
  // right boundary coefficients
  {-1.0, 0.0}
};


//------------------------------------------------------------------------------
// Divergence Cleaning Operators
//------------------------------------------------------------------------------

// coefficients for the wide stencil operator using homogeneous Dirichlet boundary conditions
#define NUM_BOUNDS_LAPLACE_WS_D0 3
#define STENCIL_WIDTH_LAPLACE_WS_D0 5

REAL constant SBP_minus_laplace_WS_D0[2*NUM_BOUNDS_LAPLACE_WS_D0+1][STENCIL_WIDTH_LAPLACE_WS_D0] = {
  // left boundary coefficients
  {0.0, 0.0, 0.0, 0.0, 0.0},
  {0.0, 0.0, 3.0/4.0, 0.0, -1.0/4.0},
  {0.0, 0.0, 1.0/2.0, 0.0, -1.0/4.0},
  // central coefficients
  {-1.0/4.0, 0.0, 1.0/2.0, 0.0, -1.0/4.0},
  // right boundary coefficients
  {-1.0/4.0, 0.0, 1.0/2.0, 0.0, 0.0},
  {-1.0/4.0, 0.0, 3.0/4.0, 0.0, 0.0},
  {0.0, 0.0, 0.0, 0.0, 0.0}
};

// coefficients for the narrow stencil operator using homogeneous Dirichlet boundary conditions
#define NUM_BOUNDS_LAPLACE_NS_D0 2
#define STENCIL_WIDTH_LAPLACE_NS_D0 3

REAL constant SBP_minus_laplace_NS_D0[2*NUM_BOUNDS_LAPLACE_NS_D0+1][STENCIL_WIDTH_LAPLACE_NS_D0] = {
  // left boundary coefficients
  {0.0, 0.0, 0.0},
  {0.0, 2.0, -1.0},
  // central coefficients
  {-1.0, 2.0, -1.0},
  // right boundary coefficients
  {-1.0, 2.0, 0.0},
  {0.0, 0.0, 0.0}
};

// coefficients for the wide stencil operator using the least norm solution
#define NUM_BOUNDS_LAPLACE_WS_LN 2
#define STENCIL_WIDTH_LAPLACE_WS_LN 5

REAL constant SBP_minus_laplace_WS_LN[2*NUM_BOUNDS_LAPLACE_WS_LN+1][STENCIL_WIDTH_LAPLACE_WS_LN] = {
  // left boundary coefficients
  {0.0, 0.0, 3.0/2.0, 1.0, -1.0/2.0},
  {0.0, 1.0/2.0, 3.0/4.0, 0.0, -1.0/4.0},
  // central coefficients
  {-1.0/4.0, 0.0, 1.0/2.0, 0.0, -1.0/4.0},
  // right boundary coefficients
  {-1.0/4.0, 0.0, 3.0/4.0, 1.0/2.0, 0.0},
  {-1.0/2.0, 1.0, 3.0/2.0, 0.0, 0.0}
};

// coefficients for the adjoint differatial operator needed for the least norm solution
REAL constant SBP_diff_adjoint[2*NUM_BOUNDS+1][STENCIL_WIDTH] = {
  // left boundary coefficients
  {0.0, -1.0, -1.0},
  // central coefficients
  {1.0/2.0, 0.0, -1.0/2.0},
  // right boundary coefficients
  {1.0, 1.0, 0.0}
};

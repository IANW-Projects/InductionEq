//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

//------------------------------------------------------------------------------
// Derivative and Dissipation Operators
//------------------------------------------------------------------------------

#define ORDER 2
#define NUM_BOUNDS 3
#define STENCIL_WIDTH 5
#define STENCIL_WIDTH_HOD 2

/*
Coefficients of the first derivative operator.

For each row, the central coefficient corresponds to the nodes at which the
derivative is computed.
*/
REAL constant SBP_diff[2*NUM_BOUNDS+1][STENCIL_WIDTH] = {
  // left boundary coefficients
  {0.0, 0.0, -1.2, 1.4, -0.2},
  {0.0, -0.5, 0.0, 0.5, 0.0},
  {0.09090909090909091, -0.6363636363636364, 0.0, 0.5454545454545454, 0.0},
  // central coefficients
  {0.0, -0.5, 0.0, 0.5, 0.0},
  // right boundary coefficients
  {0.0, -0.5454545454545454, 0.0, 0.6363636363636364, -0.09090909090909091},
  {0.0, -0.5, 0.0, 0.5, 0.0},
  {0.2, -1.4, 1.2, 0.0, 0.0}
};


/*
Coefficients of the inverse mass/norm matrix.


The central coefficient corresponds to the inner nodes of the grid and is 1.
*/
REAL constant M_INV[2*NUM_BOUNDS+1] = {
  // left boundary coefficients
  2.4,
  0.8571428571428571,
  1.0909090909090908,
  // central coefficients
  1.0,
  // right boundary coefficients
  1.0909090909090908,
  0.8571428571428571,
  2.4
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

//TODO: Implement the operators

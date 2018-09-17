//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

// Contains a number of differential operators for scalar and vector type fields (inline functions)

//--------------------------------------------------------------------------------------------------
// Operator for vector fields
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// Derivatives for vector fields
// Computes the derivative in x-direction of the vector field `d_field` at point (ix,iy,iz)
REAL4 inline diff_x(uint ix, uint iy, uint iz, global REAL4 *d_field) {

	REAL4 val = (REAL4) {0, 0, 0, 0};

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_xr(ix,i+1) - check_bound_l(ix,i+1));
	}

  for (uint i = 0; i < STENCIL_WIDTH; i++) {
    val = val + SBP_diff[NUM_BOUNDS + bound][i]*get_vector_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field);
  }

	return val / ((REAL)DX);
}

// Computes the derivative in y-direction of the vector field `d_field` at point (ix,iy,iz)
REAL4 inline diff_y(uint ix, uint iy, uint iz, global REAL4 *d_field) {

	REAL4 val = (REAL4) {0, 0, 0, 0};

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_yr(iy,i+1) - check_bound_l(iy,i+1));
	}

  for (uint i = 0; i < STENCIL_WIDTH; i++) {
    val = val + SBP_diff[NUM_BOUNDS + bound][i]*get_vector_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field);
  }

	return val / ((REAL)DY);
}

// Computes the derivative in z-direction of the vector field `d_field` at point (ix,iy,iz)
REAL4 inline diff_z(uint ix, uint iy, uint iz, global REAL4 *d_field) {

	REAL4 val = (REAL4) {0, 0, 0, 0};

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_zr(iz,i+1) - check_bound_l(iz,i+1));
	}

  for (uint i = 0; i < STENCIL_WIDTH; i++) {
    val = val + SBP_diff[NUM_BOUNDS + bound][i]*get_vector_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field);
  }

	return val / ((REAL)DZ);
}

//--------------------------------------------------------------------------------------------------
// Adjoint derivatives for vector fields
// Computes the adjoint derivative in x-direction of the vector field `d_field` at point (ix,iy,iz)
REAL4 inline diff_adj_x(uint ix, uint iy, uint iz, global REAL4 *d_field) {

	REAL4 val = (REAL4) {0, 0, 0, 0};

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_xr(ix,i+1) - check_bound_l(ix,i+1));
	}

	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		val = val + SBP_diff_adjoint[NUM_BOUNDS + bound][i]*get_vector_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field);
	}

	return val / ((REAL)DX);

}

// Computes the adjoint derivative in y-direction of the vector field `d_field` at point (ix,iy,iz)
REAL4 inline diff_adj_y(uint ix, uint iy, uint iz, global REAL4 *d_field) {

	REAL4 val = (REAL4) {0, 0, 0, 0};

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_yr(iy,i+1) - check_bound_l(iy,i+1));
	}

	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		val = val + SBP_diff_adjoint[NUM_BOUNDS + bound][i]*get_vector_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field);
	}

	return val / ((REAL)DY);
}

// Computes the adjoint derivative in z-direction of the vector field `d_field` at point (ix,iy,iz)
REAL4 inline diff_adj_z(uint ix, uint iy, uint iz, global REAL4 *d_field) {

	REAL4 val = (REAL4) {0, 0, 0, 0};

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_zr(iz,i+1) - check_bound_l(iz,i+1));
	}

	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		val = val + SBP_diff_adjoint[NUM_BOUNDS + bound][i]*get_vector_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field);
	}

	return val / ((REAL)DZ);
}

//--------------------------------------------------------------------------------------------------
// Divergence and adjoint divergence operator
// Computes the divergence of vector field `d_field` at point (ix,iy,iz)
inline REAL div(uint ix, uint iy, uint iz, global REAL4 *d_field) {

  REAL val = 0;

  val = diff_x(ix, iy, iz, d_field).x
      + diff_y(ix, iy, iz, d_field).y
      + diff_z(ix, iy, iz, d_field).z;

  return val;
}

// Computes the adjoint divergence of vector field `d_field` at point (ix,iy,iz)
inline REAL div_adj(uint ix, uint iy, uint iz, global REAL4 *d_field) {

  REAL val;

  val = diff_adj_x(ix, iy, iz, d_field).x
      + diff_adj_y(ix, iy, iz, d_field).y
      + diff_adj_z(ix, iy, iz, d_field).z;

  return val;
}

//--------------------------------------------------------------------------------------------------
// Computes the curl of vector field `d_field` at point (ix,iy,iz)
inline REAL4 curl(uint ix, uint iy, uint iz, global REAL4 *d_field) {

  REAL4 val = (REAL4) {0, 0, 0, 0};

  val.x = diff_y(ix, iy, iz, d_field).z - diff_z(ix, iy, iz, d_field).y;
  val.y = diff_z(ix, iy, iz, d_field).x - diff_x(ix, iy, iz, d_field).z;
  val.z = diff_x(ix, iy, iz, d_field).y - diff_y(ix, iy, iz, d_field).x;

  return val;
}

//--------------------------------------------------------------------------------------------------
// Operators for scalar fields
//--------------------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------------------
// Derivatives for scalar fields
// Computes the derivative in x-direction of the scalar field `d_field` at point (ix,iy,iz)
REAL inline diff_S_x(uint ix, uint iy, uint iz, global REAL *d_field) {

	REAL val = 0;

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_xr(ix,i+1) - check_bound_l(ix,i+1));
	}

  for (uint i = 0; i < STENCIL_WIDTH; i++) {
    val = val + SBP_diff[NUM_BOUNDS + bound][i]*get_scalar_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field);
  }

	return val / ((REAL)DX);
}

// Computes the derivative in y-direction of the scalar field `d_field` at point (ix,iy,iz)
REAL inline diff_S_y(uint ix, uint iy, uint iz, global REAL *d_field) {

	REAL val = 0;

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_yr(iy,i+1) - check_bound_l(iy,i+1));
	}

  for (uint i = 0; i < STENCIL_WIDTH; i++) {
    val = val + SBP_diff[NUM_BOUNDS + bound][i]*get_scalar_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field);
  }

	return val / ((REAL)DY);
}

// Computes the derivative in z-direction of the scalar field `d_field` at point (ix,iy,iz)
REAL inline diff_S_z(uint ix, uint iy, uint iz, global REAL *d_field) {

	REAL val = 0;

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_zr(iz,i+1) - check_bound_l(iz,i+1));
	}

  for (uint i = 0; i < STENCIL_WIDTH; i++) {
    val = val + SBP_diff[NUM_BOUNDS + bound][i]*get_scalar_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field);
  }

	return val / ((REAL)DZ);
}

//--------------------------------------------------------------------------------------------------
// Adjoint derivative for scalar fields
// Computes the adjoint derivative in x-direction of the scalar field `d_field` at point (ix,iy,iz)
REAL inline diff_adj_S_x(uint ix, uint iy, uint iz, global REAL *d_field) {

	REAL val = 0;

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_xr(ix,i+1) - check_bound_l(ix,i+1));
	}

	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		val = val + SBP_diff_adjoint[NUM_BOUNDS + bound][i]*get_scalar_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field);
	}

	return val / ((REAL)DX);
}

// Computes the adjoint derivative in y-direction of the scalar field `d_field` at point (ix,iy,iz)
REAL inline diff_adj_S_y(uint ix, uint iy, uint iz, global REAL *d_field) {

	REAL val = 0;

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_yr(iy,i+1) - check_bound_l(iy,i+1));
	}

	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		val = val + SBP_diff_adjoint[NUM_BOUNDS + bound][i]*get_scalar_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field);
	}

	return val / ((REAL)DY);
}

// Computes the adjoint derivative in z-direction of the scalar field `d_field` at point (ix,iy,iz)
REAL inline diff_adj_S_z(uint ix, uint iy, uint iz, global REAL *d_field) {

	REAL val = 0;

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		bound = bound + (NUM_BOUNDS - i)*(check_bound_zr(iz,i+1) - check_bound_l(iz,i+1));
	}

	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		val = val + SBP_diff_adjoint[NUM_BOUNDS + bound][i]*get_scalar_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field);
	}

	return val / ((REAL)DZ);
}


//--------------------------------------------------------------------------------------------------
// Gradient and adjoint gradient
// Computes the gradient of scalar field `d_field` at point (ix,iy,iz)
inline REAL4 grad(uint ix, uint iy, uint iz, global REAL *d_field) {

  REAL4 val;

  val.x = diff_S_x(ix, iy, iz, d_field);
  val.y = diff_S_y(ix, iy, iz, d_field);
  val.z = diff_S_z(ix, iy, iz, d_field);
  val.w = 0;

  return val;
}

// Computes the adjoint gradient of scalar field `d_field` at point (ix,iy,iz)
inline REAL4 grad_adj(uint ix, uint iy, uint iz, global REAL *d_field) {

  REAL4 val;

  val.x = diff_adj_S_x(ix, iy, iz, d_field);
  val.y = diff_adj_S_y(ix, iy, iz, d_field);
  val.z = diff_adj_S_z(ix, iy, iz, d_field);
  val.w = 0;

  return val;
}

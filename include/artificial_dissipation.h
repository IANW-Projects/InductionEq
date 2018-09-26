//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

// Contains auxillary function for artificial dissipation

//--------------------------------------------------------------------------------------------------
// Artifical Dissipation
//--------------------------------------------------------------------------------------------------

#ifdef USE_ADAPTIVE_DISSIPATION

// Needed to compute modified weights for adaptive dissipation
inline REAL kappa(uint k, REAL dx, REAL field_b1, REAL field_b2, REAL field_u1, REAL field_u2) {

  REAL s = fmax(fabs(field_u1), fabs(field_u2)) * (field_b1 - field_b2) * (field_b1 - field_b2);

  // Use fabs in sqrt because the sqrt of a negative number is NaN
  REAL val = ((bool)(k*k*MP2MIN*CMIN*dx < s) && (bool)(s < k*k*MP2MAX*CMAX*dx))
               * sqrt(fabs(s - k*k*MP2MIN*CMIN*dx)/fabs(k*k*MP2MAX*CMAX*dx - k*k*MP2MIN*CMIN*dx))
             + (bool)(s >= k*k*MP2MAX*CMAX*dx);
  return val;
}

inline REAL4 kappa4(uint k, REAL dx, REAL4 field_b1, REAL4 field_b2, REAL4 field_u1, REAL4 field_u2) {

  REAL4 s = fmax(fabs(field_u1), fabs(field_u2)) * (field_b1 - field_b2) * (field_b1 - field_b2);
  // k must not equal 0!
  k = k + (k == 0);

  REAL4 val = (REAL4) {0, 0, 0, 0};

  // Use fabs in sqrt because the sqrt of a negative number is NaN
  val.x = ((bool)(k*k*MP2MIN*CMIN*dx < s.x) && (bool)(s.x < k*k*MP2MAX*CMAX*dx))
            * sqrt(fabs(s.x - k*k*MP2MIN*CMIN*dx)/fabs(k*k*MP2MAX*CMAX*dx - k*k*MP2MIN*CMIN*dx))
          + (bool)(s.x >= k*k*MP2MAX*CMAX*dx);

  val.y = ((bool)(k*k*MP2MIN*CMIN*dx < s.y) && (bool)(s.y < k*k*MP2MAX*CMAX*dx))
            * sqrt(fabs(s.y - k*k*MP2MIN*CMIN*dx)/fabs(k*k*MP2MAX*CMAX*dx - k*k*MP2MIN*CMIN*dx))
          + (bool)(s.y >= k*k*MP2MAX*CMAX*dx);

  val.z = ((bool)(k*k*MP2MIN*CMIN*dx < s.z) && (bool)(s.z < k*k*MP2MAX*CMAX*dx))
            * sqrt(fabs(s.z - k*k*MP2MIN*CMIN*dx)/fabs(k*k*MP2MAX*CMAX*dx - k*k*MP2MIN*CMIN*dx))
          + (bool)(s.z >= k*k*MP2MAX*CMAX*dx);

  return val;
}

#endif // USE_ADAPTIVE_DISSIPATION

// Returns the componentwise absolut maximum of two vectors
// Used to computes the weights for first order and adaptive dissipation
REAL4 inline sigma(REAL4 d_field1, REAL4 d_field2) {
	return fmax(fabs(d_field1), fabs(d_field2));
}


// Auxillary function for high order dissipation
// Used to compute the intermediate result D*a of equation (?)
REAL4 inline diff_HOD_x(uint ix, uint iy, uint iz, int bound_x, int start, global REAL4 *d_field) {

  REAL4 HOD = (REAL4) {0, 0, 0, 0};

  for (int i = 0; i < STENCIL_WIDTH_HOD; i++) {
    HOD = HOD
      + M_INV[NUM_BOUNDS + bound_x]
        * D_HO[NUM_BOUNDS + bound_x][start + STENCIL_WIDTH_HOD - 1]
        * D_HO[NUM_BOUNDS][i] * pown(-1.0, (int)(ORDER / 2 + 1))
        * get_vector_field(ix, iy, iz, i + start, 0, 0, d_field);
  }

  return HOD;
}

REAL4 inline diff_HOD_y(uint ix, uint iy, uint iz, int bound_y, int start, global REAL4 *d_field) {

  REAL4 HOD = (REAL4) {0, 0, 0, 0};

  for (int i = 0; i < STENCIL_WIDTH_HOD; i++) {
    HOD = HOD
      + M_INV[NUM_BOUNDS + bound_y]
        * D_HO[NUM_BOUNDS + bound_y][start + STENCIL_WIDTH_HOD - 1]
        * D_HO[NUM_BOUNDS][i] * pown(-1.0, (int)(ORDER / 2 + 1))
        * get_vector_field(ix, iy, iz, 0, i + start, 0, d_field);
  }

  return HOD;
}

REAL4 inline diff_HOD_z(uint ix, uint iy, uint iz, int bound_z, int start, global REAL4 *d_field) {

  REAL4 HOD = (REAL4) {0, 0, 0, 0};

  for (int i = 0; i < STENCIL_WIDTH_HOD; i++) {
    HOD = HOD
      + M_INV[NUM_BOUNDS + bound_z]
        * D_HO[NUM_BOUNDS + bound_z][start + STENCIL_WIDTH_HOD - 1]
        * D_HO[NUM_BOUNDS][i] * pown(-1.0, (int)(ORDER / 2 + 1))
        * get_vector_field(ix, iy, iz, 0, 0, i + start, d_field);
  }

  return HOD;
}

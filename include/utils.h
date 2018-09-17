//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

#ifndef UTILS_H
#define UTILS_H

// Contains functions for boundary check, index calculation element access of scalar and vector type fields
//--------------------------------------------------------------------------------------------------
// Boundary checker
//--------------------------------------------------------------------------------------------------

// Checks if the node at the dimension specific index (ix, iy or iz) is at least b_nodes nodes away
// from the right box boundary.
inline bool check_interior_xr(uint idx, uint b_nodes) {
  return idx < ((uint)NODES_X - b_nodes);
}

inline bool check_interior_yr(uint idx, uint b_nodes) {
  return idx < ((uint)NODES_Y - b_nodes);
}

inline bool check_interior_zr(uint idx, uint b_nodes) {
  return idx < ((uint)NODES_Z - b_nodes);
}

// Checks if the node at the dimension specific index (ix, iy or iz) is at least b_nodes nodes away
// from the left box boundary.
inline bool check_interior_l(uint idx, uint b_nodes) {
  return idx > (0 + b_nodes - 1);
}

// Checks if node at the 3D index (ix iy, iz) is located inside the box
// (bx:NODES_X-bx, by:NODES_Y-by, bz:NODES_Z-bz).
inline bool check_interior_all(uint ix, uint iy, uint iz, uint bx, uint by, uint bz) {

  return check_interior_xr(ix,bx) && check_interior_l(ix,bx)
      && check_interior_yr(iy,by) && check_interior_l(iy,by)
      && check_interior_zr(iz,bz) && check_interior_l(iz,bz);
}


// Checks if the node at the dimension specific index (ix, iy or iz) is exactly b_nodes nodes away
// from the right box boundary.
inline bool check_bound_xr(uint idx, uint b_nodes) {
  return idx == ((uint)NODES_X - b_nodes);
}

inline bool check_bound_yr(uint idx, uint b_nodes) {
  return idx == ((uint)NODES_Y - b_nodes);
}

inline bool check_bound_zr(uint idx, uint b_nodes) {
  return idx == ((uint)NODES_Z - b_nodes);
}

// Checks if the node at the dimension specific index (ix, iy or iz) is exactly b_nodes nodes away
// from the left box boundary.
inline bool check_bound_l(uint idx, uint b_nodes) {
  return idx == (0 + b_nodes - 1);
}


//--------------------------------------------------------------------------------------------------
// Index related
//--------------------------------------------------------------------------------------------------

inline uint calc_idx(uint ix, uint iy, uint iz) {
  return (ix  + (iy*(uint)NODES_X)+(iz*(uint)NODES_Y*(uint)NODES_X));
}

// Calculated the three dimensional index (ix,iy,iz) from the one dimensional index idx needed for
// array access.
inline uint4 calc_sub_idx(uint idx){

  uint4 s_idx = (uint4) {0,0,0,0};

  s_idx.z = idx / (NODES_X*NODES_Y);

  s_idx.y = (idx - s_idx.z*NODES_X*NODES_Y) / NODES_X;

  s_idx.x = idx - s_idx.y*NODES_X - s_idx.z*NODES_Y*NODES_X;

  return s_idx;
}

// Calculate the boundary index (row of the derivative etc. operator) in x,y,z direction.
inline int get_bound_x(uint ix, uint num_bounds) {

  int bound_x = 0;

  for (uint i = 0; i < num_bounds; i++)
  {
    bound_x = bound_x + (num_bounds - i) * (check_bound_xr(ix, i + 1) - check_bound_l(ix, i + 1));
  }

  return bound_x;
}
inline int get_bound_y(uint iy, uint num_bounds) {

  int bound_y = 0;

  for (uint i = 0; i < num_bounds; i++)
  {
    bound_y = bound_y + (num_bounds - i) * (check_bound_yr(iy, i + 1) - check_bound_l(iy, i + 1));
  }

  return bound_y;
}
inline int get_bound_z(uint iz, uint num_bounds) {

  int bound_z = 0;

  for (uint i = 0; i < num_bounds; i++)
  {
    bound_z = bound_z + (num_bounds - i) * (check_bound_zr(iz, i + 1) - check_bound_l(iz, i + 1));
  }

  return bound_z;
}

#ifdef USE_PERIODIC

// Get the value of the vector field d_field at the 3D index (ix+bx, iy+by, iz+bz) with
// periodic boundary conditions.
inline REAL4 get_vector_field(uint ix, uint iy, uint iz, int bx, int by, int bz, global REAL4 *d_field) {

  uint n_ix = ix + bx + (bx < 0) * !check_interior_l(ix, abs(bx)) * NODES_X
                      - (bx > 0) * !check_interior_xr(ix, abs(bx)) * NODES_X;

  uint n_iy = iy + by + (by < 0) * !check_interior_l(iy, abs(by)) * NODES_Y
                      - (by > 0) * !check_interior_yr(iy, abs(by)) * NODES_Y;

  uint n_iz = iz + bz + (bz < 0) * !check_interior_l(iz, abs(bz)) * NODES_Z
                      - (bz > 0) * !check_interior_zr(ix, abs(bz)) * NODES_Z;

  uint idx = calc_idx(n_ix, n_iy, n_iz);

  return d_field[idx];
}

// Get the value of the scalar field d_field at the 3D index (ix+bx, iy+by, iz+bz) with
// periodic boundary conditions.
inline REAL get_scalar_field(uint ix, uint iy, uint iz, int bx, int by, int bz, global REAL *d_field) {

  uint n_ix = ix + bx + (bx < 0) * !check_interior_l(ix, abs(bx)) * NODES_X
                      - (bx > 0) * !check_interior_xr(ix, abs(bx)) * NODES_X;

  uint n_iy = iy + by + (by < 0) * !check_interior_l(iy, abs(by)) * NODES_Y
                      - (by > 0) * !check_interior_yr(iy, abs(by)) * NODES_Y;

  uint n_iz = iz + bz + (bz < 0) * !check_interior_l(iz, abs(bz)) * NODES_Z
                      - (bz > 0) * !check_interior_zr(ix, abs(bz)) * NODES_Z;

  uint idx = calc_idx(n_ix, n_iy, n_iz);

  return d_field[idx];
}

// Calculate the boundary index (row of the derivative etc. operator) in x,y,z direction. Since
// periodic boundaries are used, the central row is applied.
inline int get_bound_x(uint ix, uint num_bounds) {

  int bound_x = num_bounds + 1;

  return bound_x;
}
inline int get_bound_y(uint iy, uint num_bounds) {

  int bound_y = num_bounds + 1;

  return bound_y;
}
inline int get_bound_z(uint iz, uint num_bounds) {

  int bound_z = num_bounds + 1;

  return bound_z;
}

#else // nonperiodic boundaries

// Get the value of the vector field d_field at the 3D index (ix+bx, iy+by, iz+bz) if this index
// is inside the bounds. Indices out of bounds are set to ix, iy, or iz, respectively.
inline REAL4 get_vector_field(uint ix, uint iy, uint iz, int bx, int by, int bz, global REAL4 *d_field){

  uint n_ix = ix + ((bx < 0)*check_interior_l(ix,abs(bx)) + (bx > 0)*check_interior_xr(ix,abs(bx)))*bx;

  uint n_iy = iy + ((by < 0)*check_interior_l(iy,abs(by)) + (by > 0)*check_interior_yr(iy,abs(by)))*by;

  uint n_iz = iz + ((bz < 0)*check_interior_l(iz,abs(bz)) + (bz > 0)*check_interior_zr(iz,abs(bz)))*bz;

  uint idx = calc_idx(n_ix, n_iy, n_iz);

  return d_field[idx];
}

// Get the value of the scalar field d_field at the 3D index (ix+bx, iy+by, iz+bz) if this index
// is inside the bounds. Indices out of bounds are set to ix, iy, or iz, respectively.
inline REAL get_scalar_field(uint ix, uint iy, uint iz, int bx, int by, int bz, global REAL *d_field){

  uint n_ix = ix + ((bx < 0)*check_interior_l(ix,abs(bx)) + (bx > 0)*check_interior_xr(ix,abs(bx)))*bx;

  uint n_iy = iy + ((by < 0)*check_interior_l(iy,abs(by)) + (by > 0)*check_interior_yr(iy,abs(by)))*by;

  uint n_iz = iz + ((bz < 0)*check_interior_l(iz,abs(bz)) + (bz > 0)*check_interior_zr(iz,abs(bz)))*bz;

  uint idx = calc_idx(n_ix, n_iy, n_iz);

  return d_field[idx];
}

#endif // USE_PERIODIC

#endif //UTILS_H

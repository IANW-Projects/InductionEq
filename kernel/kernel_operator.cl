//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

// Contains kernel for norm and derivative operators for scalar and vetory type fields
//--------------------------------------------------------------------------------------------------
// Norm for scalar and vector type fields
//--------------------------------------------------------------------------------------------------

/*
Computes the squared L^2 norm of the vector field `d_field` and stores the result in `output`.
*/
kernel void norm2(global REAL4 *d_field, global REAL *output) {

  const int gid = get_global_id(0);
  const int lid = get_local_id(0);
  const int group_size = get_local_size(0);
  const int wgid = get_group_id(0);

  local REAL partial_n[(uint)W_SIZE];

  uint4 s_idx = calc_sub_idx(gid);

  int bound_x = 0;
  int bound_y = 0;
  int bound_z = 0;

  for (uint i = 0; i < NUM_BOUNDS; i++) {
    bound_x = bound_x + (NUM_BOUNDS - i)*(check_bound_xr(s_idx.x,i+1) - check_bound_l(s_idx.x,i+1));
    bound_y = bound_y + (NUM_BOUNDS - i)*(check_bound_yr(s_idx.y,i+1) - check_bound_l(s_idx.y,i+1));
    bound_z = bound_z + (NUM_BOUNDS - i)*(check_bound_zr(s_idx.z,i+1) - check_bound_l(s_idx.z,i+1));
  }

  REAL fac = ((REAL)DX / M_INV[NUM_BOUNDS + bound_x])
           * ((REAL)DY / M_INV[NUM_BOUNDS + bound_y])
           * ((REAL)DZ / M_INV[NUM_BOUNDS + bound_z]);


  partial_n[lid] = fac * dot(d_field[gid], d_field[gid]);

  barrier(CLK_LOCAL_MEM_FENCE);

  for (uint i = group_size/2; i > 0; i /= 2) {
    if(lid < i) {
      partial_n[lid] += partial_n[lid + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if (lid == 0) {
    #ifdef REDUCE
      atomic_add_global(&(output[0]), partial_n[0]);
    #else
      output[wgid] = partial_n[0];
    #endif
  }
}

/*
Computes the squared L^2 norm of the difference of the vector fields `d_field_1`, `d_field_2` and
stores the result in `output`.
*/
kernel void norm2_diff(global REAL4 *d_field_1, global REAL4 *d_field_2, global REAL *output) {

  const int gid = get_global_id(0);
  const int lid = get_local_id(0);
  const int group_size = get_local_size(0);
  const int wgid = get_group_id(0);

  local REAL partial_n[(uint)W_SIZE];

  uint4 s_idx = calc_sub_idx(gid);

  int bound_x = 0;
  int bound_y = 0;
  int bound_z = 0;

  for (uint i = 0; i < NUM_BOUNDS; i++) {
    bound_x = bound_x + (NUM_BOUNDS - i)*(check_bound_xr(s_idx.x,i+1) - check_bound_l(s_idx.x,i+1));
    bound_y = bound_y + (NUM_BOUNDS - i)*(check_bound_yr(s_idx.y,i+1) - check_bound_l(s_idx.y,i+1));
    bound_z = bound_z + (NUM_BOUNDS - i)*(check_bound_zr(s_idx.z,i+1) - check_bound_l(s_idx.z,i+1));
  }

  REAL fac = ((REAL)DX / M_INV[NUM_BOUNDS + bound_x])
           * ((REAL)DY / M_INV[NUM_BOUNDS + bound_y])
           * ((REAL)DZ / M_INV[NUM_BOUNDS + bound_z]);

  REAL4 diff = d_field_1[gid] - d_field_2[gid];
  partial_n[lid] = fac * dot(diff, diff);

  barrier(CLK_LOCAL_MEM_FENCE);

  for (uint i = group_size/2; i > 0; i /= 2) {
    if(lid < i) {
      partial_n[lid] += partial_n[lid + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if (lid == 0) {
    #ifdef REDUCE
      atomic_add_global(&(output[0]), partial_n[0]);
    #else
      output[wgid] = partial_n[0];
    #endif
  }
}


/*
Computes the squared L^2 norm of the scalar field `d_field` and stores the result in `output`.
*/
kernel void norm2_S(global REAL *d_field, global REAL *output) {

  const int gid = get_global_id(0);
  const int lid = get_local_id(0);
  const int group_size = get_local_size(0);
  const int wgid = get_group_id(0);

  local REAL partial_n[(uint)W_SIZE];

  uint4 s_idx = calc_sub_idx(gid);

  int bound_x = 0;
  int bound_y = 0;
  int bound_z = 0;

  for (uint i = 0; i < NUM_BOUNDS; i++) {
    bound_x = bound_x + (NUM_BOUNDS - i)*(check_bound_xr(s_idx.x,i+1) - check_bound_l(s_idx.x,i+1));
    bound_y = bound_y + (NUM_BOUNDS - i)*(check_bound_yr(s_idx.y,i+1) - check_bound_l(s_idx.y,i+1));
    bound_z = bound_z + (NUM_BOUNDS - i)*(check_bound_zr(s_idx.z,i+1) - check_bound_l(s_idx.z,i+1));
  }

  REAL fac = ((REAL)DX / M_INV[NUM_BOUNDS + bound_x])
           * ((REAL)DY / M_INV[NUM_BOUNDS + bound_y])
           * ((REAL)DZ / M_INV[NUM_BOUNDS + bound_z]);


  partial_n[lid] = fac * d_field[gid]*d_field[gid];

  barrier(CLK_LOCAL_MEM_FENCE);

  for (uint i = group_size/2; i > 0; i /= 2) {
    if(lid < i) {
      partial_n[lid] += partial_n[lid + i];
    }
    barrier(CLK_LOCAL_MEM_FENCE);
  }

  if (lid == 0) {
    #ifdef REDUCE
      atomic_add_global(&(output[0]), partial_n[0]);
    #else
      output[wgid] = partial_n[0];
    #endif
  }
}


//--------------------------------------------------------------------------------------------------
// Derivative Operators
//--------------------------------------------------------------------------------------------------

// Computes the divergence of the vector field `d_field` and stores it in `d_field_div`
kernel void calc_div(global REAL4 *d_field, global REAL *d_field_div) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  // Calc divergence
  d_field_div[idx] = div(ix, iy, iz, d_field);
}

// Computes the curl of the vector field `d_field` and stores it in `d_field_curl`
kernel void calc_curl(global REAL4 *d_field, global REAL4 *d_field_curl) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_curl[idx] = curl(ix, iy, iz, d_field);
}

// Computes the gradient of the scalar field `d_field` and stores it in `d_field_grad`
kernel void calc_grad(global REAL *d_field, global REAL4 *d_field_grad) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_grad[idx] = grad(ix, iy, iz, d_field);
}

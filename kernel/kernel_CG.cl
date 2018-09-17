//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

// Contains kernel specifically needed for the conjugate gradient method


// Compute x = a x + y
// Note: There are some problems (with MatCL, OpenCL) if a is given as REAL instead of gloobal REAL*
kernel void calc_axpy(global REAL *d_x, global REAL *d_y, global REAL *a) {

   uint ix = get_global_id(0) + BNODES;
   uint iy = get_global_id(1) + BNODES;
   uint iz = get_global_id(2) + BNODES;

   uint idx = calc_idx(ix,iy,iz);

   d_x[idx] = (*a) * d_x[idx] + d_y[idx];
}

// Compute x = x + a y
// Note: There are some problems (with MatCL, OpenCL) if a is given as REAL instead of gloobal REAL*
kernel void calc_xpay(global REAL *d_x, global REAL *d_y, global REAL *a) {

   uint ix = get_global_id(0) + BNODES;
   uint iy = get_global_id(1) + BNODES;
   uint iz = get_global_id(2) + BNODES;

   uint idx = calc_idx(ix,iy,iz);

   d_x[idx] = d_x[idx] + (*a) * d_y[idx];
}

//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

kernel void init_b(global REAL4 *d_field_b) {

   uint ix = get_global_id(0);
   uint iy = get_global_id(1);
   uint iz = get_global_id(2);

   uint idx = calc_idx(ix,iy,iz);

   d_field_b[idx] = b_init(ix, iy, iz);
}


kernel void init_u(global REAL4 *d_field_u) {

   uint ix = get_global_id(0);
   uint iy = get_global_id(1);
   uint iz = get_global_id(2);

   uint idx = calc_idx(ix,iy,iz);

   d_field_u[idx] = u_init(ix, iy, iz);
}


kernel void init_rho(global REAL *d_field_rho) {

   uint ix = get_global_id(0);
   uint iy = get_global_id(1);
   uint iz = get_global_id(2);

   uint idx = calc_idx(ix,iy,iz);

   d_field_rho[idx] = rho_init(ix, iy, iz);
}


kernel void analytical_b(global REAL4 *d_field_b, global REAL *time) {

   uint ix = get_global_id(0);
   uint iy = get_global_id(1);
   uint iz = get_global_id(2);

   uint idx = calc_idx(ix,iy,iz);

   d_field_b[idx] = b_analytical(ix, iy, iz, *time);
}

//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

//Contains kernel and auxillary kernel for time integration. Currently implemented are:
//SSPRK33, SSPRK43, SSPRK93, SSPRK104, KennedyCarpenterLewis2R54C, CalvoFrancoRandez2R64, CarpenterKennedy2N54, ToulorgeDesmet2N84F


// Auxillary function to calculate the approximate time derivative of `d_field_b`
// Adaptive dissipation and the Hall-Term can be enabled by setting the defines USE_ARTIFICIAL_DISSIPATION and USE_HALL respectively
inline REAL4 DB(uint ix, uint iy, uint iz, global REAL4 *d_field_b, global REAL4 *d_field_curlB_rho,
                global REAL4 *d_field_u, REAL time) {

  REAL4 db_dt = (REAL4) {0, 0, 0, 0};

  // Computes the approximate time derivative for the linear induction equation at interior nodes
  db_dt = convec_volume(ix, iy, iz, d_field_b, d_field_u);

  //TODO: merge volume terms of Hall and convec terms to improve performance?
  #ifdef USE_HALL
    db_dt = db_dt + convec_Hall_surface(ix, iy, iz, d_field_b, d_field_curlB_rho, d_field_u, time)
                  + Hall_volume(ix, iy, iz, d_field_b, d_field_curlB_rho);
  #else
    // Computes the approximate time derivative for the linear induction equation at boundary nodes
    db_dt = db_dt + convec_surface(ix, iy, iz, d_field_b, d_field_u, time);
  #endif

  #ifdef USE_ARTIFICIAL_DISSIPATION
    db_dt = db_dt + artificial_dissipation(ix, iy, iz, d_field_b, d_field_u);
  #endif

  return db_dt;
}


//--------------------------------------------------------------------------------------------------
// Time update kernels
//--------------------------------------------------------------------------------------------------

/*
Update the global time variable as time -> time + DT. This kernel takes two arguments d_field_bi
and can be used with CarpenterKennedy2N54.
*/
kernel void calc_time_2_args(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  if (idx == 0) {
	  // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of the Runge-Kutta method
	  time[0] = time[0] + (REAL)(DT);
    time[1] = time[0];
	}
}

/*
Update the global time variable as time -> time + DT. This kernel takes three arguments d_field_bi
and can be used with SSPRK33.
*/
kernel void calc_time_3_args(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
	  time[0] = time[0] + (REAL)(DT);
    time[1] = time[0];
	}
}


//--------------------------------------------------------------------------------------------------
// Velocity update kernels
//--------------------------------------------------------------------------------------------------

/*
Update the ion velocity `d_field_u`. This kernel takes three arguments d_field_bi and can be used
with CarpenterKennedy2N54.
*/
kernel void calc_u_2_args(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  // time[1] is the local time at one stage of a Runge-Kutta method
  d_field_u[idx] = u_analytical(ix, iy, iz, time[1]);
}

/*
Update the ion velocity `d_field_u`. This kernel takes three arguments d_field_bi and can be used
with SSPRK33.
*/
kernel void calc_u_3_args(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  // time[1] is the local time at one stage of a Runge-Kutta method
  d_field_u[idx] = u_analytical(ix, iy, iz, time[1]);
}


//--------------------------------------------------------------------------------------------------
// curl B / rho update kernels
//--------------------------------------------------------------------------------------------------

/*
Compute `curl B / rho` using `d_field_bj`. This kernel takes two arguments d_field_bi and can be
useed with CarpenterKennedy2N54.
*/
kernel void calc_curlB_rho_1_2args(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_curlB_rho[idx] = curl(ix, iy, iz, d_field_b1) / d_field_rho[idx];
}

kernel void calc_curlB_rho_2_2args(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_curlB_rho[idx] = curl(ix, iy, iz, d_field_b2) / d_field_rho[idx];
}

/*
Compute `curl B / rho` using `d_field_bj`. This kernel takes three arguments d_field_bi and can be
useed with SSPRK33.
*/
kernel void calc_curlB_rho_1_3args(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_curlB_rho[idx] = curl(ix, iy, iz, d_field_b1) / d_field_rho[idx];
}

kernel void calc_curlB_rho_2_3args(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_curlB_rho[idx] = curl(ix, iy, iz, d_field_b2) / d_field_rho[idx];
}

kernel void calc_curlB_rho_3_3args(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_curlB_rho[idx] = curl(ix, iy, iz, d_field_b3) / d_field_rho[idx];
}



//--------------------------------------------------------------------------------------------------
// Time integrator kernels using 3 fields
//--------------------------------------------------------------------------------------------------


/*
  SSPRK33_X

Perform one time step using the three stage, third order, strong stability preserving explicit
Runge-Kutta method SSPRK33. Starting with the old value `b1`, the new value after one time step is
obtained using the temporary arrays `b2` and `b3` via
*/
kernel void SSPRK33_1(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b2[idx] =
    d_field_b1[idx] + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}

kernel void SSPRK33_2a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + (REAL)(DT);
	}
}
kernel void SSPRK33_2b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b3[idx] =
    (REAL)(0.75) * d_field_b1[idx]
    + (REAL)(0.25) * (d_field_b2[idx] + (REAL)(DT) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]));
}

kernel void SSPRK33_3a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

	if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + (REAL)(0.5*DT);
	}
}
kernel void SSPRK33_3b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b1[idx] =
    (REAL)(1.0/3.0) * d_field_b1[idx]
    + (REAL)(2.0/3.0) * (d_field_b3[idx] + (REAL)(DT) * DB(ix, iy, iz, d_field_b3, d_field_curlB_rho, d_field_u, time[1]));
}


/*
  SSPRK104_X

Perform one time step using the ten stage, fourth order, strong stability preserving explicit
Runge-Kutta method SSPRK104.

Note: SSPRK104_06 is just a linear combination step and does not require the computation of a time
derivative.
*/
kernel void SSPRK104_01(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b2[idx] =
    d_field_b1[idx] + (REAL)(DT/6.0) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}

kernel void SSPRK104_02a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + (REAL)(DT/6.0);
	}
}
kernel void SSPRK104_02b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b3[idx] =
    d_field_b2[idx] + (REAL)(DT/6.0) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]);
}

kernel void SSPRK104_03a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

	if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + (REAL)(DT/3.0);
	}
}
kernel void SSPRK104_03b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b2[idx] =
    d_field_b3[idx] + (REAL)(DT/6.0) * DB(ix, iy, iz, d_field_b3, d_field_curlB_rho, d_field_u, time[1]);
}

kernel void SSPRK104_04a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

	if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + (REAL)(0.5*DT);
	}
}
kernel void SSPRK104_04b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b3[idx] =
    d_field_b2[idx] + (REAL)(DT/6.0) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]);
}

kernel void SSPRK104_05a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

	if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + (REAL)(DT*2.0/3.0);
	}
}
kernel void SSPRK104_05b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b2[idx] =
    d_field_b3[idx] + (REAL)(DT/6.0) * DB(ix, iy, iz, d_field_b3, d_field_curlB_rho, d_field_u, time[1]);
}

kernel void SSPRK104_06(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b3[idx] = (REAL)(0.04) * d_field_b1[idx] + (REAL)(0.36) * d_field_b2[idx];
  d_field_b2[idx] = (REAL)(0.60) * d_field_b1[idx] + (REAL)(0.40) * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + (REAL)(DT/3.0);
	}
}

kernel void SSPRK104_07(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b1[idx] =
    d_field_b2[idx] + (REAL)(DT/6.0) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]);
}

kernel void SSPRK104_08a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

	if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + (REAL)(0.5*DT);
	}
}
kernel void SSPRK104_08b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b2[idx] =
    d_field_b1[idx] + (REAL)(DT/6.0) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}

kernel void SSPRK104_09a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

	if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + (REAL)(DT*2.0/3.0);
	}
}
kernel void SSPRK104_09b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b1[idx] =
    d_field_b2[idx] + (REAL)(DT/6.0) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]);
}

kernel void SSPRK104_10a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

	if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + (REAL)(DT*5.0/6.0);
	}
}
kernel void SSPRK104_10b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b2[idx] =
    d_field_b1[idx] + (REAL)(DT/6.0) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}

kernel void SSPRK104_11a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

	if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + (REAL)(DT);
	}
}
kernel void SSPRK104_11b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_b3,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b1[idx] =
    d_field_b3[idx]
    + (REAL)(0.6) * (d_field_b2[idx] + (REAL)(DT/6.0) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]));
}


/*
  KennedyCarpenterLewis2R54C_X

Perform one time step using the five stage, fourth order, explicit Runge-Kutta method
KennedyCarpenterLewis2R54C.
*/
kernel void KennedyCarpenterLewis2R54C_1a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_db[idx] = (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void KennedyCarpenterLewis2R54C_1b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(0.22502245872571303);
  REAL coef_B = (REAL)(-0.17379315208537388);
  REAL coef_C = (REAL)(0.22502245872571303);

  d_field_b1[idx] = d_field_b1[idx] + coef_A * d_field_db[idx];
  d_field_b2[idx] = d_field_b1[idx] + coef_B * d_field_db[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void KennedyCarpenterLewis2R54C_2a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_db[idx] = (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void KennedyCarpenterLewis2R54C_2b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(0.5440433129514047);
  REAL coef_B = (REAL)(-0.1630884872250029);
  REAL coef_C = (REAL)(0.5952726195917439);

  d_field_b2[idx] = d_field_b2[idx] + coef_A * d_field_db[idx];
  d_field_b1[idx] = d_field_b2[idx] + coef_B * d_field_db[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void KennedyCarpenterLewis2R54C_3a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_db[idx] = (REAL)(DT) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void KennedyCarpenterLewis2R54C_3b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(0.14456824349399464);
  REAL coef_B = (REAL)(-0.5179208398863779);
  REAL coef_C = (REAL)(0.5767523758607357);

  d_field_b1[idx] = d_field_b1[idx] + coef_A * d_field_db[idx];
  d_field_b2[idx] = d_field_b1[idx] + coef_B * d_field_db[idx];
}

kernel void KennedyCarpenterLewis2R54C_4a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_db[idx] = (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void KennedyCarpenterLewis2R54C_4b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(0.7866643421983568);
  REAL coef_B = (REAL)(-0.19416305717199442);
  REAL coef_C = (REAL)(0.8454958781727144);

  d_field_b2[idx] = d_field_b2[idx] + coef_A * d_field_db[idx];
  d_field_b1[idx] = d_field_b2[idx] + coef_B * d_field_db[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void KennedyCarpenterLewis2R54C_5(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(0.34866717899927996);

  d_field_b1[idx] =
    d_field_b1[idx]
    + coef_B * (REAL)(DT) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]);
}


/*
    CalvoFrancoRandez2R64_X

Perform one time step using the six stage, fourth order, explicit Runge-Kutta method
CalvoFrancoRandez2R64.
*/
kernel void CalvoFrancoRandez2R64_1a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_db[idx] = (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void CalvoFrancoRandez2R64_1b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(0.17985400977138);
  REAL coef_B = (REAL)(0.10893125722541);
  REAL coef_C = (REAL)(0.28878526699679);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_db[idx];
  d_field_b2[idx] = d_field_b1[idx] + coef_A * d_field_db[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CalvoFrancoRandez2R64_2a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_db[idx] = (REAL)(DT) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void CalvoFrancoRandez2R64_2b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(0.14081893152111);
  REAL coef_B = (REAL)(0.13201701492152);
  REAL coef_C = (REAL)(0.38176720366804);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_db[idx];
  d_field_b2[idx] = d_field_b1[idx] + coef_A * d_field_db[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CalvoFrancoRandez2R64_3a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_db[idx] = (REAL)(DT) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void CalvoFrancoRandez2R64_3b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(0.08255631629428);
  REAL coef_B = (REAL)(0.38911623225517);
  REAL coef_C = (REAL)(0.71262082069639);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_db[idx];
  d_field_b2[idx] = d_field_b1[idx] + coef_A * d_field_db[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CalvoFrancoRandez2R64_4a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_db[idx] = (REAL)(DT) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void CalvoFrancoRandez2R64_4b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(0.65804425034331);
  REAL coef_B = (REAL)(-0.59203884581148);
  REAL coef_C = (REAL)(0.69606990893393);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_db[idx];
  d_field_b2[idx] = d_field_b1[idx] + coef_A * d_field_db[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CalvoFrancoRandez2R64_5a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_db[idx] = (REAL)(DT) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void CalvoFrancoRandez2R64_5b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(0.31862993413251);
  REAL coef_B = (REAL)(0.47385028714844);
  REAL coef_C = (REAL)(0.83050587987157);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_db[idx];
  d_field_b2[idx] = d_field_b1[idx] + coef_A * d_field_db[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CalvoFrancoRandez2R64_6(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2, global REAL4 *d_field_db,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(0.48812405426094);

  d_field_b1[idx] =
    d_field_b1[idx]
    + coef_B * (REAL)(DT) * DB(ix, iy, iz, d_field_b2, d_field_curlB_rho, d_field_u, time[1]);
}



//--------------------------------------------------------------------------------------------------
// Time integrator kernels using 2 fields
//--------------------------------------------------------------------------------------------------


/*
    CarpenterKennedy2N54_X

Perform one time step using the five stage, fourth order, explicit Runge-Kutta method
CarpenterKennedy2N54. Starting with the old value `b1`, the new value after one time step is
obtained using the temporary array `b2`.
*/
kernel void CarpenterKennedy2N54_1a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b2[idx] = (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void CarpenterKennedy2N54_1b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(1432997174477.0 / 9575080441755.0);
  REAL coef_C = (REAL)(1432997174477.0 / 9575080441755.0);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CarpenterKennedy2N54_2a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(-567301805773.0 / 1357537059087.0);

  d_field_b2[idx] = coef_A * d_field_b2[idx]
                    + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void CarpenterKennedy2N54_2b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(5161836677717.0 / 13612068292357.0);
  REAL coef_C = (REAL)(2526269341429.0 / 6820363962896.0);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CarpenterKennedy2N54_3a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(-2404267990393.0 / 2016746695238.0);

  d_field_b2[idx] = coef_A * d_field_b2[idx]
                    + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void CarpenterKennedy2N54_3b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(1720146321549.0 / 2090206949498.0);
  REAL coef_C = (REAL)(2006345519317.0 / 3224310063776.0);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CarpenterKennedy2N54_4a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(-3550918686646.0 / 2091501179385.0);

  d_field_b2[idx] = coef_A * d_field_b2[idx]
                    + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void CarpenterKennedy2N54_4b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(3134564353537.0 / 4481467310338.0);
  REAL coef_C = (REAL)(2802321613138.0 / 2924317926251.0);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void CarpenterKennedy2N54_5a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(-1275806237668.0 / 842570457699.0);

  d_field_b2[idx] = coef_A * d_field_b2[idx]
                    + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void CarpenterKennedy2N54_5b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(2277821191437.0 / 14882151754819.0);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];
}



/*
    ToulorgeDesmet2N84F

Perform one time step using the five stage, fourth order, explicit Runge-Kutta method
ToulorgeDesmet2N84F.
*/
kernel void ToulorgeDesmet2N84F_1a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  d_field_b2[idx] = (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void ToulorgeDesmet2N84F_1b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(0.08037936882736950);
  REAL coef_C = (REAL)(0.08037936882736950);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_2a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(-0.5534431294501569);

  d_field_b2[idx] = coef_A * d_field_b2[idx]
                    + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void ToulorgeDesmet2N84F_2b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(0.5388497458569843);
  REAL coef_C = (REAL)(0.3210064250338430);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_3a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(0.01065987570203490);

  d_field_b2[idx] = coef_A * d_field_b2[idx]
                    + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void ToulorgeDesmet2N84F_3b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(0.01974974409031960);
  REAL coef_C = (REAL)(0.3408501826604660);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_4a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(-0.5515812888932000);

  d_field_b2[idx] = coef_A * d_field_b2[idx]
                    + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void ToulorgeDesmet2N84F_4b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(0.09911841297339970);
  REAL coef_C = (REAL)(0.3850364824285470);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_5a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(-1.885790377558741);

  d_field_b2[idx] = coef_A * d_field_b2[idx]
                    + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void ToulorgeDesmet2N84F_5b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(0.7466920411064123);
  REAL coef_C = (REAL)(0.5040052477534100);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_6a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(-5.701295742793264);

  d_field_b2[idx] = coef_A * d_field_b2[idx]
                    + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void ToulorgeDesmet2N84F_6b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(1.679584245618894);
  REAL coef_C = (REAL)(0.6578977561168540);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_7a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(2.113903965664793);

  d_field_b2[idx] = coef_A * d_field_b2[idx]
                    + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void ToulorgeDesmet2N84F_7b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(0.2433728067008188);
  REAL coef_C = (REAL)(0.9484087623348481);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];

  if (idx == 0) {
    // time[0] is the time at the beginning of a time step
    // time[1] is the local time at one stage of a Runge-Kutta method
    time[1] = time[0] + coef_C * (REAL)(DT);
	}
}

kernel void ToulorgeDesmet2N84F_8a(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_A = (REAL)(-0.5339578826675280);

  d_field_b2[idx] = coef_A * d_field_b2[idx]
                    + (REAL)(DT) * DB(ix, iy, iz, d_field_b1, d_field_curlB_rho, d_field_u, time[1]);
}
kernel void ToulorgeDesmet2N84F_8b(
    global REAL4 *d_field_b1, global REAL4 *d_field_b2,
    global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u, global REAL *d_field_rho,
    global REAL *time) {

  uint ix = get_global_id(0);
  uint iy = get_global_id(1);
  uint iz = get_global_id(2);

  uint idx = calc_idx(ix,iy,iz);

  REAL coef_B = (REAL)(0.1422730459001373);

  d_field_b1[idx] = d_field_b1[idx] + coef_B * d_field_b2[idx];
}

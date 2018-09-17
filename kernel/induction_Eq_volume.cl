//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

//--------------------------------------------------------------------------------------------------
// Volume terms of the semidiscretisation
//--------------------------------------------------------------------------------------------------

/*
Extended numerical fluxes for '\partial_j(u_i B_j)' in space directions x, y, z.

Note: An aditional argument `uint dir` to access the vector via `Bk[dir]` instead
of `Bk.x` etc. can be used only for OpenCL version 2 and newer.
*/
REAL4 inline ext_num_flux_x_uiBj(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	#ifdef USE_UIBJ_CENTRAL
		return (REAL)(0.5) * (um*Bm.x + uk*Bk.x);
	#elif defined USE_UIBJ_SPLIT
		return (REAL)(0.25) * (Bm.x + Bk.x) * (um + uk);
	#elif defined USE_UIBJ_PRODUCT
		return (REAL)(0.5) * (um*Bk.x + uk*Bm.x);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_i B_j)' specified!"
  #endif
}

REAL4 inline ext_num_flux_y_uiBj(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	#ifdef USE_UIBJ_CENTRAL
		return (REAL)(0.5) * (um*Bm.y + uk*Bk.y);
	#elif defined USE_UIBJ_SPLIT
		return (REAL)(0.25) * (Bm.y + Bk.y) * (um + uk);
	#elif defined USE_UIBJ_PRODUCT
		return (REAL)(0.5) * (um*Bk.y + uk*Bm.y);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_i B_j)' specified!"
  #endif
}

REAL4 inline ext_num_flux_z_uiBj(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	#ifdef USE_UIBJ_CENTRAL
		return (REAL)(0.5) * (um*Bm.z + uk*Bk.z);
	#elif defined USE_UIBJ_SPLIT
		return (REAL)(0.25) * (Bm.z + Bk.z) * (um + uk);
	#elif defined USE_UIBJ_PRODUCT
		return (REAL)(0.5) * (um*Bk.z + uk*Bm.z);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_i B_j)' specified!"
  #endif
}


/*
Extended numerical fluxes for 'u_i \partial_j B_j' in space directions x, y, z.

Note: An aditional argument `uint dir` to access the vector via `Bk[dir]` instead
of `Bk.x` etc. can be used only for OpenCL version 2 and newer.
*/
REAL4 inline ext_num_flux_x_source(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	#ifdef USE_SOURCE_ZERO
		return (REAL4) {0, 0, 0, 0};
	#elif defined USE_SOURCE_CENTRAL
		return (REAL)(0.5) * um * (Bk.x - Bm.x);
	#elif defined USE_SOURCE_SPLIT
		return (REAL)(0.25) * (um + uk) * (Bk.x - Bm.x);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'u_i D_j B_j' specified!"
  #endif
}

REAL4 inline ext_num_flux_y_source(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	#ifdef USE_SOURCE_ZERO
		return (REAL4) {0, 0, 0, 0};
	#elif defined USE_SOURCE_CENTRAL
		return (REAL)(0.5) * um * (Bk.y - Bm.y);
	#elif defined USE_SOURCE_SPLIT
		return (REAL)(0.25) * (um + uk) * (Bk.y - Bm.y);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'u_i D_j B_j' specified!"
  #endif
}

REAL4 inline ext_num_flux_z_source(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	#ifdef USE_SOURCE_ZERO
		return (REAL4) {0, 0, 0, 0};
	#elif defined USE_SOURCE_CENTRAL
		return (REAL)(0.5) * um * (Bk.z - Bm.z);
	#elif defined USE_SOURCE_SPLIT
		return (REAL)(0.25) * (um + uk) * (Bk.z - Bm.z);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'u_i D_j B_j' specified!"
  #endif
}


/*
Extended numerical fluxes for '\partial_j(u_j B_i)' in space directions x, y, z.

Note: An aditional argument `uint dir` to access the vector via `Bk[dir]` instead
of `Bk.x` etc. can be used only for OpenCL version 2 and newer.
*/
REAL4 inline ext_num_flux_x_ujBi(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	#ifdef USE_UJBI_CENTRAL
		return (REAL)(0.5) * (um.x*Bm + uk.x*Bk);
	#elif defined USE_UJBI_SPLIT
		return (REAL)(0.25) * (um.x + uk.x) * (Bm + Bk);
	#elif defined USE_UJBI_PRODUCT
		return (REAL)(0.5) * (um.x*Bk + uk.x*Bm);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_j B_i)' specified!"
  #endif
}

REAL4 inline ext_num_flux_y_ujBi(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	#ifdef USE_UJBI_CENTRAL
		return (REAL)(0.5) * (um.y*Bm + uk.y*Bk);
	#elif defined USE_UJBI_SPLIT
		return (REAL)(0.25) * (um.y + uk.y) * (Bm + Bk);
	#elif defined USE_UJBI_PRODUCT
		return (REAL)(0.5) * (um.y*Bk + uk.y*Bm);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_j B_i)' specified!"
  #endif
}

REAL4 inline ext_num_flux_z_ujBi(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	#ifdef USE_UJBI_CENTRAL
		return (REAL)(0.5) * (um.z*Bm + uk.z*Bk);
	#elif defined USE_UJBI_SPLIT
		return (REAL)(0.25) * (um.z + uk.z) * (Bm + Bk);
	#elif defined USE_UJBI_PRODUCT
		return (REAL)(0.5) * (um.z*Bk + uk.z*Bm);
	#else
    #error "Error in induction_Eq_volume.cl: No discretization of 'D_j(u_j B_i)' specified!"
  #endif
}


/*
Extended numerical fluxes for the linear magnetic induction equation
'\partial_t B_i = \partial_j(u_i B_j - u_j B_i)'.

Note: An aditional argument `uint dir` to access the vector via `Bk[dir]` instead
of `Bk.x` etc. can be used only for OpenCL version 2 and newer.
*/
REAL4 inline ext_num_flux_x(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	return ext_num_flux_x_uiBj(Bm, um, Bk, uk) - ext_num_flux_x_source(Bm, um, Bk, uk) - ext_num_flux_x_ujBi(Bm, um, Bk, uk);
}

REAL4 inline ext_num_flux_y(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	return ext_num_flux_y_uiBj(Bm, um, Bk, uk) - ext_num_flux_y_source(Bm, um, Bk, uk) - ext_num_flux_y_ujBi(Bm, um, Bk, uk);
}

REAL4 inline ext_num_flux_z(REAL4 Bm, REAL4 um, REAL4 Bk, REAL4 uk) {

	return ext_num_flux_z_uiBj(Bm, um, Bk, uk) - ext_num_flux_z_source(Bm, um, Bk, uk) - ext_num_flux_z_ujBi(Bm, um, Bk, uk);
}


REAL4 inline convec_volume(uint ix, uint iy, uint iz, global REAL4 *d_field_b, global REAL4 *d_field_u) {

	uint idx = calc_idx(ix,iy,iz);

	int bound_x = 0;
	int bound_y = 0;
	int bound_z = 0;

	for (uint i = 0; i < NUM_BOUNDS; i++) {
		 bound_x = bound_x + (NUM_BOUNDS - i)*(check_bound_xr(ix,i+1) - check_bound_l(ix,i+1));
		 bound_y = bound_y + (NUM_BOUNDS - i)*(check_bound_yr(iy,i+1) - check_bound_l(iy,i+1));
		 bound_z = bound_z + (NUM_BOUNDS - i)*(check_bound_zr(iz,i+1) - check_bound_l(iz,i+1));
	}

	REAL4 val = (REAL4) {0, 0, 0, 0};
	REAL4 dB;

	REAL4 Bm = d_field_b[idx];
	REAL4 um = d_field_u[idx];
	REAL4 Bk;
	REAL4 uk;

	dB = (REAL4) {0, 0, 0, 0};
	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		Bk = get_vector_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field_b);
		uk = get_vector_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field_u);

		dB = dB + SBP_diff[NUM_BOUNDS + bound_x][i] * ext_num_flux_x(Bm, um, Bk, uk);
	}
	val = val + (REAL)(2.0/DX) * dB;

	dB = (REAL4) {0, 0, 0, 0};
	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		Bk = get_vector_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field_b);
		uk = get_vector_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field_u);

		dB = dB + SBP_diff[NUM_BOUNDS + bound_y][i] * ext_num_flux_y(Bm, um, Bk, uk);
	}
	val = val + (REAL)(2.0/DY) * dB;

	dB = (REAL4) {0, 0, 0, 0};
	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		Bk = get_vector_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field_b);
		uk = get_vector_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field_u);

		dB = dB + SBP_diff[NUM_BOUNDS + bound_z][i] * ext_num_flux_z(Bm, um, Bk, uk);
	}
	val = val + (REAL)(2.0/DZ) * dB;

	return val;
}


REAL4 inline Hall_volume(uint ix, uint iy, uint iz, global REAL4* d_field_b, global REAL4* d_field_curlB_rho) {

	uint idx = calc_idx(ix,iy,iz);

	int bound_x = 0;
	int bound_y = 0;
	int bound_z = 0;

	for (uint i = 0; i < NUM_BOUNDS; i++) {
		 bound_x = bound_x + (NUM_BOUNDS - i)*(check_bound_xr(ix,i+1) - check_bound_l(ix,i+1));
		 bound_y = bound_y + (NUM_BOUNDS - i)*(check_bound_yr(iy,i+1) - check_bound_l(iy,i+1));
		 bound_z = bound_z + (NUM_BOUNDS - i)*(check_bound_zr(iz,i+1) - check_bound_l(iz,i+1));
	}

	REAL4 val = (REAL4) {0, 0, 0, 0};
	REAL4 dB;
	REAL4 B;
	REAL4 curlB_rho;

	dB = (REAL4) {0, 0, 0, 0};
	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		B = get_vector_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field_b);
		curlB_rho = get_vector_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field_curlB_rho);

		dB = dB + SBP_diff[NUM_BOUNDS + bound_x][i] * (curlB_rho.x*B - curlB_rho*B.x);
	}
	val = val + (REAL)(1.0/DX) * dB;

	dB = (REAL4) {0, 0, 0, 0};
	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		B = get_vector_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field_b);
		curlB_rho = get_vector_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field_curlB_rho);

		dB = dB + SBP_diff[NUM_BOUNDS + bound_y][i] * (curlB_rho.y*B - curlB_rho*B.y);
	}
	val = val + (REAL)(1.0/DY) * dB;

	dB = (REAL4) {0, 0, 0, 0};
	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		B = get_vector_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field_b);
		curlB_rho = get_vector_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field_curlB_rho);

		dB = dB + SBP_diff[NUM_BOUNDS + bound_z][i] * (curlB_rho.z*B - curlB_rho*B.z);
	}
	val = val + (REAL)(1.0/DZ) * dB;

	return val;
}

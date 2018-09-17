//--------------------------------------------------------------------------------------------------
// Volume terms of the semidiscretisation
//--------------------------------------------------------------------------------------------------

// TODO: Delete?
/*
REAL4 inline diff_x(uint ix, uint iy, uint iz, global REAL4 *d_field) {

	REAL4 val = (REAL4) {0.0,0.0,0.0,0.0};

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		 bound = bound + (NUM_BOUNDS - i)*(check_bound_xr(ix,i+1) - check_bound_l(ix,i+1));
	}

	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		val = val + SBP_diff[NUM_BOUNDS + bound][i]*get_vector_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field);
	}

	return val/((REAL)DX);

}


REAL4 inline diff_y(uint ix, uint iy, uint iz, global REAL4 *d_field) {

	REAL4 val = (REAL4) {0.0,0.0,0.0,0.0};

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		 bound = bound + (NUM_BOUNDS - i)*(check_bound_yr(iy,i+1) - check_bound_l(iy,i+1));
	}

	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		val = val + SBP_diff[NUM_BOUNDS + bound][i]*get_vector_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field);
	}

	return val/((REAL)DY);

}


REAL4 inline diff_z(uint ix, uint iy, uint iz, global REAL4 *d_field) {

	REAL4 val = (REAL4) {0.0,0.0,0.0,0.0};

	int bound = 0;
	for (uint i = 0; i < NUM_BOUNDS; i++) {
		 bound = bound + (NUM_BOUNDS - i)*(check_bound_zr(iz,i+1) - check_bound_l(iz,i+1));
	}

	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		val = val + SBP_diff[NUM_BOUNDS + bound][i]*get_vector_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field);
	}

	return val/((REAL)DZ);

}

REAL4 inline flux_difference_split(uint ix, uint iy, uint iz, global REAL4 *d_field_b, global REAL4 *d_field_u) {

	REAL4 val = (REAL4) {0.0,0.0,0.0,0.0};

    int bound_x = 0;
	int bound_y = 0;
	int bound_z = 0;

	for (uint i = 0; i < NUM_BOUNDS; i++) {
		 bound_x = bound_x + (NUM_BOUNDS - i)*(check_bound_xr(ix,i+1) - check_bound_l(ix,i+1));
		 bound_y = bound_y + (NUM_BOUNDS - i)*(check_bound_yr(iy,i+1) - check_bound_l(iy,i+1));
		 bound_z = bound_z + (NUM_BOUNDS - i)*(check_bound_zr(iz,i+1) - check_bound_l(iz,i+1));
	}
	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		val = val + SBP_diff[NUM_BOUNDS + bound_x][i]*((get_vector_field(ix,iy,iz,0,0,0,d_field_u).x + get_vector_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field_u).x)*(get_vector_field(ix,iy,iz,0,0,0,d_field_b) + get_vector_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field_b)))/2.0f/((REAL)DX)
		 		  + SBP_diff[NUM_BOUNDS + bound_y][i]*((get_vector_field(ix,iy,iz,0,0,0,d_field_u).y + get_vector_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field_u).y)*(get_vector_field(ix,iy,iz,0,0,0,d_field_b) + get_vector_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field_b)))/2.0f/((REAL)DY)
				  + SBP_diff[NUM_BOUNDS + bound_z][i]*((get_vector_field(ix,iy,iz,0,0,0,d_field_u).z + get_vector_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field_u).z)*(get_vector_field(ix,iy,iz,0,0,0,d_field_b) + get_vector_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field_b)))/2.0f/((REAL)DZ);
	}
	return val;

}

REAL4 inline flux_difference_conservative(uint ix, uint iy, uint iz, global REAL4 *d_field_b, global REAL4 *d_field_u) {

	REAL4 val = (REAL4) {0.0,0.0,0.0,0.0};

    int bound_x = 0;
	int bound_y = 0;
	int bound_z = 0;

	for (uint i = 0; i < NUM_BOUNDS; i++) {
		 bound_x = bound_x + (NUM_BOUNDS - i)*(check_bound_xr(ix,i+1) - check_bound_l(ix,i+1));
		 bound_y = bound_y + (NUM_BOUNDS - i)*(check_bound_yr(iy,i+1) - check_bound_l(iy,i+1));
		 bound_z = bound_z + (NUM_BOUNDS - i)*(check_bound_zr(iz,i+1) - check_bound_l(iz,i+1));
	}
	for (uint i = 0; i < STENCIL_WIDTH; i++) {
		val = val + SBP_diff[NUM_BOUNDS + bound_x][i]*(get_vector_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field_u).x*get_vector_field(ix,iy,iz,(i - (STENCIL_WIDTH - 1)/2),0,0,d_field_b) + get_vector_field(ix,iy,iz,0,0,0,d_field_u).x*get_vector_field(ix,iy,iz,0,0,0,d_field_b))/((REAL)DX)
		 		  + SBP_diff[NUM_BOUNDS + bound_y][i]*(get_vector_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field_u).y*get_vector_field(ix,iy,iz,0,(i - (STENCIL_WIDTH - 1)/2),0,d_field_b) + get_vector_field(ix,iy,iz,0,0,0,d_field_u).y*get_vector_field(ix,iy,iz,0,0,0,d_field_b))/((REAL)DY)
				  + SBP_diff[NUM_BOUNDS + bound_z][i]*(get_vector_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field_u).z*get_vector_field(ix,iy,iz,0,0,(i - (STENCIL_WIDTH - 1)/2),d_field_b) + get_vector_field(ix,iy,iz,0,0,0,d_field_u).z*get_vector_field(ix,iy,iz,0,0,0,d_field_b))/((REAL)DZ);
	}

	return val;
}

REAL4 inline flux_difference_nonconservative(uint ix, uint iy, uint iz, global REAL4 *d_field_b, global REAL4 *d_field_u) {

	REAL4 val = (REAL4) {0.0,0.0,0.0,0.0};

    val = (REAL)(2)*flux_difference_split(ix,iy,iz,d_field_b,d_field_u) - flux_difference_conservative(ix,iy,iz,d_field_b,d_field_u);

	return val;
}

REAL4 inline convec(uint ix, uint iy, uint iz, global REAL4 *d_field_b, global REAL4 *d_field_u) {
	uint idx = calc_idx(ix,iy,iz);
	REAL4 val = (REAL4) {0.0,0.0,0.0,0.0};

	val.xyz = d_field_b[idx].x*diff_x(ix,iy,iz,d_field_u).xyz
	    	+ d_field_b[idx].y*diff_y(ix,iy,iz,d_field_u).xyz
			+ d_field_b[idx].z*diff_z(ix,iy,iz,d_field_u).xyz;

	#ifdef USE_SPLIT

		val.xyz = val.xyz - flux_difference_split(ix,iy,iz,d_field_b,d_field_u).xyz;

    #elif defined USE_CONSERVATIVE

		val.xyz = val.xyz - flux_difference_conservative(ix,iy,iz,d_field_b,d_field_u).xyz;

    #elif defined USE_NONCONSERVATIVE

		val.xyz = val.xyz - flux_difference_nonconservative(ix,iy,iz,d_field_b,d_field_u).xyz;

    #else
        #error "Error in kernel_time_integrator.cl: NO DISCRETIZATION SPECIFIED!"
    #endif

	return val;
}
*/

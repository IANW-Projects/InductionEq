//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

//--------------------------------------------------------------------------------------------------
// Surface terms of the semidiscretisation: Simultaneous spproximation terms (SATs)
//--------------------------------------------------------------------------------------------------

inline REAL4 convec_surface(uint ix, uint iy, uint iz, global REAL4 *d_field_b,
                            global REAL4 *d_field_u, REAL current_time) {

	uint idx =  calc_idx(ix,iy,iz);

	REAL4 val;
	REAL4 b_bound = b_boundary(ix, iy, iz, current_time);

	val = (REAL)(M_INV[0]/DX) * (
			  (check_bound_l(ix,1) * (d_field_u[idx].x > 0)) * (REAL)(-1)
			+ (check_bound_xr(ix,1) * (d_field_u[idx].x < 0)) * (REAL)(1)
		) * d_field_u[idx].x * (d_field_b[idx] - b_bound)
		+ (REAL)(M_INV[0]/DY) * (
			  (check_bound_l(iy,1) * (d_field_u[idx].y > 0)) * (REAL)(-1)
			+ (check_bound_yr(iy,1) * (d_field_u[idx].y < 0)) * (REAL)(1)
		) * d_field_u[idx].y * (d_field_b[idx] - b_bound)
		+ (REAL)(M_INV[0]/DZ) * (
			  (check_bound_l(iz,1) * (d_field_u[idx].z > 0)) * (REAL)(-1)
			+ (check_bound_zr(iz,1) * (d_field_u[idx].z < 0)) * (REAL)(1)
		) * d_field_u[idx].z * (d_field_b[idx] - b_bound);

  	return val;
}

inline REAL4 convec_Hall_surface(uint ix, uint iy, uint iz, global REAL4 *d_field_b,
                                 global REAL4 *d_field_curlB_rho, global REAL4 *d_field_u,
                                 REAL current_time) {

	uint idx = calc_idx(ix,iy,iz);
	REAL4 val = (REAL4) {0, 0, 0, 0};

	// TODO: Depending on the boundary, we want to have either inflow or outflow BCs!
	// Up to now, there are only outflow BCs!

	// outflow BC, cf. 'directional do-nothing conditions'
	REAL4 test_velocity = (REAL)(0.5) * d_field_u[idx] - d_field_curlB_rho[idx];

	val = (REAL)(M_INV[0]/DX) * (
			( 	check_bound_l(ix,1) * (test_velocity.x > 0) * (REAL)(-1)
				+ check_bound_xr(ix,1) * (test_velocity.x < 0) * (REAL)(1)
			) * test_velocity.x * d_field_b[idx]
			+
			(	  check_bound_l(ix,1) * (REAL)(-1)
				+ check_bound_xr(ix,1) * (REAL)(1)
			) * d_field_b[idx].x * d_field_curlB_rho[idx]
		)
		+ (REAL)(M_INV[0]/DY) * (
			( 	check_bound_l(iy,1) * (test_velocity.y > 0) * (REAL)(-1)
				+ check_bound_yr(iy,1) * (test_velocity.y < 0) * (REAL)(1)
			) * test_velocity.y * d_field_b[idx]
			+
			(	  check_bound_l(iy,1) * (REAL)(-1)
				+ check_bound_yr(iy,1) * (REAL)(1)
			) * d_field_b[idx].y * d_field_curlB_rho[idx]
		)
		+ (REAL)(M_INV[0]/DZ) * (
			( 	check_bound_l(iz,1) * (test_velocity.z > 0) * (REAL)(-1)
				+ check_bound_zr(iz,1) * (test_velocity.z < 0) * (REAL)(1)
			) * test_velocity.z * d_field_b[idx]
			+
			(	  check_bound_l(iz,1) * (REAL)(-1)
				+ check_bound_zr(iz,1) * (REAL)(1)
			) * d_field_b[idx].z * d_field_curlB_rho[idx]
		);

  	return val;
}

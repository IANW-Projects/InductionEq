//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

REAL constant CONST_B0 = (REAL)(1.0);
REAL constant CONST_Bz = (REAL)(1.0);
REAL constant CONST_ux = (REAL)(1.0);
REAL constant CONST_uy = (REAL)(1.0);
REAL constant CONST_uz = (REAL)(0.1);
REAL constant CONST_k  = (REAL)(1.0);
REAL constant CONST_omega = (REAL)(1.0);

/*
Analytical solution of the magnetic field. It is used to calculate the numerical
error and related quantities.
*/
inline REAL4 b_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

	REAL xi = CONST_k * z - CONST_omega * time;
  REAL Omega = ( (CONST_k * CONST_uz - CONST_omega)*(CONST_k * CONST_uz - CONST_omega)
                  - CONST_k*CONST_k * CONST_Bz*CONST_Bz)
              / (CONST_Bz * (CONST_k*CONST_uz - CONST_omega) * CONST_k*CONST_k);

	return (REAL4) {CONST_B0*cos(Omega*xi), CONST_B0*sin(Omega*xi), CONST_Bz, 0};
}

/*
Initial condition of the magnetic field.
*/
inline REAL4 b_init(uint ix, uint iy, uint iz) {

	return b_analytical(ix, iy, iz, (REAL)0);
}

/*
Boundary condition of the magnetic field.
*/
inline REAL4 b_boundary(uint ix, uint iy, uint iz, REAL time) {

	// TODO: Something else?
	return (REAL4) {0, 0, 0, 0};
}

/*inline REAL4 curlB_rho_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

	REAL xi = CONST_k * z - CONST_omega * time;
  REAL Omega = ( (CONST_k * CONST_uz - CONST_omega)*(CONST_k * CONST_uz - CONST_omega)
                  - CONST_k*CONST_k * CONST_Bz*CONST_Bz)
              / (CONST_Bz * (CONST_k*CONST_uz - CONST_omega) * CONST_k*CONST_k);
	REAL fac = CONST_B0 * CONST_k*CONST_k * (CONST_Bz*CONST_Bz - (CONST_omega-CONST_uz)*(CONST_omega-CONST_uz))
				/ (CONST_Bz * CONST_k * (CONST_k * CONST_uz - CONST_omega));

	return (REAL4) {fac*cos(Omega*xi), fac*sin(Omega*xi), 0, 0};
}*/


inline REAL4 u_analytical(uint ix, uint iy, uint iz, REAL time)
{
	//TODO: Mind box geometry
	REAL x = -1.0 + ix*(REAL)DX;
	REAL y = -1.0 + iy*(REAL)DY;
	REAL z = -1.0 + iz*(REAL)DZ;

  REAL xi = CONST_k * z - CONST_omega * time;
  REAL Omega = ( (CONST_k * CONST_uz - CONST_omega)*(CONST_k * CONST_uz - CONST_omega)
                  - CONST_k*CONST_k * CONST_Bz*CONST_Bz)
              / (CONST_Bz * (CONST_k*CONST_uz - CONST_omega) * CONST_k*CONST_k);
  REAL fac = CONST_k * CONST_Bz / (CONST_k*CONST_uz - Omega);

	return (REAL4) {fac*CONST_B0*cos(Omega*xi) + CONST_ux, fac*CONST_B0*sin(Omega*xi) + CONST_uy, CONST_uz, 0};
}

/*
Initial condition of the ion velocity.
*/
inline REAL4 u_init(uint ix, uint iy, uint iz) {

	return u_analytical(ix, iy, iz, (REAL)(0));
}

/*
Initial condition of the ion density.
*/
inline REAL rho_init(uint ix, uint iy, uint iz) {

	return (REAL)(1);
}

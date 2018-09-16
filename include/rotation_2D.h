//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

/*
Analytical solution of the magnetic field. It is used to calculate the numerical
error and related quantities.
*/
inline REAL4 b_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

	REAL s = sin(time);
  REAL c = cos(time);

  REAL xx = c*x + s*y;
  REAL yy = -s*x + c*y;
  REAL zz = z;

  REAL fac = exp(-20*((xx-(REAL)(0.5))*(xx-(REAL)(0.5)) + yy*yy));
  REAL Bx0 = -4*yy*fac;
  REAL By0 = (4*xx-2)*fac;

	return (REAL4) {c*Bx0-s*By0, s*Bx0+c*By0, (REAL)(0), (REAL)(0)};
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

	return b_analytical(ix, iy, iz, time);
}

/*
Initial condition of the ion velocity.
*/
inline REAL4 u_init(uint ix, uint iy, uint iz) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;

	return (REAL4) {-y, x, (REAL)(0), (REAL)(0)};
}

inline REAL4 u_analytical(uint ix, uint iy, uint iz, REAL time) {

	return u_init(ix, iy, iz);
}

/*
Initial condition of the ion density.
*/
inline REAL rho_init(uint ix, uint iy, uint iz) {

	return (REAL)(1);
}

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
	REAL sqrt3 = sqrt(3.0);

  REAL xx = ( (1 + 2*c) * x + (1 - c + sqrt3*s) * y + (1 - c - sqrt3*s) * z ) / 3.0;
  REAL yy = ( (1 - c - sqrt3*s) * x + (1 + 2*c) * y + (1 - c + sqrt3*s) * z ) / 3.0;
  REAL zz = ( (1 - c + sqrt3*s) * x + (1 - c - sqrt3*s) * y + (1 + 2*c) * z ) / 3.0;

	REAL alpha = exp(-(REAL)(5.0/3.0) * (3 - 2*(3+sqrt3)*xx + 12*xx*xx
                                         - 2*(-3+sqrt3)*yy + 12*yy*yy
                                         + 4*sqrt3*zz + 12*zz*zz));

	REAL b0x = alpha * (3 - sqrt3 - 4*sqrt3*yy + 4*sqrt3*zz) / 48.0;
	REAL b0y = alpha * (-3 - sqrt3 + 4*sqrt3*xx - 4*sqrt3*zz) / 48.0;
	REAL b0z = alpha * (1 - 2*xx + 2*yy) / (8*sqrt3);

	REAL bx = ( (1 + 2*c) * b0x + (1 - c - sqrt3*s) * b0y + (1 - c + sqrt3*s) * b0z ) / 3.0;
  REAL by = ( (1 - c + sqrt3*s) * b0x + (1 + 2*c) * b0y + (1 - c - sqrt3*s) * b0z ) / 3.0;
  REAL bz = ( (1 - c - sqrt3*s) * b0x + (1 - c + sqrt3*s) * b0y + (1 + 2*c) * b0z ) / 3.0;

	return (REAL4) {bx, by, bz, (REAL)0};
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
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

	REAL sqrt3 = sqrt(3.0);

	return (REAL4) {(z-y)/sqrt3, (x-z)/sqrt3, (y-x)/sqrt3, (REAL)(0)};
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

//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

REAL constant CONST_A = 0.1;
REAL constant CONST_alpha = 0.5235987755982988; // pi/6

/*
Analytical solution of the magnetic field. It is used to calculate the numerical
error and related quantities.
*/
inline REAL4 b_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  REAL cos_alpha = cos(CONST_alpha);
  REAL sin_alpha = sin(CONST_alpha);

	REAL phi = 6.283185307179586 * (x*cos_alpha + y*sin_alpha + time);
  REAL cos_phi = cos(phi);
  REAL sin_phi = sin(phi);

  REAL B1 = cos_alpha - CONST_A * sin_alpha * sin_phi;
  REAL B2 = sin_alpha - CONST_A * cos_alpha * sin_phi;
  REAL B3 = CONST_A * cos_phi;

	return (REAL4) {B1, B2, B3, (REAL)(0)};
}

/*
Initial condition of the magnetic field.
*/
inline REAL4 b_init(uint ix, uint iy, uint iz) {

	return b_analytical(ix, iy, iz, (REAL)(0));
}

/*
Boundary condition of the magnetic field.
*/
inline REAL4 b_boundary(uint ix, uint iy, uint iz, REAL time) {

	return b_analytical(ix, iy, iz, time);
}

/*
Ion velocity.
*/
inline REAL4 u_analytical(uint ix, uint iy, uint iz, REAL time) {

  REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  REAL cos_alpha = cos(CONST_alpha);
  REAL sin_alpha = sin(CONST_alpha);

	REAL phi = 6.283185307179586 * (x*cos_alpha + y*sin_alpha + time);
  REAL cos_phi = cos(phi);
  REAL sin_phi = sin(phi);

  REAL u1 = - CONST_A * sin_alpha * sin_phi;
  REAL u2 = CONST_A * cos_alpha * sin_phi;
  REAL u3 = CONST_A * cos_phi;

	return (REAL4) {u1, u2, u3, (REAL)(0)};
}

inline REAL4 u_init(uint ix, uint iy, uint iz) {

	return u_analytical(ix, iy, iz, (REAL)(0));
}

/*
Initial condition of the ion density.
*/
inline REAL rho_init(uint ix, uint iy, uint iz) {

	return (REAL)(1);
}

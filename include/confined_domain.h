//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

/*
Analytical solution of the magnetic field. It is used to calculate the numerical
error and related quantities.
*/
inline REAL4 b_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  REAL sx = sinpi(x);
  REAL cx = cospi(x);
  REAL sy = sinpi(y);
  REAL cy = cospi(y);
  REAL sz = sinpi(z);
  REAL cz = cospi(z);

	return (REAL4) {sx*cy*cz,
                  cx*sy*cz,
                  -2*cx*cy*sz,
                  (REAL)(0)};
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
Initial condition of the ion velocity.
*/
inline REAL4 u_init(uint ix, uint iy, uint iz) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

  REAL sx = sinpi(x);
  REAL cx = cospi(x);
  REAL sy = sinpi(y);
  REAL cy = cospi(y);
  REAL sz = sinpi(z);
  REAL cz = cospi(z);

	return (REAL4) {sx*cy*cz,
                  cx*sy*cz,
                  -2*cx*cy*sz,
                  (REAL)(0)};
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

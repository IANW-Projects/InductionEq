//This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

REAL constant CONST_n1 = (REAL)(0.5773502691896258); // =1/sqrt(3)
REAL constant CONST_n2 = (REAL)(0.5773502691896258);
REAL constant CONST_n3 = (REAL)(0.5773502691896258);
REAL constant CONST_a = (REAL)(1.0);
REAL constant CONST_b = (REAL)(1.0);
REAL constant CONST_c = (REAL)(1.0);
REAL constant CONST_alpha = (REAL)(0.5);

/*
Analytical solution of the magnetic field. It is used to calculate the numerical
error and related quantities.
*/
inline REAL4 b_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

	REAL CONST_k = (REAL)(1.0 - CONST_alpha * CONST_alpha) / CONST_alpha;

	REAL u1 = CONST_a * cos(CONST_k * y + CONST_alpha * CONST_k * time * CONST_n2)
          + CONST_b * sin(CONST_k * z + CONST_alpha * CONST_k * time * CONST_n3);
	REAL u2 = CONST_b * cos(CONST_k * z + CONST_alpha * CONST_k * time * CONST_n3)
          + CONST_c * sin(CONST_k * x + CONST_alpha * CONST_k * time * CONST_n1);
	REAL u3 = CONST_c * cos(CONST_k * x + CONST_alpha * CONST_k * time * CONST_n1)
          + CONST_a * sin(CONST_k * y + CONST_alpha * CONST_k * time * CONST_n2);

	return (REAL4) {CONST_alpha * u1 + CONST_n1,
                  CONST_alpha * u2 + CONST_n2,
                  CONST_alpha * u3 + CONST_n3,
                  0};
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

  // TODO: Something else?
	return (REAL4) {0, 0, 0, 0};
}

/*inline REAL4 curlB_rho_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

	REAL CONST_k = (REAL)(1.0 - CONST_alpha * CONST_alpha) / CONST_alpha;

	REAL u1 = CONST_a * cos(CONST_k * y + CONST_alpha * CONST_k * time * CONST_n2)
          + CONST_b * sin(CONST_k * z + CONST_alpha * CONST_k * time * CONST_n3);
	REAL u2 = CONST_b * cos(CONST_k * z + CONST_alpha * CONST_k * time * CONST_n3)
          + CONST_c * sin(CONST_k * x + CONST_alpha * CONST_k * time * CONST_n1);
	REAL u3 = CONST_c * cos(CONST_k * x + CONST_alpha * CONST_k * time * CONST_n1)
          + CONST_a * sin(CONST_k * y + CONST_alpha * CONST_k * time * CONST_n2);

	return (REAL4) {CONST_alpha * CONST_k * u1,
                  CONST_alpha * CONST_k * u2,
                  CONST_alpha * CONST_k * u3,
                  0};
}*/


inline REAL4 u_analytical(uint ix, uint iy, uint iz, REAL time) {

	REAL x = (REAL)XMIN + ix*(REAL)DX;
	REAL y = (REAL)YMIN + iy*(REAL)DY;
	REAL z = (REAL)ZMIN + iz*(REAL)DZ;

	REAL CONST_k = (REAL)(1.0 - CONST_alpha * CONST_alpha) / CONST_alpha;

	REAL u1 = CONST_a * cos(CONST_k * y + CONST_alpha * CONST_k * time * CONST_n2)
          + CONST_b * sin(CONST_k * z + CONST_alpha * CONST_k * time * CONST_n3);
	REAL u2 = CONST_b * cos(CONST_k * z + CONST_alpha * CONST_k * time * CONST_n3)
          + CONST_c * sin(CONST_k * x + CONST_alpha * CONST_k * time * CONST_n1);
	REAL u3 = CONST_c * cos(CONST_k * x + CONST_alpha * CONST_k * time * CONST_n1)
          + CONST_a * sin(CONST_k * y + CONST_alpha * CONST_k * time * CONST_n2);

	return (REAL4){u1, u2, u3, (REAL)(0)};
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

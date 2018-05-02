/*
	Andrea Landella

	header file rk4solver
*/
#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>

//RK4 solver implementation for scalar systems
void rungekutta4(double		odesys(double t, double y), 
				 double		t0, 
				 double		tf, 
				 double		y0, 
				 double		h, 
				 double*	&tSOL, 
				 double*	&ySOL
				 );

//RK4 solver implementation for vector systems
void rungekutta4_vec(void		odesys_vec(double t, double* y, int DIM, double dydt[]), 
					 double		t0, 
					 double		tf, 
					 double*	y0, 
					 double		h, 
					 int		DIM, 
					 double*	&tSOL, 
					 double**	&ySOL
					 );
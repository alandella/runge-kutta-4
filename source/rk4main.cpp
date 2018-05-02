/*
	Andrea Landella

	program file rk4main

	We want to solve the 1-D and DIM-D ODE equation:

	y' = F(t,y)
	y(t0) = y0, integrate up to tf with spacing h.

	with the Runge-Kutta 4th order method. The tf and h values are specified.

	In order to do that, we need two functions: one for the system, F, and the
	other for the solver, which takes as input the system function's handle.

*/
#include"rk4solver.h"

//ODE system function
double odesysfun(double t, double y);

//ODE vector system function
void odesysfun_vec(double t, double* y, int DIM, double dydt[]);

int main() {

	/******* SCALAR SYSTEM *******/

	//initial conditions
	double y0 = 1.;

	//time-step conditions
	double t0 = 0.;
	double tf = 5.;

	//integration step
	double h = 0.001;

	//number of points
	int N = (int)((tf - t0) / h + 0.5);
	//
	printf("Integrating on %d\tpoints... ", N);

	//we declare the solution vectors
	double* tSOL = (double*)malloc(N*sizeof(double)); if (!tSOL) return NULL;
	double* ySOL = (double*)malloc(N*sizeof(double)); if (!ySOL) return NULL;
	//main function call
	rungekutta4(&odesysfun, t0, tf, y0, h, tSOL, ySOL);
	//
	free(tSOL);
	free(ySOL); 

	//exitflags
	printf("Done 1!\n");

	/******* VECTOR SYSTEM *******/

	//problem dimension
	const int DIM = 4;

	//initial conditions
	double y0_vec[DIM] = { 0.3, 1.6, 0.9, 1.3 };

	//time-step conditions
	double t0_vec = 0.;
	double tf_vec = 5.;

	//integration step (all normalized components of the h vector)
	double h_vec = h / pow(DIM, 0.5);

	//number of points
	int N_vec = (int)((tf_vec - t0_vec) / h_vec + 0.5);
	//
	printf("Integrating on %d\tpoints... ", N_vec);

	//we declare the solution vectors
	double* tSOL_vec = (double*)malloc(N_vec * sizeof(double)); if (!tSOL_vec) return NULL;
	double** ySOL_vec = (double**)malloc(N_vec * sizeof(double*)); if (!ySOL_vec) return NULL;
	//
	for (int j = 0; j < N_vec; j++){
		ySOL_vec[j] = (double*)malloc(DIM * sizeof(double)); if (!ySOL_vec[j]) return NULL;
	}
	//main function call
	rungekutta4_vec(&odesysfun_vec, t0_vec, tf_vec, y0_vec, h_vec, DIM, tSOL_vec, ySOL_vec);
	//
	free(tSOL_vec);
	for (int j = 0; j < N_vec; j++){
		free(ySOL_vec[j]);
	}
	free(ySOL_vec);

	//exitflags
	printf("Done 2!\n");
	system("pause");

	return 0;
}

//ODE system function
double odesysfun(double t, double y) {

	return t*sin(y*t);

}

//ODE vector system function
void odesysfun_vec(double t, double* y, int DIM, double dydt[]) {

	for (int j = 0; j < DIM; j++) {
		//attractor for 5, linear
		//dydt[j] = y[DIM - 1 - j] - y[j];

		//attractor for 1, nonlinear
		dydt[j] = y[DIM - 1 - j] - pow(y[j], 2);
	}
}

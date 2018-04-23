/*
	Andrea Landella

	program file flowtankdyn.cpp

We want to solve the 1-D and DIM-D ODE equation:

	y' = F(t,y)
	y(t0) = y0, integrate up to tf with spacing h.

with the Runge-Kutta 4th order method. The tf and h values are specified.

In order to do that, we need two functions: one for the system, F, and the
other for the solver, which takes as input the system function's handle.

*/

#include<stdio.h>
#include<stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include"flowtankdyn.h"

int main(void) {

	/**** 
	SCALAR SYSTEM 
	****/

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

	//main function call
	rungekutta4(y0, t0, tf, h, &odesysfun);

	//exitflags
	printf("Done 1!\n");

	/**** 
	VECTOR SYSTEM 
	****/
	
	//problem dimension
	const int DIM_sys = 4;

	//initial conditions
	double y0_vec[DIM_sys] = { 0.3, 1.6, 0.9, 1.3 };

	//time-step conditions
	t0 = 0.;
	tf = 5.;

	//integration step (all normalized components of the h vector)
	double h_vec = h / pow(DIM_sys, 0.5);

	//number of points
	int N_vec = (int)((tf - t0) / h_vec + 0.5);
	//
	printf("Integrating on %d\tpoints... ", N_vec);

	//main function call
	rungekutta4_vec(y0_vec, t0, tf, h_vec, &odesysfun_vec, DIM_sys);

	//exitflags
	printf("Done 2!\n");
	system("pause");

	return 0;
}

//ODE system function
double odesysfun (double t, double y) {

	double dydt;

	dydt = t*sin(y*t);

	return dydt;

}

//RK4 solver implementation for scalar systems
void rungekutta4(double y0, double t0, double tf, double h, double odesys(double t, double y)) {

	//constants for the RK4
	double k1;
	double k2;
	double k3;
	double k4;

	//number of points
	int N = (int) ((tf - t0) / h + 0.5);

	//mallocation of the arrays with N elements
	double* tris = (double* )malloc(N * sizeof(double));
	double* yris = (double* )malloc(N * sizeof(double));

	//setting results at zero, t = t0, y = y0:
	tris[0] = t0;
	yris[0] = y0;

	//we post process the data 
	char filename[] = "rk4_sys.txt";
	//
	FILE* fp1 = fopen(filename, "w");
	//
	fprintf(fp1, "t\t\ty\n%.5e\t%.5e\n", tris[0], yris[0]);

	//RK4 algorithm
	for (int j = 1; j < N; j++) {

		k1 = odesys(tris[j-1]		, yris[j-1]		);
		k2 = odesys(tris[j-1] + h / 2.	, yris[j-1] + k1*h / 2.	);
		k3 = odesys(tris[j-1] + h / 2.	, yris[j-1] + k2*h / 2.	);
		k4 = odesys(tris[j-1] + h	, yris[j-1] + k3*h	);

		yris[j] = yris[j-1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6 * h;
		tris[j] = tris[j-1] + h;

		fprintf(fp1, "%.5e\t%.5e\n", tris[j], yris[j]);

	}

	//close the file
	fclose(fp1);

	//free memory
	free(tris);
	free(yris);
}

//ODE vector system function
double* odesysfun_vec(double t, double* y, int DIM) {

	double* dydt = (double*)malloc(DIM * sizeof(double));

	if (!dydt) return NULL;

	for (int j = 0; j < DIM; j++) {

		//attractor for 5, linear
		//dydt[j] = y[DIM - 1 - j] - y[j];

		//attractor for 1, nonlinear
		dydt[j] = y[DIM - 1 - j] - pow(y[j],2);
	
	}

	return dydt;

}

//RK4 solver implementation for vector systems
void rungekutta4_vec(double* y0, double t0, double tf, double h, double* odesys_vec(double t, double* y, int DIM), int DIM) {

	//constants for the RK4
	double* k1;
	double* k2;
	double* k3;
	double* k4;

	double* yris_col = (double*)malloc(DIM * sizeof(double));

	//number of points
	int N = (int)((tf - t0) / h + 0.5);

	//mallocation of the arrays with N elements
	double* tris = (double*)malloc(N * sizeof(double));

	//mallocation of rows
	double** yris = (double**)malloc(N * sizeof(double*));

	//mallocation of columns at each row
	for (int j = 0; j < N; j++){
		yris[j] = (double*)malloc(DIM * sizeof(double));
	}

	//in the end yris is a N * DIM matrix

	//setting results at zero, t = t0, y = y0:
	tris[0] = t0;
	for (int k = 0; k < DIM; k++){
		yris[0][k] = y0[k];
	}

	//we post process the data 
	char filename[] = "rk4_sys_vec.txt";
	//
	FILE* fp1 = fopen(filename, "w");
	// title
	fprintf(fp1, "t\t\t");
	for (int k = 0; k < DIM; k++){
		fprintf(fp1, "y%d\t\t", k + 1);
	}
	// second row t = 0
	fprintf(fp1, "\n%.5e\t", tris[0]);
	for (int k = 0; k < DIM; k++){
		fprintf(fp1, "%.5e\t", yris[0][k]);
	}

	//RK4 algorithm
	for (int j = 1; j < N; j++) {
		
		//define the row (j,:)
		for (int k = 0; k < DIM; k++){
			yris_col[k] = yris[j-1][k];
		}
		k1 = odesys_vec(tris[j-1], yris_col, DIM);

		for (int k = 0; k < DIM; k++){
			yris_col[k] = yris[j-1][k] + k1[k] * h / 2;
		}
		k2 = odesys_vec(tris[j-1] + h / 2., yris_col, DIM);

		for (int k = 0; k < DIM; k++){
			yris_col[k] = yris[j-1][k] + k2[k] * h / 2;
		}
		k3 = odesys_vec(tris[j-1] + h / 2., yris_col, DIM);

		for (int k = 0; k < DIM; k++){
			yris_col[k] = yris[j-1][k] + k3[k] * h;
		}
		k4 = odesys_vec(tris[j-1] + h, yris_col, DIM);

		for (int k = 0; k < DIM; k++){
			yris[j][k] = yris[j-1][k] + (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]) / 6 * h;
		}
		tris[j] = tris[j-1] + h;

		fprintf(fp1, "\n%.5e\t", tris[j]);
		for (int k = 0; k < DIM; k++){
			fprintf(fp1, "%.5e\t", yris[j][k]);
		}
	}

	//close the file
	fclose(fp1);

	//free memory
	free(tris);
	//
	for (int j = 0; j < N; j++){
		free(yris[j]);
	}
	free(yris);
	free(yris_col);
}

/*
	Andrea Landella

	program file rk4solver
*/
#include"rk4solver.h"

//RK4 solver implementation for scalar systems
void rungekutta4(double odesys(double t, double y), double t0, double tf, double y0, double h, double* &tSOL, double* &ySOL) {

	//we declare the constants for the RK4 algorithm
	double k1, k2, k3, k4;

	//we define the number of integration points
	int N = (int)((tf - t0) / h + 0.5);

	//we set results at zero, t = t0, y = y0:
	tSOL[0] = t0;
	ySOL[0] = y0;

	//we setup the result txt file 
	char filename[] = "rk4_sys.txt";
	//
	FILE* fp1 = fopen(filename, "w");
	//
	fprintf(fp1, "t\t\ty\n%.5e\t%.5e\n", tSOL[0], ySOL[0]);

	//RK4 algorithm cycle
	for (int j = 1; j < N; j++) {

		//we evaluate the constants
		k1 = odesys(tSOL[j - 1], ySOL[j - 1]);
		k2 = odesys(tSOL[j - 1] + h / 2., ySOL[j - 1] + k1*h / 2.);
		k3 = odesys(tSOL[j - 1] + h / 2., ySOL[j - 1] + k2*h / 2.);
		k4 = odesys(tSOL[j - 1] + h, ySOL[j - 1] + k3*h);

		//we update the solution
		ySOL[j] = ySOL[j - 1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6 * h;
		tSOL[j] = tSOL[j - 1] + h;

		fprintf(fp1, "%.5e\t%.5e\n", tSOL[j], ySOL[j]);
	}
	//close the file
	fclose(fp1);
}
   
//RK4 solver implementation for vector systems
void rungekutta4_vec(void odesys_vec(double t, double* y, int DIM, double dydt[]), double t0, double tf, double* y0, double h, int DIM, double* &tSOL, double** &ySOL) {

	//we declare the constants for the RK4 algorithm
	double* k1 = (double*)malloc(DIM * sizeof(double));	if (!k1) return;
	double* k2 = (double*)malloc(DIM * sizeof(double));	if (!k2) return;
	double* k3 = (double*)malloc(DIM * sizeof(double));	if (!k3) return;
	double* k4 = (double*)malloc(DIM * sizeof(double));	if (!k4) return;

	//we define the number of integration points
	int N = (int)((tf - t0) / h + 0.5);

	//we define the temporary vector for evaluating kj
	double* yris_col = (double*)malloc(DIM * sizeof(double)); if (!yris_col) return;

	//we set results at zero, t = t0, y = y0
	tSOL[0] = t0;
	for (int k = 0; k < DIM; k++){
		ySOL[0][k] = y0[k];
	}

	//we setup the result txt file 
	char filename[] = "rk4_sys_vec.txt";
	//
	FILE* fp1 = fopen(filename, "w");
	// title
	fprintf(fp1, "t\t\t");
	for (int k = 0; k < DIM; k++){
		fprintf(fp1, "y%d\t\t", k + 1);
	}
	//
	fprintf(fp1, "\n%.5e\t", tSOL[0]);
	for (int k = 0; k < DIM; k++){
		fprintf(fp1, "%.5e\t", ySOL[0][k]);
	}

	//RK4 algorithm cycle
	for (int j = 1; j < N; j++) {

		//we define the row (j,:)
		for (int k = 0; k < DIM; k++){
			yris_col[k] = ySOL[j - 1][k];
		}
		odesys_vec(tSOL[j - 1], yris_col, DIM, k1);

		for (int k = 0; k < DIM; k++){
			yris_col[k] = ySOL[j - 1][k] + k1[k] * h / 2;
		}
		odesys_vec(tSOL[j - 1] + h / 2., yris_col, DIM, k2);

		for (int k = 0; k < DIM; k++){
			yris_col[k] = ySOL[j - 1][k] + k2[k] * h / 2;
		}
		odesys_vec(tSOL[j - 1] + h / 2., yris_col, DIM, k3);

		for (int k = 0; k < DIM; k++){
			yris_col[k] = ySOL[j - 1][k] + k3[k] * h;
		}
		odesys_vec(tSOL[j - 1] + h, yris_col, DIM, k4);

		//we update the solution
		for (int k = 0; k < DIM; k++){
			ySOL[j][k] = ySOL[j - 1][k] + (k1[k] + 2 * k2[k] + 2 * k3[k] + k4[k]) / 6 * h;
		}
		tSOL[j] = tSOL[j - 1] + h;

		fprintf(fp1, "\n%.5e\t", tSOL[j]);
		for (int k = 0; k < DIM; k++){
			fprintf(fp1, "%.5e\t", ySOL[j][k]);
		}
	}
	//close the file
	fclose(fp1);

	//free memory
	free(k1); free(k2); free(k3); free(k4);
	//
	free(yris_col);
}
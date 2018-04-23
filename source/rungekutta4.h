/*
Andrea Landella

header file rungekutta4.cpp
*/

double odesysfun(double t, double y);

void rungekutta4(double y0, double t0, double tf, double h, double odesys(double t, double y));

double* odesysfun_vec(double t, double* y, int DIM);

void rungekutta4_vec(double* y0, double t0, double tf, double h, double* odesys_vec(double t, double* y, int DIM), int DIM);

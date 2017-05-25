/*
* rk45.h: Header file for Runge Kutta 4
*/

#ifndef RK45_H
#define RK45_H

typedef void (*func3arg) (double, double*, double*);

typedef struct {
    double *k1;
    double *k2;
    double *k3;
    double *k4;
    double *k5;
    double *k6;
} kvector;

void rk45 ( int neqn, int nparams, double *t0, double *x0, double *h, func3arg f, double TOL, double hmin, double hmax, double tfinal);

kvector rk45_kvector ( int neqn, int nparams, double t0, double *x0, double h, func3arg f);

void optimal_step_size(double TOL, double R, double *h, double hmin, double hmax, double t, double tfinal);

#endif

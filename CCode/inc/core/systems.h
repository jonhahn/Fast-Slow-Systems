#ifndef SYSTEMS_H
#define SYSTEMS_H

void dgauss_elim ( double **a, double *b, int n, int ps );

void invert(   double ** A,    //input square matrix A
                double ** B,
                    int n);          //dimension of A


#endif // SYSTEMS_H

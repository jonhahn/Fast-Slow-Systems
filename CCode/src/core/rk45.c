/*
* Runge-Kutta 4
*/

#include <stdio.h>
#include "core/rk45.h"
#include <math.h>

kvector rk45_kvector( int neqn,     //number of equations
                      int nparams,  //number of parameters
                      double t0,    //starting time
                      double *x0,   //state and parameter vector //TODO: separate state and parameters
                      double h,     //step size
                      func3arg f)   //function callback
/*
*Returns the values of k1, k2, k3, k4, k5, k6 used in the Runge-Kutta 4/5 method
*/

{
    int i;
    double xp[neqn], k1[neqn], k2[neqn], k3[neqn], k4[neqn], k5[neqn], k6[neqn];
    double  xtilde[neqn + nparams];

    for (i = 0; i < neqn + nparams; i++){
         xtilde[i] = x0[i];
     }

    //k1
    f(t0, x0, xp);
    for ( i = 0; i < neqn; i++ ) {
      k1[i] = h*xp[i];
      xtilde[i] = x0[i] + k1[i]/4.0;
    }

    //k2
    f(t0 + h/4.0, xtilde, xp);
    for ( i = 0; i < neqn; i++ ) {
      k2[i] = h*xp[i];
      xtilde[i] = x0[i] + 3.0*k1[i]/32.0 + 9.0*k2[i]/32.0;
    }

    //k3
    f(t0 + 3.0*h/8.0, xtilde, xp);
    for ( i = 0; i < neqn; i++ ) {
      k3[i] = h*xp[i];
      xtilde[i] = x0[i] + 1932.0*k1[i]/2197.0 - 7200.0*k2[i]/2197.0 + 7296.0*k3[i]/2197.0 ;
    }

    //k4
    f(t0 + 12.0*h/13.0, xtilde, xp);
    for ( i = 0; i < neqn; i++ ) {
      k4[i] = h*xp[i];
      xtilde[i] = x0[i] + 439.0*k1[i]/216.0 - 8.0*k2[i] + 3680.0*k3[i]/513.0 - 845.0*k4[i]/4104.0;
    }

    //k5
    f(t0 + h, xtilde, xp);
    for ( i = 0; i < neqn; i++ ) {
      k5[i] = h*xp[i];
      xtilde[i] = x0[i] - 8.0*k1[i]/27.0 + 2.0*k2[i] - 3544.0*k3[i]/2565.0 + 1859.0*k4[i]/4104.0 - 11.0*k5[i]/40.0;
    }

    //k6
    f(t0 + h/2.0, xtilde, xp);
    for ( i = 0; i < neqn; i++ ){
      k6[i] = h*xp[i];
    }

    kvector kvec;
    kvec.k1 = k1;
    kvec.k2 = k2;
    kvec.k3 = k3;
    kvec.k4 = k4;
    kvec.k5 = k5;
    kvec.k6 = k6;

    return kvec;
}

void rk45 ( int neqn,
            int nparams,
            double *t0,
            double *x0,
            double *h,
            func3arg f,
            double TOL,
            double hmin,
            double hmax,
            double tfinal)

/*
     PURPOSE:
          perform single time step of the classical fourth-order Runge-Kutta 
          method to approximate the solution of an initial value problem 
          (this routine works for a single equation as well as for a system 
          of equations)
          

     CALLING SEQUENCE:
          rk4 ( neqn, t0, x0, h, f );

     INPUTS:
          neqn          number of equations in initial value problem
                        whose solution is being approximated
                        type:  int
          t0        initial value for independent variable
                        type:  double
          x0        initial value for dependent variable
                        type:  *double
          h     length of time step
                        type:  double
          f             function of three arguments which defines the
                        right-hand side of the differential equation;
                        this function must be of the form
                        
                           void f ( double t, double *x, double *xp )
                           {
                              xp[0] = ...;
                              xp[1] = ...;
                              xp[2] = ...;
                                    .
                                    .
                                    .
                              xp[neqn-1] = ...;
                           }
                           
                        where xp[i] is the value of the right-hand side
                        of the i-th equaton in the system
                        type:  func3arg
          TOL: error tolerance over one unit of time
          hmin: smallest allowable time step
          hmax: largest allowable time step


     OUTPUT:
          x0        approximation to solution of initial value 
                        problem at t0+h
                        type:  *double
*/

{

     int i;
     double error;  //used to hold error
     kvector kv;

    if (*t0 + *h > tfinal){
        double finalh = tfinal - *t0;
        //Get vector of values of k1, k2, ...
        kv = rk45_kvector(neqn, nparams, *t0, x0, finalh, f); 
           //Get the new t value
        *t0 = *t0 + finalh;
    }
    else
    {

     //Get vector of values of k1, k2, ...
     kv = rk45_kvector(neqn, nparams, *t0, x0, *h, f);  

     //Estimate the error, R
     double R = 0.0;
     for (i = 0; i< neqn; i++){
        error = fabs(kv.k1[i]/360.0 - 128.0*kv.k3[i]/4275.0 - 2197.0*kv.k4[i]/75240.0 + kv.k5[i]/50.0 + 2.0*kv.k6[i]/55.0)/(*h);
        if (error > R){
          R = error;
        }
     }

     //If the error is larger than the tolerance, find the optimal step size h and try again
     if (R > TOL){

          optimal_step_size(TOL, R, h, hmin, hmax, *t0, tfinal);

         //Again get vector of values of k1, k2, ...
          kv = rk45_kvector(neqn, nparams, *t0, x0, *h, f);

          //Again find the error
          for (i = 0; i< neqn; i++){
            error = fabs(kv.k1[i]/360.0 - 128.0*kv.k3[i]/4275.0 - 2197.0*kv.k4[i]/75240.0 + kv.k5[i]/50.0 + 2.0*kv.k6[i]/55.0)/(*h);
            if (error > R){
              R = error;
            }
          }
      }

      //Get the new t value
      *t0 = *t0 + *h;

      optimal_step_size(TOL, R, h, hmin, hmax, *t0 - *h, tfinal);

    }

    //Get the new x value
    for ( i = 0; i < neqn; i++ ){
        x0[i] = x0[i] + 16.0*kv.k1[i]/135.0 + 6656.0*kv.k3[i]/12825.0 + 28561.0*kv.k4[i]/56430.0 - 9.0*kv.k5[i]/50.0 + 2.0*kv.k6[i]/55.0;
    }

}



void optimal_step_size(double TOL, double R, double *h, double hmin, double hmax, double t, double tfinal){
      if (*h < 0){
         *h = hmin;
      }

      //Again, find the optimal step size
      double q = 0.84 * pow(TOL / R, .25);
      double hnew = fmin ( fmax (q, 0.1 ), 4.0 ) * (*h);
      while (*h > hnew){
          *h = .5*(*h);
        }
      while (*h < .5*hnew){
          *h = 2.0*(*h);
      }

      //Get correct step size to make sure it is not too large or too small
        if (*h > hmax){
            *h = hmax;
        }
        if (t + *h > tfinal){
            *h = tfinal - t;
          }
        else if (*h < hmin){
          //This is bad.
          *h = hmin;
      }
}

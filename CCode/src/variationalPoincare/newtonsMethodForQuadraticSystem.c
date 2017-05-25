#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "core/systems.h"
#include "core/outputJson.h"
#include "core/rk45.h"

#include "newtonsMethodForQuadraticSystem.h"

void differential_for_quadratic(
    double t,
     double *x,
      double *xp)
{
    //Gives the differential and variational equations for the quadratic fast-slow system
    //x' = (z + x^2)/epsilon
    //z' = -x + alpha + 2 pi omega r Sin(2 pi omega t)

    //parameters
    double alpha = .01;
    double epsilon = .05;
    double r = .02;
    double omega = x[6];

    //state
    double xx = x[0];
    double zz = x[1];
    double dxdx0 = x[2];
    double dxdz0 = x[3];
    double dzdx0 = x[4];
    double dzdz0 = x[5];

    double dxdx0dt = 1/epsilon * dzdx0 + (2/epsilon) * xx * dxdx0;
    double dxdz0dt = 1/epsilon * dzdz0 + (2/epsilon) * xx * dxdz0;
    double dzdx0dt = -dxdx0;
    double dzdz0dt = -dxdz0;

    //
    //double dx1 = (x2 + lambda + (x1 - x1 * x1 * x1/3)) / epsilon;
    double dx = (zz + xx * xx)/epsilon;
    double dz = -xx - alpha + 2 * M_PI * omega * r * sin(2 * M_PI * omega * t);

    xp[0] = dx;
    xp[1] = dz;
    xp[2] = dxdx0dt;
    xp[3] = dxdz0dt;
    xp[4] = dzdx0dt;
    xp[5] = dzdz0dt;
}

void trace_one_period(
    double *initial_xz,
     double *end_xz,
      double **Jacobian,
       double period)
{
    //Solves the quadratic system forward by one period

    double omega = 1/period;

    double modelState[7] = {0,0,0,0,0,0,0};
    modelState[0] = initial_xz[0];
    modelState[1] = initial_xz[1];
    modelState[2] = 1;
    modelState[3] = 0;
    modelState[4] = 0;
    modelState[5] = 1;
    modelState[6] = omega;
    double currentTime = 0; //initial time

    //rk4 parameters
    int numberOfEquations = 6;  //x, y, ds
    int numberOfParameters = 1; //omega
    double timestep = .0001;

    double TOL = 1e-12;
    double hmin = .0001;
    double hmax = .01;


    double nextTime =  period;

    while (currentTime < nextTime){

        rk45 ( numberOfEquations,
                numberOfParameters,
                &currentTime,
                modelState,
                &timestep,
                differential_for_quadratic,
                TOL,
                hmin,
                hmax,
                nextTime);
    }

    //printf("The fuck %f %f \n", modelState[0], modelState[1]);

    //populate output
    end_xz[0] = modelState[0];
    end_xz[1] = modelState[1];
    //printf("The fuck %f %f \n", end_xz[0], end_xz[1]);

    if (Jacobian != NULL)
    {
        Jacobian[0][0] = modelState[2];
        Jacobian[0][1] = modelState[3];
        Jacobian[1][0] = modelState[4];
        Jacobian[1][1] = modelState[5];
    }

 
    //x_out[i] = modelState[0];
    //z_out[i] = modelState[1];

    
}

void find_poincare_fixed_point(
    double *initial_xz, //also the output fixed points
     double period)
{
    double end_xz[2];
    double Jacobian[2][2];
    double *pointerToJacobian[2];
    pointerToJacobian[0] = Jacobian[0];
    pointerToJacobian[1] = Jacobian[1];

    double J[2][2];
    double *pointerToJ[2];
    pointerToJ[0] = J[0];
    pointerToJ[1] = J[1];
    double f[2];

    double diff = 1;
    int j = 0;
    while(diff > .0000000001 && diff < 100 && j < 1000){
        j = j + 1;
        trace_one_period(
            initial_xz,
             end_xz,
              pointerToJacobian,
               period);

        //calculate difference
        diff = pow(pow(initial_xz[0] - end_xz[0],2) + pow(initial_xz[1] - end_xz[1],2),.5);
        //printf("%f\n", diff);

        J[0][0] = Jacobian[0][0] - 1;
        J[0][1] = Jacobian[0][1];
        J[1][0] = Jacobian[1][0];
        J[1][1] = Jacobian[1][1] - 1;

        f[0] = initial_xz[0] - end_xz[0];
        f[1] = initial_xz[1] - end_xz[1];

        dgauss_elim(pointerToJ, f, 2, 1);
        //printf("%f %f\n", f[0], f[1] );

        initial_xz[0] = initial_xz[0] + f[0];
        initial_xz[1] = initial_xz[1] + f[1];

    }

    printf("(%f, %f) -- %f -- period: %f\n", initial_xz[0], initial_xz[1], diff, period);

}








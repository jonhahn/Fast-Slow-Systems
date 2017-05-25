#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "core/systems.h"
#include "core/outputJson.h"
#include "core/rk45.h"

static void differential_for_van_der_pol(double t, double *x, double *xp)
{
    //parameters
    double alpha = 1.01;
    double epsilon = .1;
    double r = .062;
    double omega = x[6];

    //state
    double x1 = x[0];
    double x2 = x[1];

    double dx1 = (x2 + (x1 - x1 * x1 * x1/3))/epsilon;
    double dx2 = -x1 - alpha + 2 * M_PI * omega * r * sin(2 * M_PI * omega * t);
    double dx1dx10 = x[2];
    double dx1dx20 = x[3];
    double dx2dx10 = x[4];
    double dx2dx20 = x[5];

    double dx1dx10dt = 1/epsilon * dx2dx10 + (1/epsilon) *(1 - x1 * x1) * dx1dx10;
    double dx1dx20dt = 1/epsilon * dx2dx20 + (1/epsilon) *(1 - x1 * x1) * dx1dx20;
    double dx2dx10dt = -dx1dx10;
    double dx2dx20dt = -dx1dx20;

    xp[0] = dx1;
    xp[1] = dx2;
    xp[2] = dx1dx10dt;
    xp[3] = dx1dx20dt;
    xp[4] = dx2dx10dt;
    xp[5] = dx2dx20dt;
}



static void trace_one_period_vdp(
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
                differential_for_van_der_pol,
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

static void find_poincare_fixed_point_vdp(
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
    while(diff > .0000000001 && diff < 100 && j < 2000){
        j = j + 1;
        trace_one_period_vdp(
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

static void continuation_method_vdp(){
    double initial_xz[2];
    initial_xz[0] =  -1.033210;
    initial_xz[1] = 0.685003;
    double period = 3.845729;

    double previous_xz[2];
    previous_xz[0] = initial_xz[0];
    previous_xz[1] = initial_xz[1];

    double current_xz[2];

    while (period > 3.0 && !(isnan(initial_xz[0]))){
        find_poincare_fixed_point_vdp(initial_xz, period);
        current_xz[0] = initial_xz[0];
        current_xz[1] = initial_xz[1];
        initial_xz[0] = initial_xz[0] + current_xz[0] - previous_xz[0];
        initial_xz[1] = initial_xz[1] + current_xz[1] - previous_xz[1];
        previous_xz[0] = current_xz[0];
        previous_xz[1] = current_xz[1];
        period = period - .0000001;
    }
}

int main()
{
    double initial_xz[2];
    initial_xz[0] = -1.01068;
    initial_xz[1] = 0.669288;
    double period = 10.0;


    find_poincare_fixed_point_vdp(initial_xz, period);
    continuation_method_vdp();

    //printf("%f %f\n", end_xz[0], end_xz[1]);
    //printf("%f %f\n%f %f\n", Jacobian[0][0],Jacobian[0][1],Jacobian[1][0],Jacobian[1][1]);


    return 1;
}
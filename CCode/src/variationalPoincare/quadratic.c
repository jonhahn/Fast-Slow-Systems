
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "core/outputJson.h"
#include "core/rk45.h"

static void differential_for_quadratic(double t, double *x, double *xp)
{
    //parameters
    double alpha = .01;
    double epsilon = .1;
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

static void solve_VDP_with_rk45(
    int n,
    double x_0,
    double z_0,
    double max_increment,
    double * times_in,
    double * x_out,
    double * z_out,
    double *** grad,
    double omega)
{
    double modelState[7];
    modelState[0] = x_0;
    modelState[1] = z_0;
    modelState[2] = 1;
    modelState[3] = 0;
    modelState[4] = 0;
    modelState[5] = 1;
    modelState[6] = omega;
    double currentTime = 0; //initial time

    //rk4 parameters
    int numberOfEquations = 6;  //x, y, ds
    int numberOfParameters = 1; //omega
    double timestep;

    double TOL = 1e-12;
    double hmin = .001;
    double hmax = .01;


    for (int i=0; i < n; i++){

        double nextTime =  times_in[i] - times_in[0];

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
        //populate output
        x_out[i] = modelState[0];
        z_out[i] = modelState[1];
        printf("%f %f\n", x_out[i], z_out[i]);
    }
}


void testQuadratic()
{
    //double x1_0 = -1.01;
    //double x2_0 = .6666;
    double x_0 = 0;
    double z_0 = 0;
    double omega = 2;

    double max_increment = .001;
    int n = 20000;

    double times[n];
    for(int i = 0; i<n; i++){
        times[i] = i/100.0;
    }
    double x_out[n];
    double z_out[n];

    double grad[2][2][n];
    double *** gradpointer = NULL;

    solve_VDP_with_rk45(
        n,
        x_0,
        z_0,
        max_increment,
        times,
        x_out,
        z_out,
        gradpointer,
        omega);


    //HMMMM
    /*
    
    //Output things:
    struct DataSetForJson outputdata[2];

    outputdata[0].name = "poincare_in";
    outputdata[0].data = poincare_in;
    outputdata[0].length = m;
    char cated_string[2*phase][100*sizeof(char)];

    for (int i=0; i<phase; i++){
        sprintf(cated_string[i],"%s%d","poincare_out",i);

        outputdata[i+1].name = cated_string[i];
        outputdata[i+1].data = poincare_out[i];
        outputdata[i+1].length = m;
    }

    for (int i=0; i<phase; i++){
        sprintf(cated_string[phase + i],"%s%d","times_out",i);

        outputdata[phase + i+1].name = cated_string[phase + i];
        outputdata[phase + i+1].data = times_out[i];
        outputdata[phase + i+1].length = m;
    }

    outputJsonData("src/variationalPoincare/quadratic.json",
                "Quadratic",
                outputdata,
                2);
    */
    
}


int main()
{
    testQuadratic();
    return 1;
}


/*
Test RK4 and RK45 on differential equations with known solutions
*/

#include <stdio.h>
#include <math.h>
#include <sys/time.h>

//#define M_PI acos(-1.0)

#include "core/outputJson.h"
#include "core/rk45.h"

static void differential_for_van_der_pol(double t, double *x, double *xp)
{
    //parameters
    double alpha = 1.01;
    double epsilon = .1;
    double r = .02;
    double omega = x[3];
    double dlambda = omega * r / (2.0 * M_PI) * sin(omega * t/ (2.0 * M_PI));

    //state
    double x1 = x[0];
    double x2 = x[1];
    double lambda = x[2];

    //
    double dx1 = (x2 + lambda + (x1 - x1 * x1 * x1/3)) / epsilon;
    double dx2 = -x1 - alpha;

    xp[0] = dx1;
    xp[1] = dx2;
    xp[2] = dlambda;
}



static void solve_equation_with_rk45(
    int n,
    double x1_0,
    double x2_0,
    double max_increment,
    double * times_in,
    double * x1_out,
    double * x2_out,
    double omega)
{
    double modelState[4];
    modelState[0] = x1_0;
    modelState[1] = x2_0;
    modelState[2] = 0;
    modelState[3] = omega;

    double currentTime = 0; //initial time

    //rk4 parameters
    int numberOfEquations = 3;  //device Temperature
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
                    differential_for_van_der_pol,
                    TOL,
                    hmin,
                    hmax,
                    nextTime);
        }
        //populate output
        x1_out[i] = modelState[0];
        x2_out[i] = modelState[1];
    }
}


static void test_tipping_with_rk45(
    int n,
    double x1_0,
    double x2_0,
    double max_increment,
    double * times_in,
    double * x1_out,
    double * x2_out,
    int * tip,
    double omega)
{
    double modelState[4];
    modelState[0] = x1_0;
    modelState[1] = x2_0;
    modelState[2] = 0;
    modelState[3] = omega;

    double currentTime = 0; //initial time

    //rk4 parameters
    int numberOfEquations = 3;  //
    int numberOfParameters = 1; //
    double timestep;

    double TOL = 1e-12;
    double hmin = .001;
    double hmax = .01;

    *tip = 0;

    for (int i=0; i < n; i++){

        double nextTime =  times_in[i] - times_in[0];

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
        //populate output
        x1_out[i] = modelState[0];
        x2_out[i] = modelState[1];

        if (x1_out[i] > .5 && nextTime > 100){
            *tip = 1;
            break;
        }
    }
}


static void test_max_with_rk45(
    int n,
    double x1_0,
    double x2_0,
    double max_increment,
    double * times_in,
    double * x1_out,
    double * x2_out,
    double * max,
    double omega)
{
    double modelState[4];
    modelState[0] = x1_0;
    modelState[1] = x2_0;
    modelState[2] = 0;
    modelState[3] = omega;

    double currentTime = 0; //initial time

    //rk4 parameters
    int numberOfEquations = 3;  //
    int numberOfParameters = 1; //
    double timestep;

    double TOL = 1e-12;
    double hmin = .001;
    double hmax = .01;

    *max = -20;

    for (int i=0; i < n; i++){

        double nextTime =  times_in[i] - times_in[0];

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
        //populate output
        x1_out[i] = modelState[0];
        x2_out[i] = modelState[1];

        if (x1_out[i] > *max && nextTime > 300){
            *max = x1_out[i];
        }
    }

}



void testVanDerPol()
{
    double x1_0 = -1.01;
    double x2_0 = .6666;
    double omega = 15.21854;

    double max_increment = .1;
    int n = 200000;

    double times[n];
    for(int i = 0; i<n; i++){
        times[i] = i/100.0;
    }
    double x1_out[n];
    double x2_out[n];


    solve_equation_with_rk45(
        n,
        x1_0,
        x2_0,
        max_increment,
        times,
        x1_out,
        x2_out,
        omega);

    
    //Output things:
    struct DataSetForJson outputdata[2];

    outputdata[0].name = "x1";
    outputdata[0].data = x1_out;
    outputdata[0].length = n;

    outputdata[1].name = "x2";
    outputdata[1].data = x2_out;
    outputdata[1].length = n;

    outputJsonData("src/testResults/plot.json",
                "vanDerPol",
                outputdata,
                2);
}



void testVanDerPol2()
{
    double x1_0 = -1.01;
    double x2_0 = .6666;

    double max_increment = .1;
    int n = 50000;

    double times[n];
    for(int i = 0; i<n; i++){
        times[i] = i/100.0;
    }
    double x1_out[n];
    double x2_out[n];

    int tip;

    for (double omega = 15.0; omega < 22.5; omega = omega + .01){
        test_tipping_with_rk45(
            n,
            x1_0,
            x2_0,
            max_increment,
            times,
            x1_out,
            x2_out,
            &tip,
            omega);
        printf("%f %d\n", omega, tip);
    }
    //Output things:
    struct DataSetForJson outputdata[2];
}



void testVanDerPol3()
{
    double x1_0 = -1.01;
    double x2_0 = .6666;

    double max_increment = .1;
    int n = 50000;

    double times[n];
    for(int i = 0; i<n; i++){
        times[i] = i/100.0;
    }
    double x1_out[n];
    double x2_out[n];

    double max;

    double maxes[20000];

    double omega = 10;
    for (int i = 0; i < 20000; i++){
        test_max_with_rk45(
            n,
            x1_0,
            x2_0,
            max_increment,
            times,
            x1_out,
            x2_out,
            &max,
            omega);
        maxes[i] = max;
        printf("%f %f\n", omega, max);
        omega = omega + .001;
    }

    //Output things:
    struct DataSetForJson outputdata[1];

    outputdata[0].name = "max x1";
    outputdata[0].data = maxes;
    outputdata[0].length = 20000;


    outputJsonData("src/testResults/maxes3.json",
                "vanDerPol",
                outputdata,
                1);

}


void testVanDerPol_Zoom()
{
    double x1_0 = -1.01;
    double x2_0 = .6666;

    double max_increment = .1;
    int n = 50000;

    double times[n];
    for(int i = 0; i<n; i++){
        times[i] = i/100.0;
    }
    double x1_out[n];
    double x2_out[n];

    double max;

    double maxes[100];

    double omega = 15.2185;
    for (int i = 0; i < 100; i++){
        test_max_with_rk45(
            n,
            x1_0,
            x2_0,
            max_increment,
            times,
            x1_out,
            x2_out,
            &max,
            omega);
        maxes[i] = max;
        printf("%f %f\n", omega, max);
        omega = omega + .000001;
    }

    //Output things:
    struct DataSetForJson outputdata[1];

    outputdata[0].name = "max x1";
    outputdata[0].data = maxes;
    outputdata[0].length = 100;


    outputJsonData("src/testResults/maxes_zoom15point2185.json",
                "vanDerPol",
                outputdata,
                1);

}

int main()
{
    testVanDerPol_Zoom();
    return 1;
}

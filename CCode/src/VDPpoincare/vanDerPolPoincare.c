/*
Find the poincare map for a quadratic system
*/

#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "core/outputJson.h"
#include "core/rk45.h"

static void differential_for_van_der_pol(double t, double *x, double *xp)
{
    //parameters
    double alpha = .01;
    double epsilon = .1;
    double r = .062;
    double omega = x[3];

    //state
    double x1 = x[0];
    double x2 = x[1];
    double lambda = r * (-cos(omega * t/(2.0 * M_PI)) + 1);

    //
    //double dx1 = (x2 + lambda + (x1 - x1 * x1 * x1/3)) / epsilon;
    double dx1 = (x2 + lambda + x1 * x1)/epsilon;
    double dx2 = -x1 - alpha;

    xp[0] = dx1;
    xp[1] = dx2;
}


static void solve_VDP_with_rk45(
    int n,
    double x1_0,
    double x2_0,
    double max_increment,
    double * times_in,
    double * x1_out,
    double * x2_out,
    double omega,
    double * poincare_in,
    double ** poincare_out,
    double ** times_out,
    int phase,
    int m)
{
    double modelState[4];
    modelState[0] = x1_0;
    modelState[1] = poincare_in[0];
    modelState[2] = 0;
    modelState[3] = omega;

    double currentTime = 0; //initial time

    //rk4 parameters
    int numberOfEquations = 2;  //vanDerPol x, y, lambda
    int numberOfParameters = 2; //omega, lambda
    double timestep;

    double TOL = 1e-12;
    double hmin = .001;
    double hmax = .01;

    /*
    for(int i=0; i<m; i++){
        poincare_out_max[i] = -100;
        poincare_out_min[i] = 100;
    }
    */

    for (int k=0; k<phase; k++){
        double offset = k*2.0*M_PI/phase/omega;
        printf("%d\n", k);

        for (int j=0; j < m; j++){
            modelState[0] = x1_0;
            modelState[1] = poincare_in[j];
            modelState[2] = 0;

            poincare_out[k][j] = 0;
            currentTime = times_in[0] + offset;
            int square = 0;
            for (int i=0; i < n; i++){

                double nextTime =  times_in[i] - times_in[0] + offset;

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

                    if (x1_out[i-1] < x1_0 && modelState[1] > x1_0){
                        break;
                    }
                }
                if (x1_out[i-1] < x1_0 && modelState[0] > x1_0){
                    poincare_out[k][j] = modelState[1]; 
                    times_out[k][j] = nextTime;
                    //printf("%f, %f\n", modelState[1], nextTime);
                    break;
                }
                //populate output
                x1_out[i] = modelState[0];
                x2_out[i] = modelState[1];
            }
        }

        /*
        for(int i=0; i<m; i++){
            if (poincare_out_max[i] < poincare_out[i]){
                poincare_out_max[i] = poincare_out[i];
            }
            if (poincare_out_min[i] > poincare_out[i]){
                poincare_out_min[i] = poincare_out[i];
            }
        }
        */
    }
}

void testVanDerPol()
{
    //double x1_0 = -1.01;
    //double x2_0 = .6666;
    double x1_0 = 0;
    double x2_0 = 0;
    double omega = 25;

    double max_increment = .001;
    int n = 20000;

    double times[n];
    for(int i = 0; i<n; i++){
        times[i] = i/100.0;
    }
    double x1_out[n];
    double x2_out[n];

    int m = 100;
    int phase = 5;
    double poincare_in[m];
    double poincare_out[phase][m];
    double * pointerToPoincareOut[phase];
    for(int i = 0; i <phase; i++){
        pointerToPoincareOut[i] = poincare_out[i];
    }
    double times_out[phase][m];
    double * pointerToTimesOut[phase];
    for(int i = 0; i <phase; i++){
        pointerToTimesOut[i] = times_out[i];
    }
    //double poincare_out_max[m];
    //double poincare_out_min[m];
    for (int i = 0; i<m; i++){
        poincare_in[i] = -.05 + .1*(i)/100.0;
    }

    solve_VDP_with_rk45(
        n,
        x1_0,
        x2_0,
        max_increment,
        times,
        x1_out,
        x2_out,
        omega,
        poincare_in,
        pointerToPoincareOut,
        pointerToTimesOut,
        phase,
        m);

    //for (int i = 0; i<m; i++){
    //    printf("%f %f\n", poincare_in[i], poincare_out[i]);
    //}
    
    //Output things:
    struct DataSetForJson outputdata[2*phase + 1];

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

    outputJsonData("src/poincare/poincare.json",
                "vanDerPol",
                outputdata,
                2*phase + 1);

}


int main()
{
    testVanDerPol();
    return 1;
}
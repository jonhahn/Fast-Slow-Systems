#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "core/systems.h"
#include "core/outputJson.h"
#include "core/rk45.h"

#include "newtonsMethodForQuadraticSystem.h"



static void continuation_method(){
    double initial_xz[2];
    initial_xz[0] = -0.031617;
    initial_xz[1] =  0.028226;
    double period = 2.5516;

    while (period > 2.5 && !(isnan(initial_xz[0]))){
        find_poincare_fixed_point(initial_xz, period);
        period = period - .0000001;
    }
}

int main()
{
    double initial_xz[2];
    initial_xz[0] = -0.023925;
    initial_xz[1] =  0.022379;
    double period = 2.6;


    find_poincare_fixed_point(initial_xz, period);
    continuation_method();

    //printf("%f %f\n", end_xz[0], end_xz[1]);
    //printf("%f %f\n%f %f\n", Jacobian[0][0],Jacobian[0][1],Jacobian[1][0],Jacobian[1][1]);


    return 1;
}
//onePeriod.c

#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "core/systems.h"
#include "core/outputJson.h"
#include "core/rk45.h"
#include "core/parson/parson.h"

#include "newtonsMethodForQuadraticSystem.h"


static void follow_square_one_period(
    double xmin,
     double xmax,
      double zmin,
       double zmax)
{
    double initial_xz[2];
    double end_xz[2];

    double **Jacobian = NULL;
    double period = 3.0;

    int II = 20;
    int JJ = 20;

    double xs[II*JJ];
    double zs[II*JJ];

    initial_xz[1] = zmin;

    for(int i=0; i<II; i++){
        initial_xz[0] = xmin;
        for(int j=0; j<JJ; j++){
            trace_one_period(
                initial_xz,
                 end_xz,
                  Jacobian,
                   period);
            xs[10*i + j] = 0.0;
            zs[10*i + j] = 0.0;
            if (!(isnan(end_xz[0]))){
                xs[JJ*i + j] = end_xz[0];
                zs[JJ*i + j] = end_xz[1];
            }
            printf("(%f, %f) ", end_xz[0], end_xz[1]);


            initial_xz[0] = initial_xz[0] + (xmax - xmin)/JJ; 
        }
        initial_xz[1] = initial_xz[1] + (zmax - zmin)/II;
        printf("\n");
    }

    int numberOfParameters = 4;
    int numberOfDataSets = 2;
    struct DataSetForJson outputdata[numberOfDataSets];
    struct ParameterList parameterList[numberOfParameters];

    outputdata[0].name = "xs";
    outputdata[0].data = xs;
    outputdata[0].length = II * JJ;
    outputdata[1].name = "zs";
    outputdata[1].data = zs;
    outputdata[1].length = II * JJ;


    char outputPath[100] = "src/variationalPoincare/output/square.json";
    parameterList[0].name = "xmin";
    parameterList[1].name = "xmax";
    parameterList[2].name = "zmin";
    parameterList[3].name = "zmax";
    parameterList[0].value = xmin;
    parameterList[1].value = xmax;
    parameterList[2].value = zmin;
    parameterList[3].value = zmax;
    

    outputJsonData_and_Parameters(
        outputPath,
        "Square",
        outputdata,
        numberOfDataSets,
        parameterList,
        numberOfParameters);
    

}

static void follow_one(){
    double initial_xz[2];
    initial_xz[0] = -0.014213;
    initial_xz[1] = .009918;
    double end_xz[2];
    double period = 3.0;
    double **Jacobian = NULL;

    trace_one_period(
        initial_xz,
        end_xz,
        Jacobian,
        period);
    printf("(%f, %f) %f--%f ", initial_xz[0], initial_xz[1], end_xz[0], end_xz[1]);

}

int main()
{
    //follow_one();
    follow_square_one_period(-0.2, .4, -0.05,  0.07);
}


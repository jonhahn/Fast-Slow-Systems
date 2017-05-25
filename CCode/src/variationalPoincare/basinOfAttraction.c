//basinOfAttraction.h

#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include "core/systems.h"
#include "core/outputJson.h"
#include "core/rk45.h"
#include "newtonsMethodForQuadraticSystem.h"


static int is_it_in_basin(
    double *initial_xz,
    double period)
{

    double **Jacobian = NULL;

    for(int i=0; i<50; i++){
        double end_xz[2];
        trace_one_period(        
            initial_xz,
             end_xz,
              Jacobian,
               period);
        //printf("%f %f\n", end_xz[0], end_xz[1]);
        if (end_xz[0] > 2.0){
            return 0;
        }

        if (isnan(end_xz[0])){
            return 0;
        }

        if (fabs(end_xz[0] - initial_xz[0])< .00001 && fabs(end_xz[1] - initial_xz[1]) < .00001){
            return 1;
        }

        initial_xz[0] = end_xz[0];
        initial_xz[1] = end_xz[1];
        //printf("New: %f %f\n", initial_xz[0], initial_xz[1]);
    }
    return 1;
}





static void find_basin(){
    double initial_xz[2];
    initial_xz[0] = -0.014392;
    initial_xz[1] =  0.009918;
    double period = 3.0;

    //int what = is_it_in_basin(initial_xz, period);
    //printf("%d\n", what);

    find_poincare_fixed_point(initial_xz, period);


    double end_xz[2];
    double **Jacobian = NULL;
    trace_one_period(
        initial_xz,
         end_xz,
          Jacobian,
           period);
    //printf("end_xz: %f %f\n", end_xz[0], end_xz[1]);

    double xmin = -.2;
    double xmax = .4;
    double zmin = -.05;
    double zmax = .07;

    double xs_inbasin[2500];
    double zs_inbasin[2500];

    int numberInBasin = 0;

    double z = zmin;
    int II = 50;
    int JJ = 50;
    int inBasin = 0;
    for(int i=0; i<II; i++){
        double x = xmin;
        for(int j=0; j<JJ; j++){
            initial_xz[0] = x;
            initial_xz[1] = z;

            //printf("%f\n", initial_xz[0]);
            inBasin = is_it_in_basin(initial_xz, period);
            if(inBasin == 1){
                xs_inbasin[numberInBasin] = x;
                zs_inbasin[numberInBasin] = z;
                numberInBasin = numberInBasin + 1;
            }
            printf("%d ", inBasin);
        x = x + (xmax - xmin)/JJ;
        }
        printf("\n");
        z = z + (zmax - zmin)/II;
    }
    printf("%d\n", numberInBasin);

    double newXsInBasin[numberInBasin];
    double newZsInBasin[numberInBasin];

    for(int i =0; i<numberInBasin; i++){
        newXsInBasin[i] = xs_inbasin[i];
        newZsInBasin[i] = zs_inbasin[i];
    }
    

    //continuation_method();

    //printf("%f %f\n", end_xz[0], end_xz[1]);
    //printf("%f %f\n%f %f\n", Jacobian[0][0],Jacobian[0][1],Jacobian[1][0],Jacobian[1][1]);
    
    int numberOfParameters = 4;
    int numberOfDataSets = 2;
    struct DataSetForJson outputdata[numberOfDataSets];
    struct ParameterList parameterList[numberOfParameters];

    outputdata[0].name = "xs";
    outputdata[0].data = newXsInBasin;
    outputdata[0].length = numberInBasin;
    outputdata[1].name = "zs";
    outputdata[1].data = newZsInBasin;
    outputdata[1].length = numberInBasin;

    char title[5] = "Basin";
    char outputPath[100] = "src/variationalPoincare/output/basin.json";
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
        "Basin",
        outputdata,
        numberOfDataSets,
        parameterList,
        numberOfParameters);
    

}

int main()
{
    find_basin();

    return 1;
}
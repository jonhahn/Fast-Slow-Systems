#include "core/systems.h"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void dgauss_elim ( double **a, double *b, int n, int ps )

/*
PURPOSE:   calculate the solution of system of linear algebraic equations
           using Gaussian elimination with back substitution -- all
           calculations are performed in double precision arithmetic

CALLING SEQUENCE:
           dgauss_elim ( a, b, n, ps );

INPUTS:
           a        coefficient matrix
                    type:  **double
           b        right-hand side vector
                    type:  *double
           n        number of equations in system
                    type:  int
           ps       flag indicating which pivoting strategy to use
                         ps == 0:   no pivoting
                         ps == 1;   partial pivoting
                         ps == 2;   scaled partial pivoting
                    type:  int

OUTPUT:
           b        solution vector  (originl contents overwritten)
                    type:  *double
                    
NOTE:
           This routine overwrites the contents of the a matrix.

*/

{
     int pass, row, col, j, temp;
     double rmax,
           ftmp,
           mult,
           sum;

/*
   initialize row pointer array
*/
       
     int pvt[n];
     for ( row = 0; row < n; row++ )
         pvt[row] = row;

/*
   if scaled partial pivoting option was selected,
   initialize scale vector
*/
     double s[n];
     if ( ps == 2 ) {
        for ( row = 0; row < n; row++ ) {
            s[row] = fabs( a[row][0] );
            for ( col = 1; col < n; col++ )
                if ( fabs( a[row][col] ) > s[row] ) 
                   s[row] = fabs( a[row][col] );
        }
     }

/*
   elimination phase
*/
           
     for ( pass = 0; pass < n; pass++ ) {

/*
   perform requested pivoting strategy
   
   even if no pivoting option is requested, still must check for
   zero pivot
*/
     
         if ( ps != 0 ) {
            rmax = ( ps == 1 ? fabs( a[pvt[pass]][pass] ) : 
                               fabs( a[pvt[pass]][pass] ) / s[pvt[pass]] );
            j = pass;
            for ( row = pass+1; row < n; row++ ) { 
                ftmp = ( ps == 1 ? fabs( a[pvt[row]][pass] ) : 
                                   fabs( a[pvt[row]][pass] ) / s[pvt[row]] );
                if ( ftmp > rmax ) {
                   rmax = ftmp;
                   j = row;
                }
            }
            
            if ( j != pass ) {
               temp = pvt[j];
               pvt[j] = pvt[pass];
               pvt[pass] = temp;
            }
         }
         else {
            if ( a[pvt[pass]][pass] == 0.0 ) {
               for ( row = pass+1; row < n; row++ )
                   if ( a[pvt[row]][pass] != 0.0 ) break;
               temp = pvt[row];
               pvt[row] = pvt[pass];
               pvt[pass] = temp;  
            }
         }
     
         for ( row = pass + 1; row < n; row++ ) {
             mult = - a[pvt[row]][pass] / a[pvt[pass]][pass];
             a[pvt[row]][pass] = 0;
             for ( col = pass+1; col < n; col++ )
                 a[pvt[row]][col] += mult * a[pvt[pass]][col];
             b[pvt[row]] += mult * b[pvt[pass]];
         }
     }
     
/*
   solve step
*/

     double x[n];
     x[n-1] = b[pvt[n-1]] / a[pvt[n-1]][n-1];
     for ( row = n-2; row >= 0; row-- ) {
         sum = b[pvt[row]];
         for ( col = row+1; col < n; col++ )
             sum -= x[col] * a[pvt[row]][col];
         x[row] = sum / a[pvt[row]][row];
     } 
     
     for ( row = 0; row < n; row++ )
         b[row] = x[row];
         
     //free(pvt);
     //if ( ps == 2 ) free(s);
     //free(x);
     
}




void invert(   double ** A,
          double ** B,    //input square matrix A
          int n)          //dimension of A
{
    //TODO: there's probably a faster function in systems.c

    int i,j,k;
    double tempVector[n];
    double tempA[n][n];
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            tempA[i][j] = A[i][j];
        }
    }

    double * pointerToTempA[n];
    for(i=0; i<n; i++){
        pointerToTempA[i] = A[i];
    }

    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            tempVector[j] = 0;
        }
        tempVector[i] = 1;

        for(j=0; j<n; j++){
            for(k=0; k<n; k++){
                tempA[j][k] = A[j][k];
            }
        }
        dgauss_elim(pointerToTempA, tempVector, n, 1);

        for(j=0; j<n; j++){
            B[j][i] = tempVector[j];
        }
    }

}


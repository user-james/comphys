#include <stdio.h>
#include <math.h>

/************************************************************************************************/
int ranposn( double x[], double y[], double z[], double L, int N )
{
/* This function places N particles with a uniform random distribution in a cube of side L,
   with -L/2 < x,y,z < L/2
  
Input variables:
  L    =  size of box 
  N    =  number of particles
  

Output Variables:
  (x,y,z) =  positions of particles (cartesian coords)

*/
    int  i;
    static long int iseed = 0;

    double dran1_(long int *);    /*  uniform random number generator function (FORTRAN)  */
    
    
    for(i=0; i< N; i++)
    {
        x[i] = (2*dran1_(&iseed)-1.0)*L/2;
        y[i] = (2*dran1_(&iseed)-1.0)*L/2;
        z[i] = (2*dran1_(&iseed)-1.0)*L/2;
    }
    return 0;
}


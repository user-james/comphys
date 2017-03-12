#include <math.h>

double scalevel(double vx[], double vy[], double vz[], 
                int N, double scale)
{
/*  
   This function scales the velocities of N particles 
   by the factor scale.


Input variables:
  (vx,vy,vz)  =  particle velocities (cartesian coords)
  N       =  number of particles
  scale   =  average energy increment per particle

Output Variables:
  (vx,vy,vz)  =  scaled particle velocities (cartesian coords)

The function returns the change in the total kinetic energy.

*/

  double  KEold, KEnew;
  int i;

  KEold = 0.0;  KEnew = 0.0;

  for (i=0; i<N; i++)
  {
    KEold += (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]) / 2 ;

    vx[i] *= scale; 
    vy[i] *= scale; 
    vz[i] *= scale;

    KEnew += (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]) / 2 ;
  }
  return KEnew - KEold;
}

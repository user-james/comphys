#include <math.h>

double ranveladd(double vx[], double vy[], double vz[], 
                 int N, double dE, long int *iseed)
{
/*  
   This function adds random increments to the velocities of N particles 
   so that the average energy per particle increases by dE.
   (The mean square of velocity increment = 2*dE.)

   This function uses the fortran subroutine dran1 to generate
   uniform pseudorandom numbers on the interval [0,1].

Input variables:
  N       =  number of particles
  dE      =  average energy increment per particle
  iseed   =  random number seed (usually 0, except possibly for first call)

Output Variables:
  (vx,vy,vz)  =  particle velocities (cartesian coords)

The function returns the change in the total kinetic energy of the particles.

*/
  double dran1_(long int *);

  double r, t, KEold, KEnew;
  double pi2, kT2;
  int i;

  pi2 = 8 * atan(1.0);
  kT2 = 4.0/3.0 * dE;
  KEold = 0.0;  KEnew = 0.0;

  for (i=0; i<N; i++)
  {
    KEold += (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]) / 2 ;
/* 
            Box-Muller Transformation 
*/
    r = sqrt(- kT2 * log( dran1_(iseed) ));  
    t = pi2 * dran1_(iseed);
    vx[i] += r * sin(t);
    vy[i] += r * cos(t);

    r = sqrt(- kT2 * log( dran1_(iseed) ));  
    t = pi2 * dran1_(iseed);
    vz[i] += r * sin(t);

    KEnew += (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]) / 2 ;
  }
  return KEnew - KEold;
}

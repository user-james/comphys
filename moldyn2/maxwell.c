#include <math.h>

double maxwell(double vx[], double vy[], double vz[], 
               int N, double kT, int *iseed)
{
/*  
   This function assigns random velocities to N particles 
   from a Maxwell distribution at temperature T. 
   (The mean square velocity = 3*kT.)

   This function uses the fortran subroutine dran1 to generate
   uniform pseudorandom numbers on the interval [0,1].

Input variables:
  N       =  number of particles
  kT      =  temperature times Boltzmann's constant.
  iseed   =  random number seed (usually 0, except for first call)

Output Variables:
  (vx,vy,vz)  =  particle velocities (cartesian coords)

The function returns the total kinetic energy of the particles.

*/
  double dran1_(int *);

  double r, t, KE;
  double pi2, kT2;
  int i;

  pi2 = 8 * atan(1.0);
  kT2 = 2 * kT;
  KE = 0.0;

  for (i=0; i<N; i++)
  {
/* 
            Box-Muller Transformation 
*/
    r = sqrt(- kT2 * log( dran1_(iseed) ));  
    t = pi2 * dran1_(iseed);
    vx[i] = r * sin(t);
    vy[i] = r * cos(t);

    r = sqrt(- kT2 * log( dran1_(iseed) ));  
    t = pi2 * dran1_(iseed);
    vz[i] = r * sin(t);

    KE += (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]) / 2 ;
  }
  return KE;
}

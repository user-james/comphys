#include <stdio.h>
#include <math.h>


int  leapfrog(double  x, double v, int nstep, double dt, double gamma)
{
/*  
      This function implements the damped leap-frog algorithm
      for a particle in one dimension.

      The external force on the particle depends on position only and
      is calculated in function "force".
      The damping force equals -gamma * v and is incorporated directly
      in the integration of the equations of motion.


Input variables:
  x    =  initial position of particle
  v    =  initial velocity of particle
  nstep   =  number of time-steps for integration
  dt      =  length of each time-step
  gamma   =  velocity damping factor (damping force = -gamma * v)

Output Variables:
  x    =  final position of particle
  v    =  final velocity of particle

The function returns 0 for normal completion, 1 for error.

*/
  double step2, gs2, velfac, forfac, F;
  int istep, i;

  double force( double );

/*
   velfac and forfac incorporate the damping into the leap-frog method. 
   When the damping (gamma) is zero, velfac = 1 and forfac = step/2
*/
  step2 = dt/2.0;
  gs2 = gamma*step2;
  velfac = exp( -gs2 );
  if ( abs(gs2) < 0.00001 ) 
    forfac = step2*(1.0 + 0.5*gs2*(1.0-gs2/3.0));
  else
    forfac = (1.0-velfac)/gamma ;

/*
          integrate equations of motion
*/
 
  F = force(x);
  for(istep=0; istep < nstep; istep++)
  {
    v = velfac * v + forfac * F;
    x = x + dt * v;
    F = force(x);
    v = velfac * v + forfac * F;
  }
  return 0;
}

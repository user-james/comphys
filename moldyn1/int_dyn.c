#include <stdio.h>
#include <math.h>

#define NMAX 100
    



int  intdyn(double  x[], double  y[], double  z[],
            double vx[], double vy[], double vz[], 
            int N, int nstep, double dt, double gamma,
            double alpha, double beta, 
            double Fxi[], double Fyi[], double Fzi[], 
            double L, double V0,
            double Fxe[], double Fye[], double Fze[], 
            double *Vint, double *Vext,
            double rij[][NMAX], double *KE_average )
{
/*  
   This function integrates the equations of motion for N particles, 
   interacting via the internal (Lennard-Jones) and external (confining
   box) forces, as calculated in function force. 

   The forces are found by calling function forces.

Input variables:
  (x,y,z)    =  initial positions of particles (cartesian coords)
  (vx,vy,vz) =  initial velocities of particles (cartesian coords)
  N       =  number of particles
  nstep   =  number of time-steps for integration
  dt      =  length of each time-step
  gamma   =  velocity damping factor (damping force = -gamma * v)

  alpha   =  coefficient of 1/r^12
  beta    =  coefficient of -1/r^6
  L       =  size of confining "box"
  V0      =  strength of confining potential  - not used if we are using forces_m

Output Variables:
  (x,y,z)    =  final positions of particles (cartesian coords)
  (vx,vy,vz) =  final velocities of particles (cartesian coords)
  (Fxi,Fyi,Fzi)  =  final internal forces on particles 
  (Fxe,Fye,Fze)  =  final external forces on particles
  Vint  =  final total internal potential energy of system 
  Vext  =  final total external potential energy of system 
  rij   =  final interatomic distances
  KE_average = time-average of system kinetic energy

The function returns 0 for normal completion, 1 for error.

*/
  double step2, gs2, velfac, forfac;
  int istep, i;

  int forces( double *,double *,double *, int, 
              double, double, double, double,
              double *, double *, double *, 
              double *, double *, double *,
              double *, double *, double rij[][NMAX]);
  int bin_rij(double rij[][NMAX], int, int);

  if( N > NMAX )
  {
    printf("%d is too many particles for forces. NMAX = %d\n", N, NMAX);
    return 1;
  }
/*
   velfac and forfac incorporate the damping into the leap-frog method. 
   When the damping (gamma) is zero, velfac = 1 and forfac = step/2
*/
  step2 = dt/2.0;
  gs2 = gamma*step2;
  velfac = exp( -gs2 );
  if ( fabs(gs2) < 1.0e-4 ) 
        forfac = step2*(1.0 - 0.5*gs2*(1.0-gs2/3.0));
  else
        forfac = (1.0-velfac)/gamma ;
    
  *KE_average = 0.0;

/*
          integrate equations of motion
*/
 
  if( forces(x,y,z, N, alpha, beta, L, V0,
             Fxi, Fyi, Fzi, Fxe, Fye, Fze, Vint,Vext, rij) != 0)
    return 1;    /*
    if( forces_m(x,y,z, N, vx,vy,vz, alpha, beta, L, V0,
               Fxi, Fyi, Fzi, Fxe, Fye, Fze, Vint,Vext, rij) != 0)
        return 1;*/
  for(istep=0; istep < nstep; istep++)
  {
    for(i=0; i < N; i++)
    {
      vx[i] = velfac * vx[i] + forfac * (Fxi[i]+Fxe[i]);
      vy[i] = velfac * vy[i] + forfac * (Fyi[i]+Fye[i]);
      vz[i] = velfac * vz[i] + forfac * (Fzi[i]+Fze[i]);
      x[i] = x[i] + dt * vx[i];
      y[i] = y[i] + dt * vy[i];
      z[i] = z[i] + dt * vz[i];
    }
    if( forces(x,y,z, N, alpha, beta, L, V0,
       Fxi, Fyi, Fzi, Fxe, Fye, Fze, Vint,Vext, rij) != 0)
       return 1;     
/*      if( forces_m(x,y,z, N, vx,vy,vz, alpha, beta, L, V0,
                   Fxi, Fyi, Fzi, Fxe, Fye, Fze, Vint,Vext, rij) != 0)
          return 1;*/
      for(i=0; i < N; i++)
      {
          vx[i] = velfac * vx[i] + forfac * (Fxi[i]+Fxe[i]);
          vy[i] = velfac * vy[i] + forfac * (Fyi[i]+Fye[i]);
          vz[i] = velfac * vz[i] + forfac * (Fzi[i]+Fze[i]);
          *KE_average += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
      }
      bin_rij(rij, N, 1);     /* accumulate interatomic distance statistics */
  }
    *KE_average = *KE_average / (2*nstep);
      return 0;
}

 
           
           
           

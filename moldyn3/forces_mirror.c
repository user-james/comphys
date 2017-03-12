#include <stdio.h>
#include <math.h>

#define NMAX 100

int  forces_m(double x[], double y[], double z[], int N, 
              double vx[], double vy[], double vz[],
              double alpha, double beta, double L, double V0,
              double Fxi[], double Fyi[], double Fzi[], 
              double Fxe[], double Fye[], double Fze[], 
              double *Vint, double *Vext, double rij[][NMAX] )
{
/*  
   This function calculates the forces on N particles at positions
   (x,y,z) for the Lennard-Jones interaction,
                   V(r) = alpha/r^12 - beta/r^6 . 

   It also reverses the velocities and mirrors the position if the particle 
   coordinates have moved outside the box, -L/2 < x,y,z, < L/2.
 
Input variables:
  (x,y,z) =  positions of particles (cartesian coords)
  N       =  number of particles
  (vx,vy,vz) =  velocities of particles (cartesian coords)
  alpha   =  coefficient of 1/r^12
  beta    =  coefficient of -1/r^6
  L       =  size of confining "box"
  V0      =  strength of confining potential

Output Variables:
  (Fxi,Fyi,Fzi)  =  internal force on particle (cartesian coords)
  (Fxe,Fye,Fze)  =  external force on particle (cartesian coords)
  (x,y,z) =  positions of particles (cartesian coords)  - mirrored if coordinate has exceeded bounds of box
  (vx,vy,vz) =  velocities of particles (cartesian coords) - reversed if coordinate has exceeded bounds of box
  Vint  =  total internal potential energy of system
  Vext  =  total external potential energy of system
  rij   =  interatomic distances

function returns 0 for normal completion, 1 for error

*/
  double xij, yij, zij;
  double d2, d6, Fr;
  int i,j;

  double mirror(double, double *);

  if( N > NMAX )
  {
    printf("%d is too many particles for forces. NMAX = %d\n", N, NMAX);
    return 1;
  }

  *Vext = 0.0;  *Vint = 0.0;

  for(i=0; i<N; i++)
  {
/* reverse velocity if coordinate has reached wall  */

      vx[i] = vx[i] * mirror(L, x+i);
      vy[i] = vy[i] * mirror(L, y+i);
      vz[i] = vz[i] * mirror(L, z+i);
   
/* zero the internal and external forces  */
    
    Fxi[i] = 0.0;   Fyi[i] = 0.0;   Fzi[i] = 0.0;
    Fxe[i] = 0.0;   Fye[i] = 0.0;   Fze[i] = 0.0;
  }

/* accumulate the internal forces  */

  for(i=0; i<N; i++)
  {
    rij[i][i] = 0.0;

    for(j=0; j<i; j++)
    {
      xij = x[i] - x[j];     /*  relative position vectors and distances */
      yij = y[i] - y[j];
      zij = z[i] - z[j];
      d2 = xij*xij + yij*yij + zij*zij; 
      rij[j][i] = ( rij[i][j] = sqrt(d2) );

      d6 = d2*d2*d2; 
      Fr = (6*beta - 12*alpha/d6)/(d2*d6);       /* Lennard-Jones */
      *Vint += (alpha/d6 - beta)/d6;

      Fxi[i] += -xij * Fr;
      Fyi[i] += -yij * Fr;
      Fzi[i] += -zij * Fr;

      Fxi[j] +=  xij * Fr;
      Fyi[j] +=  yij * Fr;
      Fzi[j] +=  zij * Fr;
    }
  }
  return 0;
}
double mirror(double L, double *x )
{
/*
This function "bounces" the coordinate x off walls at +-L/2
 i.e. if |x| > L/2, x is mirrored in the wall to where it would be if the velocity had reversed when it reached x = +- L/2
 
The function returns +1.0 if the particle has not collided with wall, -1.0 if it has

*/
    if( *x > L/2 ){ 
        *x = L - *x;
        return -1.0;
    }
    if( *x < -L/2 ) {
        *x = -L - *x;
        return -1.0;
    }
    return 1.0;
}

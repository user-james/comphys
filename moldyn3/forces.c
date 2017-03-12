#include <stdio.h>
#include <math.h>

#define NMAX 100

int  forces(double x[], double y[], double z[], int N, 
            double alpha, double beta, double L, double V0,
            double Fxi[], double Fyi[], double Fzi[], 
            double Fxe[], double Fye[], double Fze[], 
            double *Vint, double *Vext, double rij[][NMAX] )
{
/*  
   This function calculates the forces on N particles at positions
   (x,y,z) for the Lennard-Jones interaction,
                   V(r) = alpha/r^12 - beta/r^6 . 

   It also calculates the forces on each particle from an external 
   confining "box" potential,
   Vext(x,y,z) = V0 * [exp{(2x/L)^4} + exp{(2x/L)^4} + exp{(2x/L)^4} -3].
   The external forces are found by calling function boxpot.

Input variables:
  (x,y,z) =  positions of particles (cartesian coords)
  N       =  number of particles
  alpha   =  coefficient of 1/r^12
  beta    =  coefficient of -1/r^6
  L       =  size of confining "box"
  V0      =  strength of confining potential

Output Variables:
  (Fxi,Fyi,Fzi)  =  internal force on particle (cartesian coords)
  (Fxe,Fye,Fze)  =  external force on particle (cartesian coords)
  Vint  =  total internal potential energy of system
  Vext  =  total external potential energy of system
  rij   =  interatomic distances

function returns 0 for normal completion, 1 for error

*/
  double xij, yij, zij;
  double d2, d6, Fr;
  int i,j;

  double boxpot(double, double, double, double *);

  if( N > NMAX )
  {
    printf("%d is too many particles for forces. NMAX = %d\n", N, NMAX);
    return 1;
  }

  *Vext = 0.0;  *Vint = 0.0;

  for(i=0; i<N; i++)
  {
/* find the external forces  */

    *Vext += boxpot(L,V0, x[i], Fxe+i );
    *Vext += boxpot(L,V0, y[i], Fye+i );
    *Vext += boxpot(L,V0, z[i], Fze+i );

/* zero the internal forces  */
    
    Fxi[i] = 0.0;   Fyi[i] = 0.0;   Fzi[i] = 0.0;
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
double boxpot(double L, double V0, double x, double *Force )
{
/*
   Calculate the x component of external force on particle at x 
   due to a confining "box" potential,
                  V(x) = V0 * [exp{(2x/L)^4} - 1].
return variables:
  Force = component of the force

function returns the value of the potential at x

*/
  double x3, ex;  

  x = 2*x/L; 
  x3 = x*x*x;
  ex = exp(x3*x);

  *Force =  -8 * V0 * x3 * ex / L;
  return (V0*(ex - 1.0));
}

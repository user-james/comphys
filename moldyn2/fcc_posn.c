#include <stdio.h>
#include <math.h>

#define NSHMAX 5

int  fccposn(double x[], double y[], double z[], 
            double alpha, double beta, int Nshell) 
{
/*  
   This function places particles on a fragment of an fcc lattice so 
   that the nearest neighbour distance is 
                   r_nn = (2*alpha/beta)^(1/6),
   i.e. at the minimum of the Lennard-Jones potential,
                   V(r) = alpha/r^12 - beta/r^6 . 

   The atoms are placed in concentric "shells" about a central atom 
   which is at (0,0,0). For the moment, only up to 5 shells are allowed, 
   so that the maximum number of atoms is 79.

Input variables:
  alpha   =  coefficient of 1/r^12
  beta    =  coefficient of -1/r^6
  Nshell  =  number of shells (must be < 6)

Output Variables:
  (x,y,z) =  positions of particles (cartesian coords)

The function returns the total number of particles in the cluster.

*/
  double a;
  int i, atom;

  if( Nshell > NSHMAX )
  {
    printf("%d is too many shells for fccposn. NSHMAX = %d\n", Nshell, NSHMAX);
    return 0;
  }

  a = sqrt(2.0) * pow((2*alpha/beta),(1.0/6.0));   /* fcc lattice constant */

  atom = 0;
  for(i=0; i<Nshell+1; i++)
  {
    switch (i)
    {
      case 0:          /* central atom */

        x[atom] = 0.0;  y[atom] = 0.0;   z[atom++] = 0.0; 
        break;

      case 1:          /* first shell */

        x[atom] =  0.0;  y[atom] =  a/2;   z[atom++] =  a/2; 
        x[atom] =  0.0;  y[atom] = -a/2;   z[atom++] = -a/2; 
        x[atom] =  a/2;  y[atom] =  0.0;   z[atom++] =  a/2; 
        x[atom] = -a/2;  y[atom] =  0.0;   z[atom++] = -a/2; 
        x[atom] =  a/2;  y[atom] =  a/2;   z[atom++] =  0.0; 
        x[atom] = -a/2;  y[atom] = -a/2;   z[atom++] =  0.0; 
        x[atom] =  0.0;  y[atom] = -a/2;   z[atom++] =  a/2; 
        x[atom] =  0.0;  y[atom] =  a/2;   z[atom++] = -a/2; 
        x[atom] = -a/2;  y[atom] =  0.0;   z[atom++] =  a/2; 
        x[atom] =  a/2;  y[atom] =  0.0;   z[atom++] = -a/2; 
        x[atom] = -a/2;  y[atom] =  a/2;   z[atom++] =  0.0; 
        x[atom] =  a/2;  y[atom] = -a/2;   z[atom++] =  0.0; 
        break;

      case 2:          /* second shell */

        x[atom] =   a ;  y[atom] =  0.0;   z[atom++] =  0.0; 
        x[atom] =  -a ;  y[atom] =  0.0;   z[atom++] =  0.0; 
        x[atom] =  0.0;  y[atom] =   a ;   z[atom++] =  0.0; 
        x[atom] =  0.0;  y[atom] =  -a ;   z[atom++] =  0.0; 
        x[atom] =  0.0;  y[atom] =  0.0;   z[atom++] =   a ; 
        x[atom] =  0.0;  y[atom] =  0.0;   z[atom++] =  -a ; 
        break;

      case 3:          /* third shell */

        x[atom] =   a ;  y[atom] =  a/2;   z[atom++] =  a/2; 
        x[atom] =  -a ;  y[atom] = -a/2;   z[atom++] = -a/2; 
        x[atom] =   a ;  y[atom] = -a/2;   z[atom++] =  a/2; 
        x[atom] =  -a ;  y[atom] =  a/2;   z[atom++] = -a/2; 
        x[atom] =   a ;  y[atom] =  a/2;   z[atom++] = -a/2; 
        x[atom] =  -a ;  y[atom] = -a/2;   z[atom++] =  a/2; 
        x[atom] =   a ;  y[atom] = -a/2;   z[atom++] = -a/2; 
        x[atom] =  -a ;  y[atom] =  a/2;   z[atom++] =  a/2; 

        x[atom] =  a/2;  y[atom] =   a ;   z[atom++] =  a/2; 
        x[atom] = -a/2;  y[atom] =  -a ;   z[atom++] = -a/2; 
        x[atom] = -a/2;  y[atom] =   a ;   z[atom++] =  a/2; 
        x[atom] =  a/2;  y[atom] =  -a ;   z[atom++] = -a/2; 
        x[atom] =  a/2;  y[atom] =   a ;   z[atom++] = -a/2; 
        x[atom] = -a/2;  y[atom] =  -a ;   z[atom++] =  a/2; 
        x[atom] = -a/2;  y[atom] =   a ;   z[atom++] = -a/2; 
        x[atom] =  a/2;  y[atom] =  -a ;   z[atom++] =  a/2; 
    
        x[atom] =  a/2;  y[atom] =  a/2;   z[atom++] =   a ; 
        x[atom] = -a/2;  y[atom] = -a/2;   z[atom++] =  -a ; 
        x[atom] = -a/2;  y[atom] =  a/2;   z[atom++] =   a ; 
        x[atom] =  a/2;  y[atom] = -a/2;   z[atom++] =  -a ; 
        x[atom] =  a/2;  y[atom] = -a/2;   z[atom++] =   a ; 
        x[atom] = -a/2;  y[atom] =  a/2;   z[atom++] =  -a ; 
        x[atom] = -a/2;  y[atom] = -a/2;   z[atom++] =   a ; 
        x[atom] =  a/2;  y[atom] =  a/2;   z[atom++] =  -a ; 
        break;

      case 4:          /* fourth shell */

        x[atom] =  0.0;  y[atom] =   a ;   z[atom++] =   a ; 
        x[atom] =  0.0;  y[atom] =  -a ;   z[atom++] =  -a ; 
        x[atom] =   a ;  y[atom] =  0.0;   z[atom++] =   a ; 
        x[atom] =  -a ;  y[atom] =  0.0;   z[atom++] =  -a ; 
        x[atom] =   a ;  y[atom] =   a ;   z[atom++] =  0.0;  
        x[atom] =  -a ;  y[atom] =  -a ;   z[atom++] =  0.0;  
        x[atom] =  0.0;  y[atom] =  -a ;   z[atom++] =   a ; 
        x[atom] =  0.0;  y[atom] =   a ;   z[atom++] =  -a ; 
        x[atom] =  -a ;  y[atom] =  0.0;   z[atom++] =   a ; 
        x[atom] =   a ;  y[atom] =  0.0;   z[atom++] =  -a ; 
        x[atom] =  -a ;  y[atom] =   a ;   z[atom++] =  0.0;  
        x[atom] =   a ;  y[atom] =  -a ;   z[atom++] =  0.0;  
        break;

      case 5:          /* fifth shell */

        x[atom] =  3*a/2;  y[atom] =  0.0;   z[atom++] =  a/2; 
        x[atom] = -3*a/2;  y[atom] =  0.0;   z[atom++] = -a/2; 
        x[atom] =  3*a/2;  y[atom] =  0.0;   z[atom++] = -a/2; 
        x[atom] = -3*a/2;  y[atom] =  0.0;   z[atom++] =  a/2; 
        x[atom] =  3*a/2;  y[atom] =  a/2;   z[atom++] =  0.0; 
        x[atom] = -3*a/2;  y[atom] = -a/2;   z[atom++] =  0.0; 
        x[atom] =  3*a/2;  y[atom] = -a/2;   z[atom++] =  0.0; 
        x[atom] = -3*a/2;  y[atom] =  a/2;   z[atom++] =  0.0; 

        x[atom] =  0.0;  y[atom] =  3*a/2;   z[atom++] =  a/2; 
        x[atom] =  0.0;  y[atom] = -3*a/2;   z[atom++] = -a/2; 
        x[atom] =  0.0;  y[atom] =  3*a/2;   z[atom++] = -a/2; 
        x[atom] =  0.0;  y[atom] = -3*a/2;   z[atom++] =  a/2; 
        x[atom] =  a/2;  y[atom] =  3*a/2;   z[atom++] =  0.0; 
        x[atom] = -a/2;  y[atom] = -3*a/2;   z[atom++] =  0.0; 
        x[atom] = -a/2;  y[atom] =  3*a/2;   z[atom++] =  0.0; 
        x[atom] =  a/2;  y[atom] = -3*a/2;   z[atom++] =  0.0; 

        x[atom] =  0.0;  y[atom] =  a/2;   z[atom++] =  3*a/2; 
        x[atom] =  0.0;  y[atom] = -a/2;   z[atom++] = -3*a/2; 
        x[atom] =  0.0;  y[atom] = -a/2;   z[atom++] =  3*a/2; 
        x[atom] =  0.0;  y[atom] =  a/2;   z[atom++] = -3*a/2; 
        x[atom] =  a/2;  y[atom] =  0.0;   z[atom++] =  3*a/2; 
        x[atom] = -a/2;  y[atom] =  0.0;   z[atom++] = -3*a/2; 
        x[atom] = -a/2;  y[atom] =  0.0;   z[atom++] =  3*a/2; 
        x[atom] =  a/2;  y[atom] =  0.0;   z[atom++] = -3*a/2; 
        break;
    }   
  }
  return atom;
}

#include <stdio.h>
#include <math.h>

#define NRMAX 500   /* maximum number of grid points used for finite difference solution of Schrodinger equation */

FILE *wavefn;

int main( )
{
    /*  This program finds the energy eigenvalues and wave functions of a particle in a 1-dimensional harmonic potential
     The dimensionless Schrodinger equation for this system is:
     
         - d^2 psi / dr^2  + psi * r^2/2 =  E psi
     
     where psi = wavefunction * r  and E = energy eigenvalue.
     The boundary conditions are: 
          psi(r) = 0 at r = rmin and r = rmax
     
             */
    
    int n_wf = 6;
    int nrmax = NRMAX;     
    
    double rmax = 10.0, rmin= -10.0;   /*  default parameters for grid */
    
    
    double r[NRMAX], rab[NRMAX], u[NRMAX];
    int nr=NRMAX;
    
    double sd[NRMAX], e[NRMAX], psi[NRMAX][NRMAX];
    
    int ierr, i, j;
    double rad, r_av, r2_av, w_av;
    char wavefile[30];
    
    

    
    for(j=0; j < nr; j++)                   /* generate uniform grid */
    {
        r[j] = rmin + (rmax-rmin)*(double)j/(double)nr; 
        rab[j] = (rmax-rmin)/(double)nr;
    }

    printf("uniform grid: min %lf, max %lf, with %d grid points \n", rmin, rmax, nr);
    
    
    for( i=1; i< nr ; i++ )                  /* write out the potential on the grid  */
    {
        u[i] = r[i]*r[i] ;
    }
    
    ierr = findiff( r, rab, u, nr, nrmax, sd, e, psi );   /* find the wave functions and energies */
    

    for(i=0; i < n_wf; i++)                 /* write out the first n_wf wave functions and energies */
    {
        printf( "energy %d = %lf \n", i, e[i]  );
        sprintf(wavefile, "SHO_wavefn_%d", i);
        wavefn = fopen(wavefile, "w");
        r_av = 0.0;
        r2_av = 0.0;
        w_av = 0.0;
        for(j=0; j < nr; j++)
        {
            fprintf(wavefn, " %lf \t %lf \n", r[j], psi[i][j]);
            r_av += psi[i][j]*psi[i][j] * r[j] * rab[j];
            r2_av += psi[i][j]*psi[i][j] * r[j] * r[j] * rab[j];
            w_av += psi[i][j]*psi[i][j] * rab[j] ;
        }
        printf( " average coord for state %d is %lf, coord^2 is %lf and weight is %lf \n", i, r_av/w_av, r2_av/w_av, w_av );
    }
    
}



int findiff( double r[], double rab[], double u[], int nr, int nrmax, 
             double sd[], double e[], double psi[][NRMAX] )
{
/********************************************************************************  
    This function finds the solutions of the radial schrodinger equation 
    on a non-uniform radial grid, using the finite difference method and 
    finding the eigenvalues and eigenvectors of the appropriate tridiagonal 
    matrix using the standard (Fortran) LAPACK subroutine, tql2.

    The boundary condition on the wavefunctions is that 
    psi(r) = 0 at r = r[0] and at r = r[nr-1].

Input:

r[i]   =  radii      (usually generated in subroutine "grid")
rab[i] =  dr/di      (also from "grid")
u[i]   =  potential, 2*V(r) + l(l+1)/r^2, at r[i]
nr     =  number of grid points used
nrmax  =  leading dimension of array psi (must equal parameter NRMAX)

sd[ ]   =   a work space for the matrix diagonalization routine 
            (of length at least nr)


Output:

e[j]      =  jth eigenvalue  (for j=0,...,nr-3)
psi[j][i] =  jth radial function * r at r[i]

Note: the radial wavefunctions satisfy the orthonormalization condition

   sum_{i=0}^{nr-1} psi[j][i]*psi[j'][i] rab[i]  =   1  for j=j'
                                                 =   0  for j not = j'

This function returns the error code from tql2 (zero for successful completion).

********************************************************************************/

  int n, i, j, ierr;
/*
    check consistency of dimensions of the array psi
*/
  if(nrmax != NRMAX)
  {
    ierr = 1;
    printf( "dimensions of psi do not match in call to findiff \n" );
    return ierr;
   }

/*  
    construct the diagonal and subdiagonal elements of the matrix  
    and construct the unit matrix in psi
*/

  for(i=1; i < nr-1; i++)
  {
    e[i-1] = u[i] + ( 2.0/(rab[i]+rab[i+1]) + 2.0/(rab[i]+rab[i-1]) )/rab[i];

    sd[i] = -( 2.0/(rab[i]+rab[i+1]) ) / sqrt(rab[i]*rab[i+1]);

    for(j=0; j < nr-2; j++) 
    {
      psi[j][i-1] = 0.0;
    }
    psi[i-1][i-1] = 1.0;
  }
/*
             find nr-2 eigenvalues and eigenvectors 
*/
  n = nr-2;
  tql2_( &nrmax, &n, e, sd, psi, &ierr );    
  if( ierr != 0 )
  {
    printf(" diagonalization did not succeed. error code = %d \n", ierr);
    return ierr;
  }

  for(i=0; i < nr-2; i++)
  {
    for(j=nr-3; j > 0; j--) 
    {
      psi[i][j+1] = psi[i][j] / sqrt(rab[j+1]) ;
    }
    psi[i][0] = 0.0;
    psi[i][nr-1] = 0.0;
  } 
  return ierr;
}

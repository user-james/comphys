#include <stdio.h>
#include <math.h>

#define NRMAX 500

FILE *wavefn;
FILE *energies;

int main( )
{
    /*  This program finds the energy eigenvalues and wave functions of the diatomic vibrations
        for a Lennard-Jones interaction with dimensionless interaction parameter, A. 
        The dimensionless radial Schroedinger equation for this system is:
     
         - d^2 psi / dr^2  + 4 A ( 1/r^12 - 1/r^6 ) psi + l(l+1)/ r^2 psi  =  E psi
     
        where psi = wavefunction * r  and E = energy eigenvalue.
     
        Typical values of A are 5.556 for He, 115.48 for Ne, and 1205.43 for Ar


        OUTPUT:
          n_wf state energies are printed to standard output (usually the terminal) -> this has now been modified
                                                                                        to print to the file 
                                                                                        energies_<atom symbol>
          wave functions are output (one per file) in ASCII format to files, named "wavefn_0", "wavefn_1", etc.
     */
    
   
    double A = 115.48;  /*   Helium = 5.556, Neon = 115.48, Argon = 1205.43   */
    
    int l = 0;          /* angular momentum quantum number l */
    
    int n_wf = 10;
    int nrmax = NRMAX;     
    
    double rmax = 0.0, aa= 0.0, bb = 0.0;    /*  default parameters for log-radial grid */
    double rmin = 0.5;                         /*  minimum radius of radial grid  */
    
    
    double r[NRMAX], rab[NRMAX], u[NRMAX],  a, b;
    int nr;
    
    double sd[NRMAX], e[NRMAX], psi[NRMAX][NRMAX];
    
    int ierr, i, j;
    double rad;
    char wavefile[50];
    char efile[50];
    char atom[10];

    printf("Enter Atomic Symbol (he, ar, ne): ");
    fgets(atom, 3, stdin);

    if(strcmp(atom,  "he") == 0){
        A = 5.556;
    }
    else if(strcmp(atom,  "ne") == 0){
        A = 115.48;
    }
    else if(strcmp(atom,  "ar") == 0){
        A = 1205.43;
    }
    else{
        printf("Atom symbol unknown, defaulted to Ne\n");
    }

    
    
    
    grid_(&rmax, &aa, &bb, &nrmax, r,rab, &nr, &a, &b);   /* generate the standard log-radial grid */
    for(i=0; i<nr; i++)
    {
        r[i] = r[i] + rmin;      /* add offset to all radii, to avoid divergence of LJ potnl for r<<1 */
    }
    rmax = rmax + rmin;
    printf("radial grid: max radius %lf, with %d grid points \n", rmax, nr);
    
    
    for( i=0; i< nr ; i++ )                  /* write out the potential on the grid  */
    {
        rad = r[i];
        u[i] = 4*A * ( pow(rad,-12) - pow(rad,-6) ) + l*(l+1)/pow(rad,2);
    }
    ierr = findiff( r, rab, u, nr, nrmax, sd, e, psi );   /* find the wave functions and energies */
    
    sprintf(efile, "../my_files/ex2/energy_%s", atom);
    energies = fopen(efile, "w");

    /* write out the first n_wf wave functions and energies */
    for(i=0; i < n_wf; i++)
    {
        fprintf(energies, "%d \t %lf \n", i, e[i]/A  );
        sprintf(wavefile, "../my_files/ex2/%s_wavefn_%d",atom,  i);
        wavefn = fopen(wavefile, "w");
        for(j=0; j < nr; j++)
        {
            fprintf(wavefn, " %lf \t %lf \n", r[j], psi[i][j]);
        }
    }
    fclose(energies);
    
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

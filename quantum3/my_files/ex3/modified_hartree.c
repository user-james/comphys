#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NRMAX 500  /*  maximum number of radial grid points allowed  */

FILE *charge;
FILE *energy;

int main(int argc, char** argv)
{
    /*  This program finds the energy eigenvalues and wave functions of the electrons in
        the helium atom using the Hartree approximation for the mean-field of the electronic charge. 
        The dimensionless radial Schroedinger equation for this system is:
     
         - d^2 psi / dr^2  - 4/r psi + l(l+1)/ r^2 psi + V_h(r) psi + V_xc(r) psi =  E psi
     
        where psi = wavefunction * r  and E = energy eigenvalue.
     
        V_h is the "Hartree potential" 
        - i.e. the classical electrostatic potential due to the charge of the lowest energy state, psi_0,
        which is doubly-occupied (with spin-up and spin-down electrons).
     
        V_xc is the "exchange-correlation potential" of density functional theory. 
        [W. Kohn & L. J. Sham, Physical Review 140, A1133–A1138 (1965). 
         See also http://en.wikipedia.org/wiki/Kohn–Sham_equations]
     
    */
    
    double Z_nuc = 2;   /* atomic number  */
    int l = 0;          /*  angular momentum quantum number l for occupied state */
    int nt = 10;        /*  number of iterations of self-consistent field  */
    
    int n_wf = 2;
    int nrmax = NRMAX;     
    
    double rmax = 0.0, aa= 0.0, bb = 0.0;    /*  default parameters for log-radial grid */
    
    
    double r[NRMAX], rab[NRMAX], u[NRMAX],  a, b;
    int nr;
    
    double sd[NRMAX], e[NRMAX], psi[NRMAX][NRMAX], prev[NRMAX];                 // PREV WILL STORE PREVIOUS ITERATIONS CHARGE DENSITY CONTRIBUTION
    double Q[NRMAX], v_hartree[NRMAX], v_exc[NRMAX];                 
    
    int it, ierr, i, j;
    double external_energy, v_average, vxc_average, exc_average;
    double cdd, cdu, rad, vxcd, vxcu, exct;
    char rho_file[50];
    char e_file[40];

    float ratio = 1;

    if(argc == 2){
        ratio = strtof(argv[1], NULL);          /* Note: if the argument entered isn't a number ratio will default to zero */
        if(!(ratio<=1 && ratio >=0)){
            printf("Ratio entered outside range (0-1)\nDefaulted to ratio=1\n");                /* Handles error if ratio is outside [0, 1]  */
            ratio =1;
        }
        printf("Ratio = %.2f\n", ratio);
    }
    else{
        printf("Number of arguments entered incorrect\n");
    }    
    
    grid_(&rmax, &aa, &bb, &nrmax, r,rab, &nr, &a, &b);   /* generate the standard log-radial grid */

    printf("radial grid: max radius %lf, with %d grid points \n", rmax, nr);
    
    for( i=nr-1; i>0 ; i-- )               /* set the Hartree potential initially to zero  */
    {
        v_hartree[i] = 0.0;     v_exc[i] = 0.0;
    }
   

    sprintf(e_file, "./energies_%.2f", ratio );
    energy = fopen(e_file, "w");
    for( it=0; it<nt; it++ )               /* iterate the self-consistent field equations nt times */
    {
        for( i=1; i< nr ; i++ )                                /* calculate the potential on the grid          */
        {
            u[i] = -2*Z_nuc/r[i] + l*(l+1)/(r[i]*r[i]) + v_hartree[i] + v_exc[i];
        }
    
        ierr = findiff( r, rab, u, nr, nrmax, sd, e, psi );   /* find the wave functions and energies          */
        
        v_average = 0.0;                                      /* calculate the average of the Hartree interaction  */
        vxc_average = 0.0;
        for( i=1; i< nr ; i++ )
        {
            v_average = v_average + v_hartree[i]*2*psi[0][i]*psi[0][i]*rab[i]; 
            vxc_average = vxc_average + v_exc[i]*2*psi[0][i]*psi[0][i]*rab[i];
        }
        external_energy = 2*e[0] - v_average - vxc_average ;                /* kinetic energy + electron-ion energy  */

    
        Q[0] = 0.0;                                           /* find the updated electronic charge Q inside each radius */



        /* MAIN MODIFICATION*/

        for( i=1; i< nr ; i++ )
        {
            if(it==0){
                Q[i] = Q[i-1]+ 2*psi[0][i]*psi[0][i]*rab[i];      /* IF IT'S THE FIRST ITERATION THE CHARGE DENSITY IS CALCULATED NORMALLY  */
            }
            else{
                /* FOR OTHER ITERATIONS THE CHARGE IS CALCULATED AS A RATIO OF THE PREVIOUS ITERATION (STORED AT END OF PROGRAM) AND THE CURRENT ONE  */
                Q[i] = Q[i-1]+ (ratio*2*psi[0][i]*psi[0][i] + (1-ratio)*prev[i])*rab[i];        
            }
        }
    
        
        
        
        v_hartree[nr-1] = Q[nr-1]/r[nr-1];                    /* find the Hartree potential, starting from maximum radius  */
        for( i=nr-1; i>0 ; i-- )  
        {
            v_hartree[i-1] = v_hartree[i] + Q[i-1]/(r[i]*r[i])*rab[i];
        }
        v_average = 0.0;                                     /* calculate the average of the Hartree interaction  */
        for( i=1; i< nr ; i++ )
        {
            v_average = v_average + v_hartree[i]*2*psi[0][i]*psi[0][i]*rab[i];     
        }
              
        exc_average = 0.0;
        for( i=1; i< nr ; i++ )                                /* calculate the exchange-correlation potential on the grid          */
        {
            cdd = ( cdu = psi[0][i]*psi[0][i] ) ;    rad = r[i];
            vexcorr_(&cdd, &cdu, &rad, &vxcd, &vxcu, &exct);
            //vxcd = ( vxcu = ( exct = 0.0 ) ) ;                 /*  !!!! TEST without DFT exchange-correlation terms !!!!  */
            v_exc[i] = vxcd;
            exc_average = exc_average + exct*(cdd+cdu)*rab[i];
        }
        fprintf(energy, "%d \t %lf \n", it, external_energy + v_average/2 + exc_average) ;

        
                                                             /* write out the electronic charge density for this iteration */
        sprintf(rho_file, "./charge_density_iter_%d", it);
        charge = fopen(rho_file, "w");
        for(j=0; j < nr; j++)
        {
            fprintf(charge, " %lf \t %lf \n", r[j], 2*psi[0][j]*psi[0][j] );
            prev[j] = 2*psi[0][j]*psi[0][j];                            // STORE THE PREVIOUS CHARGE DENSITIES FOR NEXT ITERATION
        }

    }

    fclose(energy);
    fclose(charge);
        
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

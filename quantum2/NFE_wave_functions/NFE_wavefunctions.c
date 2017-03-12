#include <stdio.h>
#include <math.h>

#define NGMAX 50      /* maximum number of plane-waves used in representation of wave functions */
#define Nk    40      /* number of values of k in Brillouin zone, for which the bands are calculated */

FILE *bands;

int main( )
{
    /*  This program finds the energy eigenvalues and wave functions of a particle in a 1-dimensional periodic potential
     
                 V(x) = V1 cos(2 pi x/a) + V2 cos(4 pi x/a),   where a is the "lattice constant",
     
        using the momentum (or "plane-wave") representation (equivalent to expanding the wave function in a Fourier series). 
        From Bloch's theorem, the wavefunctions are of the form:
               psi_k(x)  =  exp[ikx] u_k(x),     where u_k(x) is a periodic function, with u(x+a) = u(x).
        Thus,
               psi_k(x)  =  exp[ikx] sum_j c_j exp[i G_j x],    (1)  
        where 
               G_j = j 2 pi/a,  for j = +-1, +-2, +-3, ... are the reciprocal lattice vectors.
     
        The matrix elements of the Hamiltonian for each value of k are:
     
               H_{j,j} = (k+G_j)^2/2;  H_{j,j+1} = V/2 = H_{j+1,j}              [Units: particle mass = 1, hbar = 1].
        
     
        We vary k from -pi/a to pi/a  (i.e. k is in the first Brillouin zone). 
        For each k, the eigenvectors of H = {c_j}, define the wave function in Eq. (1) above, 
        different eigenvectors of H corresponding to different energy bands.
     
     INPUT:
     
       Input parameters (lattice constant, number of k-points, number of Fourier components in wavefunctions, V1 and V2)
       are currently set at complilation of the code - see below and 'define' statements above.
     
     OUTPUT:
     
       The band energies are output (in ASCII text) to file "NFE_bands" for momenta k between -pi/a and +pi/a 
       in the following format:
          ka/(2 pi), E_1(k), E_2(k), ... , E_{n_wf}(k)
       
       The function 'plot' plots the first four bands to the screen and to a postscript file, 'NFE_band_structure.ps', using gnuplot.
     
        */
    
    int n_wf = 6;         /*  number of bands  */
    long int ngmax = NGMAX;     
    long int ng = NGMAX;
    
    double a = 1.0;       /* lattice constant */
    double V1 = 1.5, V2 = 3.5;       /* amplitude of sinusoidal potential components */
    
    double G[NGMAX];
    double H[NGMAX][NGMAX], e[NGMAX], psiR[NGMAX][NGMAX], psiI[NGMAX][NGMAX];   
    
    double ES[NGMAX], TAU[2*NGMAX], GRA[NGMAX], GIA[NGMAX];   /* scratch space arrays for matrix diagonalization routines */

    int ierr, i, j;
    double pi = acos(-1.0e0);
    double k;
    
    for(i=0; i < ng ; i++)               /*  define reciprocal lattice vectors */
    {
        G[i] = (2*i-ng)*pi/a; 
    }
    
    bands = fopen("NFE_bands", "w");     /*  open file for output of band energies */
    
    for(k=-pi/a; k <= pi/a; k = k + 2*pi/(Nk*a))
    {
    /*  
     construct the diagonal and subdiagonal (real-part of H) elements of the Hamiltonian matrix  
     */
    
        for(i=0; i < ng; i++)
        {
            for(j=0; j < ng; j++) 
            {
              H[i][j] = 0.0;
            }
            H[i][i] = (k+G[i])*(k+G[i])/2 ;         /*  kinetic energy of free particle of momentum k+G      */
            if( i+1 < ng) H[i][i+1] = V1/2;        /*  matrix element <k+G_i| V |k+G_{i+1}>                 */
            if( i+2 < ng) H[i][i+2] = V2/2;        /*  matrix element <k+G_i| V |k+G_{i+2}>                 */
        }
    
    /*     find ng eigenvalues and eigenvectors of H_k   */
    
         diagonalize_(&ngmax, &ng, H, e, psiR, psiI, ES, TAU, GRA, GIA);      /* call FORTRAN subroutine "diagonalize" */
        
    /*     write out the lowest n_wf band energies for each k    */
        
        fprintf(bands, " %lf \t", k*a/(2*pi) );  
        for(i=0; i < n_wf; i++)
        {
            fprintf(bands, " %lf \t", e[i]);         /*  ith band energy  */
        }
        fprintf(bands, " \n");
    }
    plot(V1,V2);                      /*  plot out the band structure  */
}
/******************************************************************/
int plot(double V1, double V2) 
{
    /*  this function calls gnuplot to plot the results */
    
	FILE *pipe = popen("gnuplot -persist","w");
    
    /* x (momentum) plot range */    
    fprintf(pipe, "set xrange [-0.5:0.5]\n");
    
    /* graph title and labels for axes */    
    fprintf(pipe, "set title 'Nearly Free Particle Band Structure, V1=%f, V2=%f'\n", V1, V2);
	fprintf(pipe, "set xlabel 'Momentum (pi/a)'\n");
	fprintf(pipe, "set ylabel 'Energy'\n");
    
    /* output to postcript file */    
    fprintf(pipe, "set terminal postscript \n");
    fprintf(pipe, "set output 'NFE_band_structure.ps'\n"); 
	fprintf(pipe, "plot 'NFE_bands' using 1:2 title \"band 1\" with lines, 'NFE_bands' using 1:3 title \"band 1\" with lines, 'NFE_bands' using 1:4 title \"band 1\" with lines, 'NFE_bands' using 1:5 title \"band 1\" with lines \n");
    
    /* output to screen */    
    fprintf(pipe, "set terminal x11 \n"); 
    fprintf(pipe, "set output \n"); 
    fprintf(pipe, "replot\n");
    
    fflush(pipe); 
    
    close(pipe);
    return 0;
}



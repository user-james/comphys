#include <stdio.h>
#include <math.h>

/*   
 
 This program investigates the quantum ground state of the Lennard-Jones liquid/gas, using the Metropolis algorithm 
 to sample the ground state probability distribution. 
 The variational ground state wave function is a product of pair-wise positive terms, with two parameters, 
 which may be varied to minimize the energy, as discussed by McMillan (see reference, below). 
 [Refs: W. L. McMillan, Phys. Rev. 138, A442 (1965) and 
        N. Metropolis, A. W. Rosenbluth, M. N. Rosenbluth, A. H. Teller, E. Teller, J. Chem. Phys.21, 1087 (1953)].
 
 The system of N particles is confined in a cube of side L, with periodic boundary conditions. 
 The potential energy of the system is a sum of pair-wise Lennard-Jones interaction terms,
                              V_LJ(r)  =  alpha/r^12  - beta/r^6   =  4 epsilon [ (r_0/r)^12 - (r_0/r)^6 ].

 The Metropolis method uses a random walk in the configuration space {(x_i,y_i,z_i), i=1,N} of the particles, which samples 
 the variational ground-state wavefunction distribution,
                 P({x_i,y_i,z_i})  =  |psi({r_i})|^2  =   P(R)  =  exp[-2V(R)] / Z
 where 
                 V(R)  =  sum_{i=1,N} sum_{j<i} (a_1/r_ij)^a_2
 and 
                  Z  =  int exp[-V(R)] d^{3N}R   is the normalization factor for the wave function.
                
 In this code: 
      length is measured in units of the Lennard-Jones radius, r_0 (where V_LJ(r_0) = 0)
      energy is measured in units of the quantum confinement energy, hbar^2/(mr_0^2), where m is the mass of a Helium atom
 
      in these units,
                      V_LJ(r)  =   4A [1/r^12 - 1/r^6]
      and the kinetic energy operator of the system is 
                            T  =   sum_i { - grad_i^2 / 2 }

 
 Here we vary the parameter a_1, keeping a_2 fixed at 5, and calculate the average energy (expectation value of the Hamiltonian)
 as a function of a_1.
 The best estimate of the ground state wave function and properties is for the value of a_1 that minimizes the average energy.
 
 */


#define NMAX 100              /* maximum number of atoms (dimension of atomic arrays) */


FILE *output;


/***********************************************************************************************************/
int main( )
 {
     
   /*  parameters of the Lennard-Jones interaction potential  */
     
     double A = 5.55;                                /*   quantum dimensionless constant  */
                                                     /*   A equals 5.556, 115.48 and 1205.43 for He, Ne and Ar, respectively. */  
     
     double alpha = 4*A, beta = 4*A;                 /*   parameters for Lennard-Jones potential */
     
     
   /*  particle number and density parameters   */
     
     double L_eq = 6.0069, L;                                 /* side length of box */
     int N = 80;                                     /* number of particles in box */
     /* For He, the experimental equiibrium liquid density is 2.20 x 10^{22} per cm^3, corresponding to L = 6.0069 for N = 80 */
     
     
   /*  wave function parameters - values that McMillan found to minimize the energy  */
    
     double a_1_min = 2.6/2.56,   a_1;
     double a_2_min = 5.0,   a_2;
     
     
   /*  Metropolis random walk parameters  */
     
     long int iseed = 56398765;                      /* random number seed  */        
     
     int nstep_settle = 10000;                        /* number of random walk steps to equilibrate distribution after each change of parameters */
     int nstep_average = 20000;                      /* number of time-steps to average Hamiltonian after each change of parameters */
     double dx = 0.3;                                /* max coordinate-step for random walk */
     
     
     double  x[NMAX],  y[NMAX],  z[NMAX], rij[NMAX][NMAX][4];   /*  particle coords, interatomic distances  */

     double energy_av;
     int i;

     
    /* functions used   - see descriptions below  */   
     
     int  ranposn(double *, double *, double *, double, int);
     double random_walk_quantum( double *, double *, double *, double, int, double, double, double, double, double (*)[NMAX][4], double, int);
     double  dran1_(long int *);
     int bin_rij_PBC_quantum(double (*)[NMAX][4], int, int );
     int rho1_PBC(double *, double *, double *, int, double, double, double, int);
     
 
     output=fopen("./data/energy_vs_a1","w");             /* name and open output file */
     
     
/*   -----------------------------------------------------------------------------  */ 
 
     dran1_(&iseed); iseed = 0;        /* initialize random number generator and set iseed to zero for subsequent calls  */  
     
     printf("Lennard-Jones interaction: alpha = %lf, beta = %lf \n", alpha, beta);
                            
     L = L_eq;
     printf("simulation of %d particles at density = %lf \n", N, (double)N/(L*L*L) );
     a_2 = a_2_min;
     /* a_1 = a_1_min; */  
     
     for(a_1 = a_1_min*0.95; a_1 < a_1_min*1.05; a_1 += a_1_min*0.01)         /* vary the parameter a_1 in 1% steps about McMillan's energy-minimum value  */
     { 
         ranposn(x,y,z, L, N);             /* start with particles distributed uniformly at random in cube  */

         printf("\n equilibration stage: %d steps  \n", nstep_settle);
         energy_av = random_walk_quantum( x,y,z,L,N, alpha,beta, a_1,a_2, rij, dx, nstep_settle );    /* run random walk for nstep_settle steps to equilibrate distribution */
         
         printf("averaging stage: %d steps  \n", nstep_average);
         energy_av = random_walk_quantum( x,y,z,L,N, alpha,beta, a_1,a_2, rij, dx, nstep_average );   /* run random walk for nstep_average steps to estimate average energy */
         
         fprintf(output, " %lf \t %lf \n", a_1, energy_av);                                              /* output energy versus a_1  */
         bin_rij_PBC_quantum(rij, N, 2);                                                                 /* output pair-distribution histogram */
         rho1_PBC(x,y,z, N, a_1,a_2, L, 2);                                                              /* output single-particle density matrix  */
     }
     
     
     plot( a_1_min*0.95, a_1_min*1.05 );     /*  plot out the results of energy vs temperature  */
 }


 /***********************************************************************************************************/
int plot(double a_min, double a_max) 
{
    /*  this function calls gnuplot to plot the results 
     
     INPUT:
       
       max_temp = the maximum temperature (on the x-axis of the plot) 
     
     */
    
	FILE *pipe = popen("gnuplot -persist","w");
    
    /* x and y plot ranges */    
    /* fprintf(pipe, "set xrange [ %lf : %lf ]\n", Vint0, Energy); */
	fprintf(pipe, "set xrange [%lf : %lf ]\n", a_min, a_max);
    
    /* graph title and labels for axes */    
    fprintf(pipe, "set title 'Lennard-Jones Quantum Ground State: Energy vs variational parameter, a_1'\n");
	fprintf(pipe, "set ylabel 'Energy'\n");
	fprintf(pipe, "set xlabel 'a_1'\n");
    
    /* output to postcript file */    
    fprintf(pipe, "set terminal postscript \n");
    fprintf(pipe, "set output './plots/energy_vs_a1.ps'\n"); 
	fprintf(pipe, "plot './data/energy_vs_a1' using 1:2 notitle  with lines \n");
    
    /* output to screen */    
    fprintf(pipe, "set terminal x11 \n"); 
    fprintf(pipe, "set output \n"); 
    fprintf(pipe, "replot\n");
    
    fflush(pipe); 
    
    close(pipe);
    return 0;
}

/************************************************************************************************************
 
NOTES:
 
 
Function files used here...

random_walk_quantum.c:  (contains function random_walk)
Integrates the equations of motion over a specified number of time steps using the Metropolis algorithm.
Calls functions pot_energy_PBC, kin_energy_PBC and wfn_log_ratio_PBC.

energy_PBC_quantum.c:  (contains functions pot_energy_PBC, kin_energy_PBC and wfn_log_ratio_PBC, called from random_walk)
   pot_energy_PBC calculates Lennard-Jones interatomic energy for system 
   kin_energy_PBC calculates the kinetic energy for the pair-wise product wave function and 
   wfn_log_ratio_PBC calculates the ratio of old to new wave function after moving one of the particles
all are calculated with periodic boundary conditions.

random_positions.c: (contains function ranposn)
Places N particles uniformly at random in a cube of side L.
Calls function dran1.
 

dran1.f:  Fortran function to generate uniform pseudorandom numbers.
 
distance_histrogram_PBC_quantum.c: (contains function bin_rij_PBC_quantum)
 Atomic pair-distribution histogram accumulation, initialization or output to file
 
 
 */



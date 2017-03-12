#include <stdio.h>
#include <math.h>

/*   
 
 This program investigates the heat-capacity of a noble gas cluster, confined in a large box,
 using molecular dynamics (MD) with the Lennard-Jones interatomic potential.
 
 A fragment of several shells of atoms from a fcc crystal is set up near equilibrium and then annealed to equilibrium
 using damped leap-frog dynamics.
 The system is then gradually heated by intermittently adding small random increments to the particle velocities.
 Between each such random increment of the velocities, 
    (a) the system is allowed to re-equilibrate for a number of time-steps at the new energy level,
    (b) the system kinetic energy is then averaged over a larger number of time-steps to estimate the temperature.
    (c) a histogram of the pair-distribution function is also calculated during stage (b)
 
 */


#define NMAX 100              /* maximum number of atoms (dimension of atomic arrays) */
#define NSHMAX 5              /* largest number of shells allowed in current version of code is 5 */


FILE *output;


/***********************************************************************************************************/
int main( )
 {
     long int iseed = 56398765;                      /* random number seed  */        
     
     double alpha = 1.0, beta = 1.0;                 /* parameters for Lennard-Jones potential */
     double L = 30.0, V0 = 0.1;                      /* parameters for confining box potential */
     
     double dt = 0.01;                               /* time-step for leap-frog dynamics */
      
     double dE = 0.01;                                /* increment of the energy per particle for heating system */
     double max_temp = 0.05;                          /* maximum temperature / k_B to be reached  */
     
     int nstep_anneal = 10000;                       /* number of time-steps to anneal system to equilibrium structure */
     int nstep_settle = 1000;                        /* number of time-steps to allow settling after each velocity increment */
     int nstep_average = 100000;                     /* number of time-steps to average temperature after each velocity increment */

     int Nshell = 5;                                 /* number of shells of atoms in cluster (fcc fragment) */
     
     
     int N, i;                                                  /*  number of particles and particle index */
     double  x[NMAX],  y[NMAX],  z[NMAX], rij[NMAX][NMAX];      /*  particle coords, interatomic distances  */
     double vx[NMAX], vy[NMAX], vz[NMAX];                       /*  particle velocities  */
     double Fxi[NMAX], Fyi[NMAX], Fzi[NMAX];                    /*  internal particle forces      */
     double Fxe[NMAX], Fye[NMAX], Fze[NMAX];                    /*  external particle forces      */

     double gamma;                                              /* velocity damping constant */
     double Energy, Vint, Vint0, Vext, KE_average, temp;
     
/* functions used   - see descriptions below  */   
     
     int  fccposn(double *, double *, double *, double, double, int);
     double ranveladd(double *, double *, double *, int, double, long int *);
     int  intdyn(double *, double *, double *, double *, double *, double *, 
                 int, int, double, double, double, double, double *, double *, double *, 
                 double, double, double *, double *, double *, double *, double *,
                 double rij[][NMAX], double * );
     double  dran1_(long int *);
     int bin_ij(double rij[][NMAX], int, int );

 
     output=fopen("output.txt","w");             /* name and open output file */

/*   -----------------------------------------------------------------------------
      set up a cluster of Lennard-Jones atoms, centred at (0,0,0),
      with the atomic positions arranged in up to 5 shells of an fcc lattice, 
      with the inter-atomic distance close to equilibrium,  */
     
     
     Nshell = 5;                        /* largest number of shells allowed is 5 */

     N = fccposn(x,y,z, alpha,beta, Nshell );
     printf("%d\n", N);
     for(i=0; i < N; i++) 
     {
         vx[i]=0.0;   vy[i]=0.0;   vz[i]=0.0; 
     }
 
     
/*   -----------------------------------------------------------------------------
      anneal the positions of the atoms to equilibrium using damped leap-frog MD */ 
 
  
     gamma = 1.0;                                          /* damped motion allows system to settle at equilibrium  */
     intdyn(x,y,z, vx,vy,vz, N, 
             nstep_anneal, dt, gamma,   alpha,beta, 
             Fxi,Fyi,Fzi, L,V0, Fxe,Fye,Fze, 
             &Vint,&Vext, rij, &KE_average );
 
/*  run one time-step to find annealed potential and kinetic energies and pair-distribution  */
     bin_rij(rij, N, 0);                                /* initialize pair-distribution histogram to zero */
     intdyn(x,y,z, vx,vy,vz, N, 
            1, dt, gamma,   alpha,beta, 
            Fxi,Fyi,Fzi, L,V0, Fxe,Fye,Fze, 
            &Vint,&Vext, rij, &KE_average );
     printf("Initial Potential Energy: \t %lf \t %lf \t Kinetic Energy: \t %lf \n", Vint, Vext, KE_average);
     Vint0 = Vint;
     bin_rij(rij, N, 2);  /* output pair-distribution histogram */
/*
         printf(" \n POSITIONS \n " );
     for(i=0; i < N; i++) 
     {
         printf("%lf \t %lf \t %lf \t %lf \n", x[i],y[i],z[i], rij[0][i]); 
     }

 */
 
     
/*   -----------------------------------------------------------------------------
       simulate heating process...         */
 
  
     gamma = 0.0;                     /* remove damping to allow conservative dynamics  */
     Energy = Vint0;                  /* take equilibrium structure potential energy   */
        
     
     dran1_(&iseed); iseed = 0;        /* initialize random number generator and set iseed to zero for subsequent calls  */   
     do {  
         Energy += ranveladd( vx,vy,vz, N, dE, &iseed );    /*  add small random Gaussian velocity increments  */
         intdyn(x,y,z, vx,vy,vz, N, 
                 nstep_settle, dt, gamma,   alpha,beta, 
                 Fxi,Fyi,Fzi, L,V0, Fxe,Fye,Fze, 
                 &Vint,&Vext, rij, &KE_average );           /*  allow system to settle after each velocity increment */
         
         bin_rij(rij, N, 0);                                /* initialize pair-distribution histogram to zero */

         intdyn(x,y,z, vx,vy,vz, N, 
                 nstep_average, dt, gamma,  alpha,beta, 
                 Fxi,Fyi,Fzi, L,V0, Fxe,Fye,Fze, 
                 &Vint,&Vext, rij, &KE_average );           /*  integrate eqns of motion, taking average of kinetic energy  */
         temp = 2 * KE_average / (3*N);                     /*  average kinetic energy per particle = 3/2 k_B T  */
         
         printf("Energy \t %lf  \t Temperature \t %lf \n", Energy, temp );
         
         bin_rij(rij, N, 2);  /* output pair-distribution histogram */


       
         fprintf(output,"%lf\t%lf \n", Energy, temp  );       /* plot temperature vs. energy added  */
     }     
     while (temp < max_temp);
     
/*
The heat capacity is the inverse of the derivative of the "temperature" with respect to system energy. 
To avoid large statistical noise, the temperature vs. energy should be fitted with a smooth function before 
taking the derivative.  */
     
     plot(Vint0, Energy, max_temp );     /*  plot out the results  */
 }


 /***********************************************************************************************************/
int plot(double Vint0, double Energy, double max_temp) 
{
    /*  this function calls gnuplot to plot the results 
     
     INPUT:
       Vint0 = the minumum energy (on the x-axis of the plot) 
       Energy = the maximum (on the x-axis of the plot)
       max_temp = the maximum temperature (on the y-axis of the plot) 
     
     */
    
	FILE *pipe = popen("gnuplot -persist","w");
    
    /* x and y plot ranges */    
    fprintf(pipe, "set xrange [ %lf : %lf ]\n", Vint0, Energy);
	fprintf(pipe, "set yrange [0.00: %lf ]\n", max_temp);
    
    /* graph title and labels for axes */    
    fprintf(pipe, "set title 'Lennard-Jones Cluster: Temperature vs Energy'\n");
	fprintf(pipe, "set xlabel 'System Energy'\n");
	fprintf(pipe, "set ylabel 'Temperature * k_B'\n");
    
    /* output to postcript file */    
    fprintf(pipe, "set terminal postscript \n");
    fprintf(pipe, "set output 'temp_vs_energy.ps'\n"); 
	fprintf(pipe, "plot 'output.txt' using 1:2 notitle  with lines \n");
    
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
 
In the low-temperature limit (harmonic oscillation regime) the heat capacity should be 3 N, 
where N is the number of particles  (Dulong-Petit Law).
 
In the high-temperature limit (ideal gas) the heat capacity should be 3N/2

Note: we are using units, where the Boltzmann constant is unity
 
 
Function files used here...

int_dyn.c:  (contains function intdyn)
Integrates the equations of motion over a specified number of 
time steps using the leap-frog methods (allows for damping).
Calls function forces.

forces.c:  (called from intdyn)
Calculates Lennard-Jones interatomic forces and the external force from a confining "box" potential.
forces_mirror.c: calculates internal forces only and reverses the particle perpendicular velocity component when it hits a wall of the box.

fcc_posn.c: (contains function fccposn)
Constructs a fragment of an fcc crystal with the bond length
equal to the equilibrium separation of two atoms in the 
Lennard-Jones potential.

ran_vel_add.c: (contains function ranveladd)
Adds random increment to particle velocities from a Maxwell 
distribution at specified energy increment.  This subroutine 
calls the fortran subroutine dran1 (seen in C program as dran1_) 
which takes a non-zero integer seed on its first call, 0 on 
subsequent calls.

dran1.f:  Fortran function to generate uniform pseudorandom numbers.
 
distance_histrogram.c: (contains function bin_rij)
 Atomic pair-distribution histogram accumulation, initialization or output to file
 
 
 */



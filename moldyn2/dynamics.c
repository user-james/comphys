#include <stdio.h>
#include <math.h>

/*   
 
 This program investigates the velocity-velocity correlation function of a noble gas cluster, confined in a large box,
 using molecular dynamics (MD) with the Lennard-Jones interatomic potential.
 
 A fragment of several shells of atoms from a fcc crystal is set up near equilibrium and then annealed to equilibrium
 using damped leap-frog dynamics.
 The system is then gradually heated by intermittently adding random increments to the particle velocities.
 Between each such random increment of the velocities, 
    (a) the system is allowed to re-equilibrate for a number of time-steps at the new energy level,
    (b) the velocity-velocity correlation function C_v then averaged over a larger number of time-steps and
    (c) we find the Fourier Transform of C_v to get the vibration spectrum.
 
 */


#define NMAX 100              /* maximum number of atoms (dimension of atomic arrays) */
#define NSHMAX 5              /* largest number of shells allowed in current version of code is 5 */

#define MAXTAU 1050           /* maximum number of time-steps for velocity-velocity correlation function */


FILE *spectrum;


/***********************************************************************************************************/
int main( )
 {
     long int iseed = 56398765;                      /* random number seed  */        
     
     double alpha = 1.0, beta = 1.0;                 /* parameters for Lennard-Jones potential */
     double L = 30.0,    V0 = 0.1;                   /* parameters for confining box potential */
     
     int Nshell = 5;                                 /* number of shells of atoms in cluster (fcc fragment) */
     
     double dt = 0.01;                               /* time-step for leap-frog dynamics */
      
     double dE = 0.01;                                /* increment of the energy per particle for heating system */
     double max_temp = 0.05;                          /* maximum temperature / k_B to be reached  */
     
     int nstep_anneal = 10000;                       /* number of time-steps to anneal system to equilibrium structure */
     int nstep_settle = 10000;                       /* number of time-steps to allow settling after each velocity increment */
     int nstep = 20;                                 /* nstep = # MD time-steps between averaging times for C_v  */
     int Nav = 10000;                                /* Nav is the number of values of t used in averaging C_v at each temperature */
     int maxtau = 1024;                              /* T = maxtau*nstep*dt  */
     double max_freq = 2.5;                          /* max_freq = maximum frequency (cycles per unit time) output for velocity spectrum */

     
     int N, i;                                                  /*  number of particles and particle index */
     double  x[NMAX],  y[NMAX],  z[NMAX], rij[NMAX][NMAX];      /*  particle coords, interatomic distances  */
     double vx[NMAX], vy[NMAX], vz[NMAX];                       /*  particle velocities  */
     double Fxi[NMAX], Fyi[NMAX], Fzi[NMAX];                    /*  internal particle forces      */
     double Fxe[NMAX], Fye[NMAX], Fze[NMAX];                    /*  external particle forces      */

     double gamma;                                              /* velocity damping constant */
     double Energy, Vint, Vint0, Vext, KE_average, temp;
     
     double C_v[MAXTAU];                                                    /* accumulated velocity-velocity correlation function at time-difference tau = itau * nstep * dt   */
     double vxold[NMAX][MAXTAU], vyold[NMAX][MAXTAU], vzold[NMAX][MAXTAU];  /* saved velocities for recent MD steps, used in calculating C_v   */
     int itau, iav, i_file;
     
     double A[4*MAXTAU+2], work[4*MAXTAU+2], trig[4*MAXTAU+2];              /* arrays and variables for use in Fourier Transform  (see FORTRAN function fft1d) */
     long int n_ft, ft_dimension;
     long int isign;
     double freq, signal;
     char  filename[30];
     
     
/* functions used   - see descriptions below  */   
     
     int  fccposn(double *, double *, double *, double, double, int);
     double ranveladd(double *, double *, double *, int, double, long int *);
     int  intdyn(double *, double *, double *, double *, double *, double *, 
                 int, int, double, double, double, double, double *, double *, double *, 
                 double, double, double *, double *, double *, double *, double *,
                 double rij[][NMAX], double * );
     double  dran1_(long int *);
     int bin_ij(double rij[][NMAX], int, int );
     void fft1d_(double *, double *, long int *, long int *, long int *, double *);

 
     /* output=fopen("output.txt","w");             /* name and open output file */

/*   -----------------------------------------------------------------------------
      set up a cluster of Lennard-Jones atoms, centred at (0,0,0),
      with the atomic positions arranged in up to 5 shells of an fcc lattice, 
      with the inter-atomic distance close to equilibrium,  */
     
     
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
 
/*  run one time-step to find annealed potential and kinetic energies */
     
     intdyn(x,y,z, vx,vy,vz, N, 
            1, dt, gamma,   alpha,beta, 
            Fxi,Fyi,Fzi, L,V0, Fxe,Fye,Fze, 
            &Vint,&Vext, rij, &KE_average );
     printf("Initial Potential Energy: \t %lf \t %lf \t Kinetic Energy: \t %lf \n", Vint, Vext, KE_average);
     Vint0 = Vint;

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
     i_file = 0;                      /* index of current output spectrum file */
        
     
     dran1_(&iseed); iseed = 0;        /* initialize random number generator and set iseed to zero for subsequent calls  */   
     do {  
         Energy += ranveladd( vx,vy,vz, N, dE, &iseed );    /*  add small random Gaussian velocity increments  */
         intdyn(x,y,z, vx,vy,vz, N, 
                nstep_settle, dt, gamma,   alpha,beta, 
                Fxi,Fyi,Fzi, L,V0, Fxe,Fye,Fze, 
                &Vint,&Vext, rij, &KE_average );           /*  allow system to settle after each velocity increment */

         bin_rij(rij, N, 0);                                /* initialize pair-distribution histogram to zero */
         
/*
    ---------------------------------------------------------------------------------------------------------------------
        calculate the velocity-velocity correlation function C_v(tau) for tau = itau * nstep * dt, 0 < itau < maxtau
          
          initial set-up of old velocities, etc.
          
          NOTE:  during averaging vxold[i][itau] is the x velocity at 
          a time itau*nstep*dt earlier than the present time.
          */
         for(itau=maxtau-1; itau >= 0; itau--)
         {
             intdyn(x,y,z, vx,vy,vz, N, 
                    nstep, dt, gamma,  alpha,beta, 
                    Fxi,Fyi,Fzi, L,V0, Fxe,Fye,Fze, 
                    &Vint,&Vext, rij, &KE_average ); 
             
             C_v[itau] =  0.0;
             for(i=0; i < N; i++)
             {
                 vxold[i][itau] = vx[i]; 
                 vyold[i][itau] = vy[i]; 
                 vzold[i][itau] = vz[i];
             }
         }
         
/*
    -----------------------------------------------------------------------
    do averaging of C_v at Nav times.
*/
         for(iav=0; iav < Nav ; iav++)
         {
             intdyn(x,y,z, vx,vy,vz, N, 
                    nstep, dt, gamma,  alpha,beta, 
                    Fxi,Fyi,Fzi, L,V0, Fxe,Fye,Fze, 
                    &Vint,&Vext, rij, &KE_average );         
             
             for(itau=maxtau-1; itau>0 ; itau--)
             {
                 for(i=0; i < N; i++)
                 {
                     vxold[i][itau] = vxold[i][itau-1];     /*  shift old velocities */
                     vyold[i][itau] = vyold[i][itau-1];     /*  on by one time step  */
                     vzold[i][itau] = vzold[i][itau-1];
                     C_v[itau] += vx[i]*vxold[i][itau] + vy[i]*vyold[i][itau] + vz[i]*vzold[i][itau] ;  /* accumulate average */
                 }
             }
             for(i=0; i < N; i++)      /* equal-time average */
             {
                 vxold[i][0] = vx[i];
                 vyold[i][0] = vy[i];
                 vzold[i][0] = vz[i];
                 C_v[0] += vx[i]*vxold[i][0] + vy[i]*vyold[i][0] + vz[i]*vzold[i][0] ;
             }
         }
         
/*
    ------------------------------------
    normalize C_v
*/
         
         for(itau=0; itau < maxtau ; itau++)
         {
             C_v[itau] = C_v[itau]/(N*Nav);
         }
         temp = C_v[0] / 3.0 ;                   /* average particle speed squared = 3 * temperature    */
         
/*
    ------------------------------------
    Fourier Transform C_v to get spectrum
*/
         
         for(i=0 ; i < maxtau; i++)
         {
             A[4*maxtau-2*i]   = (A[2*i] = C_v[i]);    /* mirror real part of C_v and double the repeat time to avoid discontinuities */
             A[4*maxtau-2*i+1] = (A[2*i+1] = 0.0);     /* imaginary part of C_v = 0 */
         }
         
         isign = 1;   n_ft = 2*maxtau;   ft_dimension = 2*MAXTAU;
         
         fft1d_(A, work, &n_ft, &isign, &ft_dimension, trig);
         
         
         i_file += 1;                                  /* open a new spectrum file for each temperature */
         sprintf(filename, "spectrum_%d", i_file);
         spectrum = fopen(filename, "w");
         
         freq = 0.0;
         i = 0;
         while ( freq < max_freq )                     /* output the spectrum for frequencies up to maximum frequency */ 
         {
             signal = n_ft * sqrt( A[2*i]*A[2*i]+A[2*i+1]*A[2*i+1] );
             fprintf(spectrum, " %lf \t %lf \n", freq, signal );
             i += 1;
             freq = i/(2*maxtau*nstep*dt);
         }
         fclose(spectrum);
         printf(" system energy %lf, \t system temperature %lf \n", Energy, temp );
     }     
     while (temp < max_temp);
          
    /* plot(Vint0, Energy, max_temp );     /*  plot out the results  */
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
	fprintf(pipe, "set ylabel 'Temperature / k_B'\n");
    
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
Calculates Lennard-Jones interatomic forces and the external
force from a confining "box" potential.

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



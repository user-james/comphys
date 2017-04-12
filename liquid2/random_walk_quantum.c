#include <stdio.h>
#include <math.h>

#define NMAX 100



/************************************************************************************************/
double random_walk_quantum( double x[], double y[], double z[], double L, int N, 
                            double alpha, double beta, double a_1, double a_2, 
                            double rij[][NMAX][4], double dx, int nstep)
{
/*    
   Metropolis random walk, based on ground-state Boson wave function = product of positive pair-wise terms:
   [Ref: W. L. McMillan, Phys. Rev. 138, A442 (1965)]
 
       psi({r_i})  =  exp[ - sum_{i=1,N} sum_{j < i} (a_1/|r_i-r_j|)^a_2 ]     [see Eqs. (3) & (5) of McMillan]
 
   This function takes nstep (trial) steps of each of the N particles in a box of side L (periodic boundary conditions)
   and returns the average potential energy (using Lennard-Jones pair-wise interaction). 
   It also accumulates a histogram of the pair-distribution function.
      

  
Input variables:
  (x,y,z) =  positions of particles (cartesian coords)
  L       =  size of "box" for periodic boundary conditions
  N       =  number of particles
  alpha   =  coefficient of 1/r^12  in LJ potential
  beta    =  coefficient of -1/r^6  in LJ potential
  a_1     =  first coefficient in pair-wise product wave function 
  a_2     =  second coefficient in pair-wise product wave function
  rij     =  interparticle distances
  dx      =  maximum step size
  nstep   =  number of (trial) random steps for each particle


Output Variables:
 
  rij   =  interatomic distances (including periodic images, indexed by 3rd index)

function returns V_average  =  average potential energy of system over the nstep steps

*/
    double x_new, y_new, z_new;
    double H_average, dV, V, T;
    int istep, i, n_accept;
    static long int iseed = 0;

 
    /* functions used  */
    
    int bin_rij_PBC_quantum(double (*)[NMAX][4], int, int );
    int rho1_PBC(double *, double *, double *, int, double, double, double, int);
    double dran1_(long int *);
    double  wfn_log_ratio_PBC(double *, double *, double *, int, 
                              double, double, double,
                              int, double, double, double);
    double  pot_energy_PBC(double *, double *, double *, int,  
                           double, double, double,
                           double (*)[NMAX][4] );
    double  kin_energy_PBC(double *, double *, double *, int,  
                           double, double, double,
                           double (*)[NMAX][4] );
    

  if( N > NMAX )
  {
      printf("%d is too many particles for random_walk. NMAX = %d\n", N, NMAX);
      return 1;
  }

         /* accumulate averages over nstep configurations of system  */
    
    bin_rij_PBC_quantum(rij, N, 0);        /*  initialize pair-distribution histogram to zero */
    rho1_PBC(x,y,z, N, a_1,a_2, L, 0);     /* initialize single-particle density matrix to zero */
    H_average = 0.0;                       /*  average energy  */
    n_accept = 0;                          /*  trial move acceptance rate for random walk  */
    
    for(istep=0; istep<nstep; istep++)    
    {
        for(i=0; i< N; i++)               /*  loop through all particles, attempting move of each one  */
        {
            /* define new trial position */
            
            x_new = x[i] + dx * (2*dran1_(&iseed)-1.0);  
            if(x_new > L/2) x_new += -L; if(x_new < -L/2) x_new += L;
            
            y_new = y[i] + dx * (2*dran1_(&iseed)-1.0);  
            if(y_new > L/2) y_new += -L; if(y_new < -L/2) y_new += L;
            
            z_new = z[i] + dx * (2*dran1_(&iseed)-1.0);  
            if(z_new > L/2) z_new += -L; if(z_new < -L/2) z_new += L;
      
            /* find change in V = -log of wave function, and accept or reject move at random 
                based on Metropolis method for sampling quantum distribution, exp[-2V]  */
            
            dV  =  wfn_log_ratio_PBC(x,y,z, N, a_1,a_2, L, i, x_new,y_new,z_new);
            
            if( dV < 0.0 ) dV = -1.0;  /* handle large number overflow  in exponential */
            if( dV > 20.0) dV = 20;
            
            if( exp(-2*dV) > dran1_(&iseed) )  /* accept or reject trial move  */
            {
                x[i] = x_new;   y[i] = y_new;   z[i] = z_new;
                n_accept += 1;
            }
        }
        
        V = pot_energy_PBC(x,y,z, N, alpha,beta, L, rij );    /*  V is the potential energy of the new configuration  */
        T = kin_energy_PBC(x,y,z, N, a_1,a_2, L, rij );       /*  T is the kinetic energy of the new configuration  */
        H_average += V + T;
        
        bin_rij_PBC_quantum(rij, N, 1);   /* accumulate statistics for pair-distribution function  */
        rho1_PBC(x,y,z, N, a_1,a_2, L, 1);     /* accumulate statistics for single-particle density matrix  */
    }
    
    H_average = H_average/nstep;
    
    printf(" step size = %lf,  acceptance rate = %lf \n", 
           dx, (double)n_accept/(double)nstep/(double)N );

    
    return H_average;
}


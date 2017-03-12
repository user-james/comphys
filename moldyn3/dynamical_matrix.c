#include <stdio.h>
#include <math.h>

/*   
 
 This program calculates the binding energy and dynamical matrix of a noble gas cluster with the Lennard-Jones interatomic potential 
 and finds its normal mode frequencies for a series of scaled atomic configurations about the equilibrium geometry.
 
 A fragment of several shells of atoms from a fcc crystal is set up near equilibrium and then annealed to equilibrium
 using damped leap-frog dynamics. 
 The binding energy and normal mode frequencies are then found for a series of cluster configurations, 
 in which the equilibrium atomic coordinates are scaled by a factor s, where 0.98 < s < 1.02. 
 
 For each (scaled) atomic geometry, with atoms at positions (xs[i],ys[i],zs[i]), the dynamical matrix is found as follows:
 Each atom is then displaced slightly from equilibrium in each of the three (x,y,z) directions and the resulting forces on all atoms are calculated.
 The dynamical matrix is the derivative of the forces with respect to displacement. (The atomic mass is assumed to be m=1.)
 The eigenvectors of the dynamical matrix are normal modes and the eigevalues are the normal mode frequencies squared.
 
 */


#define NMAX 100              /* maximum number of atoms (dimension of atomic arrays) */
#define NSHMAX 5              /* largest number of shells allowed in current version of code is 5 */
#define MAXSCALE 20           /* maximum number of scaled atomic distances on either side of equilibrium */



FILE *spectrum;
FILE *matrix_out;
FILE *evals;
FILE *grun_params;

/***********************************************************************************************************/
int main( )
 {

     double alpha = 1.0, beta = 1.0;                 /* parameters for Lennard-Jones potential */
     double L = 3000.0,    V0 = 0.1;                 /* parameters for confining box potential */
     
     int Nshell = 5;                                 /* number of shells of atoms in cluster (fcc fragment) */
     
     double dt = 0.01;                               /* time-step for leap-frog dynamics */
      
     double dE = 0.01;                               /* increment of the energy per particle for heating system */
     double max_temp = 0.05;                         /* maximum temperature / k_B to be reached  */
     
     int nstep_anneal = 10000;                       /* number of time-steps to anneal system to equilibrium structure */
     
     double dx = 0.001;                              /* displacement of atoms from equilibrium in finding derivative of forces */
     int nscale = 20;                                 /* number of scaled geometries calculated on either side of equilibrium  */
     double fac = 0.02;                              /* scale atomic coordinates from (1-fac) to (1+fac)  */


     
     int N, i,j;                                                /*  number of particles and particle indices */
     double  x[NMAX],  y[NMAX],  z[NMAX], rij[NMAX][NMAX];      /*  particle coords, interatomic distances  */
     double vx[NMAX], vy[NMAX], vz[NMAX];                       /*  particle velocities  */
     double Fxi[NMAX], Fyi[NMAX], Fzi[NMAX];                    /*  internal particle forces      */
     double Fxe[NMAX], Fye[NMAX], Fze[NMAX];                    /*  external particle forces      */

     double gamma;                                              /* velocity damping constant */
     double Vint, Vint0, Vext, KE_average, temp;
     double D[3*NMAX];                                          /* eigenvalues of the dynamical matrix */
     
     double xs[NMAX], ys[NMAX], zs[NMAX];
     double Energy[2*MAXSCALE+1];                               /* energy of scaled system */
     double omega[2*MAXSCALE+1][3*NMAX];                        /* normal mode frequencies for scaled system  */
     double pressure, gruneisen, d_omega;
     int iscale;
     
/* functions used    */   
     
     int  fccposn(double *, double *, double *, double, double, int);
     double ranveladd(double *, double *, double *, int, double, long int *);
     int  intdyn(double *, double *, double *, double *, double *, double *, 
                 int, int, double, double, double, double, double *, double *, double *, 
                 double, double, double *, double *, double *, double *, double *,
                 double rij[][NMAX], double * );
     int forces_m(double *, double *, double *, int, double *, double *, double *,
              double, double, double, double,
              double *, double *, double *, double *, double *, double *, double *, double *, double rij[][NMAX] );

     int bin_ij(double rij[][NMAX], int, int );
     int dynmat(double *, double *, double *, int, double *, double *, double *, double,
                double, double, double, double, double rij[][NMAX], double *, int);


 
     /* output=fopen("output.txt","w");             /* name and open output file */

/*   -----------------------------------------------------------------------------
      set up a cluster of Lennard-Jones atoms, centred at (0,0,0),
      with the atomic positions arranged in up to 5 shells of an fcc lattice, 
      with the inter-atomic distance close to equilibrium,  */
     
     
     N = fccposn(x,y,z, alpha,beta, Nshell );
     printf(" \n NUMBER OF ATOMS IN CLUSTER = %d\n", N);
     for(i=0; i < N; i++) 
     {
         vx[i]=0.0;   vy[i]=0.0;   vz[i]=0.0; 
     } 
 
     bin_rij(rij, N, 0);
/*   -----------------------------------------------------------------------------
      anneal the positions of the atoms to equilibrium using damped leap-frog MD 
*/     
  
     gamma = 1.0;                                          /* damped motion allows system to settle at equilibrium  */
     intdyn(x,y,z, vx,vy,vz, N, 
             nstep_anneal, dt, gamma,   alpha,beta, 
             Fxi,Fyi,Fzi, L,V0, Fxe,Fye,Fze, 
             &Vint,&Vext, rij, &KE_average );
 

     printf(" \n EQUILIBRIUM ATOMIC POSITIONS AND DISTANCES FROM CENTRAL ATOM \n " );
     for(i=0; i < N; i++) 
     {
         printf("%lf \t %lf \t %lf \t %lf \n", x[i],y[i],z[i], rij[0][i]); 
     }

/*     -----------------------------------------------------------------------------
        scale the atomic positions and find the potential energy and the mode frequencies  */
     
     
     for(iscale = -nscale; iscale < nscale+1; iscale++)
     {
         for(i=0; i<N; i++)
         {
             xs[i] = x[i] * (1.0 + iscale*fac/nscale);
             ys[i] = y[i] * (1.0 + iscale*fac/nscale);
             zs[i] = z[i] * (1.0 + iscale*fac/nscale);
         }
         
         forces_m(xs, ys, zs, N, vx, vy, vz,
                  alpha, beta, L, V0,
                  Fxi, Fyi, Fzi, Fxe, Fye, Fze, &Vint, &Vext, rij );
         Energy[iscale+nscale] = Vint;                                     /* potential energy for scaled cluster */
         
         dynmat( xs,ys,zs, N, vx,vy,vz, dx, alpha,beta,L,V0, rij, D, iscale);      /*  find the dynamical matrix, diagonalize it and return its eigenvalues  */
           
         /*  save the normal mode frequencies 
               -  we set to zero the lowest 6 frequencies, which correspond to translations and rotations of the cluster  
          */
         
         for(j=0; j<6; j++)               
         {
             omega[iscale+nscale][j] = 0.0;
         }

         for(j=6; j<3*N; j++)               
         {
             omega[iscale+nscale][j] = sqrt( D[j]);
         }
     }
     
/* ----------------------------------------------------------------------------
   write out the spectrum of frequencies for the equilibrium geometry and 
   write out the pressure and sum of mode Gruneisen parameters for scaling near equilibrium
 */
    
 
     grun_params = fopen("gruneisen.txt", "w");
     spectrum = fopen("equil_mode_frequencies", "w");
     printf(" nscale = %d \n", nscale);
     for(iscale=1; iscale < 2*nscale; iscale++) 
     {
         gruneisen = 0.0;
         for(j=6; j<3*N; j++)
         {
             d_omega =  ( omega[iscale-1][j] - omega[iscale+1][j] ) * nscale/fac;
             gruneisen += d_omega/omega[iscale][j];   
             if(iscale == nscale) 
                 fprintf(spectrum, " %d \t %lf \t %lf  \n", j, omega[iscale][j]/(2*3.1415926), d_omega/omega[iscale][j] ); 
         }
         pressure = (Energy[iscale-1] - Energy[iscale+1]) * nscale/fac;
         fprintf(grun_params, "%lf \t %lf \t %lf \n ", 1.0+(iscale-nscale)*fac/nscale, pressure, gruneisen );
     }
     fclose(grun_params);
     fclose(spectrum); 
     
 }
 

/***********************************************************************************************************/
int dynmat(double x[], double y[], double z[], int N, double vx[], double vy[], double vz[], double dx,
           double alpha, double beta, double L, double V0, double rij[][NMAX], double D[], int iscale)
{
    /*   -----------------------------------------------------------------------------
     calculate dynamical matrix by finite-difference calculation of derivatives of the forces w.r.t. atomic displacements   
     and find its eigenvalues and eigenvectors (normal modes)
     
     NOTE: maximum number of atoms allowed is NMAX  (defined parameter)
    
     INPUT:
         x,y,z         atomic coordinates
         N             number of atoms
         vx,vy,vz      atomic velocities (values not significant in this function)
         dx            small displacement of atoms used to find numerical derivatives of forces
         alpha, beta   parameters in Lennard-Jones forces
         L, V0         parameters defining the confining box (not used in this function)
         rij           interatomic distances array
     OUTPUT:
         D             eigenvalues of the dynamical matrix
     
     */
    
    
    double Vint, Vext;
    
    double DM[3*NMAX][3*NMAX];                                    /*  dynamical matrix */
    double Fxim[NMAX], Fyim[NMAX], Fzim[NMAX];                    /*  internal particle forces at x-dx or y-dx or z-dx   */
    double Fxip[NMAX], Fyip[NMAX], Fzip[NMAX];                    /*  internal particle forces at x+dx or y+dx or z+dx   */
    double Fxe[NMAX], Fye[NMAX], Fze[NMAX];                    /*  external particle forces   */
    
    long int n_dim, n_def;
    double A[3*NMAX][3*NMAX], ZR[3*NMAX][3*NMAX], ZI[3*NMAX][3*NMAX];    /* scratch space arrays for matrix diagonalization routines */
    double E[3*NMAX], TAU[6*NMAX], GRA[3*NMAX], GIA[3*NMAX];
    
    int i,j;
    
    void diagonalize_(long int *, long int *, 
                      double A[][3*NMAX], double *, double ZR[][3*NMAX], double ZI[][3*NMAX], 
                      double *, double *, double *, double *); 
    
    
    for(i=0; i < N; i++)
    {
        
        /*  displace x-coordinate of each atom by +dx and -dx and get central difference of forces on all atoms */
        
        x[i] = x[i] + dx;                         
        forces_m(x, y, z, N, vx, vy, vz,
                 alpha, beta, L, V0,
                 Fxip, Fyip, Fzip, Fxe, Fye, Fze, &Vint, &Vext, rij );
        x[i] = x[i] - 2*dx;
        forces_m(x, y, z, N, vx, vy, vz,
                 alpha, beta, L, V0,
                 Fxim, Fyim, Fzim, Fxe, Fye, Fze, &Vint, &Vext, rij );
        x[i] = x[i] + dx;
        for(j=0; j < N; j++)
        {
            DM[3*i][3*j]   = -(Fxip[j]-Fxim[j]) / (2*dx); 
            DM[3*i][3*j+1] = -(Fyip[j]-Fyim[j]) / (2*dx); 
            DM[3*i][3*j+2] = -(Fzip[j]-Fzim[j]) / (2*dx); 
        }
        
        /*  displace y-coordinate of each atom by +dx and -dx and get central difference of forces on all atoms */
        
        y[i] = y[i] + dx;
        forces_m(x, y, z, N, vx, vy, vz,
                 alpha, beta, L, V0,
                 Fxip, Fyip, Fzip, Fxe, Fye, Fze, &Vint, &Vext, rij );
        y[i] = y[i] - 2*dx;
        forces_m(x, y, z, N, vx, vy, vz,
                 alpha, beta, L, V0,
                 Fxim, Fyim, Fzim, Fxe, Fye, Fze, &Vint, &Vext, rij );
        y[i] = y[i] + dx;
        for(j=0; j < N; j++)
        {
            DM[3*i+1][3*j]   = -(Fxip[j]-Fxim[j]) / (2*dx); 
            DM[3*i+1][3*j+1] = -(Fyip[j]-Fyim[j]) / (2*dx); 
            DM[3*i+1][3*j+2] = -(Fzip[j]-Fzim[j]) / (2*dx); 
        }
        
        /*  displace z-coordinate of each atom by +dx and -dx and get central difference of forces on all atoms */
        
        z[i] = z[i] + dx;
        forces_m(x, y, z, N, vx, vy, vz,
                 alpha, beta, L, V0,
                 Fxip, Fyip, Fzip, Fxe, Fye, Fze, &Vint, &Vext, rij );
        z[i] = z[i] - 2*dx;
        forces_m(x, y, z, N, vx, vy, vz,
                 alpha, beta, L, V0,
                 Fxim, Fyim, Fzim, Fxe, Fye, Fze, &Vint, &Vext, rij );
        z[i] = z[i] + dx;
        for(j=0; j < N; j++)
        {
            DM[3*i+2][3*j]   = -(Fxip[j]-Fxim[j]) / (2*dx); 
            DM[3*i+2][3*j+1] = -(Fyip[j]-Fyim[j]) / (2*dx); 
            DM[3*i+2][3*j+2] = -(Fzip[j]-Fzim[j]) / (2*dx); 
        }
    }
 
    
    /*   TEST:
          write out all the elements of the dynamical matrix  */    
    /*    
     matrix_out = fopen("dyn_mat", "w");
     for(j=0; j<3*N; j++)
     {
     for(i=0; i<3*N; i++)
     {
     fprintf(matrix_out, "%lf \n", (DM[j][i]+ DM[i][j])/2 );
     }
     }
     fclose(matrix_out);
     */    
    /*   END OF TEST  */
 
    
    /*  diagonalize the dynamical matrix DM  */
    
    for(j=0; j < 3*N; j++)
    {
        for(i=0; i < j; i++)
        {
            A[i][j] = (DM[j][i]+ DM[i][j])/2 ;   A[j][i] = 0.0;
        }
        A[j][j] = DM[j][j];
    }
    
    n_dim = 3*NMAX;
    n_def = 3*N;
    diagonalize_(&n_dim, &n_def, A, D, ZR, ZI, E, TAU, GRA, GIA);      /* call FORTRAN subroutine "diagonalize" */
    
    if(iscale==0){
        evals = fopen("e-vals.txt", "w");
        for(i=0; i<3*N; i++){
            fprintf(evals, "%lf\n", D[i]);
        }
        fclose(evals);
    }
}



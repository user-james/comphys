#include <stdio.h>
#include <math.h>

#define NMAX 100



/************************************************************************************************/
double  pot_energy_PBC(double x[], double y[], double z[], int N, 
            double alpha, double beta, double L,
            double rij[][NMAX][4] )
{
/*  
   This function calculates the potential energy of N particles at positions
   (x,y,z) for the Lennard-Jones interaction,
                   V(r) = alpha/r^12 - beta/r^6 . 

   with periodic boundary conditions on a cube with
     -L/2 < x,y,z < L/2 
   including the six periodic images of each particle at (x+-L, y, z), (x, y+-L, z), (x, y, z+-L)

Input variables:
  (x,y,z) =  positions of particles (cartesian coords)
  N       =  number of particles
  alpha   =  coefficient of 1/r^12
  beta    =  coefficient of -1/r^6 
  L       =  size of confining "box"

Output Variables:
 
  rij   =  interatomic distances and vectors

function returns V  =  total potential energy of system 

*/
  double xij, yij, zij;
  double d2, d6, V;
  int i,j,k;

  if( N > NMAX )
  {
      printf("%d is too many particles for energy_PBC. NMAX = %d\n", N, NMAX);
      return 1;
  }
    
  V = 0.0; 
  for(i=0; i<N; i++)
  {
      rij[i][i][0] = 0.0;
      rij[i][i][1] = 0.0;
      rij[i][i][2] = 0.0;
      rij[i][i][3] = 0.0;
      
      for(j=0; j<i; j++)
      {
          xij = x[j] - x[i];     /*  relative position vectors  */
          yij = y[j] - y[i];
          zij = z[j] - z[i];
          if( fabs(xij) > L/2 ) xij = xij - L*xij/fabs(xij);    /* find closest image in PBC  */
          if( fabs(yij) > L/2 ) yij = yij - L*yij/fabs(yij);
          if( fabs(zij) > L/2 ) zij = zij - L*zij/fabs(zij);
          
          d2 = xij*xij + yij*yij + zij*zij; 
          rij[j][i][0] =  ( rij[i][j][0] = sqrt(d2) );
          rij[j][i][1] = -( rij[i][j][1] = xij );
          rij[j][i][2] = -( rij[i][j][2] = yij );
          rij[j][i][3] = -( rij[i][j][3] = zij );
          d6 = d2*d2*d2; 
          V += (alpha/d6 - beta)/d6;    
      }
  }
    
    return V;
}

/************************************************************************************************/
double  kin_energy_PBC(double x[], double y[], double z[], int N, 
                       double a_1, double a_2, double L,
                       double rij[][NMAX][4] )
{
    /*  
     This function calculates the kinetic energy term for the McMillan Boson wave function psi,
     for which
                      log(psi) = - sum_{i<j}  (a_1/r_ij)^a_2,
     
     for N particles at positions (x,y,z) with periodic boundary conditions on a cube, where
     -L/2 < x,y,z < L/2. 
     The kinetic energy is of the form:
     
                  T  =  sum_i { 2*T_i - |F_i|^2 },
     where
          T_i  =  -grad_i^2 log(psi)/4    =  sum_{j != i}  -1/(2r_ij) d^2/dr_ij^2 { r (a_1/r_ij)^a_2 } 
          F_i  = grad_i log(psi)/sqrt(2)  =  sum_{j != i}  grad_i { (a_1/r_ij)^a_2 } / sqrt(2) .
     
     
 Input variables:
     (x,y,z) =  positions of particles (cartesian coords)
     N       =  number of particles
     a_1     =  1st coefficient of pair-wise term
     a_2     =  2nd coefficient (exponent) of pair-wise term
     L       =  size of confining "box"
     rij     =  interatomic distances and vectors (calculated in potential energy function, pot_energy_PBC)
     
 Output Variables:
     
     function returns T  =  total kinetic energy of system 
     
     */
    double xij, yij, zij, r;
    double T, Fij, Fxi, Fyi, Fzi, Ti;
    int i,j;
    
    if( N > NMAX )
    {
        printf("%d is too many particles for energy_PBC. NMAX = %d\n", N, NMAX);
        return 1;
    }
    
    T = 0.0; 
    for(i=0; i< N; i++)
    {
        Ti=0.0;   Fxi = 0.0; Fyi = 0.0; Fzi = 0.0;
        
        for(j=0; j<N; j++)
        {
            if(j != i)
            {
                r = rij[i][j][0];
                xij = rij[i][j][1] / r;    /*   direction cosines   */
                yij = rij[i][j][2] / r;
                zij = rij[i][j][3] / r;
                
                Fij = pow((a_1/r), a_2) * (a_2/r);
                Fxi += -xij * Fij;
                Fyi += -yij * Fij;
                Fzi += -zij * Fij;
                
                Ti += (a_2-1.0) * Fij / (4*r);
            }
        } 
        T += 2*Ti - (Fxi*Fxi + Fyi*Fyi + Fzi*Fzi) / 2;
    }    
    return T;
}



/************************************************************************************************/
double  wfn_log_ratio_PBC(double x[], double y[], double z[], int N, 
                          double a_1, double a_2, double L,
                          int i, double x_new, double y_new, double z_new)
{
    /*  
     This function calculates the change of the log of the McMillan pair-wise product Boson wave 
     function, 
                   log(psi) = - sum_{i<j}  (a_1/r_ij)^a_2

     for N particles at positions (x,y,z) with periodic boundary conditions on a cube, where
                    -L/2 < x,y,z < L/2, 
     if particle i were to be moved to position (x_new,y_new,z_new).
           
     
Input variables:
     (x,y,z) =  current positions of particles (cartesian coords)
     N       =  number of particles
     a_1     =  first parameter of McMillan function (see defn above)
     a_2     =  second parameter of McMillan function (see defn above)
     L       =  size of confining "box"
     i       =  index of particle to be moved
     (x_new, y_new, z_new) = new (trial) position of particle i

     
     Output Variables:
     
     function returns dV  =  change in  sum_{i<j} (a_1/r_ij)^a_2    
                             when particle i is moved from (x[i],y[i],z[i]) to (x_new,y_new,z_new) 
     
     */
    double xij, yij, zij;
    double r, V, V_old;
    int j;
    
    if( N > NMAX )
    {
        printf("%d is too many particles for energy_change_PBC. NMAX = %d\n", N, NMAX);
        return 1;
    }    
    
      /*  sum the pair-wise energy terms for particle i at its current position   */
    
    V = 0.0; 
    for(j=0; j<N; j++)
    {
        if(j != i)
        {
            xij = fabs(x[i] - x[j]);     /*  relative position vectors  */
            yij = fabs(y[i] - y[j]);
            zij = fabs(z[i] - z[j]);
            if( xij > L/2 ) xij = L - xij; /* find closest image in PBC  */
            if( yij > L/2 ) yij = L - yij;
            if( zij > L/2 ) zij = L - zij;
            
            r = sqrt(xij*xij + yij*yij + zij*zij); 
            V +=  pow( (a_1/r), a_2 ) ;   
        }    
    }
    V_old = V; 
    
    /*  sum the pair-wise energy terms for particle i at its new position   */
        
    V = 0.0;     
    for(j=0; j<N; j++)
    {
        if(j != i)
        {
            xij = fabs(x_new - x[j]);     /*  relative position vectors  */
            yij = fabs(y_new - y[j]);
            zij = fabs(z_new - z[j]);
            if( xij > L/2 ) xij = L - xij; /* find closest image in PBC  */
            if( yij > L/2 ) yij = L - yij;
            if( zij > L/2 ) zij = L - zij;
            
            r = sqrt(xij*xij + yij*yij + zij*zij); 
            V +=  pow( (a_1/r), a_2 ) ; 
       }    
    }
        return V-V_old;
}


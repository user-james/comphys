#include <stdio.h>
#include <math.h>

#define NMAX 100         /* dimension of interatomic distance array */
#define RMAX 3.0         /* largest interatomic distance that will be binned */
#define NHIST 100        /* number of bins in histogram */


FILE *outhist; 

/*************************************************************************/
int rho1_PBC(double x[], double y[], double z[], int N, double a_1, double a_2, double L, int key)
{
    /* this function accumulates statistics of the single-particle density function 
     for use with the variational quantum Monte Carlo code
     
     Statistics are accumulated for the density function when rho1_PBC is called with key=1
          
     Several histograms can be output during a simulation run: 
        each time rho1_PBC is called with key=2, a new output file "rho1_XXX" is generated and XXX is incremented by 1.
     
     The histogram is printed (for key=2) to the output file: 
           each line of output contains "distance, single-particle density matrix value".
           output file name = "rho1_XXX", where XXX = i_file. 
           i_file is incremented each time bin_rij is called with key=2
     
     The histogram can be reset to zero by calling rho1_PBC with key=0
     
     
     INPUT:
         x[i],y[i],z[i]  = position of particle i
         N    =  number of particles
         a_1  =  first parameter in variational wave function
         a_2  =  second parameter in variational wave function
         L    =  size of box for periodic boundary conditions
     
         key = 1 if called to accumulate statistics in histogram
             = 0 if called to reset histogram to zero
             = 2 if called to print out histogram
     
     
     OUTPUT:    (note these are static variables, saved between calls to this function) 
     
         n_sample  =  number of cluster samples included in constructing histogram 
                      set = 0 for key=0, incremented by 1 for key=1, unchanged for key=2
     
         rho1_hist[i] = sum of values of psi(R')/psi(R) for one particle in R' between i*RMAX/NHIST and (i+1)*RMAX/NHIST from particle in R
                       set = 0 for key=0, incremented for key=1, output (per sample) to file for key=2
     
         N1_hist[i] = number of alternate positions at distance between i*RMAX/NHIST and (i+1)*RMAX/NHIST from original particle position
         large_dist = number of pairs occurring with distance beyond histogram range
         
     */

    
    static double rho1_hist[NHIST];     /* accumulation of psi(R')/psi(R) for |r_i-r_i'| in given histogram bin */
    static long int N1_hist[NHIST];     /* number of samples with |r_i-r_i'| in given histogram bin */
    static long int large_dist;     /* number of atom pairs occurring at distance larger than histogram range */
    static int n_sample;            /* number of samples already accumulated in current histogram */
    static long int iseed = 0;

    static int i_file=0;            /* number of histograms already printed */
    
    int i, i_hist;
    char filename[25];              /* name of histogram output file */
    
    double x_new, y_new, z_new;
    double dx, dy, dz, r;
    
    double dran1_(long int *);
    double  wfn_log_ratio_PBC(double *, double *, double *, int, 
                              double, double, double,
                              int, double, double, double);
     
    switch (key)    
    {
        case 1:                    /* if key = 1,  accumulate statistics */
            
            for(i=0; i<N; i++)
            {
                /* define new position of particle i  */
                
                x_new = L*(dran1_(&iseed)-0.5); 
                y_new = L*(dran1_(&iseed)-0.5); 
                z_new = L*(dran1_(&iseed)-0.5); 
            
                dx = fabs(x_new - x[i]);     /*  relative position vectors  */
                dy = fabs(y_new - y[i]);
                dz = fabs(z_new - z[i]);
                if( dx > L/2 ) dx = L - dx; /* find closest image in PBC  */
                if( dy > L/2 ) dy = L - dx;
                if( dz > L/2 ) dz = L - dx;
                r = sqrt(dx*dx + dy*dy + dz*dz);
 
                i_hist = NHIST * r / RMAX;
                if(i_hist < NHIST) 
                {
                    rho1_hist[i_hist] += exp(-wfn_log_ratio_PBC(x,y,z, N, a_1,a_2, L, i, x_new,y_new,z_new));  
                    N1_hist[i_hist] += 1;
                }
                else large_dist += 1;
            }
            n_sample += 1;
            return 1;


        case 0:                   /* if key = 0, reset histogram to zero */

            for(i_hist=0; i_hist < NHIST; i_hist++)
            {
                rho1_hist[i_hist] = 0.0;
                N1_hist[i_hist] = 0;
            }
            large_dist = 0;
            n_sample = 0;
            return 0;
   
        case 2:                   /* if key = 2, output histogram to file  */
    
            sprintf(filename, "./data/rho1_%d", i_file);
            i_file += 1;
            outhist = fopen(filename, "w");
        
            for(i_hist=0; i_hist<NHIST; i_hist++)
            {
                fprintf(outhist, "%lf\t%lf\n", (i_hist*RMAX)/NHIST, rho1_hist[i_hist]/N1_hist[i_hist]) ;
            }
        
            fclose(outhist);
            return 2;
    }
    return 3;
}
           
           
           
           
           
           
           
           
           

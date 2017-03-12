#include <stdio.h>
#include <math.h>

#define NMAX 100         /* dimension of interatomic distance array */
#define RMAX 7.0         /* largest interatomic distance that will be binned */
#define NHIST 700        /* number of bins in histogram */


FILE *outhist; 

/*************************************************************************/
int bin_rij(double rij[][NMAX], int N, int key)
{
    /* this function accumulates statistics of the atomic pair-correlation function 
     for use with the Lennard-Jones molecular dynamics code
     
     Statistics are accumulated for the pair-distribution when bin_rij is called with key=1
          
     Several histograms can be output during a simulation run: 
        each time bin_rij is called with key=2, a new output file "pair_distributionXXX" is generated and XXX is incremented by 1.
     
     The histogram is printed (for key=2) to the output file: 
           each line of output contains "distance, pair distribution".
           output file name = "pair_distributionXXX", where XXX = i_file. 
           i_file is incremented each time bin_rij is called with key=2
     
     The histogram can be reset to zero by calling bin_rij with key=0
     
     
     INPUT:
         rij = current interatomic distances
         N  =  number of particles
     
         key = 1 if called to accumulate statistics in histogram
             = 0 if called to reset histogram to zero
             = 2 if called to print out histogram
     
     
     OUTPUT:    (note these are static variables, saved between calls to this function) 
     
         n_sample  =  number of cluster samples included in constructing histogram 
                      set = 0 for key=0, incremented by 1 for key=1, unchanged for key=2
     
         rij_hist[i] = number of atom pairs at distance between i*RMAX/NHIST and (i+1)*RMAX/NHIST
                       set = 0 for key=0, incremented for key=1, output (per sample) to file for key=2
     
         large_dist = number of pairs occurring with distance beyond histogram range
         
     */

    
    static long int rij_hist[NHIST];     /* number of atom pairs occurring in given histogram bin */
    static long int large_dist;     /* number of atom pairs occurring at distance larger than histogram range */
    static int n_sample;            /* number of samples already accumulated in current histogram */

    static int i_file=0;            /* number of histograms already printed */
    
    int i,j, i_hist;
    char filename[25];              /* name of histogram output file */
    
     
    switch (key)    
    {
        case 1:                    /* if key = 1,  accumulate statistics */
            
            for(i=0; i<N; i++)
            {
                for(j=0; j<i; j++)
                {
                    i_hist = NHIST * rij[i][j] / RMAX;
                    if(i_hist < NHIST) rij_hist[i_hist] += 1;
                    else large_dist += 1;
                }
            }
            n_sample += 1;
            return 1;


        case 0:                   /* if key = 0, reset histogram to zero */

            for(i_hist=0; i_hist < NHIST; i_hist++)
            {
                rij_hist[i_hist] = 0;
            }
            large_dist = 0;
            n_sample = 0;
            return 0;
   
        case 2:                   /* if key = 2, output histogram to file  */
    
            sprintf(filename, "pair_distribution_%d", i_file);
            i_file += 1;
            outhist = fopen(filename, "w");
        
            for(i_hist=0; i_hist<NHIST; i_hist++)
            {
                fprintf(outhist, "%lf\t%lf\n", (i_hist*RMAX)/NHIST, (double)rij_hist[i_hist]/(double)n_sample );
            }
            printf("pair fraction outside histogram: %lf\n",(double)large_dist/((double)n_sample*N*(N-1)) );
            fclose(outhist);
            return 2;
    }
    return 3;
}
           
           
           
           
           
           
           
           
           

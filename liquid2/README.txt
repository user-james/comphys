This directory contains programs related to the Monte Carlo simulation of the quantum ground state of 4He, using the Lennard-Jones interaction. The ground state wave function is a product of positive, pair-wise terms, as used by McMillan.

PROGRAMS

  distance_histogram_PBC_quantum.c
    C function for accumulating statistics of the pair-distribution
    function, with periodic boundary conditions on a cube
    similar to function used in classical liquid simulation, but including full vector in terms r_ij

  dran1.f
    FORTRAN function for generating pseudo-random numbers

  energy_PBC_quantum.c
    contains C functions for calculating the potential energy (pot_energy_PBC) and kinetic energy (kin_energy_PBC), 
    with periodic boundary conditions on a cube;
    also contains C function to find the log of the ratio of the wave functions 
    for new and old positions of one of the particles (wfn_log_ratio_PBC)
    

  metropolis_quantum.c
    Main C code for doing simulation of classical liquid

  random_positions.c
    C function to set up initial random positions of particles

  random_walk_quantum.c
    C function to run a number of steps of the random walk
    in the Metropolis method. Returns the average potential energy.

  rho1_histogram_PBC.c
    C function to calculate a histogram of the single-particle density matrix rho(|r-r'|).
    with period boundary conditions on a cube,
    similar to distance_histogram_PBC_quantum.c
    
COMPILATION

use gfortran on est1 systems:
	gfortran metropolis_quantum.c random_positions.c random_walk_quantum.c energy_PBC_quantum.c distance_histogram_PBC_quantum.c dran1.f rho1_histogram_PBC.c -o metropolis

OUTPUT
execution  (./metropolis_quantum) will produce 
    a plot of energy vs variational parameter a1 (and associated text output file)
    pair-distribution function histograms for each value of a1
    single-particle density matrix histograms for each value of a1

METHOD:
Program (metropolis_quantum.c) generates for each value of the variational parameter a1 
(after an initial warm-up stage) a sequence of configurations, R = {r1,r2,r3,…, rN}, 
of the system that are distributed according to the wave function, 
   P(R )  =  exp[- 2 sum_{i<j} (a1/r_ij)^a2 ] / Z    
where Z is the normalization factor. (We take a2 = 5.)
The pair-distribution function histogram is generated during the walk and the average potential energy is calculated.
The single-particle density matrix
rho(r,r') = int Psi(r,r2,r3,r4,…,rN) Psi(r, r2,r3,r4,…,rN) dr2 dr3…drN
is also averaged as a function of |r-r'|. 
The large distance limit of rho(|r-r'|) equals the fraction of particles in the Bose condensate.

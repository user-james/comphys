This directory contains programs related to the Monte Carlo simulation of the classical thermal state of 4He, using the Lennard-Jones interaction
with periodic boundary conditions.

PROGRAMS

  distance_histogram_PBC.c
    C function for accumulating statistics of the pair-distribution
    function, with periodic boundary conditions on a cube

  dran1.f
    FORTRAN function for generating pseudo-random numbers

  energy_PBC.c
    C function for calculating the potential energy, 
    with periodic boundary conditions on a cube

  metropolis.c
    Main C code for doing simulation of classical liquid

  random_positions.c
    C function to set up initial random positions of particles

  random_walk.c
    C function to run a number of steps of the random walk
    in the Metropolis method. Returns the average potential energy.
    
COMPILATION

use gfortran on est1 systems:
  gfortran metropolis.c random_positions.c random_walk.c energy_PBC.c distance_histogram_PBC.c dran1.f -o metropolis


OUTPUT
execution  (./metropolis) will produce 
    a plot of energy vs temperature (and associated text output file)
    pair-distribution function histograms for each temperature in simulation

Method:
Program (metropolis.c) generates for each temperature (after an initial warm-up stage) a sequence 
of configurations, R = {r1,r2,r3…, rN}, of the system that are distributed according to the 
Boltzmann distribution, 
   P(R ) = exp[-U(R)/kT] / Z    
where Z is the parition function.
The pair-distribution function histogram is generated during the walk and the average potential energy is calculated.

This directory and its sub-directories contains programs for solving Assignments 4 and 5 in the 
PY4109 Computational Physics course at UCC.

PROGRAMS:

 diagonalize.f:
    Fortran code to diagonalize a Hermitian matrix

 eispc.f
    Standard linear algebra routines used by diagonalize.f

 eispt.f
     TQL2 routine for use in finite difference method 

 grid.f
     Fortran code to generate a non-uniform radial grid for use in finite-difference solver

SUBDIRECTORIES:

  ./box_wave_functions/
     code for a 1-d particle in a square well  [assignment 4(a)]

  ./NFE_wave_functions
      code for a 1-d particle in a periodic potential  [assignment 4(b)]
      (NFE stands for Nearly Free Electron method - as standard method in band structure theory)

  ./H_wave_functions/
     code to solve the radial Schrodinger equation for the Hydrogen atom  [assignment 5(a & b)]

  ./LJ_wave_functions/
     code for quantum vibrations of a Lennard-Jones dimer  [assignment 5c]

   ./harmonic_wave_functions/
     code for a 1-d simple harmonic oscillator




Course Syllabus extract  -------------------------------------------------------

2. NUMERICAL QUANTUM MECHANICS 

2a. Single-particle:  (weeks 5 & 6)
   1-D (quantum well, etc) - grid (coordinate) representation
   periodic potentials:  - Fourier (momentum) representation
       
  Assignments: 
      (4) 1-D problems - coordinate and momentum representations
          (a) 1-D particle in a square well - compare numerical with analytical soln 
          (b) 1-D periodic potential - find band gaps and compare with 2-wave approximation ]    
 
      (5) atomic wave-functions and molecular vibrations He2, Ne2, Ar2  -
          (a) H atom - compare with analytical soln for n=1 (l=0) 
          (b) plot excited states n=2, l=0,1; n=3, l=0,1,2
          (c) vibration states of Ar2, Ne2, He2 with Lennard-Jones interaction  
      

     
This directory contains code required to numerically solve the 1-D Schrodinger equation, using the coordinate
representation, for a particle confined in a box.

The boundary conditions are that the wavefunction is zero at the walls of the box, r=r_min and r=r_max.
The potential energy function is zero within the box.


PROGRAM:

box_wavefunctions.c
   This code is the main driver and should be compiled with ../eispt.f in the parent directory
   with the command:
          gfortran box_wavefunctions.c  ../eispt.f  -o  box_wavefunctions
   to produce an executable file, 'box_wavefunctions'
   The Schrodinger equation is solved using the finite-difference expansion for the grad^2 operator on a uniform grid.
   See header comments in this file for a more precise specification of the method.


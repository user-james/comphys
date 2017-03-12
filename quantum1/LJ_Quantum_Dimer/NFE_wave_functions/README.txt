This directory contains code required to numerically solve the 1-D Schrodinger equation, using the momentum
representation, for a particle in a periodic potential,

    V(x) =  V1 cos(2 pi x/a) + V2 cos(4 pi x/a),

where a is the period of the potential (corresponding to the lattice constant in a crystal).
V1 and V2 are adjustable parameters, which give the strength of the first and second Fourier components of the potential.

Codes:

NFE_wavefunctions.c
   This code is the main driver and should be compiled with ../diagonalize.c and ../eispc.f in the parent directory
   with the command:
          gfortran NFE_wavefunctions.c ../diagonalize.f ../eispc.f  -o  NFE_wavefunctions
   to produce an executable file, 'NFE_wavefunctions'
   See header comments in this file for a more precise specification of the method.

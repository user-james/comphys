# CALCULATES THE FRACTION OF PARTICLES IN THE BOSE-EINSTEIN CONDENSATE
#
# PLOTS THE PAIR DISTRIBUTION FUNCTION GIVEN THE PAIR DISTRIBUTION DATA
# ONCE THAT PLOT IS CLOSED IT PLOT THE SINGLE-PARTICLE DENSITY MATRIX 


import numpy as np
import matplotlib.pyplot as plt



if __name__ == '__main__':

    # we know the index where the energy is minimised is 6
    # this was found from the program, a1_plot.py
    ind = 6 
    
    # grabs data from relevant file
    data = np.genfromtxt("./data/pair_distribution_{:d}".format(ind))
    single_part_mat = np.genfromtxt("./data/rho1_{:d}".format(ind))

    freq = data[:, 1]
    sep1 = data[:, 0]

    particle_dens = single_part_mat[:, 1]
    sep2 = single_part_mat[:, 0]

    # estimates fraction of particles in Bose-Einstein condesate by taking an average
    # of the single-particle density once it flattens out (20% away from the end)
    start_pt = int(len(particle_dens)*.8)       # start_pt is where we take an average from, here it's chosen as 0.8 the entire plot range
    BE_frac = sum(particle_dens[start_pt:])/len(particle_dens[start_pt:])       #average taken to find fraction

    print("Fraction of atoms in Bose-Einstein condesate =", BE_frac)

    x = np.linspace(0, 3, 100)

    # plots pair distribution function and saves it
    plt.title("Pair-Distribution Function")
    plt.ylabel("Fraction of Particles")
    plt.xlabel("Separation")
    plt.plot(sep1, freq, 'k')
    #plt.savefig("./plots/pair_dist")
    plt.show()

    # plots single-particle density matrix and saves it
    plt.title("Single-Partcle Density Matrix")
    plt.ylabel("Fraction of Particles")
    plt.xlabel("Separation")
    plt.plot(sep2, particle_dens, 'k.')
    plt.plot(x, np.ones(len(x))*BE_frac, 'r--')
    #plt.savefig("./plots/density_matrix")
    plt.show()

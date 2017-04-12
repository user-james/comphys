# CALCULATES BULK MODULUS
# SEE SCRIPT 'THERMAL_EXP' AS IT WORKS IN A VERY SIMILAR WAY


import matplotlib.pyplot as plt
import numpy as np


if __name__=='__main__':


    # sets up arrays for densities, lengths and temperatures used
    lengths = np.arange(4.47, 4.15,-0.05)
    densities = 80*lengths**(-3)
    temperatures = np.genfromtxt("./data/length_4.17/energy_vs_temp")[:, 0]
    data = []
    
    for i in range(len(temperatures)):
        energies = []

        # grabs energies for different densities given the temperature
        for l in lengths:
            energies.append(np.genfromtxt("./data/length_{:.2f}/energy_vs_temp".format(l))[i])
        energies = np.array(energies)
        energies = energies[:, 1]
        min_e = energies[0]
        ind = 0

        # calculates the minimum energy for this temperature
        # records minimum energy and index where it occurs
        for j in range(len(energies)):
            if(energies[j] < min_e):
                min_e = energies[j]
                ind = j

        # stores minimum energy and the density where it occurs for the given temperature
        data.append([densities[ind], min_e])

    
    data = np.array(data)
    # converts densities to volumes
    data[:, 0] = 80/data[:, 0]
    
    # calculates bulk modulus
    bulk_mods = data[:, 1] / data[:, 0]
    avg_bulk_mod = abs(sum(bulk_mods) / len(bulk_mods))
    print("Bulk Modulus =", avg_bulk_mod)
    
    

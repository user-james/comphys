# CALCULATES THE THERMAL EXPANSIVITY AND PLOTS VOLUME VS ENERGY


import matplotlib.pyplot as plt
import numpy as np


if __name__=='__main__':

    # creating arrays for lengths, densities and temperatures used
    lengths = np.arange(4.47, 4.15,-0.05)
    densities = 80*lengths**(-3)
    temperatures = np.genfromtxt("./data/length_4.17/energy_vs_temp")[:, 0]
    thermal_exp = []
    
    for i in range(len(temperatures)):
        energies = []

        # grabs energies for each denisty given the temperature
        for l in lengths:
            energies.append(np.genfromtxt("./data/length_{:.2f}/energy_vs_temp".format(l))[i])
        energies = np.array(energies)
        energies = energies[:, 1]
        min_e = energies[0]
        ind = 0

        # records minimum energy and index at which it occurs
        # this is the minimum energy across different densities for a single temperature
        for j in range(len(energies)):
            if(energies[j] < min_e):
                min_e = energies[j]
                ind = j
        
        # stores density and temperature for point where density vs energy graph is minimum
        thermal_exp.append([temperatures[ind], densities[ind]])

    
    thermal_exp = np.array(thermal_exp)
    
    # calculates the slope of the density vs energy plot and thus finds the thermal expansivity
    slope = -(thermal_exp[-1, 1] - thermal_exp[0, 1]) / (thermal_exp[-1, 0] - thermal_exp[0, 0])
    expansivity = 2*slope / (thermal_exp[-1, 1] + thermal_exp[0, 1])
    print("Thermal Expansitivty =", expansivity)

    # converts stored densities to volumes
    thermal_exp[:, 1] = 80/thermal_exp[:, 1]


    # formats and plots volume vs temperature
    plt.xlabel("Temperature")
    plt.ylabel("Volume")
    plt.plot(thermal_exp[:, 0], thermal_exp[:, 1])
    #plt.savefig("./plots/vol_vs_temp.jpg")
    plt.show()

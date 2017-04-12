# PLOTS ENERGY VS DENSITY FOR DIFFERENT TEMPERATURES


import matplotlib.pyplot as plt
import numpy as np


if __name__=='__main__':

    lengths = np.arange(4.47, 4.15, -0.05)      # creates array filled with lengths used in sim
    densities = 80*lengths**(-3)                # create array of densities used
    temperatures = np.genfromtxt("./data/length_4.17/energy_vs_temp")[:, 0]    # grabs temperature data for one density (temperatures unchanged for each denisty)

    
    for i in range(len(temperatures)):
        energies = []

        # for each temperature grabs different energies for each density
        for l in lengths:
            energies.append(np.genfromtxt("./data/length_{:.2f}/energy_vs_temp".format(l))[i])


        # formats and saves graphs for each temperature
        energies = np.array(energies)
        plt.xlabel("Density")
        plt.ylabel("Energy")
        plt.title("Energy vs Density at T = {:.2f}".format(temperatures[i]))
        plt.plot(densities, energies[:, 1])
        #plt.savefig("./plots/energy_vs_density_{:.2f}.jpg".format(temperatures[i]))       # uncomment to save graphs
        plt.close()

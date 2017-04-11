import matplotlib.pyplot as plt
import numpy as np


if __name__=='__main__':

    lengths = np.arange(4.47, 4.15, -0.05)
    densities = 80*lengths**(-3)
    temperatures = np.genfromtxt("length_4.17/energy_vs_temp")[:, 0]

    
    for i in range(len(temperatures)):
        energies = []
        for l in lengths:
            energies.append(np.genfromtxt("length_{:.2f}/energy_vs_temp".format(l))[i])
        energies = np.array(energies)
        plt.xlabel("Density")
        plt.ylabel("Energy")
        plt.title("Energy vs Density at T = {:.2f}".format(temperatures[i]))
        plt.plot(densities, energies[:, 1])
        #plt.savefig("../plots/energy_vs_density_{:.2f}.jpg".format(temperatures[i]))
        plt.close()

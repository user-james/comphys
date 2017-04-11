import matplotlib.pyplot as plt
import numpy as np


if __name__=='__main__':

    lengths = np.arange(4.47, 4.15,-0.05)
    densities = 80*lengths**(-3)
    temperatures = np.genfromtxt("length_4.17/energy_vs_temp")[:, 0]
    thermal_exp = []
    
    for i in range(len(temperatures)):
        energies = []
        for l in lengths:
            energies.append(np.genfromtxt("length_{:.2f}/energy_vs_temp".format(l))[i])
        energies = np.array(energies)
        energies = energies[:, 1]
        min_e = energies[0]
        ind = 0
        for j in range(len(energies)):
            if(energies[j] < min_e):
                min_e = energies[j]
                ind = j

        thermal_exp.append([temperatures[ind], densities[ind]])

    
    thermal_exp = np.array(thermal_exp)
    
    slope = -(thermal_exp[-1, 1] - thermal_exp[0, 1]) / (thermal_exp[-1, 0] - thermal_exp[0, 0])
    expansivity = 2*slope / (thermal_exp[-1, 1] + thermal_exp[0, 1])

    thermal_exp[:, 1] = 1/thermal_exp[:, 1]

    print("Thermal Expansitivty =", expansivity)
    plt.xlabel("Temperature")
    plt.ylabel("Volume")
    plt.plot(thermal_exp[:, 0], thermal_exp[:, 1])
    #plt.savefig("../plots/vol_vs_temp.jpg")
    plt.show()

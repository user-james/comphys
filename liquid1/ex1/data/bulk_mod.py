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

        thermal_exp.append([densities[ind], min_e])

    
    thermal_exp = np.array(thermal_exp)
    thermal_exp[:, 0] = 80/thermal_exp[:, 0]
    
    bulk_mods = thermal_exp[:, 1] / thermal_exp[:, 0]
    avg_bulk_mod = abs(sum(bulk_mods) / len(bulk_mods))
    print(avg_bulk_mod)
    print(thermal_exp)
    
    #print("Bulk Modulus =", expansivity)
    #plt.xlabel("Temperature")
    #plt.ylabel("Volume")
    #plt.plot(thermal_exp[:, 0], thermal_exp[:, 1])
    #plt.savefig("../plots/vol_vs_temp.jpg")
    #plt.show()

import matplotlib.pyplot as plt
import numpy as np


if __name__=='__main__':

    lengths = np.arange(4.47, 4.15, -0.05)
    total = 0
    for l in lengths:
        data = np.genfromtxt("length_{:.2f}/energy_vs_temp".format(l))
        
        temp = data[:, 0]
        energy = data[:, 1]
        
        heat_cap = []

        for i in range(len(temp) - 1):
            heat_cap.append((energy[i+1] - energy[i]) / (temp[i+1] - temp[i]))

        heat_cap = (energy[3] - energy[0]) / (temp[3] - temp[0])
        total += heat_cap
        print("heat capacity =", heat_cap)
        plt.xlabel("Temperature")
        plt.ylabel("Energy")
        plt.plot(temp, energy, label = "Density = {:.2f}".format(80/(l*l*l)))
        
    avg_heat_cap = total / len(lengths)
    print("\nHeat Capacity per Particle =", avg_heat_cap/80)
    plt.legend(loc=2)
    #plt.savefig("../plots/heat_cap.jpg")
    plt.show()

    

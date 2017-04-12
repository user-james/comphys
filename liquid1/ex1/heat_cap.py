# CALCULATES HEAT CAPACITY PER PARTICLE

import matplotlib.pyplot as plt
import numpy as np


if __name__=='__main__':

    # constructs array of different lengths used
    lengths = np.arange(4.47, 4.15, -0.05)
    total = 0

    for l in lengths:
        data = np.genfromtxt("./data/length_{:.2f}/energy_vs_temp".format(l))
        
        temp = data[:, 0]
        energy = data[:, 1]
        
        # calculates heat capacity for each density as the slope of the beginning of the E vs T curve
        heat_cap = (energy[3] - energy[0]) / (temp[3] - temp[0])
        total += heat_cap

        # prints different heat capacities 
        print("heat capacity for density {:.2f} =".format(80/(l)**(3)), heat_cap)
        plt.xlabel("Temperature")
        plt.ylabel("Energy")
        plt.plot(temp, energy, label = "Density = {:.2f}".format(80/(l*l*l)))
        
    # calculates average of heat capacities as final value, divides by number of particles
    avg_heat_cap = total / len(lengths)
    print("\nHeat Capacity per Particle =", avg_heat_cap/80)
    plt.legend(loc=2)
    #plt.savefig("./plots/heat_cap.jpg")
    plt.show()

    

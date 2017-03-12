import numpy as np
import matplotlib.pyplot as plt

if __name__=='__main__':

    energies = np.genfromtxt("energies")

    plt.plot(energies[:, 0], energies[:, 1])
    
    plt.xlabel("Iteration")
    plt.ylabel("Energy")
    plt.savefig("energy_func.jpg")
    plt.show()

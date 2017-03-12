import numpy as np
import matplotlib.pyplot as plt

if __name__=='__main__':

    he = np.genfromtxt("charge_density_iter_8")
    e = np.genfromtxt("../ex1/charge_density_iter_8")


    plt.plot(e[:, 0], e[:, 1], label="DFT")
    plt.plot(he[:, 0], he[:, 1], label="Hartree")
    
    plt.xlim(0, 10)
    plt.xlabel("Radius")
    plt.ylabel("Charge Density")
    plt.legend()
    plt.savefig("charge_density_comp.jpg")
    plt.show()

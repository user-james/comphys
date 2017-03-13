import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

    energies = np.genfromtxt("energies")
    e_list = list(energies[:, 1])
    min_index = e_list.index(min(e_list))

    plt.plot(energies[:, 0], energies[:, 1], label = r"$<E>$")
    plt.plot(energies[min_index, 0], energies[min_index, 1], 'r.', label = r"$<E>_{min}$")
    plt.xlabel(r"$x_0$")
    plt.ylabel(r"$<E(x_0)>$")
    plt.xlim(0.5, 1.5)
    plt.legend(numpoints = 1)
    plt.savefig("new_energy_x0.jpg")
    plt.show()


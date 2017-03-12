import numpy as np
import matplotlib.pyplot as plt


if __name__=='__main__':

    ne = np.genfromtxt("energy_ne")
    ar = np.genfromtxt("energy_ar")



    for i in range(len(ne[:, 0])):
        if ne[i, 1] > 0:
            ne = ne[:i]
            break
    for i in range(len(ar[:, 0])):
        if ar[i, 1] > 0:
            ar = ar[:i]
            break

    plt.plot(ne[:, 0] + 1, ne[:, 1], label = "Ne")
    plt.plot(ar[:, 0]+1, ar[:, 1], label = "Ar")

    plt.legend(loc=4)
    plt.savefig("lj_energies.jpg")
    plt.show()

import numpy as np
import matplotlib.pyplot as plt


if __name__=='__main__':

    files = [0, 1, 2, 6, 7, 8]

    for num in files:

        chg_density = np.genfromtxt("charge_density_iter_{:d}".format(num))

        plt.plot(chg_density[:, 0], chg_density[:, 1], label="{:d}".format(num))

    plt.xlabel("Radius")
    plt.ylabel("Charge Density")
    plt.xlim(0, 10)
    plt.legend()
    plt.savefig("charge_density.jpg")
    plt.show()

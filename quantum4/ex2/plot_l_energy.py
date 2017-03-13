import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    
    x0 = [0.5, 1.0, 1.5, 1.32]
    
    for x in x0:
        local_e = np.genfromtxt("localenergy_x0_{:.2f}".format(x))

        plt.plot(local_e[:, 0], local_e[:, 1], label = r"$x0 = {:.2f}$".format(x))
    
    plt.xlabel("x")
    plt.ylabel("Local Energy")
    plt.xlim(-5, 5)
    plt.ylim(-100, 100)
    x = 1.324431

    plt.legend(numpoints = 1)
    plt.savefig("localenergy.jpg")
    plt.show()


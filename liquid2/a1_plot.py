# PLOTS CURVE OF EXPECTED ENERGY VS A1
# FINDS VALUE OF ENERGY THAT MINIMISES A1
# RETURNS THE MIN ENERGY, THE CORRESPONDING A1 AND THE INDEX AT WHICH IT OCCURS

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    # grabs data and creates variables to be used
    data = np.genfromtxt("./data/energy_vs_a1")
    a1 = data[:, 0]
    energy = data[:, 1]
    min_e = energy[0]
    min_ind = 0

    # finds min energy and correspondin index
    for i in range(1, len(energy)):
        if energy[i] < min_e:
            min_ind = i
            min_e = energy[i]


    # prints results
    print("miniumum expected energy =", min_e)
    print("index at which min energy occurs =", min_ind)
    print("a1 corresponding to minimum energy =", a1[min_ind])

    # plots <E> vs a1 and saves the graph
    plt.ylabel(r'$<E>$')
    plt.xlabel(r'$a_1$')
    plt.plot(a1, energy, 'k')
    plt.plot(a1[min_ind], energy[min_ind], 'ro')
    #plt.savefig("./plots/energy_vs_a1")
    plt.show()

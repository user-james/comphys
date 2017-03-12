import numpy as np
import matplotlib.pyplot as plt


if __name__=='__main__':

    bands = np.genfromtxt("../LJ_Quantum_Dimer/NFE_wave_functions/NFE_bands")
    n_columns = len(bands[0, :])

    for i in range(1, n_columns):

        plt.plot(bands[:, 0], bands[:, i])
    
       
    plt.xlabel(r'$k$')
    plt.ylabel(r'$E(k)$')
    plt.xlim(-0.5, 0.5)
    plt.ylim(-10, 78)
    plt.savefig("NFE_bandgaps.jpg")
    plt.show()

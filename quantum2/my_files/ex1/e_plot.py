import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

    energies = np.genfromtxt("energies")
    
    # separate data into n and E_n
    #n = energies[:, 0]+4        # n = 3, therefore the principal qunatum number is (l + n + 1) = l+4
    n = energies[:, 0]+1        # l = 0, the principal quantum number is (l + n + 1) = n + 1
    e = energies[:, 1]

    # create 2 axes
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    # create xspace over which we plot expected results
    x = np.linspace(1, 6, 100)

    # plot simulated data and expected form 
    ax1.plot(n, e, 'r-', label='Simulated')
    ax2.plot(x, -1/(x*x), 'b-', label=r'$-\frac{1}{n^2}$')
   

    # formatting and plotting
    ax1.set_xlabel("$n$")
    ax2.set_ylabel(r'$\frac{1}{n^2}$')
    ax1.set_ylabel("$E_{n}$")
    ax1.legend(loc=4)
    ax2.legend(loc=3)
    plt.tight_layout()
    plt.savefig("energy_comp.jpg")         # uncomment to save figure
    plt.show()

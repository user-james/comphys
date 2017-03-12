import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

    # simulated energies
    e_data = np.genfromtxt("e_100points")#[9.889365, 39.557363, 89.003699, 158.227886, 247.229238, 356.006875]
    #e_data_x = [1, 2, 3, 4, 5, 6]

    # analytical energy in our units
    real_energy = lambda n: n**2 * np.pi**2
    x = np.linspace(0, 7, 100)


    # plotting 
    plt.plot(e_data[:, 0], e_data[:, 1], 'r.', label="Simulated Results")
    plt.plot(x, real_energy(x), label="Expected Results")

    # formatting and final plot
    plt.title("Simulated vs Actual for 1000 grid Points")
    plt.ylabel(r'$E_{n}$')
    plt.xlabel(r'$n$')
    plt.legend(loc=2, fontsize=11)
    plt.savefig("energy_comp_1000.jpg")
    plt.show()

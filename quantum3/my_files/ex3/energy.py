import numpy as np
import matplotlib.pyplot as plt

if __name__=='__main__':


    ratios = [0.25, 0.50, 0.75, 0.9]
    
    for r in ratios:
        e = np.genfromtxt("./energies_{:.2f}".format(r))


        plt.plot(e[:, 0], e[:, 1], label="Ratio: {:.2f}".format(r))
        
        plt.xlabel("Iteration")
        plt.ylabel("Energy")
    plt.legend()
    plt.savefig("convergence_comp.jpg")
    plt.show()

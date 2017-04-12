import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

       
    density = input("What density are you interested in (high = h/ low = l): ")

    if density == 'h':
        length = 4.17
        figname = "high"
        
    elif density == 'l':
        length = 4.47
        figname = "low"
    else:
        print("Invalid density chosen, defaulting to high density")
        length = 4.17
        figname = "high"



    temperatures = np.arange(0.30, 0.029, -0.03)

    for i in range(len(temperatures)):
        pair_low = np.genfromtxt("./data/length_{:.2f}/pair_distribution_{:d}".format(length, i))
        plt.plot(pair_low[:, 0], pair_low[:, 1], label="temp = {:.2f}".format(temperatures[i]))
    
    
    plt.ylim(0, 160)
    plt.ylabel("Number of pairs")
    plt.xlabel("Atomic Distance")
    plt.title("Pair Distribution for  Density = {:.2f}".format(80*length**(-3)))
    plt.legend(loc = 2, fontsize = 12)
    plt.savefig("./plots/pair_distribution_{:s}".format(figname))
    plt.show()
